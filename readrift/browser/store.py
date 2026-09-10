"""The browser cache: one analysed run, stored columnar and queried by region.

Why this exists
---------------

``olesdata/DRR325755.btop`` is 985 MB.  Reading it, classifying every read and
accumulating the statistics takes minutes.  An interactive view cannot pay that
per pan, so the analysis is written once into a single ``.npz`` archive and every
later session loads it in well under a second.

Layout
------

Columnar, and CSR-style where a row owns a variable number of children:

* **hits** are stored sorted by ``(contig, slow)``.  A region query is then two
  :func:`numpy.searchsorted` calls plus one boolean mask -- no window table, no
  R-tree, no per-query allocation proportional to the genome.
* **reads** index into that sorted array through ``read_hit_idx``, a permutation
  grouped by read, with ``read_hit_off`` / ``read_hit_len`` as the offsets.
* **strings** (read names, labels, annotation names) are one UTF-8 blob plus an
  offset array.  Nothing is decoded until something asks for it, so opening a
  cache of half a million reads does not build half a million Python strings.

Coverage in the depth arrays is *matched bases only* -- decision D3, the same
definition :mod:`readrift.stats` uses -- so the track and the reported scalar
always agree.
"""

from __future__ import annotations

import contextlib
import json
import math
import time
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from readrift.btop_trace import parse as parse_btop
from readrift.inputs.reference import Reference
from readrift.labels import gene_label
from readrift.models import Junction, ReadClass, ReadGroup
from readrift.params import Params
from readrift.stats import Stats

#: Bumped whenever the array layout changes.  An older cache is rebuilt rather
#: than misread.
#:
#: 2 -- annotations carry ``product``, ``locus_tag`` and ``note``, for the
#:      browser's feature detail panel.
FORMAT_VERSION = 2

#: Class and junction orderings.  Stored as small integers; the names live in
#: the metadata so a reader never has to import :mod:`readrift.models`.
CLASS_ORDER: tuple[ReadClass, ...] = (
    ReadClass.UNDIVIDED,
    ReadClass.SHORT_DIVIDED,
    ReadClass.LONG_DIVIDED,
)
CLASS_NAMES: tuple[str, ...] = tuple(c.value for c in CLASS_ORDER)
_CLASS_INDEX = {c: i for i, c in enumerate(CLASS_ORDER)}

JUNCTION_ORDER: tuple[Junction, ...] = (
    Junction.PLAIN,
    Junction.CIRCLE,
    Junction.CONTIG_JOIN,
)
JUNCTION_NAMES: tuple[str, ...] = tuple(j.value for j in JUNCTION_ORDER)
_JUNCTION_INDEX = {j: i for i, j in enumerate(JUNCTION_ORDER)}

EVENT_KINDS: tuple[str, ...] = ("junction", "inversion", "circle", "contig_join")
_EVENT_INDEX = {k: i for i, k in enumerate(EVENT_KINDS)}

#: Finest depth resolution worth storing, and the cap on bins per contig.  A
#: 5 Mb contig lands on 25 bp bins; a 3 Gb one on ~7.5 kb bins and still only
#: 1.6 MB of array.
MIN_DEPTH_BIN = 25
MAX_DEPTH_BINS = 400_000

#: Parameters that change *what the cache contains*.  A run that differs in any
#: of these produces different reads, so the cache must be rebuilt.  Cosmetic
#: options (page scale, line width, label limits) are deliberately absent: they
#: only affect the printed map.
CLASSIFYING_PARAMS: tuple[str, ...] = (
    "min_read_length",
    "min_match_length",
    "min_hsp_length",
    "micro_match_length",
    "division_cut",
    "cov_max",
    "keep_fold_back",
    "compat",
    "allow_unknown_contigs",
    "identity",
)


# --------------------------------------------------------------------------
# String blobs
# --------------------------------------------------------------------------


def _pack_strings(values: Sequence[str]) -> tuple[np.ndarray, np.ndarray]:
    """Pack strings into one UTF-8 buffer plus ``len(values) + 1`` offsets."""
    encoded = [v.encode("utf-8", "replace") for v in values]
    offsets = np.zeros(len(encoded) + 1, dtype=np.int64)
    if encoded:
        offsets[1:] = np.cumsum(np.fromiter((len(e) for e in encoded), dtype=np.int64))
    blob = np.frombuffer(b"".join(encoded), dtype=np.uint8)
    # frombuffer gives a read-only view of a temporary; copy so savez owns it.
    return blob.copy(), offsets


class _StringColumn:
    """Lazy random access into a packed string blob."""

    __slots__ = ("_blob", "_offsets")

    def __init__(self, blob: np.ndarray, offsets: np.ndarray) -> None:
        self._blob = blob
        self._offsets = offsets

    def __len__(self) -> int:
        return max(0, len(self._offsets) - 1)

    def __getitem__(self, index: int) -> str:
        start = int(self._offsets[index])
        end = int(self._offsets[index + 1])
        return self._blob[start:end].tobytes().decode("utf-8", "replace")

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]


# --------------------------------------------------------------------------
# Building
# --------------------------------------------------------------------------


def _depth_bin(length: int) -> int:
    if length <= 0:
        return MIN_DEPTH_BIN
    return max(MIN_DEPTH_BIN, math.ceil(length / MAX_DEPTH_BINS))


def _classifying_snapshot(params: Params) -> dict[str, object]:
    out: dict[str, object] = {}
    for key in CLASSIFYING_PARAMS:
        value = getattr(params, key)
        out[key] = bool(value) if isinstance(value, bool) else value
    return out


def _source_stamp(path: str) -> dict[str, object]:
    """Path, size and mtime -- enough to notice an input has been replaced.

    The path is resolved, so running from a different working directory does
    not look like a different input and trigger a pointless multi-minute
    rebuild.
    """
    p = Path(path)
    # Unreachable on a normal filesystem, but resolve() can raise on a broken
    # mount or a path that vanished between calls.  An unresolved path is still
    # a usable stamp, so there is nothing to do about it.
    with contextlib.suppress(OSError):
        p = p.resolve()
    try:
        stat = p.stat()
    except OSError:
        return {"path": str(p), "size": -1, "mtime": -1.0}
    return {"path": str(p), "size": stat.st_size, "mtime": round(stat.st_mtime, 3)}


def build_store(
    path: str | Path,
    reference: Reference,
    reads: Iterable[ReadGroup],
    stats: Stats,
    params: Params,
    notes: Sequence[str] = (),
) -> Path:
    """Write the cache for one analysed run.  Returns the path written."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)

    contigs = list(reference.contigs)
    contig_index = {c.name: i for i, c in enumerate(contigs)}
    n_contigs = len(contigs)

    read_list = list(reads)
    n_reads = len(read_list)
    n_hits = sum(len(g.hits) for g in read_list)

    # ---- flat columns ----------------------------------------------------

    read_contig = np.zeros(n_reads, dtype=np.int32)
    read_cls = np.zeros(n_reads, dtype=np.int8)
    read_junction = np.zeros(n_reads, dtype=np.int8)
    read_inverted = np.zeros(n_reads, dtype=np.uint8)
    read_qlen = np.zeros(n_reads, dtype=np.int64)
    read_distance = np.zeros(n_reads, dtype=np.int64)

    # The strand of the read's *first* HSP in BTOP order.  Stored explicitly
    # because the hit arrays are re-sorted by position below, which throws that
    # order away -- and it is the reference every inversion is measured against
    # (`ReadGroup.inverted` is `any(h.strand != hits[0].strand)`).  Recovering it
    # from the leftmost hit instead would mark the wrong segment as inverted on
    # any read whose first piece is not its leftmost.
    read_strand = np.zeros(n_reads, dtype=np.int8)

    hit_read = np.zeros(n_hits, dtype=np.int32)
    hit_contig = np.zeros(n_hits, dtype=np.int32)
    hit_slow = np.zeros(n_hits, dtype=np.int64)
    hit_shigh = np.zeros(n_hits, dtype=np.int64)
    hit_qstart = np.zeros(n_hits, dtype=np.int64)
    hit_qend = np.zeros(n_hits, dtype=np.int64)
    hit_strand = np.zeros(n_hits, dtype=np.int8)
    hit_identity = np.full(n_hits, np.nan, dtype=np.float32)

    names: list[str] = []
    labels: list[str] = []

    want_identity = bool(params.identity)
    cursor = 0
    for i, group in enumerate(read_list):
        names.append(group.read)
        labels.append(group.label)
        read_contig[i] = contig_index.get(group.contig, 0)
        read_cls[i] = _CLASS_INDEX[group.cls]
        read_junction[i] = _JUNCTION_INDEX[group.junction]
        read_inverted[i] = 1 if group.inverted else 0
        read_qlen[i] = group.qlen
        read_distance[i] = group.distance
        read_strand[i] = group.strand

        for hit in group.hits:
            hit_read[cursor] = i
            hit_contig[cursor] = contig_index.get(hit.contig, 0)
            hit_slow[cursor] = hit.slow
            hit_shigh[cursor] = hit.shigh
            hit_qstart[cursor] = hit.qstart
            hit_qend[cursor] = hit.qend
            hit_strand[cursor] = hit.strand
            if want_identity and hit.btop:
                trace = parse_btop(hit.btop)
                if trace.aligned_columns:
                    hit_identity[cursor] = trace.identity
            cursor += 1

    # ---- sort hits by (contig, slow) ------------------------------------
    #
    # This ordering is what makes a region query two searchsorted calls.  The
    # read -> hit mapping is rebuilt afterwards as a permutation grouped by
    # read, so both access patterns stay O(result size).

    order = np.lexsort((hit_slow, hit_contig))
    hit_read = hit_read[order]
    hit_contig = hit_contig[order]
    hit_slow = hit_slow[order]
    hit_shigh = hit_shigh[order]
    hit_qstart = hit_qstart[order]
    hit_qend = hit_qend[order]
    hit_strand = hit_strand[order]
    hit_identity = hit_identity[order]

    read_hit_idx = np.argsort(hit_read, kind="stable").astype(np.int64)
    read_hit_len = np.bincount(hit_read, minlength=n_reads).astype(np.int32)
    read_hit_off = np.zeros(n_reads + 1, dtype=np.int64)
    if n_reads:
        read_hit_off[1:] = np.cumsum(read_hit_len, dtype=np.int64)

    contig_hit_start = np.searchsorted(
        hit_contig, np.arange(n_contigs + 1, dtype=np.int32)
    ).astype(np.int64)

    contig_max_span = np.zeros(max(n_contigs, 1), dtype=np.int64)
    for c in range(n_contigs):
        lo, hi = int(contig_hit_start[c]), int(contig_hit_start[c + 1])
        if hi > lo:
            contig_max_span[c] = int((hit_shigh[lo:hi] - hit_slow[lo:hi]).max())

    # ---- annotations -----------------------------------------------------

    annotations = [a for a in reference.annotations if a.contig in contig_index]
    annotations.sort(key=lambda a: (contig_index[a.contig], a.start, a.end))

    group_vocab: list[str] = []
    feature_vocab: list[str] = []
    group_ids: dict[str, int] = {}
    feature_ids: dict[str, int] = {}

    n_anno = len(annotations)
    anno_contig = np.zeros(n_anno, dtype=np.int32)
    anno_start = np.zeros(n_anno, dtype=np.int64)
    anno_end = np.zeros(n_anno, dtype=np.int64)
    anno_strand = np.zeros(n_anno, dtype=np.int8)
    anno_group = np.zeros(n_anno, dtype=np.int8)
    anno_feature = np.zeros(n_anno, dtype=np.int16)
    anno_names: list[str] = []
    anno_products: list[str] = []
    anno_locus_tags: list[str] = []
    anno_notes: list[str] = []

    for i, annotation in enumerate(annotations):
        anno_contig[i] = contig_index[annotation.contig]
        anno_start[i] = annotation.start
        anno_end[i] = annotation.end
        anno_strand[i] = 1 if annotation.strand >= 0 else -1
        gid = group_ids.get(annotation.group)
        if gid is None:
            gid = len(group_vocab)
            group_ids[annotation.group] = gid
            group_vocab.append(annotation.group)
        anno_group[i] = gid
        fid = feature_ids.get(annotation.feature)
        if fid is None:
            fid = len(feature_vocab)
            feature_ids[annotation.feature] = fid
            feature_vocab.append(annotation.feature)
        anno_feature[i] = fid
        anno_names.append(annotation.name)
        anno_products.append(annotation.product)
        anno_locus_tags.append(annotation.locus_tag)
        anno_notes.append(annotation.note)

    contig_anno_start = np.searchsorted(
        anno_contig, np.arange(n_contigs + 1, dtype=np.int32)
    ).astype(np.int64)
    contig_anno_max_span = np.zeros(max(n_contigs, 1), dtype=np.int64)
    for c in range(n_contigs):
        lo, hi = int(contig_anno_start[c]), int(contig_anno_start[c + 1])
        if hi > lo:
            contig_anno_max_span[c] = int((anno_end[lo:hi] - anno_start[lo:hi]).max())

    # ---- events ----------------------------------------------------------

    name_to_read = {group.read: i for i, group in enumerate(read_list)}
    events = [e for e in stats.events() if e.contig in contig_index]
    events.sort(key=lambda e: (contig_index[e.contig], e.position))

    n_events = len(events)
    event_contig = np.zeros(n_events, dtype=np.int32)
    event_pos = np.zeros(n_events, dtype=np.int64)
    event_kind = np.zeros(n_events, dtype=np.int8)
    event_read = np.full(n_events, -1, dtype=np.int32)
    for i, event in enumerate(events):
        event_contig[i] = contig_index[event.contig]
        event_pos[i] = event.position
        event_kind[i] = _EVENT_INDEX.get(event.kind, 0)
        event_read[i] = name_to_read.get(event.read, -1)

    contig_event_start = np.searchsorted(
        event_contig, np.arange(n_contigs + 1, dtype=np.int32)
    ).astype(np.int64)

    # ---- depth -----------------------------------------------------------

    arrays: dict[str, np.ndarray] = {}
    contig_meta: list[dict[str, object]] = []
    for c, contig in enumerate(contigs):
        bin_size = _depth_bin(contig.length)
        n_bins = max(1, math.ceil(contig.length / bin_size)) if contig.length else 1
        lo, hi = int(contig_hit_start[c]), int(contig_hit_start[c + 1])
        if hi > lo:
            b0 = np.minimum(hit_slow[lo:hi] // bin_size, n_bins - 1)
            b1 = np.minimum(hit_shigh[lo:hi] // bin_size, n_bins - 1)
            diff = np.bincount(b0, minlength=n_bins + 1).astype(np.int64)
            diff -= np.bincount(b1 + 1, minlength=n_bins + 1).astype(np.int64)
            depth = np.cumsum(diff[: n_bins + 1])[:n_bins]
        else:
            depth = np.zeros(n_bins, dtype=np.int64)
        arrays[f"depth_{c}"] = depth.astype(np.int32)

        cs = stats.per_contig.get(contig.name)
        contig_meta.append(
            {
                "name": contig.name,
                "length": contig.length,
                "depth_bin": bin_size,
                "reads": cs.reads if cs else 0,
                "matched_bases": cs.matched_bases if cs else 0,
                "mean_coverage": round(cs.mean_coverage, 4) if cs else 0.0,
                "inverted": cs.inverted if cs else 0,
                "counts": (
                    {c2.value: cs.counts[c2] for c2 in CLASS_ORDER} if cs else
                    {c2.value: 0 for c2 in CLASS_ORDER}
                ),
                "events": len(cs.events) if cs else 0,
            }
        )

    # ---- metadata --------------------------------------------------------

    from readrift import __version__

    meta = {
        "format": FORMAT_VERSION,
        "readrift_version": __version__,
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
        "sample": Path(params.out).name,
        "contigs": contig_meta,
        "classes": list(CLASS_NAMES),
        "junctions": list(JUNCTION_NAMES),
        "event_kinds": list(EVENT_KINDS),
        "anno_groups": group_vocab,
        "anno_features": feature_vocab,
        "has_identity": want_identity,
        "is_genbank": reference.is_genbank,
        "params": _classifying_snapshot(params),
        "sources": {
            "reference": _source_stamp(params.seq_file),
            "btop": _source_stamp(params.btop_file),
        },
        "project": {
            "organism": reference.metadata.organism,
            "bio_project": reference.metadata.bio_project,
            "bio_sample": reference.metadata.bio_sample,
            "sra": reference.metadata.sra,
            "assembly_method": reference.metadata.assembly_method,
            "sequencing_technology": reference.metadata.sequencing_technology,
        },
        "summary": {
            "total_reads": stats.total_reads,
            "counts": {c.value: stats.counts[c] for c in CLASS_ORDER},
            "inverted_reads": stats.inverted_reads,
            "matched_bases": stats.matched_bases,
            "mean_coverage": round(stats.mean_coverage, 4),
            "n50": stats.n50,
            "reference_length": stats.reference_length,
            "events": n_events,
        },
        "notes": list(notes),
    }

    name_blob, name_off = _pack_strings(names)
    label_blob, label_off = _pack_strings(labels)
    anno_name_blob, anno_name_off = _pack_strings(anno_names)
    anno_product_blob, anno_product_off = _pack_strings(anno_products)
    anno_locus_blob, anno_locus_off = _pack_strings(anno_locus_tags)
    anno_note_blob, anno_note_off = _pack_strings(anno_notes)

    arrays.update(
        {
            "read_contig": read_contig,
            "read_cls": read_cls,
            "read_junction": read_junction,
            "read_inverted": read_inverted,
            "read_qlen": read_qlen,
            "read_distance": read_distance,
            "read_strand": read_strand,
            "read_hit_idx": read_hit_idx,
            "read_hit_off": read_hit_off,
            "read_hit_len": read_hit_len,
            "hit_read": hit_read,
            "hit_contig": hit_contig,
            "hit_slow": hit_slow,
            "hit_shigh": hit_shigh,
            "hit_qstart": hit_qstart,
            "hit_qend": hit_qend,
            "hit_strand": hit_strand,
            "hit_identity": hit_identity,
            "contig_hit_start": contig_hit_start,
            "contig_max_span": contig_max_span,
            "anno_contig": anno_contig,
            "anno_start": anno_start,
            "anno_end": anno_end,
            "anno_strand": anno_strand,
            "anno_group": anno_group,
            "anno_feature": anno_feature,
            "contig_anno_start": contig_anno_start,
            "contig_anno_max_span": contig_anno_max_span,
            "event_contig": event_contig,
            "event_pos": event_pos,
            "event_kind": event_kind,
            "event_read": event_read,
            "contig_event_start": contig_event_start,
            "name_blob": name_blob,
            "name_off": name_off,
            "label_blob": label_blob,
            "label_off": label_off,
            "anno_name_blob": anno_name_blob,
            "anno_name_off": anno_name_off,
            "anno_product_blob": anno_product_blob,
            "anno_product_off": anno_product_off,
            "anno_locus_blob": anno_locus_blob,
            "anno_locus_off": anno_locus_off,
            "anno_note_blob": anno_note_blob,
            "anno_note_off": anno_note_off,
            "meta": np.array(json.dumps(meta, ensure_ascii=False)),
        }
    )

    # Uncompressed: the archive is tens of MB and loads in well under a second,
    # which is the whole point.  Compression would trade that away for disk.
    temporary = target.with_suffix(target.suffix + ".tmp")
    with open(temporary, "wb") as handle:
        np.savez(handle, **arrays)
    temporary.replace(target)
    return target


# --------------------------------------------------------------------------
# Reading
# --------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class ContigInfo:
    index: int
    name: str
    length: int
    depth_bin: int


class StoreFormatError(Exception):
    """The file is not a readrift cache, or is one this build cannot read."""


class BrowserStore:
    """Read-only access to one cached run."""

    def __init__(self, path: Path, data: dict[str, np.ndarray], meta: dict) -> None:
        self.path = path
        self.meta = meta
        self._d = data

        self.contigs: list[ContigInfo] = [
            ContigInfo(i, c["name"], int(c["length"]), int(c["depth_bin"]))
            for i, c in enumerate(meta["contigs"])
        ]
        self._contig_by_name = {c.name: c for c in self.contigs}

        self.names = _StringColumn(data["name_blob"], data["name_off"])
        self.labels = _StringColumn(data["label_blob"], data["label_off"])
        self.anno_names = _StringColumn(data["anno_name_blob"], data["anno_name_off"])
        self.anno_products = _StringColumn(
            data["anno_product_blob"], data["anno_product_off"]
        )
        self.anno_locus_tags = _StringColumn(
            data["anno_locus_blob"], data["anno_locus_off"]
        )
        self.anno_notes = _StringColumn(data["anno_note_blob"], data["anno_note_off"])

        self._read_lookup: dict[str, int] | None = None
        self._anno_lookup: list[tuple[str, int]] | None = None

    # ---- construction ----------------------------------------------------

    @classmethod
    def open(cls, path: str | Path) -> BrowserStore:
        """Load a cache.  Raises :class:`StoreFormatError` if it cannot be read."""
        target = Path(path)
        try:
            with np.load(target, allow_pickle=False) as archive:
                data = {key: archive[key] for key in archive.files}
        except (OSError, ValueError) as exc:
            raise StoreFormatError(f"cannot read {target}: {exc}") from exc

        if "meta" not in data:
            raise StoreFormatError(f"{target} is not a readrift cache")
        try:
            meta = json.loads(str(data["meta"].item()))
        except (ValueError, AttributeError) as exc:
            raise StoreFormatError(f"{target} has unreadable metadata: {exc}") from exc

        if meta.get("format") != FORMAT_VERSION:
            raise StoreFormatError(
                f"{target} was written by cache format {meta.get('format')}, "
                f"this build reads {FORMAT_VERSION}"
            )
        return cls(target, data, meta)

    def is_stale(self, params: Params) -> tuple[bool, str]:
        """Whether *params* would produce different contents.

        Returns ``(stale, reason)`` so the caller can say *why* it is rebuilding
        rather than silently spending minutes.
        """
        wanted = _classifying_snapshot(params)
        stored = self.meta.get("params", {})
        for key, value in wanted.items():
            if stored.get(key) != value:
                return True, f"--{key.replace('_', '-')} changed ({stored.get(key)} -> {value})"

        for label, attribute in (("reference", "seq_file"), ("BTOP", "btop_file")):
            path = getattr(params, attribute, "")
            if not path:
                continue
            now = _source_stamp(path)
            before = self.meta.get("sources", {}).get(
                "reference" if label == "reference" else "btop", {}
            )
            if Path(str(before.get("path", ""))) != Path(now["path"]):
                return True, f"{label} file is a different path"
            if before.get("size") != now["size"] or before.get("mtime") != now["mtime"]:
                return True, f"{label} file has changed on disk"

        return False, ""

    # ---- contigs ---------------------------------------------------------

    def contig(self, name: str) -> ContigInfo | None:
        return self._contig_by_name.get(name)

    def require_contig(self, name: str) -> ContigInfo:
        info = self._contig_by_name.get(name)
        if info is None:
            raise KeyError(name)
        return info

    # ---- region queries --------------------------------------------------

    def hits_in(self, contig: str, start: int, end: int) -> np.ndarray:
        """Indices of every HSP on *contig* overlapping ``[start, end]``.

        Sorted by start position, because the hit arrays are.
        """
        info = self._contig_by_name.get(contig)
        if info is None:
            return np.empty(0, dtype=np.int64)

        c = info.index
        lo_all = int(self._d["contig_hit_start"][c])
        hi_all = int(self._d["contig_hit_start"][c + 1])
        if hi_all <= lo_all:
            return np.empty(0, dtype=np.int64)

        slow = self._d["hit_slow"][lo_all:hi_all]
        shigh = self._d["hit_shigh"][lo_all:hi_all]
        span = int(self._d["contig_max_span"][c])

        lo = int(np.searchsorted(slow, start - span, side="left"))
        hi = int(np.searchsorted(slow, end, side="right"))
        if hi <= lo:
            return np.empty(0, dtype=np.int64)

        keep = np.flatnonzero(shigh[lo:hi] >= start)
        return (keep + lo + lo_all).astype(np.int64)

    def annotations_in(self, contig: str, start: int, end: int) -> np.ndarray:
        """Indices of every annotation on *contig* overlapping ``[start, end]``."""
        info = self._contig_by_name.get(contig)
        if info is None:
            return np.empty(0, dtype=np.int64)

        c = info.index
        lo_all = int(self._d["contig_anno_start"][c])
        hi_all = int(self._d["contig_anno_start"][c + 1])
        if hi_all <= lo_all:
            return np.empty(0, dtype=np.int64)

        starts = self._d["anno_start"][lo_all:hi_all]
        ends = self._d["anno_end"][lo_all:hi_all]
        span = int(self._d["contig_anno_max_span"][c])

        lo = int(np.searchsorted(starts, start - span, side="left"))
        hi = int(np.searchsorted(starts, end, side="right"))
        if hi <= lo:
            return np.empty(0, dtype=np.int64)

        keep = np.flatnonzero(ends[lo:hi] >= start)
        return (keep + lo + lo_all).astype(np.int64)

    def read_hits(self, read_index: int) -> np.ndarray:
        """Indices of every HSP belonging to one read."""
        off = int(self._d["read_hit_off"][read_index])
        length = int(self._d["read_hit_len"][read_index])
        return self._d["read_hit_idx"][off : off + length]

    def read_hits_on(self, read_index: int, contig_index: int) -> np.ndarray:
        """Indices of one read's HSPs that lie on one contig."""
        hits = self.read_hits(read_index)
        if hits.size == 0:
            return hits
        return hits[self._d["hit_contig"][hits] == contig_index]

    # ---- column access ---------------------------------------------------

    def column(self, name: str) -> np.ndarray:
        return self._d[name]

    @property
    def read_count(self) -> int:
        return int(self._d["read_contig"].size)

    @property
    def hit_count(self) -> int:
        return int(self._d["hit_slow"].size)

    def class_name(self, index: int) -> str:
        return CLASS_NAMES[int(index)]

    def junction_name(self, index: int) -> str:
        return JUNCTION_NAMES[int(index)]

    # ---- depth -----------------------------------------------------------

    def depth(
        self, contig: str, start: int, end: int, bins: int = 800
    ) -> dict[str, list]:
        """Binned depth over ``[start, end]``, downsampled to at most *bins*.

        Returns the minimum, mean and maximum within each output bin.  The band
        between min and max is what makes a narrow coverage trough visible when
        a whole contig is on screen -- a mean alone smooths it away, and a
        trough is exactly where structural variation lives.
        """
        info = self._contig_by_name.get(contig)
        empty = {"positions": [], "min": [], "mean": [], "max": [], "bin": 0}
        if info is None:
            return empty

        depth = self._d.get(f"depth_{info.index}")
        if depth is None or depth.size == 0:
            return empty

        bin_size = info.depth_bin
        lo = max(0, int(start) // bin_size)
        hi = min(depth.size, int(end) // bin_size + 1)
        if hi <= lo:
            return empty

        segment = depth[lo:hi].astype(np.float64)
        bins = max(1, int(bins))

        if segment.size <= bins:
            positions = (np.arange(lo, hi) * bin_size).astype(np.int64)
            values = segment
            return {
                "positions": positions.tolist(),
                "min": values.tolist(),
                "mean": values.tolist(),
                "max": values.tolist(),
                "bin": bin_size,
            }

        edges = np.linspace(0, segment.size, bins + 1).astype(np.int64)
        edges[-1] = segment.size
        starts = edges[:-1]
        widths = np.diff(edges).astype(np.float64)

        means = np.add.reduceat(segment, starts) / widths
        maxima = np.maximum.reduceat(segment, starts)
        minima = np.minimum.reduceat(segment, starts)
        positions = ((lo + starts) * bin_size).astype(np.int64)

        return {
            "positions": positions.tolist(),
            "min": minima.tolist(),
            "mean": [round(v, 3) for v in means.tolist()],
            "max": maxima.tolist(),
            "bin": int(bin_size * widths[0]),
        }

    # ---- read detail -----------------------------------------------------

    def annotation_detail(self, anno_index: int) -> dict:
        """Everything known about one GenBank feature, for the detail panel.

        Fetched on demand rather than carried in every region payload: a whole
        contig on screen is thousands of features, and the products and notes
        that make this worth reading are far larger than the coordinates that
        get drawn.
        """
        d = self._d
        if not 0 <= anno_index < int(d["anno_start"].size):
            raise IndexError(anno_index)

        groups = self.meta.get("anno_groups", [])
        features = self.meta.get("anno_features", [])
        gid = int(d["anno_group"][anno_index])
        fid = int(d["anno_feature"][anno_index])
        start = int(d["anno_start"][anno_index])
        end = int(d["anno_end"][anno_index])

        return {
            "id": int(anno_index),
            "contig": self.contigs[int(d["anno_contig"][anno_index])].name,
            "start": start,
            "end": end,
            "length": end - start,
            "strand": int(d["anno_strand"][anno_index]),
            "feature": features[fid] if fid < len(features) else "",
            "group": groups[gid] if gid < len(groups) else "",
            "name": self.anno_names[anno_index],
            "label": gene_label(self.anno_names[anno_index]),
            "locus_tag": self.anno_locus_tags[anno_index],
            "product": self.anno_products[anno_index],
            "note": self.anno_notes[anno_index],
        }

    def read_detail(self, read_index: int) -> dict:
        """Everything known about one read, for the detail panel."""
        if not 0 <= read_index < self.read_count:
            raise IndexError(read_index)

        d = self._d
        hits = self.read_hits(read_index)
        rows = []
        for h in hits.tolist():
            identity = float(d["hit_identity"][h])
            # ref_start / ref_end are low / high, not BLAST's sstart / send:
            # the cache stores the ordered pair and carries the direction in
            # `strand`, so a minus-strand hit is not a backwards-looking range.
            rows.append(
                {
                    "contig": self.contigs[int(d["hit_contig"][h])].name,
                    "strand": int(d["hit_strand"][h]),
                    "ref_start": int(d["hit_slow"][h]),
                    "ref_end": int(d["hit_shigh"][h]),
                    "ref_span": int(d["hit_shigh"][h] - d["hit_slow"][h]),
                    "qstart": int(d["hit_qstart"][h]),
                    "qend": int(d["hit_qend"][h]),
                    "identity": None if math.isnan(identity) else round(identity, 5),
                }
            )
        rows.sort(key=lambda r: (r["contig"], r["ref_start"]))

        return {
            "id": int(read_index),
            "name": self.names[read_index],
            "label": self.labels[read_index],
            "cls": self.class_name(d["read_cls"][read_index]),
            "junction": self.junction_name(d["read_junction"][read_index]),
            "inverted": bool(d["read_inverted"][read_index]),
            "strand": int(d["read_strand"][read_index]),
            "qlen": int(d["read_qlen"][read_index]),
            "distance": int(d["read_distance"][read_index]),
            "contig": self.contigs[int(d["read_contig"][read_index])].name,
            "matched_bases": sum(r["ref_span"] for r in rows),
            "hits": rows,
            "has_identity": bool(self.meta.get("has_identity")),
        }

    # ---- search ----------------------------------------------------------

    def _ensure_read_lookup(self) -> dict[str, int]:
        """Read name and label -> index.  Built on first use, not on open.

        Half a million names is ~100 ms and ~60 MB of Python strings; a session
        that never searches should not pay it.
        """
        if self._read_lookup is None:
            lookup: dict[str, int] = {}
            for i in range(self.read_count):
                lookup.setdefault(self.names[i], i)
                lookup.setdefault(self.labels[i], i)
            self._read_lookup = lookup
        return self._read_lookup

    def _ensure_anno_lookup(self) -> list[tuple[str, int]]:
        if self._anno_lookup is None:
            self._anno_lookup = [
                (self.anno_names[i].lower(), i) for i in range(len(self.anno_names))
            ]
        return self._anno_lookup

    def read_extent(self, read_index: int, contig_index: int | None = None) -> tuple[str, int, int]:
        """Contig and reference span of one read, for jumping to it."""
        d = self._d
        if contig_index is None:
            contig_index = int(d["read_contig"][read_index])
        hits = self.read_hits_on(read_index, contig_index)
        if hits.size == 0:
            hits = self.read_hits(read_index)
            if hits.size == 0:
                return self.contigs[contig_index].name, 0, 0
            contig_index = int(d["hit_contig"][hits[0]])
            hits = self.read_hits_on(read_index, contig_index)
        return (
            self.contigs[contig_index].name,
            int(d["hit_slow"][hits].min()),
            int(d["hit_shigh"][hits].max()),
        )

    def search(self, query: str, limit: int = 25) -> list[dict]:
        """Resolve a coordinate, a contig, a gene name or a read name."""
        text = query.strip()
        if not text:
            return []

        results: list[dict] = []

        # contig:start-end / contig:start..end / contig:position
        if ":" in text:
            contig_part, _, span = text.rpartition(":")
            info = self._contig_by_name.get(contig_part.strip())
            if info is not None:
                bounds = _parse_span(span)
                if bounds is not None:
                    start, end = bounds
                    results.append(
                        {
                            "kind": "region",
                            "label": f"{info.name}:{start:,}-{end:,}",
                            "contig": info.name,
                            "start": max(0, start),
                            "end": min(info.length, end),
                        }
                    )
                    return results

        # a bare contig name -> the whole contig
        info = self._contig_by_name.get(text)
        if info is not None:
            return [
                {
                    "kind": "contig",
                    "label": f"{info.name} ({info.length:,} bp)",
                    "contig": info.name,
                    "start": 0,
                    "end": info.length,
                }
            ]

        # a read name or its short label
        read_index = self._ensure_read_lookup().get(text)
        if read_index is not None:
            contig, start, end = self.read_extent(read_index)
            pad = max(500, (end - start) // 5)
            results.append(
                {
                    "kind": "read",
                    "label": f"{self.labels[read_index]}  ({self.names[read_index]})",
                    "contig": contig,
                    "start": max(0, start - pad),
                    "end": end + pad,
                    "read": int(read_index),
                }
            )

        # gene names: exact, then prefix, then substring
        needle = text.lower()
        d = self._d
        exact: list[int] = []
        prefix: list[int] = []
        loose: list[int] = []
        for name, index in self._ensure_anno_lookup():
            if name == needle:
                exact.append(index)
            elif name.startswith(needle):
                prefix.append(index)
            elif needle in name:
                loose.append(index)
            if len(exact) + len(prefix) + len(loose) >= limit * 3:
                break

        for index in (exact + prefix + loose)[: max(0, limit - len(results))]:
            contig = self.contigs[int(d["anno_contig"][index])]
            start = int(d["anno_start"][index])
            end = int(d["anno_end"][index])
            pad = max(500, (end - start))
            results.append(
                {
                    "kind": "gene",
                    "label": f"{self.anno_names[index]}  {contig.name}:{start:,}-{end:,}",
                    "contig": contig.name,
                    "start": max(0, start - pad),
                    "end": min(contig.length, end + pad),
                    "feature_start": start,
                    "feature_end": end,
                }
            )

        return results

    # ---- events ----------------------------------------------------------

    def events_in(self, contig: str, start: int, end: int, kinds: Sequence[str] | None = None) -> list[dict]:
        info = self._contig_by_name.get(contig)
        if info is None:
            return []
        d = self._d
        lo_all = int(d["contig_event_start"][info.index])
        hi_all = int(d["contig_event_start"][info.index + 1])
        if hi_all <= lo_all:
            return []
        positions = d["event_pos"][lo_all:hi_all]
        lo = int(np.searchsorted(positions, start, side="left"))
        hi = int(np.searchsorted(positions, end, side="right"))
        wanted = self._kind_mask(kinds)
        out: list[dict] = []
        for i in range(lo_all + lo, lo_all + hi):
            kind = EVENT_KINDS[int(d["event_kind"][i])]
            if kind not in wanted:
                continue
            out.append(self._event_dict(i))
        return out

    def next_event(
        self,
        contig: str,
        position: int,
        direction: int = 1,
        kinds: Sequence[str] | None = None,
        wrap: bool = True,
    ) -> dict | None:
        """The nearest structural event after (or before) *position*.

        Walks off the end of a contig onto the next one, and wraps at the end of
        the reference, so pressing the key repeatedly tours every event in the
        run without the user having to change contig by hand.
        """
        if not self.contigs:
            return None
        info = self._contig_by_name.get(contig)
        start_index = info.index if info else 0
        wanted = self._kind_mask(kinds)
        d = self._d
        n = len(self.contigs)

        for step in range(n + 1):
            c = (start_index + direction * step) % n
            lo_all = int(d["contig_event_start"][c])
            hi_all = int(d["contig_event_start"][c + 1])
            if hi_all <= lo_all:
                continue
            positions = d["event_pos"][lo_all:hi_all]

            if step == 0:
                if direction > 0:
                    order = range(
                        lo_all + int(np.searchsorted(positions, position, side="right")),
                        hi_all,
                    )
                else:
                    order = range(
                        lo_all + int(np.searchsorted(positions, position, side="left")) - 1,
                        lo_all - 1,
                        -1,
                    )
            elif direction > 0:
                order = range(lo_all, hi_all)
            else:
                order = range(hi_all - 1, lo_all - 1, -1)

            for i in order:
                if EVENT_KINDS[int(d["event_kind"][i])] in wanted:
                    return self._event_dict(i)

            if not wrap and step > 0:
                return None

        return None

    def _kind_mask(self, kinds: Sequence[str] | None) -> set[str]:
        if not kinds:
            return set(EVENT_KINDS)
        return {k for k in kinds if k in _EVENT_INDEX}

    def _event_dict(self, index: int) -> dict:
        d = self._d
        read_index = int(d["event_read"][index])
        contig = self.contigs[int(d["event_contig"][index])]
        return {
            "contig": contig.name,
            "position": int(d["event_pos"][index]),
            "kind": EVENT_KINDS[int(d["event_kind"][index])],
            "read": read_index if read_index >= 0 else None,
            "read_label": self.labels[read_index] if read_index >= 0 else "",
        }

    @property
    def event_count(self) -> int:
        return int(self._d["event_pos"].size)


def _parse_span(text: str) -> tuple[int, int] | None:
    """``120000-160000``, ``120000..160000`` or ``120000`` -> a range."""
    cleaned = text.strip().replace(",", "").replace("_", "").replace("..", "-")
    if not cleaned:
        return None
    if "-" in cleaned[1:]:
        left, _, right = cleaned[1:].rpartition("-")
        left = cleaned[0] + left
        try:
            start, end = int(left), int(right)
        except ValueError:
            return None
        return (start, end) if start <= end else (end, start)
    try:
        centre = int(cleaned)
    except ValueError:
        return None
    return max(0, centre - 5_000), centre + 5_000
