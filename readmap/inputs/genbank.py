"""GenBank reference parser: contig lengths, project metadata, and features.

Differences from the Perl (see ``port_plan.md`` §4):

* Metadata fields default to ``""`` instead of being left undefined, which is
  what killed the Perl on any record lacking ``BioProject:`` (finding B13).
* Feature blocks may begin at any standalone feature key, not only ``gene``,
  so ``repeat_region`` / ``mobile_element`` / ``misc_feature`` are no longer
  folded into whichever gene happened to precede them (finding B24).
* The whole RNA family is recognised, not just ``tRNA|rRNA|RNA``.
* Coordinates come from feature *locations* only, never from qualifier text.
* ``ACCESSION`` and ``VERSION`` are recorded as aliases of the ``LOCUS`` name,
  because a BLAST database built from the same record usually carries the
  versioned accession while the Perl matched on the ``LOCUS`` name alone --
  a silent, total mismatch (finding B11).
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO

from readmap.models import Annotation, Contig, Metadata

# --------------------------------------------------------------------------
# Line patterns
# --------------------------------------------------------------------------

_LOCUS_RE = re.compile(r"^LOCUS\s+(\S+)(?:\s+(\d+)\s+bp)?", re.IGNORECASE)
_FEATURE_RE = re.compile(r"^ {5}(\S+)\s+(\S.*)$")
_CONTINUATION_RE = re.compile(r"^ {6,}(\S.*)$")
_QUALIFIER_RE = re.compile(r"^/([A-Za-z_]\w*)(?:=(.*))?$", re.DOTALL)
_RANGE_RE = re.compile(r"[<>]?(\d+)\s*\.\.\s*[<>]?(\d+)")
_SINGLE_POS_RE = re.compile(r"^\s*[<>]?(\d+)\s*$")

_META_PATTERNS: dict[str, re.Pattern[str]] = {
    "organism": re.compile(r"^\s*ORGANISM\s+(.+?)\s*$"),
    "bio_project": re.compile(r"BioProject:\s*(\S+)"),
    "bio_sample": re.compile(r"BioSample:\s*(\S+)"),
    "sra": re.compile(r"Sequence Read Archive:\s*(\S+)"),
    "assembly_method": re.compile(r"Assembly Method\s*::\s*(.+?)\s*$"),
    "sequencing_technology": re.compile(r"Sequencing Technology\s*::\s*(.+?)\s*$"),
}

# Feature keys that stand on their own rather than qualifying a gene.
_BLOCK_START = frozenset(
    {
        "gene",
        "repeat_region",
        "mobile_element",
        "misc_feature",
        "misc_binding",
        "misc_structure",
        "rep_origin",
        "oriT",
        "operon",
        "regulatory",
        "STS",
        "protein_bind",
        "stem_loop",
    }
)

# Never interesting for the map.
_SKIP_KEYS = frozenset({"source"})

_PHAGE_RE = re.compile(r"\b(phage|tail)\b", re.IGNORECASE)
_HYP_RE = re.compile(r"\bhypothetical\b", re.IGNORECASE)
_IS_RE = re.compile(r"\b(transposase|recombinase|integrase|insertion sequence)\b", re.IGNORECASE)


def _open(path: str | Path) -> IO[str]:
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, "rt", encoding="utf-8", errors="replace")
    return open(p, "r", encoding="utf-8", errors="replace")


# --------------------------------------------------------------------------
# Feature assembly
# --------------------------------------------------------------------------


@dataclass(slots=True)
class _Feature:
    key: str
    location: str
    qualifiers: dict[str, list[str]] = field(default_factory=dict)

    def add_qualifier(self, raw: str) -> None:
        match = _QUALIFIER_RE.match(raw)
        if not match:
            return
        name = match.group(1)
        value = match.group(2) or ""
        value = value.strip()
        if value.startswith('"') and value.endswith('"') and len(value) >= 2:
            value = value[1:-1]
        elif value.startswith('"'):
            value = value[1:]
        self.qualifiers.setdefault(name, []).append(value)

    def first(self, name: str) -> str:
        values = self.qualifiers.get(name)
        return values[0] if values else ""


def _coords(locations: list[str]) -> tuple[int, int] | None:
    """Smallest and largest coordinate mentioned in a set of locations."""
    values: list[int] = []
    for loc in locations:
        for a, b in _RANGE_RE.findall(loc):
            values.append(int(a))
            values.append(int(b))
        if not _RANGE_RE.search(loc):
            single = _SINGLE_POS_RE.match(loc)
            if single:
                values.append(int(single.group(1)))
    if not values:
        return None
    return min(values), max(values)


def _classify_block(contig: str, block: list[_Feature]) -> Annotation | None:
    """Reduce one feature block to a single drawable annotation."""
    locations = [f.location for f in block if f.location]
    span = _coords(locations)
    if span is None:
        return None
    start, end = span

    strand = -1 if any("complement" in f.location for f in block) else 1

    # Prefer a descriptive sub-feature (CDS, tRNA, ...) over the bare `gene`.
    descriptive = next((f for f in block if f.key != "gene"), block[0])
    feature_key = descriptive.key

    def first_in_block(qualifier: str) -> str:
        for f in block:
            value = f.first(qualifier)
            if value:
                return value
        return ""

    locus_tag = first_in_block("locus_tag")
    name = first_in_block("gene") or locus_tag
    product = first_in_block("product")
    note = first_in_block("note")

    group = ""
    if feature_key.endswith("RNA"):
        group = "RNA"
    # Order below matches the Perl: each test overwrites the previous, so the
    # effective precedence is pseudo > IS > hyp > phage > RNA.
    if product and _PHAGE_RE.search(product):
        group = "phage"
    if product and _HYP_RE.search(product):
        group = "hyp"
    if product and _IS_RE.search(product):
        group = "IS"
    if any("pseudo" in f.qualifiers or "pseudogene" in f.qualifiers for f in block):
        group = "pseudo"

    return Annotation(
        contig=contig,
        strand=strand,
        start=start,
        end=end,
        feature=feature_key,
        name=name,
        group=group,
        product=product,
        locus_tag=locus_tag,
        note=note,
    )


# --------------------------------------------------------------------------
# Record parser
# --------------------------------------------------------------------------


def parse_genbank(
    path: str | Path,
) -> tuple[list[Contig], list[Annotation], Metadata, dict[str, str]]:
    """Parse a (possibly multi-record) GenBank file."""
    contigs: list[Contig] = []
    annotations: list[Annotation] = []
    aliases: dict[str, str] = {}
    meta = Metadata(source_format="gb")

    # per-record state
    name: str | None = None
    declared_length = 0
    counted_length = 0
    in_features = False
    in_origin = False
    record_aliases: list[str] = []

    block: list[_Feature] = []
    current: _Feature | None = None
    pending_qualifier: list[str] = []

    def flush_qualifier() -> None:
        nonlocal pending_qualifier
        if current is not None and pending_qualifier:
            current.add_qualifier(" ".join(pending_qualifier))
        pending_qualifier = []

    def flush_feature() -> None:
        nonlocal current
        flush_qualifier()
        if current is not None:
            block.append(current)
        current = None

    def flush_block() -> None:
        nonlocal block
        flush_feature()
        if block and name:
            annotation = _classify_block(name, block)
            if annotation is not None:
                annotations.append(annotation)
        block = []

    def flush_record() -> None:
        nonlocal name, declared_length, counted_length, record_aliases
        nonlocal in_features, in_origin
        flush_block()
        if name:
            length = counted_length or declared_length
            contigs.append(Contig(name=name, length=length))
            for alias in [name, *record_aliases]:
                if alias:
                    aliases.setdefault(alias, name)
                    if "." in alias:
                        aliases.setdefault(alias.rsplit(".", 1)[0], name)
        name = None
        declared_length = 0
        counted_length = 0
        record_aliases = []
        in_features = False
        in_origin = False

    with _open(path) as handle:
        for raw in handle:
            line = raw.rstrip("\r\n").expandtabs(8)

            if line.startswith("//"):
                flush_record()
                continue

            locus = _LOCUS_RE.match(line)
            if locus:
                flush_record()
                name = locus.group(1)
                declared_length = int(locus.group(2)) if locus.group(2) else 0
                continue

            if name is None:
                continue

            if in_origin:
                counted_length += sum(1 for ch in line if ch.isalpha())
                continue

            if line.upper().startswith("ORIGIN"):
                flush_block()
                in_origin = True
                in_features = False
                continue

            if line.startswith("FEATURES"):
                in_features = True
                continue

            if not in_features:
                # header section: metadata and accession aliases
                if line.startswith("ACCESSION"):
                    record_aliases.extend(line[len("ACCESSION") :].split())
                elif line.startswith("VERSION"):
                    parts = line[len("VERSION") :].split()
                    if parts:
                        record_aliases.append(parts[0])
                for key, pattern in _META_PATTERNS.items():
                    if getattr(meta, key):
                        continue
                    found = pattern.search(line)
                    if found:
                        setattr(meta, key, found.group(1).strip())
                if not meta.first_line:
                    meta.first_line = line.strip()
                continue

            # ---- inside FEATURES ----
            feature = _FEATURE_RE.match(line)
            if feature:
                key, location = feature.group(1), feature.group(2).strip()
                if key in _SKIP_KEYS:
                    flush_feature()
                    continue
                if key in _BLOCK_START:
                    flush_block()
                else:
                    flush_feature()
                current = _Feature(key=key, location=location)
                continue

            cont = _CONTINUATION_RE.match(line)
            if cont and current is not None:
                text = cont.group(1)
                if text.startswith("/"):
                    flush_qualifier()
                    pending_qualifier = [text]
                elif pending_qualifier:
                    pending_qualifier.append(text)
                else:
                    # continuation of the location line
                    current.location += text
                continue

    flush_record()

    if not contigs:
        raise ValueError(f"no LOCUS records found in {path}")

    return contigs, annotations, meta, aliases
