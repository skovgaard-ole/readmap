"""Summary statistics: counts, coverage, and the inputs to the plots.

Coverage here is **matched bases only** (decision D3): the sum of aligned
reference spans.  Reference that a divided read jumps *over* is not covered by
that read.  The Perl used this definition for undivided and long-divided reads
but summed the whole first-to-last span for short-divided ones, so the number
it reported depended on how rearranged the sample was (finding B10).
"""

from __future__ import annotations

import math
from array import array
from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass, field
from itertools import pairwise

import numpy as np

from readrift.btop_trace import Trace
from readrift.btop_trace import parse as parse_btop
from readrift.models import Contig, Junction, ReadClass, ReadGroup

#: Resolution of the coverage profile.  Enough detail for a printed plot
#: without holding a per-base array for a large reference.
MAX_DEPTH_BINS = 4_000


@dataclass(frozen=True, slots=True)
class Event:
    """A structural feature worth marking on the event map."""

    contig: str
    position: int
    kind: str
    """``junction`` | ``inversion`` | ``circle`` | ``contig_join``."""

    read: str


@dataclass(slots=True)
class ContigStats:
    contig: Contig
    bin_size: int
    _diff: np.ndarray
    matched_bases: int = 0
    counts: Counter = field(default_factory=Counter)
    inverted: int = 0
    events: list[Event] = field(default_factory=list)

    @property
    def mean_coverage(self) -> float:
        return self.matched_bases / self.contig.length if self.contig.length else 0.0

    @property
    def reads(self) -> int:
        return sum(self.counts.values())

    def depth(self) -> np.ndarray:
        """Depth per bin, as a plain array."""
        return np.cumsum(self._diff)[:-1]

    def positions(self) -> np.ndarray:
        """Reference coordinate at the centre of each bin."""
        n = len(self._diff) - 1
        return (np.arange(n) + 0.5) * self.bin_size

    def add_interval(self, low: int, high: int) -> None:
        b0 = min(low // self.bin_size, len(self._diff) - 2)
        b1 = min(high // self.bin_size, len(self._diff) - 2)
        self._diff[b0] += 1
        self._diff[b1 + 1] -= 1


@dataclass(slots=True)
class Stats:
    """Everything the front page, the plots, and the TSV report need."""

    per_contig: dict[str, ContigStats] = field(default_factory=dict)
    counts: Counter = field(default_factory=Counter)
    inverted_reads: int = 0
    junction_counts: Counter = field(default_factory=Counter)
    read_lengths: dict[ReadClass, array] = field(default_factory=dict)
    distances: array = field(default_factory=lambda: array("q"))
    reference_length: int = 0
    identities: array = field(default_factory=lambda: array("d"))
    """Per-read alignment identity, only filled when ``--identity`` is on."""

    # ---- totals ----------------------------------------------------------

    @property
    def total_reads(self) -> int:
        return sum(self.counts.values())

    @property
    def matched_bases(self) -> int:
        return sum(cs.matched_bases for cs in self.per_contig.values())

    @property
    def mean_coverage(self) -> float:
        return self.matched_bases / self.reference_length if self.reference_length else 0.0

    def fraction(self, cls: ReadClass) -> float:
        total = self.total_reads
        return 100.0 * self.counts[cls] / total if total else 0.0

    @property
    def inverted_fraction(self) -> float:
        total = self.total_reads
        return 100.0 * self.inverted_reads / total if total else 0.0

    def all_read_lengths(self) -> np.ndarray:
        chunks = [np.asarray(v, dtype=np.int64) for v in self.read_lengths.values() if len(v)]
        return np.concatenate(chunks) if chunks else np.empty(0, dtype=np.int64)

    @property
    def n50(self) -> int:
        return n50(self.all_read_lengths())

    def events(self) -> list[Event]:
        out: list[Event] = []
        for cs in self.per_contig.values():
            out.extend(cs.events)
        return out


def n50(lengths: np.ndarray) -> int:
    """Length at which half the sequenced bases sit in reads at least that long."""
    if lengths.size == 0:
        return 0
    ordered = np.sort(lengths)[::-1]
    half = ordered.sum() / 2.0
    running = np.cumsum(ordered)
    index = int(np.searchsorted(running, half))
    return int(ordered[min(index, ordered.size - 1)])


def _new_contig_stats(contig: Contig) -> ContigStats:
    bin_size = max(1, math.ceil(contig.length / MAX_DEPTH_BINS)) if contig.length else 1
    n_bins = max(1, math.ceil(contig.length / bin_size)) if contig.length else 1
    return ContigStats(
        contig=contig,
        bin_size=bin_size,
        _diff=np.zeros(n_bins + 1, dtype=np.int32),
    )


def _record_events(per_contig: dict[str, ContigStats], group: ReadGroup) -> None:
    """Record the structural events one read is evidence for.

    Each event lands on the contig it actually sits on, so a contig-joining
    read marks both sides of the join rather than one arbitrary side.
    """

    def add(contig: str, position: int, kind: str) -> None:
        cs = per_contig.get(contig)
        if cs is not None:
            cs.events.append(Event(contig, position, kind, group.read))

    if group.junction is Junction.CIRCLE:
        # The junction is at the origin, so mark both ends of the span.
        home = group.contig
        add(home, group.ref_low, "circle")
        add(home, group.ref_high, "circle")
    elif group.junction is Junction.CONTIG_JOIN:
        # One mark per contig, at the end of that contig's piece nearest the
        # join -- i.e. the outermost coordinate the read reaches there.
        for contig in {h.contig for h in group.hits}:
            pieces = [h for h in group.hits if h.contig == contig]
            add(contig, max(h.shigh for h in pieces), "contig_join")
    elif group.cls is ReadClass.LONG_DIVIDED:
        # A plain long-distance junction: the gap the read jumps over. Only
        # meaningful within a single contig, which PLAIN guarantees.
        ordered = sorted(group.hits, key=lambda h: h.slow)
        for left, right in pairwise(ordered):
            if right.slow > left.shigh:
                add(left.contig, (left.shigh + right.slow) // 2, "junction")

    if group.inverted:
        reference_strand = group.hits[0].strand
        for hit in group.hits:
            if hit.strand != reference_strand:
                add(hit.contig, hit.slow, "inversion")


def collect(
    groups: Iterable[ReadGroup],
    contigs: Iterable[Contig],
    with_identity: bool = False,
) -> tuple[Stats, list[ReadGroup]]:
    """Accumulate statistics, returning them alongside the materialised reads.

    The reads are needed again by :mod:`readrift.layout`, so they are collected
    here rather than making the caller iterate the BTOP file twice.
    """
    contig_list = list(contigs)
    stats = Stats(
        per_contig={c.name: _new_contig_stats(c) for c in contig_list},
        read_lengths={cls: array("q") for cls in ReadClass},
        reference_length=sum(c.length for c in contig_list),
    )

    kept: list[ReadGroup] = []

    for group in groups:
        kept.append(group)
        stats.counts[group.cls] += 1
        stats.junction_counts[group.junction] += 1
        stats.read_lengths[group.cls].append(group.qlen)
        if group.cls is not ReadClass.UNDIVIDED:
            stats.distances.append(group.distance)
        if group.inverted:
            stats.inverted_reads += 1

        if with_identity:
            total = Trace()
            for hit in group.hits:
                part = parse_btop(hit.btop)
                total = Trace(
                    matches=total.matches + part.matches,
                    mismatches=total.mismatches + part.mismatches,
                    query_gaps=total.query_gaps + part.query_gaps,
                    subject_gaps=total.subject_gaps + part.subject_gaps,
                )
            if total.aligned_columns:
                stats.identities.append(total.identity)

        # Coverage is attributed to the contig each HSP actually lies on, so a
        # contig-joining read contributes to both.
        for hit in group.hits:
            cs = stats.per_contig.get(hit.contig)
            if cs is None:
                continue
            cs.matched_bases += hit.sspan
            if cs.contig.length:
                cs.add_interval(hit.slow, hit.shigh)

        # The read is counted once, on the contig its first HSP lies on.
        home = stats.per_contig.get(group.contig)
        if home is not None:
            home.counts[group.cls] += 1
            if group.inverted:
                home.inverted += 1
        _record_events(stats.per_contig, group)

    return stats, kept
