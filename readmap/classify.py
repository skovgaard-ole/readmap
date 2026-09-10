"""HSP selection and read classification.

This module is a faithful reimplementation of the algorithm in the Perl
original (``Sort_reads`` / ``_process_bt_read`` / ``_process_multi_quest``),
documented step by step in ``code_structure.md`` §5.  The numbers it
produces are the scientific output of the whole program, so the *logic* here
is deliberately unchanged; only its implementation and its parameterisation
are.

Two things moved here from the drawing code, where they never belonged:

* the ``Circle`` / ``Contig_Join`` junction tag, which the Perl read from a
  stale loop variable and therefore applied to the wrong reads (finding B04);
* the inversion flag, which the Perl only noticed for segments that happened
  to fall inside one page (finding B17).
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field, replace

from readmap.inputs.reference import Reference
from readmap.labels import LabelAssigner
from readmap.models import Hit, Junction, ReadClass, ReadGroup
from readmap.params import Params

#: Two HSPs whose query spans differ by less than this are treated as
#: alternative placements of the same repeat.  Hardcoded in the Perl.
REPEAT_SPAN_TOLERANCE = 50

#: How close the reference span must come to the contig length before the read
#: is taken as evidence of circularity.  Hardcoded in the Perl.
CIRCLE_TOLERANCE = 20

#: How much of an opposite-strand HSP's *reference* span must already have been
#: mapped by the same read before the HSP is treated as a fold-back artefact.
#:
#: Not a small tolerance but a majority test, and deliberately so.  Inversions
#: in real genomes very often occur between inverted repeats, so a read
#: crossing a genuine breakpoint can legitimately share a few hundred bases of
#: reference with the piece before it.  Requiring most of the HSP to be
#: re-covered separates "read the same stretch again, backwards" -- where the
#: overlap is essentially total -- from "turned a corner at a repeat".
FOLD_BACK_OVERLAP = 0.5


@dataclass(slots=True)
class ClassifyStats:
    """Bookkeeping for the run log and the PDF stamp."""

    reads_in: int = 0
    reads_classified: int = 0
    reads_without_hits: int = 0
    hits_dropped_short: int = 0
    hits_dropped_overlap: int = 0
    hits_dropped_fold_back: int = 0
    reads_folded_back: int = 0
    unknown_contig_reads: int = 0
    unknown_contigs: set[str] = field(default_factory=set)
    coverage_used: float = 0.0
    truncated_by_cov_max: bool = False

    def warnings(self, reference: Reference) -> list[str]:
        notes: list[str] = []
        if self.unknown_contigs:
            shown = ", ".join(sorted(self.unknown_contigs)[:8])
            more = "" if len(self.unknown_contigs) <= 8 else f" (+{len(self.unknown_contigs) - 8} more)"
            notes.append(
                f"{self.unknown_contig_reads} read(s) map to "
                f"{len(self.unknown_contigs)} accession(s) that are not in the "
                f"reference and were dropped: {shown}{more}. "
                f"Reference has: {', '.join(reference.names[:8])}. "
                f"The Perl dropped these silently (finding B11)."
            )
        if self.reads_folded_back:
            notes.append(
                f"{self.reads_folded_back} read(s) fold back on themselves: "
                f"{self.hits_dropped_fold_back} HSP(s) re-read reference the "
                f"same read had already mapped, in the opposite direction. "
                f"Treated as an end-ligation artefact and dropped, so they are "
                f"not counted as junctions or inversions and their bases are "
                f"not counted twice. Pass --keep-fold-back to keep them."
            )
        if self.truncated_by_cov_max:
            notes.append(
                f"Stopped after {self.coverage_used:.2f}x coverage because of "
                f"--cov-max; the rest of the BTOP file was not read."
            )
        return notes


# --------------------------------------------------------------------------
# Step 1-3: choose which HSPs represent the read
# --------------------------------------------------------------------------


def collapse_repeats(hits: list[Hit]) -> list[Hit]:
    """Step 2 -- drop alternative placements of the same repeated segment.

    When two consecutive HSPs cover almost the same stretch of the *read*,
    they are competing placements of one repeat.  The one whose reference
    position sits further from where the read started is discarded.

    The Perl expressed the discard as ``$bt_read[$j-1] = $bt_read[$j]`` -- it
    duplicated the survivor and relied on :func:`select_non_overlapping` to
    throw the copy away.  That is reproduced here, including the comparison
    itself, which adds a read length to a reference-coordinate difference
    (finding B20).  It is preserved verbatim because changing it would change
    published numbers; it is flagged, not fixed.
    """
    if len(hits) < 3:
        return list(hits)

    out = list(hits)
    origin = hits[0].sstart
    prev_span = hits[0].qspan
    prev_hit = hits[0]

    for j in range(1, len(hits)):
        hit = hits[j]
        span = hit.qspan
        if j >= 2 and abs(prev_span - span) < REPEAT_SPAN_TOLERANCE:
            here = abs(prev_hit.sstart - origin)
            there = abs(hit.sstart - origin + hit.qlen)
            if here > there:
                out[j - 1] = hit
        prev_span = span
        prev_hit = hit

    return out


def select_non_overlapping(
    hits: Iterable[Hit], tolerance: float
) -> tuple[list[Hit], int]:
    """Step 3 -- greedily keep HSPs that claim fresh parts of the read.

    Walk the HSPs in file order and keep one only if none of its query bases
    are already claimed.  On keeping it, claim ``[qstart + tol, qend - tol]``:
    the margin lets a short repeat shared by two junctions be used twice.

    The Perl did this with a per-base array and a loop over every base of the
    read -- O(read length) per HSP, so ~100 000 iterations per HSP on a long
    ONT read (finding B27).  This is the same result by interval arithmetic.
    """
    claimed: list[tuple[int, int]] = []
    kept: list[Hit] = []
    dropped = 0

    for hit in hits:
        lo, hi = hit.qlow, hit.qhigh
        if any(lo <= end and hi >= start for start, end in claimed):
            dropped += 1
            continue

        # math.floor matches Perl's integer range truncation.  When the
        # tolerance is wider than the HSP the range inverts and nothing is
        # claimed -- the Perl iterated zero times, so do the same (finding B38).
        start = math.floor(lo + tolerance)
        end = math.floor(hi - tolerance)
        if start <= end:
            claimed.append((start, end))
        kept.append(hit)

    return kept, dropped


def _covered_fraction(lo: int, hi: int, claimed: list[tuple[int, int]]) -> float:
    """How much of ``[lo, hi)`` the union of *claimed* already covers, 0-1.

    The claims can overlap each other -- two same-strand HSPs may both cover a
    tandem duplication -- so this merges the intersections instead of adding
    them up, which would otherwise report more than 100% coverage.
    """
    pieces = sorted(
        (max(lo, start), min(hi, end))
        for start, end in claimed
        if start < hi and end > lo
    )
    if not pieces:
        return 0.0

    span = hi - lo
    if span <= 0:
        # A degenerate HSP that something already claims. Nothing to divide by,
        # and it is inside a claim, so call it fully covered.
        return 1.0

    covered = 0
    reach = lo
    for start, end in pieces:
        if end <= reach:
            continue
        covered += end - max(start, reach)
        reach = end
    return covered / span


def drop_fold_back(
    hits: list[Hit], min_overlap: float = FOLD_BACK_OVERLAP
) -> tuple[list[Hit], int]:
    """Step 3b -- drop the fold-back tail of an end-ligation artefact.

    A long read is sometimes sequenced through its template and then straight
    back along it: adapter or end ligation joins the molecule to its own
    reverse complement, and the resulting read maps as one long HSP followed by
    a run of HSPs that re-read *the same stretch of reference* in the opposite
    direction.  The signature is unmistakable::

        ctgA  1,248,192-1,269,037  -   read      1- 21,000     <- the molecule
        ctgA  1,248,422-1,249,382  +   read 21,145- 22,076     <- reading back
        ctgA  1,249,517-1,250,021  +   read 22,129- 22,578
        ...

    The Perl's ``@taken`` rule (:func:`select_non_overlapping`) cannot see it:
    every one of those HSPs claims *fresh read bases*, so all of them are
    legitimately kept.  It is the reference that is being covered twice.

    Left in, one such read is classified ``LONG_DIVIDED``, flagged inverted,
    counted as a junction and an inversion at every fold, and its bases are
    added to coverage twice.  A handful of them can invent a structural signal
    that is not in the sample.

    So: walking in file order, an HSP is dropped when it runs **opposite to the
    read's first HSP** *and* at least *min_overlap* of its reference span is
    already covered by HSPs kept for this read on the same contig.  Both halves
    are needed -- the strand test alone would delete genuine inversions, and
    the overlap test alone would delete genuine tandem duplications.

    File order is what makes the first HSP the anchor.  BLAST writes HSPs best
    first, so the real molecule leads and the artefact follows; the same
    assumption the ``@taken`` rule has always rested on.
    """
    if len(hits) < 2:
        return list(hits), 0

    anchor = hits[0].strand
    claimed: dict[str, list[tuple[int, int]]] = {}
    kept: list[Hit] = []
    dropped = 0

    for hit in hits:
        here = claimed.setdefault(hit.contig, [])
        folded = (
            hit.strand != anchor
            and _covered_fraction(hit.slow, hit.shigh, here) >= min_overlap
        )
        if folded:
            dropped += 1
            continue
        here.append((hit.slow, hit.shigh))
        kept.append(hit)

    return kept, dropped


def select_hits(hits: list[Hit], params: Params) -> tuple[list[Hit], int, int, int]:
    """Steps 1-3b.  ``(kept, dropped_short, dropped_overlap, dropped_fold_back)``."""
    cut = params.hsp_cut
    long_enough = [h for h in hits if h.qspan > cut]
    dropped_short = len(hits) - len(long_enough)

    if not long_enough:
        return [], dropped_short, 0, 0

    collapsed = collapse_repeats(long_enough)
    kept, dropped_overlap = select_non_overlapping(collapsed, params.junction_tolerance)

    dropped_fold_back = 0
    if not params.keep_fold_back:
        kept, dropped_fold_back = drop_fold_back(kept)

    return kept, dropped_short, dropped_overlap, dropped_fold_back


# --------------------------------------------------------------------------
# Step 4-5: classify
# --------------------------------------------------------------------------


def classify_hits(
    kept: list[Hit],
    label: str,
    contig_lengths: dict[str, int],
    params: Params,
) -> ReadGroup | None:
    """Steps 4 and 5 -- turn a set of surviving HSPs into a classified read."""
    if not kept:
        return None

    read = kept[0].read

    if len(kept) == 1:
        return ReadGroup(
            read=read,
            label=label,
            hits=kept,
            cls=ReadClass.UNDIVIDED,
            junction=Junction.PLAIN,
            distance=0,
            inverted=False,
        )

    map_min = min(h.slow for h in kept)
    map_max = max(h.shigh for h in kept)
    read_min = min(h.qlow for h in kept)
    read_max = max(h.qhigh for h in kept)

    # How much reference the read skips over, relative to its own length.
    # Small for a local rearrangement, large for a structural one.
    distance = abs((map_max - map_min) - (read_max - read_min))

    first_contig = kept[0].contig
    if any(h.contig != first_contig for h in kept[1:]):
        junction = Junction.CONTIG_JOIN
    else:
        length = contig_lengths.get(first_contig)
        if length is not None and abs(length - (map_max - map_min)) < CIRCLE_TOLERANCE:
            junction = Junction.CIRCLE
        else:
            junction = Junction.PLAIN

    cls = (
        ReadClass.SHORT_DIVIDED
        if distance < params.division_cut
        else ReadClass.LONG_DIVIDED
    )

    return ReadGroup(
        read=read,
        label=label,
        hits=kept,
        cls=cls,
        junction=junction,
        distance=distance,
        inverted=any(h.strand != kept[0].strand for h in kept),
    )


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------


def classify_stream(
    groups: Iterable[list[Hit]],
    reference: Reference,
    params: Params,
    labeller: LabelAssigner | None = None,
    stats: ClassifyStats | None = None,
) -> Iterator[ReadGroup]:
    """Classify every read in *groups*, resolving contigs and honouring ``--cov-max``."""
    st = stats if stats is not None else ClassifyStats()
    lab = labeller if labeller is not None else LabelAssigner()
    contig_lengths = {c.name: c.length for c in reference.contigs}
    total_length = reference.total_length or 1

    for hits in groups:
        st.reads_in += 1

        resolved: list[Hit] = []
        unknown = False
        for hit in hits:
            canonical = reference.resolve(hit.contig)
            if canonical is None:
                st.unknown_contigs.add(hit.contig)
                unknown = True
                continue
            resolved.append(hit if canonical == hit.contig else replace(hit, contig=canonical))

        if unknown:
            st.unknown_contig_reads += 1
        if not resolved:
            continue

        kept, dropped_short, dropped_overlap, dropped_fold_back = select_hits(
            resolved, params
        )
        st.hits_dropped_short += dropped_short
        st.hits_dropped_overlap += dropped_overlap
        st.hits_dropped_fold_back += dropped_fold_back
        if dropped_fold_back:
            st.reads_folded_back += 1

        group = classify_hits(kept, lab.label(resolved[0].read), contig_lengths, params)
        if group is None:
            st.reads_without_hits += 1
            continue

        st.reads_classified += 1
        st.coverage_used += group.qlen / total_length

        # The Perl added the *incoming* read's length and then returned without
        # flushing the group it was holding (finding B16).  Here the group that
        # crosses the limit is emitted, then the scan stops.
        if params.cov_max and st.coverage_used > params.cov_max:
            group.truncated = True
            st.truncated_by_cov_max = True
            yield group
            return

        yield group


def drawable_hits(group: ReadGroup, params: Params) -> list[Hit]:
    """HSPs of *group* that are long enough to draw.

    Long-divided reads hide segments below ``--min-match-length`` so the map
    is not cluttered by the short anchors either side of a junction.  Note
    that those segments still count towards coverage -- they are matched
    bases (decision D3).
    """
    if group.cls is not ReadClass.LONG_DIVIDED:
        return group.hits
    return [h for h in group.hits if h.sspan > params.min_match_length]
