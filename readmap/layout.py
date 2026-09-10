"""Placing classified reads onto map pages.

Everything here works in *drawing units*: one unit is ``--map-scale``
reference bases, and one unit is also one PDF point, so a page holds
``axis_width = 10000 * page_scale`` units.  That is the same coordinate
system the PostScript used, which keeps the visual result comparable.

Two changes from the Perl:

* Lane packing keeps the intervals a lane holds instead of scanning a
  character string per candidate lane (finding B28), and reads are packed in
  ascending start order, which turns the vertical arrangement from file-order
  scatter into a proper pileup.
* A read may cross any number of page boundaries.  The Perl handled exactly
  one, so a read longer than a page was drawn wrong.
"""

from __future__ import annotations

import bisect
import math
from collections.abc import Iterable

from readmap.classify import drawable_hits
from readmap.models import (
    Contig,
    ContigLayout,
    Hit,
    Junction,
    Label,
    Marker,
    PlacedRead,
    ReadClass,
    ReadGroup,
    Segment,
    SegmentKind,
)
from readmap.params import Params

# Offsets copied from the PostScript so the output stays recognisable.
_ARROW_DX_PAGE_END_UP = -6.0
_ARROW_DX_PAGE_END_DOWN = -3.0
_ARROW_X_NEXT_PAGE = -28.0
_ARROW_X_CONTIG_START = -80.0
_ARROW_DX_CONTIG_END = 60.0
_LABEL_DX = 10.0
_LABEL_DX_AFTER_BREAK = 18.0
_LABEL_DY = -3.0
_EDGE_TOLERANCE = 10.0
"""How close to a contig end a read must reach to earn a terminal arrow."""

_INVERSION_ARROW_DX = -90.0
_INVERSION_ARROW_DY = 55.0


def page_count(contig: Contig, params: Params) -> int:
    """How many pages this contig needs.

    The Perl added ``map_scale`` -- a scale factor -- to a scaled length
    (finding B37).  The quantity actually needed is the horizontal room the
    read labels overhang by, which is ``--space``.
    """
    units = contig.length / params.map_scale + params.space
    return max(1, math.ceil(units / params.axis_width))


def split_across_pages(
    x0: float, x1: float, axis_width: float
) -> list[tuple[int, float, float]]:
    """Cut an absolute x-range into per-page local ranges."""
    if x1 < x0:
        x0, x1 = x1, x0
    first = int(x0 // axis_width)
    last = int(x1 // axis_width)
    out: list[tuple[int, float, float]] = []
    for page in range(first, last + 1):
        lo = max(x0, page * axis_width) - page * axis_width
        hi = min(x1, (page + 1) * axis_width) - page * axis_width
        out.append((page, lo, hi))
    return out


class LanePacker:
    """First-fit lane assignment, one packer per contig and strand.

    Public because the interactive browser packs the visible window with the
    same rule (``readmap/browser/region.py``); there must be exactly one
    implementation or a region would not look like the matching slice of the
    printed map.

    A lane remembers **every** interval it holds, not just its rightmost end.
    While a caller's requests ascend those are the same thing, and that is how
    both callers feed most of them.  The difference is the *class boundary*:
    both place all long-divided reads first, then short-divided, then undivided
    (:func:`_requests`, ``region.CLASS_LANE_ORDER``), so the sweep restarts at
    the left edge twice.  Knowing only a rightmost end, one long-divided read
    near the end of a contig makes its entire lane unusable to every later
    class -- a band of lanes reserved across the whole reference for a handful
    of reads.  Reads bridging the origin of a circular replicon are the worst
    case: by definition they put a piece at *each* end of the contig, so a
    lane they touch is closed from base 1 to the last base.  Remembering the
    intervals lets the later classes fill those holes, which is where most of
    the lanes were going.

    Touching is not overlapping: a request beginning exactly where another ends
    fits beside it.
    """

    __slots__ = ("_claims", "_ends", "_k_max", "overflow")

    def __init__(self, k_max: int) -> None:
        self._k_max = k_max
        self._ends: list[float] = [float("-inf")] * (k_max + 1)
        # Lazily, so a packer offered thousands of lanes costs nothing for the
        # ones it never reaches.  A lane with an end has a claim list.
        self._claims: dict[int, list[tuple[float, float]]] = {}
        self.overflow = 0

    def place(self, begin: float, end: float) -> int:
        for lane in range(1, self._k_max + 1):
            # Past everything this lane holds: the common case, and the only
            # one that can arise while the caller's requests ascend.
            if self._ends[lane] <= begin:
                self._ends[lane] = end
                self._claims.setdefault(lane, []).append((begin, end))
                return lane
            # Behind this lane's right edge -- but the hole may still be free.
            at = self._free_at(lane, begin, end)
            if at is not None:
                self._claims[lane].insert(at, (begin, end))
                return lane
        # No free lane.  The Perl walked off the end of its lane array and died
        # with an uninitialized-value error (finding B03); overplot and count.
        self.overflow += 1
        return self._k_max

    def _free_at(self, lane: int, begin: float, end: float) -> int | None:
        """Index ``[begin, end)`` would take in *lane*, or ``None`` if it collides."""
        claims = self._claims[lane]
        at = bisect.bisect_left(claims, (begin,))
        if at < len(claims) and claims[at][0] < end:
            return None
        if at and claims[at - 1][1] > begin:
            return None
        return at


# --------------------------------------------------------------------------
# Placement requests
# --------------------------------------------------------------------------


class _Request:
    """One thing that needs a lane: a whole read, or one segment of one."""

    __slots__ = ("group", "hits", "begin", "end", "strand")

    def __init__(self, group: ReadGroup, hits: list[Hit], params: Params) -> None:
        self.group = group
        self.hits = hits
        self.strand = hits[0].strand
        lo = min(h.slow for h in hits)
        hi = max(h.shigh for h in hits)
        self.begin = lo / params.map_scale
        self.end = hi / params.map_scale


def _requests(
    groups: Iterable[ReadGroup], contig_name: str, params: Params
) -> list[_Request]:
    """Build placement requests for one contig, in the Perl's drawing order.

    Order matters: long-divided reads are placed first so they land in the
    lanes nearest the axis, where they are easiest to see.  Short-divided
    reads follow, then undivided ones.

    Only the HSPs that lie on *this* contig are placed.  A contig-joining read
    is therefore drawn on both contigs, each showing the piece that belongs
    there -- which is what makes the join visible from either side.
    """
    by_class: dict[ReadClass, list[_Request]] = {
        ReadClass.LONG_DIVIDED: [],
        ReadClass.SHORT_DIVIDED: [],
        ReadClass.UNDIVIDED: [],
    }

    for group in groups:
        hits = [h for h in drawable_hits(group, params) if h.contig == contig_name]
        if not hits:
            continue
        if group.cls is ReadClass.SHORT_DIVIDED:
            # One placement for the whole read: its pieces share a lane.
            by_class[ReadClass.SHORT_DIVIDED].append(_Request(group, hits, params))
        else:
            # Undivided reads have one HSP anyway; long-divided ones are drawn
            # as independent lines carrying the same label, which is how the
            # reader spots both ends of a structural event.
            for hit in hits:
                by_class[group.cls].append(_Request(group, [hit], params))

    ordered: list[_Request] = []
    for cls in (ReadClass.LONG_DIVIDED, ReadClass.SHORT_DIVIDED, ReadClass.UNDIVIDED):
        ordered.extend(sorted(by_class[cls], key=lambda r: (r.begin, r.end)))
    return ordered


# --------------------------------------------------------------------------
# Drawing decisions
# --------------------------------------------------------------------------


def _continuation_style(group: ReadGroup, strand: int) -> str:
    """Which arrow marks a read running off the edge of the drawn region."""
    if group.junction is Junction.CIRCLE:
        return "arrow_circular"
    if group.junction is Junction.CONTIG_JOIN:
        return "arrow_ct"
    return "arrow_fwd" if strand > 0 else "arrow_rev"


def _place_one(
    request: _Request,
    lane: int,
    contig: Contig,
    params: Params,
    label_reads: bool,
) -> PlacedRead:
    group = request.group
    axis_width = params.axis_width
    strand = request.strand
    y = (params.space_to_reads + lane * params.line_space) * strand
    arrow_dx = _ARROW_DX_PAGE_END_UP if y > 0 else _ARROW_DX_PAGE_END_DOWN
    style = _continuation_style(group, strand)

    placed = PlacedRead(group=group, lane=lane)
    ref_strand = request.hits[0].strand

    for hit in request.hits:
        inverted = hit.strand != ref_strand
        kind = SegmentKind.INVERSION if inverted else SegmentKind.READ
        x0 = hit.slow / params.map_scale
        x1 = hit.shigh / params.map_scale

        pieces = split_across_pages(x0, x1, axis_width)
        for index, (page, lo, hi) in enumerate(pieces):
            placed.segments.append(
                Segment(page=page, x0=lo, x1=hi, y=y, kind=kind, cls=group.cls)
            )
            if index < len(pieces) - 1:
                # The line continues on the next page: arrow out, arrow in.
                placed.markers.append(
                    Marker(page=page, x=axis_width + arrow_dx, y=y, style=style, cls=group.cls)
                )
                placed.markers.append(
                    Marker(
                        page=page + 1,
                        x=_ARROW_X_NEXT_PAGE,
                        y=y,
                        style=style,
                        cls=group.cls,
                    )
                )

        if inverted:
            page, lo, _hi = pieces[0]
            placed.markers.append(
                Marker(
                    page=page,
                    x=lo + _INVERSION_ARROW_DX,
                    y=y + _INVERSION_ARROW_DY,
                    style="inversion_arrow",
                    cls=group.cls,
                )
            )

    # Junction tick marks, for reads drawn as several pieces on one lane.
    if group.cls is ReadClass.SHORT_DIVIDED and len(request.hits) > 1:
        for hit in request.hits:
            for coordinate in (hit.sstart, hit.send):
                x = coordinate / params.map_scale
                page = int(x // axis_width)
                placed.markers.append(
                    Marker(
                        page=page,
                        x=x - page * axis_width,
                        y=y,
                        style="tick",
                        cls=group.cls,
                    )
                )

    # Terminal arrows: the read reaches the very start or the very end of the
    # contig.  Restricted to the first page, which is what the Perl meant --
    # it tested a page-local coordinate and so drew the start arrow on every
    # page whose first read began near the fold.
    start_page = int(request.begin // axis_width)
    if start_page == 0 and request.begin < _EDGE_TOLERANCE:
        placed.markers.append(
            Marker(page=0, x=_ARROW_X_CONTIG_START, y=y, style=style, cls=group.cls)
        )

    end_page, _lo, end_local = split_across_pages(request.begin, request.end, axis_width)[-1]
    if contig.length / params.map_scale - request.end < _EDGE_TOLERANCE:
        placed.markers.append(
            Marker(
                page=end_page,
                x=end_local + _ARROW_DX_CONTIG_END,
                y=y,
                style=style,
                cls=group.cls,
            )
        )

    if label_reads:
        crossed = end_page != start_page
        dx = _LABEL_DX_AFTER_BREAK if crossed else _LABEL_DX
        placed.labels.append(
            Label(
                page=end_page,
                x=end_local + dx,
                y=y + _LABEL_DY,
                text=group.label,
                cls=group.cls,
            )
        )

    return placed


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------


def build_layout(
    groups: Iterable[ReadGroup],
    contigs: Iterable[Contig],
    params: Params,
) -> dict[str, ContigLayout]:
    """Place every read.  Returns one :class:`ContigLayout` per contig."""
    contig_list = list(contigs)
    layouts = {
        c.name: ContigLayout(contig=c, pages=page_count(c, params)) for c in contig_list
    }

    # A read is bucketed under every contig it has an HSP on, so that a
    # contig-joining read is drawn on both.
    per_contig: dict[str, list[ReadGroup]] = {c.name: [] for c in contig_list}
    for group in groups:
        for name in {h.contig for h in group.hits}:
            bucket = per_contig.get(name)
            if bucket is not None:
                bucket.append(group)

    for name, layout in layouts.items():
        requests = _requests(per_contig[name], name, params)
        label_reads = len(requests) <= params.max_labels
        packers = {1: LanePacker(params.k_max), -1: LanePacker(params.k_max)}

        for request in requests:
            packer = packers[1 if request.strand > 0 else -1]
            lane = packer.place(request.begin, request.end + params.space)
            layout.placed.append(
                _place_one(request, lane, layout.contig, params, label_reads)
            )

        layout.overflow = packers[1].overflow + packers[-1].overflow

    return layouts
