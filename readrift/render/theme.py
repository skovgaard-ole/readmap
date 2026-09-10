"""Colours, geometry and arrow shapes.

Replaces the Perl's ``define_PS_settings``.

**Palette provenance.** The read-class colours are taken from the data-viz
reference palette's pre-validated slots rather than carried over verbatim from
the PostScript.  The original used ``0.5 0.0 0.0`` for long-divided reads and
``0.8 0.0 0.0`` for inversions -- two dark reds that are near-indistinguishable
under any colour-vision deficiency, and both below the lightness band against
white.  The three *class* colours here are the reference palette's first three
slots, which are documented as clearing the all-pairs gate in light mode; the
inversion colour is the palette's fixed ``critical`` status step, which is the
right role for it -- an inversion is a *state* of a divided read, not a fourth
peer class.

Two consequences follow from that choice and are honoured throughout:

* Status colour never carries meaning alone: inversions are always drawn with
  the rotated inversion arrow beside them, and always labelled in the legend.
* Slot 3 (the green used for undivided reads) sits below 3:1 contrast on a
  white surface, so the relief rule applies -- every plot that uses it also
  carries direct labels or a value table.

Annotation colours are kept close to the original, because they are already
redundant with a symbol suffix (``*`` RNA, ``#`` phage, ``+`` IS, ``^`` pseudo,
``?`` hypothetical) and the user reads them by that symbol.
"""

from __future__ import annotations

import math

from readrift.models import ReadClass

# --------------------------------------------------------------------------
# Surface and ink
# --------------------------------------------------------------------------

SURFACE = "#ffffff"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#333333"
FRAME = "#e1e0d9"
TITLE = "#1a0099"

# --------------------------------------------------------------------------
# Read classes
# --------------------------------------------------------------------------

CLASS_COLORS: dict[ReadClass, str] = {
    ReadClass.UNDIVIDED: "#1baf7a",
    ReadClass.SHORT_DIVIDED: "#2a78d6",
    ReadClass.LONG_DIVIDED: "#eb6834",
}

INVERSION = "#d03b3b"

CLASS_LABELS: dict[ReadClass, str] = {
    ReadClass.UNDIVIDED: "Undivided reads",
    ReadClass.SHORT_DIVIDED: "Divided reads, short distance",
    ReadClass.LONG_DIVIDED: "Divided reads, long distance",
}

CLASS_ORDER: tuple[ReadClass, ...] = (
    ReadClass.UNDIVIDED,
    ReadClass.SHORT_DIVIDED,
    ReadClass.LONG_DIVIDED,
)

EVENT_COLORS: dict[str, str] = {
    "junction": CLASS_COLORS[ReadClass.LONG_DIVIDED],
    "inversion": INVERSION,
    "circle": "#4a3aa7",
    "contig_join": "#0b0b0b",
}

EVENT_LABELS: dict[str, str] = {
    "junction": "Long-distance junction",
    "inversion": "Inversion breakpoint",
    "circle": "Circular join",
    "contig_join": "Contig join",
}

#: One sentence saying what each mark is evidence *of*.  Here rather than in
#: the front end for the same reason the colours are: a mark drawn from
#: :data:`EVENT_COLORS` must not be explained differently depending on which
#: surface the reader is looking at.  ``stats._record_events`` decides where
#: each one is placed; these say what the reader should make of it.
EVENT_DESCRIPTIONS: dict[str, str] = {
    "junction": (
        "A read continues at a distant position: it crosses more reference "
        "than its own length accounts for. Deletion, transposition or "
        "rearrangement."
    ),
    "inversion": (
        "A piece of a read aligns against the run of the rest of it. The two "
        "orientations meet here — an inversion breakpoint."
    ),
    "circle": (
        "A read spans essentially the whole replicon, leaving one end and "
        "arriving at the other: evidence the sequence is circular."
    ),
    "contig_join": (
        "A read continues onto a different contig. The last base it reaches "
        "on this one, and evidence the two belong together."
    ),
}

# --------------------------------------------------------------------------
# Annotations
# --------------------------------------------------------------------------

ANNOTATION_COLORS: dict[str, str] = {
    "": "#000099",
    "RNA": "#cc0000",
    "pseudo": "#66004d",
    "hyp": "#663300",
    "phage": "#b32600",
    "IS": "#cc3333",
}

ANNOTATION_SYMBOLS: dict[str, str] = {
    "": "",
    "RNA": "*",
    "pseudo": "^",
    "hyp": "?",
    "phage": "#",
    "IS": "+",
}

ANNOTATION_LEGEND: tuple[tuple[str, str], ...] = (
    ("", "Gene, other"),
    ("RNA", "* RNA coding gene"),
    ("phage", "# Phage related gene"),
    ("IS", "+ Transposase related gene"),
    ("pseudo", "^ Pseudogene"),
    ("hyp", "? Hypothetical gene"),
)

# --------------------------------------------------------------------------
# Page geometry
# --------------------------------------------------------------------------

#: Front page and plot pages: A3 landscape, in points.  The map pages keep the
#: original poster format (decision D5).
A3_LANDSCAPE = (1190.55, 841.89)

MONO = "monospace"
SANS = "sans-serif"

TITLE_FONT_SIZE = 24.0
LABEL_FONT_SIZE = 12.0
ANNOTATION_FONT_SIZE = 3.0
AXIS_FONT_SIZE = 12.0

#: Vertical offsets of the two annotation tracks, relative to the axis.
ANNOTATION_Y_PLUS = -18.0
ANNOTATION_Y_MINUS = -44.0
ANNOTATION_ROW_HEIGHT = 6.0
ANNOTATION_ROWS = 3

# --------------------------------------------------------------------------
# Arrow shapes
# --------------------------------------------------------------------------
#
# Reproduced from the PostScript macros in PS_init.  Each returns a list of
# (x, y) points forming one polyline, relative to the anchor point.  The
# PostScript prefixed every arrow with `12 0 rmoveto`, which is folded in here.


def _chevron(width: float, height: float, direction: int) -> list[tuple[float, float]]:
    """One arrowhead: a short shaft, then a two-stroke chevron at its end."""
    d = 1.0 if direction >= 0 else -1.0
    return [
        (12.0 + 2.0, 0.0),
        (12.0 + 2.0 + d * (width - 1.0), 0.0),
        (12.0 + 2.0 + d * (width - 1.0) - d * width, height),
        (12.0 + 2.0 + d * (width - 1.0) - d * width + d * 1.0, 0.0),
        (12.0 + 2.0 + d * (width - 1.0) - d * width, -height),
        (12.0 + 2.0 + d * (width - 1.0), 0.0),
    ]


def arrow_points(style: str, line_width: float) -> list[tuple[float, float]]:
    """Polyline for an arrow *style*, in points relative to its anchor."""
    w = line_width * 2.5
    h = w * 0.5

    if style == "arrow_fwd":
        return _chevron(w, h, +1)
    if style == "arrow_rev":
        return _chevron(w, h, -1)

    if style in ("arrow_plain", "inversion_arrow", "arrow_circular"):
        # `arrow`: a 20-unit shaft with the chevron at its far end.
        return [
            (12.0, 0.0),
            (32.0, 0.0),
            (32.0 - w, h),
            (32.0 - w + 1.0, 0.0),
            (32.0 - w, -h),
            (32.0, 0.0),
        ]

    if style == "arrow_ct":
        # `arrow_ct`: two chevrons -- the read joins another contig's end.
        return [
            (12.0, 0.0),
            (32.0, 0.0),
            (32.0 - w, h),
            (32.0 - w + 1.0, 0.0),
            (32.0 - w, -h),
            (32.0, 0.0),
            (44.0, 0.0),
            (44.0 - w, h),
            (44.0 - w + 1.0, 0.0),
            (44.0 - w, -h),
            (44.0, 0.0),
        ]

    raise ValueError(f"unknown arrow style: {style!r}")


#: Where the ring of `arrow_circular` sits, relative to the anchor.
CIRCULAR_RING_X = 62.0

#: Rotation and scale of the inversion arrow, copied from the PostScript.
INVERSION_ARROW_ROTATION = math.radians(-30.0)
INVERSION_ARROW_SCALE = 3.0


def transform(
    points: list[tuple[float, float]],
    x: float,
    y: float,
    scale: float = 1.0,
    rotation: float = 0.0,
) -> list[tuple[float, float]]:
    """Scale, rotate and translate an arrow polyline."""
    if rotation:
        cos_r, sin_r = math.cos(rotation), math.sin(rotation)
    else:
        cos_r, sin_r = 1.0, 0.0

    out: list[tuple[float, float]] = []
    for px, py in points:
        sx, sy = px * scale, py * scale
        rx = sx * cos_r - sy * sin_r
        ry = sx * sin_r + sy * cos_r
        out.append((x + rx, y + ry))
    return out


def thousands(value: float) -> str:
    """``1234567`` -> ``1,234,567``."""
    return f"{int(round(value)):,}"
