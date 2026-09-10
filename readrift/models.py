"""Core data types shared by every stage of the pipeline.

Nothing here does I/O, computes science, or draws.  These are the contracts
between ``inputs/``, ``classify``, ``layout``, ``stats`` and ``render/``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

# --------------------------------------------------------------------------
# Alignment
# --------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Hit:
    """One line of a BLAST BTOP file -- a single HSP.

    Column order is fixed by the ``-outfmt`` string the tool documents::

        qseqid sframe qstart qend sstart send qlen sseqid btop

    ``sstart`` may be greater than ``send``; that is how BLAST reports a
    minus-strand match and the drawing code depends on it.
    """

    read: str
    """Original ``qseqid``.  Never rewritten (finding B01)."""

    strand: int
    """``sframe``: ``+1`` or ``-1``."""

    qstart: int
    qend: int
    sstart: int
    send: int
    qlen: int
    contig: str
    btop: str
    """Raw BLAST alignment trace.  Used by :mod:`readrift.btop_trace`."""

    # -- query side --------------------------------------------------------

    @property
    def qlow(self) -> int:
        return min(self.qstart, self.qend)

    @property
    def qhigh(self) -> int:
        return max(self.qstart, self.qend)

    @property
    def qspan(self) -> int:
        return self.qhigh - self.qlow

    # -- subject (reference) side -----------------------------------------

    @property
    def slow(self) -> int:
        return min(self.sstart, self.send)

    @property
    def shigh(self) -> int:
        return max(self.sstart, self.send)

    @property
    def sspan(self) -> int:
        return self.shigh - self.slow


# --------------------------------------------------------------------------
# Classification
# --------------------------------------------------------------------------


class ReadClass(Enum):
    """How a read maps onto the reference."""

    UNDIVIDED = "undivided"
    """A single surviving HSP -- the read maps contiguously."""

    SHORT_DIVIDED = "short_divided"
    """Several HSPs whose reference span matches the read span closely."""

    LONG_DIVIDED = "long_divided"
    """Several HSPs separated by more reference than the read can account for."""


class Junction(Enum):
    """What kind of discontinuity joins the pieces of a divided read."""

    PLAIN = "C"
    """An ordinary junction within one contig."""

    CIRCLE = "Circle"
    """The read spans the whole contig -- evidence the replicon is circular."""

    CONTIG_JOIN = "Contig_Join"
    """The pieces map to different contigs."""


@dataclass(slots=True)
class ReadGroup:
    """One read, after HSP selection and classification.

    ``junction`` and ``inverted`` are decided *here*, at classification time.
    The Perl decided them inside the drawing loops, which is the direct cause
    of findings B04 (wrong read decorated) and B17 (page-crossing inversions
    never counted).
    """

    read: str
    """Original read name."""

    label: str
    """Short display name drawn on the map, e.g. ``R_12345``."""

    hits: list[Hit]
    """Surviving HSPs, in the order they appeared in the BTOP file."""

    cls: ReadClass
    junction: Junction = Junction.PLAIN

    distance: int = 0
    """``|reference span - read span|`` -- the short/long discriminator."""

    inverted: bool = False
    """At least one HSP runs opposite to the read's first HSP."""

    truncated: bool = False
    """Set on the last group kept when ``--cov-max`` stopped the scan."""

    @property
    def contig(self) -> str:
        """Contig of the first HSP.  For a ``CONTIG_JOIN`` the rest differ."""
        return self.hits[0].contig

    @property
    def strand(self) -> int:
        """Strand of the first HSP; sets which side of the axis the read sits."""
        return self.hits[0].strand

    @property
    def qlen(self) -> int:
        return self.hits[0].qlen

    @property
    def matched_bases(self) -> int:
        """Reference bases actually covered by an alignment.

        This is the coverage definition used throughout (decision D3): the
        reference a divided read *skips over* is not covered by it.
        """
        return sum(h.sspan for h in self.hits)

    @property
    def ref_low(self) -> int:
        return min(h.slow for h in self.hits)

    @property
    def ref_high(self) -> int:
        return max(h.shigh for h in self.hits)


# --------------------------------------------------------------------------
# Reference
# --------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Contig:
    """A reference sequence.  Only the length is ever needed (finding B26)."""

    name: str
    length: int


@dataclass(frozen=True, slots=True)
class Annotation:
    """One GenBank feature: what the map draws, and what the browser can show.

    The first seven fields are what a map page needs.  The last three are the
    qualifiers a reader asks for once a feature is worth a second look --
    "what *is* this gene?" -- and are carried only so the browser's detail
    panel can answer that.  ``product`` was already being read to decide
    :attr:`group` and then thrown away.
    """

    contig: str
    strand: int
    start: int
    end: int
    feature: str
    """GenBank feature key, e.g. ``CDS``, ``tRNA``."""

    name: str
    """``/gene``, falling back to ``/locus_tag``."""

    group: str
    """``""`` | ``RNA`` | ``phage`` | ``IS`` | ``hyp`` | ``pseudo``."""

    product: str = ""
    """``/product`` -- what the gene is said to make."""

    locus_tag: str = ""
    """``/locus_tag``, kept separately: when ``/gene`` supplied :attr:`name`,
    this is the only way back to the assembly's own identifier."""

    note: str = ""
    """``/note``, where curators put whatever did not fit anywhere else."""


@dataclass(slots=True)
class Metadata:
    """Project information scraped from the GenBank header.

    Every field defaults to the empty string; the Perl left them undefined and
    died on the first interpolation when a record lacked them (finding B13).
    """

    organism: str = ""
    bio_project: str = ""
    bio_sample: str = ""
    sra: str = ""
    assembly_method: str = ""
    sequencing_technology: str = ""
    source_format: str = "fa"
    """``"gb"`` or ``"fa"``."""

    first_line: str = ""
    """First line of the reference file, shown for FASTA input."""

    def rows(self) -> list[tuple[str, str]]:
        """Non-empty ``(label, value)`` pairs, for the front page."""
        candidates = [
            ("Organism", self.organism),
            ("BioProject", self.bio_project),
            ("BioSample", self.bio_sample),
            ("SRA", self.sra),
            ("Assembly method", self.assembly_method),
            ("Sequencing technology", self.sequencing_technology),
        ]
        return [(k, v) for k, v in candidates if v]


# --------------------------------------------------------------------------
# Layout
# --------------------------------------------------------------------------


class SegmentKind(Enum):
    """What a drawn horizontal line represents."""

    READ = "read"
    """A normal aligned segment."""

    INVERSION = "inversion"
    """A segment running opposite to the rest of its read."""


@dataclass(frozen=True, slots=True)
class Segment:
    """A horizontal line ready to draw, in page-local drawing units."""

    page: int
    x0: float
    x1: float
    y: float
    kind: SegmentKind
    cls: ReadClass


@dataclass(frozen=True, slots=True)
class Marker:
    """A decoration attached to a placed read."""

    page: int
    x: float
    y: float
    style: str
    """``arrow_fwd`` | ``arrow_rev`` | ``arrow_plain`` | ``arrow_ct`` |
    ``arrow_circular`` | ``tick`` | ``inversion_arrow``."""

    cls: ReadClass


@dataclass(frozen=True, slots=True)
class Label:
    """A read name drawn at the end of its line."""

    page: int
    x: float
    y: float
    text: str
    cls: ReadClass


@dataclass(slots=True)
class PlacedRead:
    """Everything needed to draw one read, already split across pages."""

    group: ReadGroup
    lane: int
    segments: list[Segment] = field(default_factory=list)
    markers: list[Marker] = field(default_factory=list)
    labels: list[Label] = field(default_factory=list)


@dataclass(slots=True)
class ContigLayout:
    """All placed reads for one contig, plus how many pages they need."""

    contig: Contig
    pages: int
    placed: list[PlacedRead] = field(default_factory=list)
    overflow: int = 0
    """Reads that found no free lane and were drawn on the last one."""
