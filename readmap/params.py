"""Parameter definitions for readmap.

The :data:`OPTIONS` table is the *single* source of truth for command-line
options.  ``cli.py`` builds the ``argparse`` parser from it and
``render/frontpage.py`` builds the parameter listing from it, so the two can
never drift apart.

(The Perl original kept ``@PARAM_META`` for the help text and a separate
``GetOptions`` call for the actual parsing, which is how ``--space2reads``
came to be documented but not implemented -- finding B18.)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

Kind = Literal["int", "float", "flag", "str", "append"]


@dataclass(frozen=True, slots=True)
class Option:
    """One command-line option, as shown in help and on the PDF front page."""

    key: str
    """Attribute name on :class:`Params`."""

    long: str
    """Primary long flag, without the leading ``--``."""

    help: str

    default: Any

    kind: Kind = "int"

    short: str | None = None
    """Legacy single-letter flag from the Perl version.  Preserved exactly."""

    aliases: tuple[str, ...] = ()
    """Additional accepted long flags (without ``--``)."""

    metavar: str | None = None

    on_frontpage: bool = True
    """Whether to list this option in the PDF parameter block."""

    @property
    def flags(self) -> list[str]:
        out: list[str] = []
        if self.short:
            out.append(f"-{self.short}")
        out.append(f"--{self.long}")
        out.extend(f"--{a}" for a in self.aliases)
        return out

    @property
    def display(self) -> str:
        """Compact ``-x  --long`` form used in help and on the front page."""
        if self.short:
            return f"-{self.short}  --{self.long}"
        return f"    --{self.long}"


# --------------------------------------------------------------------------
# Option table
# --------------------------------------------------------------------------
#
# Order matters: it is the order shown by --help and on the front page.

OPTIONS: tuple[Option, ...] = (
    # ---- filtering -------------------------------------------------------
    Option(
        key="min_read_length",
        short="r",
        long="min-read-length",
        help="Minimum total read length (qlen) to consider",
        default=10_000,
    ),
    Option(
        key="min_match_length",
        short="m",
        long="min-match-length",
        help="Minimum aligned span of a drawn segment",
        default=2_000,
    ),
    Option(
        key="min_hsp_length",
        long="min-hsp-length",
        help="Minimum query span of a single BLAST HSP",
        default=100,
    ),
    Option(
        key="micro_match_length",
        short="y",
        long="micro-match-length",
        help="Repeat tolerance at junctions (tol = value/4)",
        default=100,
    ),
    Option(
        key="division_threshold",
        long="division-threshold",
        help="Short/long divided cut-off (default: min-match-length)",
        default=None,
        metavar="N",
    ),
    Option(
        key="cov_max",
        short="x",
        long="cov-max",
        help="Stop after this much coverage (0 = unlimited)",
        default=0,
    ),
    # ---- geometry --------------------------------------------------------
    Option(
        key="map_scale",
        short="u",
        long="map-scale",
        help="Reference bases per drawing unit",
        default=100.0,
        kind="float",
    ),
    Option(
        key="page_scale",
        short="p",
        long="page-scale",
        help="Page size multiplier",
        default=1.0,
        kind="float",
    ),
    Option(
        key="space",
        short="s",
        long="space",
        help="Horizontal gap reserved after each read",
        default=100,
    ),
    Option(
        key="line_space",
        short="c",
        long="line-space",
        help="Vertical distance between read lanes",
        default=15,
    ),
    Option(
        key="space_to_reads",
        short="v",
        long="space-to-reads",
        aliases=("space2reads",),
        help="Vertical gap between the axis and the first lane",
        default=70,
    ),
    Option(
        key="line_width",
        short="l",
        long="line-width",
        aliases=("l-w",),
        help="Base line width in points",
        default=2,
    ),
    # ---- output ----------------------------------------------------------
    Option(
        key="index_window",
        long="index-window",
        help="Window size (bp) of the read-extraction index",
        default=1_000,
        on_frontpage=False,
    ),
    Option(
        key="max_labels",
        long="max-labels",
        help="Suppress read labels above this many reads per page",
        default=20_000,
        on_frontpage=False,
    ),
)

#: Options above that carry a numeric value, i.e. the ones the Perl version
#: listed in ``@PARAM_META``.  Used for the front-page parameter block.
NUMERIC_OPTIONS: tuple[Option, ...] = tuple(
    o for o in OPTIONS if o.kind in ("int", "float")
)


# --------------------------------------------------------------------------
# Resolved parameters
# --------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class ExtractSpec:
    """One ``-e FILE,CONTIG,START,END`` request."""

    reads_file: str
    contig: str
    start: int
    end: int

    @classmethod
    def parse(cls, raw: str) -> ExtractSpec:
        parts = [p.strip() for p in raw.split(",")]
        if len(parts) != 4:
            raise ValueError(
                f"--extract-reads needs FILE,CONTIG,START,END (got {len(parts)} "
                f"field(s) in {raw!r})"
            )
        try:
            start, end = int(parts[2]), int(parts[3])
        except ValueError as exc:
            raise ValueError(
                f"--extract-reads START and END must be integers (in {raw!r})"
            ) from exc
        if start > end:
            start, end = end, start
        return cls(reads_file=parts[0], contig=parts[1], start=start, end=end)

    @property
    def stem(self) -> str:
        return f"extract-reads-list_{self.contig}_{self.start}_{self.end}"


@dataclass(frozen=True, slots=True)
class Params:
    """Fully resolved run parameters.

    Field names match ``Option.key`` for everything in :data:`OPTIONS`; the
    remaining fields come from positional arguments and boolean switches.
    """

    # positional
    seq_file: str
    btop_file: str

    # numeric (see OPTIONS)
    min_read_length: int = 10_000
    min_match_length: int = 2_000
    min_hsp_length: int = 100
    micro_match_length: int = 100
    division_threshold: int | None = None
    cov_max: int = 0
    map_scale: float = 100.0
    page_scale: float = 1.0
    space: int = 100
    line_space: int = 15
    space_to_reads: int = 70
    line_width: int = 2
    index_window: int = 1_000
    max_labels: int = 20_000

    # output / behaviour
    out: str = ""
    extract_reads: tuple[ExtractSpec, ...] = ()
    report_analysis: bool = False
    no_plots: bool = False
    identity: bool = False
    """Parse the BTOP trace column and plot alignment identity."""

    no_cache: bool = False
    """Skip writing ``<out>.readmapdb.npz``, the interactive browser's cache.

    The cache is written by default: reading and classifying a multi-gigabyte
    BTOP file costs minutes, writing the cache costs seconds, and without it
    every browsing session would have to repeat the expensive half.
    """

    no_pdf: bool = False
    """Write only the cache and the optional report -- skip the map document."""

    keep_fold_back: bool = False
    """Keep the fold-back tail of an end-ligation artefact instead of dropping it.

    A read sequenced through its template and back along it maps as one long
    HSP followed by a run of HSPs re-reading the same reference backwards.
    ``classify.drop_fold_back`` removes those by default, because left in they
    are counted as junctions and inversions and their bases are added to
    coverage twice.  Set this to keep every HSP the ``@taken`` rule kept, which
    is what every version before this one did.
    """

    compat: bool = False
    allow_unknown_contigs: bool = False
    debug: bool = False

    # provenance, for the PDF stamp
    argv: tuple[str, ...] = field(default_factory=tuple)

    # ---- derived ---------------------------------------------------------

    @property
    def axis_width(self) -> float:
        """Length of the drawn x-axis, in PDF points."""
        return 10_000.0 * self.page_scale

    @property
    def page_width(self) -> float:
        return 10_750.0 * self.page_scale

    @property
    def page_height(self) -> float:
        return 7_600.0 * self.page_scale

    @property
    def k_max(self) -> int:
        """Number of usable read lanes above (and below) the axis."""
        return int(3_700 / self.line_space) - 8

    @property
    def division_cut(self) -> int:
        """Resolved short/long divided threshold.

        Defaults to ``min_match_length``, which is what the Perl used --
        the same option controlled three unrelated things (finding B15).
        """
        return (
            self.min_match_length
            if self.division_threshold is None
            else self.division_threshold
        )

    @property
    def hsp_cut(self) -> int:
        """Minimum query span of an HSP.

        The Perl hardcoded ``100`` here while advertising
        ``--micro-match-length`` (finding B05).  ``--compat`` restores that.
        """
        return 100 if self.compat else self.min_hsp_length

    @property
    def junction_tolerance(self) -> float:
        """Repeat tolerance when claiming query bases, in bases."""
        return self.micro_match_length / 4.0

    # ---- validation ------------------------------------------------------

    def validate(self) -> None:
        """Raise :class:`ValueError` for parameter combinations that cannot work."""
        problems: list[str] = []

        if self.map_scale <= 0:
            problems.append("--map-scale must be > 0")
        if self.page_scale <= 0:
            problems.append("--page-scale must be > 0")
        if self.line_space <= 0:
            problems.append("--line-space must be > 0")
        elif self.k_max < 1:
            problems.append(
                f"--line-space {self.line_space} leaves no room for reads "
                f"(k_max = {self.k_max}); use a value below "
                f"{int(3_700 / 9)} or increase --page-scale"
            )
        if self.line_width <= 0:
            problems.append("--line-width must be > 0")
        if self.space < 0:
            problems.append("--space must be >= 0")
        if self.index_window <= 0:
            problems.append("--index-window must be > 0")
        if self.cov_max < 0:
            problems.append("--cov-max must be >= 0 (0 means unlimited)")
        if self.micro_match_length < 0:
            problems.append("--micro-match-length must be >= 0")

        if problems:
            raise ValueError("; ".join(problems))

    # ---- presentation ----------------------------------------------------

    def as_frontpage_rows(self) -> list[tuple[str, str, str]]:
        """``(flag, value, default)`` triples for the PDF parameter block."""
        rows: list[tuple[str, str, str]] = []
        for opt in NUMERIC_OPTIONS:
            if not opt.on_frontpage:
                continue
            value = getattr(self, opt.key)
            if opt.key == "division_threshold":
                value = self.division_cut
            rows.append((opt.display, _fmt(value), _fmt(opt.default)))
        return rows


def _fmt(value: Any) -> str:
    if value is None:
        return "auto"
    if isinstance(value, float):
        return f"{value:g}"
    if isinstance(value, bool):
        return "on" if value else "off"
    return str(value)
