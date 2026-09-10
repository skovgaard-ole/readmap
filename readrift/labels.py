"""Short display names for reads, and for genes that have no real name.

The Perl rewrote every read name in place with::

    $btop_read[0] =~ /\\w+[\\._](\\d+)/;
    $btop_read[0] = "R_" . $1;

The match result was never checked, so when a name did not fit the pattern
``$1`` still held the *previous* read's capture and the read silently
inherited the previous read's identity -- its HSPs were merged into that
read's group.  Oxford Nanopore names are UUIDs (``0a1b2c3d-e4f5-...``) and
never match, so on an ONT dataset the whole file collapses (finding B01).

Here the read name is never touched.  The ``R_<n>`` form survives only as a
*label* drawn on the map, and collisions are made visible rather than silently
merging two reads.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

#: The Perl pattern, reproduced exactly so labels look the same for the
#: SRA-style names it was written for (``DRR325755.12345`` -> ``R_12345``).
LEGACY_PATTERN = re.compile(r"\w+[._](\d+)")

#: A name shaped like ``<prefix><separator><number>``, the number optionally
#: carrying one leading letter.
#:
#: Requires a letter in the prefix and digits to the end, so ``SS37A_41660``
#: matches and ``recA`` or ``yjdC`` cannot.  Anchored, because only a whole
#: name of this shape qualifies -- ``sup_1234_like`` is a real name.
#:
#: The optional letter is for the convention that gives a genome's RNA genes
#: their own series inside the same tag prefix: ``SS37A_r00010`` alongside
#: ``SS37A_41660``.  It belongs to the number, not the prefix -- dropping it
#: would make an rRNA indistinguishable from the protein-coding gene at the
#: same index -- so it is captured separately and put back on the label, while
#: the width test below still counts digits only.
LOCUS_TAG_PATTERN = re.compile(r"^(?=.*[A-Za-z])[A-Za-z0-9]+[_-]([A-Za-z]?)(\d+)$")

#: How many digits that number needs before the name is taken for a locus tag.
#:
#: Shape alone is not enough, because a *real* gene name also takes a numeric
#: suffix when an assembly has several copies of the gene -- ``ftsH_5`` is the
#: fifth ``ftsH``, and cutting it down to ``5`` throws away the only part that
#: says what it is.  What separates the two is how the number is written: a
#: locus tag counts positions across a whole replicon and is zero-padded to a
#: fixed width, four or five places (``SS37A_41660``, ``TST_0001``), while a
#: paralog suffix counts copies of one gene and so is a small integer.
#:
#: Set high enough that no plausible paralog count reaches it.  A locus tag
#: padded to three places falls the safe way: it stays long rather than losing
#: its prefix.
LOCUS_TAG_MIN_DIGITS = 4


def gene_label(name: str) -> str:
    """The part of a gene's name worth drawing.

    A GenBank feature carries a real name (``/gene``: ``recA``, ``dnaA``) or,
    far more often, only a locus tag (``/locus_tag``: ``SS37A_41660``).  The
    tag's prefix is the same on every feature of the assembly, so on a map it
    is thousands of repetitions of one string with the only distinguishing part
    -- the number -- pushed to the right where it collides with the neighbour::

        SS37A_41660   ->  41660     a locus tag: the number is the identity
        SS37A_r00010  ->  r00010    the RNA series of the same tag prefix
        ftsH_5        ->  ftsH_5    a real name, numbered; the name is the identity
        recA          ->  recA

    Anything not shaped like a tag, and anything whose number is shorter than
    :data:`LOCUS_TAG_MIN_DIGITS`, is returned untouched.  Display only --
    ``Annotation.name`` keeps what GenBank said, and the browser's detail panel
    shows it in full.
    """
    match = LOCUS_TAG_PATTERN.match(name)
    if match is None:
        return name
    series, digits = match.groups()
    if len(digits) < LOCUS_TAG_MIN_DIGITS:
        return name
    return series + digits


@dataclass(slots=True)
class LabelAssigner:
    """Assigns a short, unique display label to every read name."""

    _by_read: dict[str, str] = field(default_factory=dict)
    _taken: dict[str, str] = field(default_factory=dict)
    _fallback_counter: int = 0

    unmatched: int = 0
    """Reads whose name did not fit :data:`LEGACY_PATTERN`."""

    collisions: int = 0
    """Distinct reads that wanted a label already used by another read."""

    def label(self, read: str) -> str:
        """Return a stable, unique label for *read*."""
        cached = self._by_read.get(read)
        if cached is not None:
            return cached

        match = LEGACY_PATTERN.search(read)
        if match:
            base = f"R_{match.group(1)}"
        else:
            self.unmatched += 1
            self._fallback_counter += 1
            base = f"R#{self._fallback_counter}"

        label = base
        if self._taken.get(base, read) != read:
            self.collisions += 1
            suffix = 1
            while label in self._taken:
                suffix += 1
                label = f"{base}~{suffix}"

        self._taken[label] = read
        self._by_read[read] = label
        return label

    @property
    def total(self) -> int:
        return len(self._by_read)

    def warnings(self) -> list[str]:
        """Human-readable notes about label quality, for stderr and the PDF."""
        notes: list[str] = []
        if self.unmatched:
            pct = 100.0 * self.unmatched / max(self.total, 1)
            notes.append(
                f"{self.unmatched} of {self.total} read names ({pct:.1f}%) do not "
                f"match the legacy 'R_<digits>' pattern. The Perl version would "
                f"have merged each of these into the preceding read (finding B01)."
            )
        if self.collisions:
            notes.append(
                f"{self.collisions} read name(s) produced a duplicate short label "
                f"and were disambiguated with a '~n' suffix."
            )
        return notes
