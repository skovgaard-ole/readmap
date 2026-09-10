"""Parsing the BLAST BTOP alignment trace.

The Perl read this column out of every line and then never looked at it
(finding B34) -- the tool is named for BTOP but derived no alignment
statistics from it at all.

BTOP grammar is a flat sequence of two kinds of token:

* a run of digits -- that many consecutive identical positions;
* a pair of characters -- one alignment column that is *not* a match:
  ``AG`` a mismatch, ``-A`` a gap in the query, ``A-`` a gap in the subject.

Example: ``120AG45-T7`` is 120 matches, an A/G mismatch, 45 matches, a
one-base query gap, 7 matches.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

_TOKEN_RE = re.compile(r"(\d+)|([A-Za-z*-]{2})")


@dataclass(frozen=True, slots=True)
class Trace:
    """Counts derived from one BTOP string."""

    matches: int = 0
    mismatches: int = 0
    query_gaps: int = 0
    subject_gaps: int = 0

    @property
    def aligned_columns(self) -> int:
        return self.matches + self.mismatches + self.query_gaps + self.subject_gaps

    @property
    def identity(self) -> float:
        """Fraction of alignment columns that are identical, 0.0 - 1.0."""
        total = self.aligned_columns
        return self.matches / total if total else 0.0

    @property
    def mismatch_rate(self) -> float:
        total = self.aligned_columns
        return (self.mismatches + self.query_gaps + self.subject_gaps) / total if total else 0.0


EMPTY = Trace()


def parse(btop: str) -> Trace:
    """Count matches, mismatches and gaps in a BTOP string.

    Unrecognised text is ignored rather than raising -- a malformed trace
    should not stop a run that is otherwise fine.
    """
    if not btop:
        return EMPTY

    matches = mismatches = query_gaps = subject_gaps = 0

    for run, pair in _TOKEN_RE.findall(btop):
        if run:
            matches += int(run)
            continue
        first, second = pair[0], pair[1]
        if first == "-":
            query_gaps += 1
        elif second == "-":
            subject_gaps += 1
        else:
            mismatches += 1

    return Trace(
        matches=matches,
        mismatches=mismatches,
        query_gaps=query_gaps,
        subject_gaps=subject_gaps,
    )
