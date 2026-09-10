"""Streaming reader for BLAST BTOP tabular output.

Expected ``-outfmt``::

    6 delim=<TAB> qseqid sframe qstart qend sstart send qlen sseqid btop

The final ``btop`` column is optional -- the tool's own usage text says so, and
the Perl tolerated it by accident (it only ever indexed columns 0..8 and never
read column 8).  Here it is optional on purpose.

Files may be gzipped; the extension decides.
"""

from __future__ import annotations

import gzip
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO

from readrift.models import Hit

#: Number of columns before the optional trace string.
_MIN_COLUMNS = 8

#: Stop tracking read names for the out-of-order check past this many, to keep
#: the check from costing more memory than the data.
_ORDER_CHECK_CAP = 2_000_000


@dataclass(slots=True)
class BtopStats:
    """What the reader saw, for the run log and the PDF stamp."""

    lines: int = 0
    hits: int = 0
    reads: int = 0
    malformed: int = 0
    short_reads: int = 0
    """Lines dropped because ``qlen < --min-read-length``."""

    contigs: set[str] = field(default_factory=set)
    malformed_examples: list[str] = field(default_factory=list)
    out_of_order: int = 0
    """Reads whose lines were not contiguous in the file."""

    order_check_complete: bool = True

    def note_malformed(self, lineno: int, why: str) -> None:
        self.malformed += 1
        if len(self.malformed_examples) < 5:
            self.malformed_examples.append(f"line {lineno}: {why}")

    def warnings(self) -> list[str]:
        notes: list[str] = []
        if self.malformed:
            detail = "; ".join(self.malformed_examples)
            notes.append(f"{self.malformed} unparsable BTOP line(s). First: {detail}")
        if self.out_of_order:
            notes.append(
                f"{self.out_of_order} read(s) appear in more than one block. BLAST "
                f"groups its output by query; a re-sorted file will be classified "
                f"per block rather than per read."
                + ("" if self.order_check_complete else " (check truncated)")
            )
        return notes


def _open(path: str | Path) -> IO[str]:
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, "rt", encoding="utf-8", errors="replace")
    return open(p, "r", encoding="utf-8", errors="replace")


def parse_line(line: str) -> Hit:
    """Parse one BTOP record.  Raises :class:`ValueError` on bad input."""
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) < _MIN_COLUMNS:
        raise ValueError(f"expected >= {_MIN_COLUMNS} tab-separated columns, got {len(fields)}")

    try:
        strand = int(fields[1])
        qstart = int(fields[2])
        qend = int(fields[3])
        sstart = int(fields[4])
        send = int(fields[5])
        qlen = int(fields[6])
    except ValueError as exc:
        raise ValueError(f"non-numeric coordinate ({exc})") from exc

    if strand not in (1, -1):
        # BLAST writes sframe as +1/-1 for blastn.  Anything else means the
        # -outfmt string did not match what this tool expects.
        raise ValueError(f"sframe must be 1 or -1, got {strand}")

    return Hit(
        read=fields[0],
        strand=strand,
        qstart=qstart,
        qend=qend,
        sstart=sstart,
        send=send,
        qlen=qlen,
        contig=fields[7],
        btop=fields[8] if len(fields) > _MIN_COLUMNS else "",
    )


def iter_read_groups(
    path: str | Path,
    min_read_length: int = 0,
    stats: BtopStats | None = None,
) -> Iterator[list[Hit]]:
    """Yield the HSPs of each read, in file order.

    Lines are grouped on consecutive equal ``qseqid``, which is how BLAST
    writes them.  Unlike the Perl, the grouping key is the *real* read name --
    see :mod:`readrift.labels` for why that matters.
    """
    st = stats if stats is not None else BtopStats()

    seen: set[int] | None = set()
    current: list[Hit] = []
    current_name: str | None = None

    with _open(path) as handle:
        for lineno, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            st.lines += 1

            try:
                hit = parse_line(line)
            except ValueError as exc:
                st.note_malformed(lineno, str(exc))
                continue

            if hit.qlen < min_read_length:
                st.short_reads += 1
                continue

            st.hits += 1
            st.contigs.add(hit.contig)

            if current_name is not None and hit.read != current_name:
                st.reads += 1
                yield current
                current = []

            if hit.read != current_name:
                current_name = hit.read
                if seen is not None:
                    key = hash(hit.read)
                    if key in seen:
                        st.out_of_order += 1
                    else:
                        seen.add(key)
                        if len(seen) >= _ORDER_CHECK_CAP:
                            st.order_check_complete = False
                            seen = None

            current.append(hit)

    if current:
        st.reads += 1
        yield current
