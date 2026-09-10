"""Reading FASTA / FASTQ read files, optionally gzipped.

The Perl sniffed the format by reading one line and then calling
``seek($fh, 0, 0)``.  For a ``.gz`` input the handle was a pipe from
``gzip -dc``, where ``seek`` fails; the return value was never checked, so
parsing resumed at line 2 -- the first record was lost and, for FASTQ, the
four-line phase was broken for the rest of the file (finding B07).

It also built that pipe with ``open($fh, "-|", "gzip -dc $file")``, which
interpolates the path into a shell command unquoted -- broken for any path
containing a space, and not available at all on stock Windows, even though a
pure-Perl gunzip had been probed for at startup and then never used
(finding B08).

Here decompression is done in-process and the first line is kept rather than
seeking back to it.
"""

from __future__ import annotations

import gzip
from collections.abc import Iterator
from pathlib import Path
from typing import IO


def open_reads(path: str | Path) -> IO[str]:
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, "rt", encoding="utf-8", errors="replace")
    return open(p, "r", encoding="utf-8", errors="replace")


def detect_format(first_line: str) -> str:
    if first_line.startswith("@"):
        return "fastq"
    if first_line.startswith(">"):
        return "fasta"
    raise ValueError(
        f"cannot tell whether the reads file is FASTA or FASTQ; "
        f"first line begins {first_line[:40]!r}"
    )


def iter_records(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield ``(read_id, verbatim_record_text)`` for every record."""
    with open_reads(path) as handle:
        first = ""
        for line in handle:
            if line.strip():
                first = line
                break
        if not first:
            return

        fmt = detect_format(first)
        if fmt == "fastq":
            yield from _iter_fastq(handle, first)
        else:
            yield from _iter_fasta(handle, first)


def _iter_fasta(handle: IO[str], first: str) -> Iterator[tuple[str, str]]:
    header = first.rstrip("\n")
    body: list[str] = []

    for line in handle:
        if line.startswith(">"):
            yield _fasta_record(header, body)
            header = line.rstrip("\n")
            body = []
        else:
            body.append(line.rstrip("\n"))

    yield _fasta_record(header, body)


def _fasta_record(header: str, body: list[str]) -> tuple[str, str]:
    read_id = header[1:].split(None, 1)[0] if len(header) > 1 else ""
    text = header + "\n" + "\n".join(body) + "\n"
    return read_id, text


def _iter_fastq(handle: IO[str], first: str) -> Iterator[tuple[str, str]]:
    """Length-matched FASTQ reader.

    Quality characters may include ``@`` and ``+``, so records are delimited by
    counting: read sequence lines until a line starting with ``+``, then read
    quality lines until as many characters have been seen as there were bases.
    """
    header = first

    while header is not None:
        read_id = header[1:].split(None, 1)[0] if len(header) > 1 else ""

        sequence: list[str] = []
        plus: str | None = None
        for line in handle:
            if line.startswith("+"):
                plus = line
                break
            sequence.append(line)

        if plus is None:
            return

        wanted = sum(len(line.strip()) for line in sequence)
        quality: list[str] = []
        seen = 0
        for line in handle:
            quality.append(line)
            seen += len(line.strip())
            if seen >= wanted:
                break

        yield read_id, "".join([header, *sequence, plus, *quality])

        header = None
        for line in handle:
            if line.startswith("@"):
                header = line
                break
