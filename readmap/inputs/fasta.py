"""FASTA reference parser.

Only sequence *lengths* are kept.  The Perl held every base of the reference
in ``%seqs`` and then used nothing but ``length()`` on it (finding B26).
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import IO

from readmap.models import Contig, Metadata


def _open(path: str | Path) -> IO[str]:
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, "rt", encoding="utf-8", errors="replace")
    return open(p, "r", encoding="utf-8", errors="replace")


def accession_from_header(header: str) -> tuple[str, list[str]]:
    """Return ``(accession, aliases)`` for a FASTA header line.

    The Perl used ``/\\|(\\S+)/`` with a greedy ``\\S+``, which turns
    ``>gi|12345|ref|NC_000913.3|`` into the accession
    ``12345|ref|NC_000913.3|`` (finding B21).
    """
    token = header[1:].split(None, 1)[0] if len(header) > 1 else ""
    if not token:
        return "", []

    aliases: list[str] = [token]

    if "|" in token:
        parts = [p for p in token.split("|") if p]
        # NCBI pipe form: the last non-empty field is the accession.
        acc = parts[-1] if parts else token
        aliases.extend(parts)
    else:
        acc = token

    # An accession with a version suffix should also match the bare form and
    # vice versa -- the commonest cause of a silent BTOP/reference mismatch
    # (finding B11).
    if "." in acc:
        aliases.append(acc.rsplit(".", 1)[0])

    return acc, aliases


def parse_fasta(path: str | Path) -> tuple[list[Contig], Metadata, dict[str, str]]:
    """Parse a FASTA reference.

    Returns the contigs in file order, the (minimal) metadata, and a map from
    every recognised alias to the canonical contig name.
    """
    contigs: list[Contig] = []
    aliases: dict[str, str] = {}
    meta = Metadata(source_format="fa")

    name: str | None = None
    length = 0
    pending_aliases: list[str] = []

    def flush() -> None:
        nonlocal name, length, pending_aliases
        if name is None:
            return
        contigs.append(Contig(name=name, length=length))
        for alias in [name, *pending_aliases]:
            aliases.setdefault(alias, name)
        name, length, pending_aliases = None, 0, []

    with _open(path) as handle:
        for lineno, line in enumerate(handle):
            line = line.rstrip("\r\n")
            if lineno == 0:
                meta.first_line = line.strip()
            if line.startswith(">"):
                flush()
                name, pending_aliases = accession_from_header(line)
                if not name:
                    name = f"unnamed_{len(contigs) + 1}"
                length = 0
                continue
            if name is None:
                continue
            # Count every letter.  The Perl stripped anything outside
            # [acgtACGTnN], silently shortening any sequence containing IUPAC
            # ambiguity codes (finding B22).
            length += sum(1 for ch in line if ch.isalpha())

    flush()

    if not contigs:
        raise ValueError(f"no FASTA records found in {path}")

    return contigs, meta, aliases
