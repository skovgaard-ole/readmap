"""Reference-sequence loading: format detection and dispatch."""

from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from pathlib import Path

from readrift.inputs.fasta import parse_fasta
from readrift.inputs.genbank import parse_genbank
from readrift.models import Annotation, Contig, Metadata


@dataclass(slots=True)
class Reference:
    """Everything the pipeline needs to know about the reference."""

    contigs: list[Contig]
    annotations: list[Annotation] = field(default_factory=list)
    metadata: Metadata = field(default_factory=Metadata)
    aliases: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        # Longest first -- the Perl's @seqs order, which determines page order.
        self.contigs.sort(key=lambda c: (-c.length, c.name))

    @property
    def total_length(self) -> int:
        return sum(c.length for c in self.contigs)

    @property
    def names(self) -> list[str]:
        return [c.name for c in self.contigs]

    @property
    def is_genbank(self) -> bool:
        return self.metadata.source_format == "gb"

    def by_name(self) -> dict[str, Contig]:
        return {c.name: c for c in self.contigs}

    def resolve(self, name: str) -> str | None:
        """Map a BTOP ``sseqid`` onto a canonical contig name, or ``None``."""
        canonical = self.aliases.get(name)
        if canonical is not None:
            return canonical
        if "." in name:
            return self.aliases.get(name.rsplit(".", 1)[0])
        return None


def sniff_format(path: str | Path) -> str:
    """Return ``"gb"`` or ``"fa"``.

    The Perl left ``$ref_seq_type`` undefined when a file matched neither and
    then died on ``$ref_seq_type eq "gb"`` with an uninitialized-value error
    rather than telling the user what was wrong (finding B12).
    """
    p = Path(path)
    opener = gzip.open if p.suffix == ".gz" else open
    with opener(p, "rt", encoding="utf-8", errors="replace") as handle:  # type: ignore[operator]
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                return "fa"
            if stripped.upper().startswith("LOCUS"):
                return "gb"
            raise SystemExit(
                f"ERROR: cannot tell what format {path} is.\n"
                f"       Expected a FASTA file (first line starts with '>') or a\n"
                f"       GenBank file (first line starts with 'LOCUS').\n"
                f"       First line was: {stripped[:120]!r}"
            )
    raise SystemExit(f"ERROR: reference file is empty: {path}")


def load_reference(path: str | Path) -> Reference:
    """Load a FASTA or GenBank reference."""
    fmt = sniff_format(path)

    if fmt == "gb":
        contigs, annotations, meta, aliases = parse_genbank(path)
        return Reference(
            contigs=contigs, annotations=annotations, metadata=meta, aliases=aliases
        )

    contigs, meta, aliases = parse_fasta(path)
    return Reference(contigs=contigs, metadata=meta, aliases=aliases)
