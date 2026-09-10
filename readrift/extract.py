"""Pulling region-specific reads back out of the original reads file (``-e``).

Two fixes over the Perl:

* The ``extract-reads-list_*.txt`` file gets the list of reads it advertises.
  The Perl opened it, wrote a two-line header, and never wrote the list
  (finding B06).
* The window index is built only when ``-e`` is actually given.  The Perl built
  it on every run -- one hash entry per kilobase of every match -- whether or
  not extraction was requested, which on a real dataset was the single largest
  memory consumer in the program and pure waste in the default case
  (finding B25).
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

from readrift.inputs.reads import iter_records
from readrift.models import Hit
from readrift.params import ExtractSpec, Params


@dataclass(frozen=True, slots=True)
class Placement:
    """Where one read matched, as recorded in the index."""

    read: str
    start: int
    end: int


class ReadIndex:
    """Window index from reference position to the reads that align there."""

    __slots__ = ("_windows", "_window_size", "_min_span")

    def __init__(self, window_size: int, min_span: int) -> None:
        self._windows: dict[str, dict[int, list[Placement]]] = {}
        self._window_size = max(1, window_size)
        self._min_span = min_span

    def add(self, hit: Hit) -> None:
        if hit.sspan < self._min_span:
            return
        placement = Placement(read=hit.read, start=hit.sstart, end=hit.send)
        contig = self._windows.setdefault(hit.contig, {})
        first = hit.slow // self._window_size
        last = hit.shigh // self._window_size
        for window in range(first, last + 1):
            contig.setdefault(window, []).append(placement)

    def query(self, contig: str, start: int, end: int) -> list[Placement]:
        windows = self._windows.get(contig)
        if not windows:
            return []
        seen: set[str] = set()
        hits: list[Placement] = []
        for window in range(start // self._window_size, end // self._window_size + 1):
            for placement in windows.get(window, ()):
                if placement.read in seen:
                    continue
                low = min(placement.start, placement.end)
                high = max(placement.start, placement.end)
                if high >= start and low <= end:
                    seen.add(placement.read)
                    hits.append(placement)
        hits.sort(key=lambda p: min(p.start, p.end))
        return hits

    @property
    def contigs(self) -> list[str]:
        return sorted(self._windows)


def index_groups(
    groups: Iterable[list[Hit]], index: ReadIndex
) -> Iterator[list[Hit]]:
    """Pass groups through, indexing each HSP on the way."""
    for hits in groups:
        for hit in hits:
            index.add(hit)
        yield hits


def run_extraction(
    spec: ExtractSpec,
    index: ReadIndex,
    params: Params,
    out_dir: Path | None = None,
) -> tuple[int, int, Path]:
    """Write the read list and the matching sequences.

    Output lands beside the PDF, next to ``--out``, rather than in whatever
    happens to be the working directory.

    Returns ``(reads_requested, reads_written, sequence_path)``.
    """
    directory = out_dir if out_dir is not None else Path(params.out).parent
    directory.mkdir(parents=True, exist_ok=True)
    placements = index.query(spec.contig, spec.start, spec.end)
    wanted = {p.read: p for p in placements}

    list_path = directory / f"{spec.stem}.txt"
    with open(list_path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"# region\t{spec.contig}\t{spec.start}\t{spec.end}\n")
        handle.write(f"# reads\t{len(placements)}\n")
        handle.write("read\tcontig\tstart\tend\n")
        for placement in placements:
            handle.write(
                f"{placement.read}\t{spec.contig}\t{placement.start}\t{placement.end}\n"
            )

    if not Path(spec.reads_file).exists():
        raise SystemExit(
            f"ERROR: reads file for --extract-reads not found: {spec.reads_file}"
        )

    remaining = dict(wanted)
    written = 0
    suffix = ".fastq"
    seen_any = False
    records: list[tuple[str, str]] = []

    for read_id, text in iter_records(spec.reads_file):
        if not seen_any:
            seen_any = True
            suffix = ".fastq" if text.startswith("@") else ".fasta"
        if read_id in remaining:
            records.append((read_id, text))
            del remaining[read_id]
            written += 1
            if not remaining:
                break

    sequence_path = directory / f"{spec.stem}{suffix}"
    with open(sequence_path, "w", encoding="utf-8", newline="\n") as handle:
        for _read_id, text in records:
            handle.write(text)

    return len(placements), written, sequence_path
