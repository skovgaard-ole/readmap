"""The batch-analysis TSV (``-a``).

The Perl emitted a 14-column summary row followed by per-contig rows that put
the accession in column 16, with no header line despite the comment in the
source describing one -- so nothing downstream could parse it without
guesswork (finding B19).  It also wrote the file on every run, because the
default value of the flag was an empty array reference, which is true in Perl
(finding B02).

Here every row has the same shape, a ``record_type`` column says what it is,
and there is a header.  The file loads directly with
``pandas.read_csv(path, sep="\\t")``.
"""

from __future__ import annotations

from pathlib import Path

from readrift.inputs.reference import Reference
from readrift.models import ReadClass
from readrift.params import Params
from readrift.stats import Stats

COLUMNS: tuple[str, ...] = (
    "record_type",
    "sample",
    "bio_project",
    "bio_sample",
    "sra",
    "assembly_method",
    "sequencing_technology",
    "organism",
    "contig",
    "length_bp",
    "mean_coverage",
    "reads_total",
    "undivided",
    "short_divided",
    "long_divided",
    "inverted",
    "events",
    # The thresholds that produced the numbers above, so a batch of runs
    # aggregates without having to remember how each one was invoked.
    "min_read_length",
    "min_match_length",
    "division_threshold",
    "cov_max",
)


def write_report(
    path: str | Path,
    stats: Stats,
    reference: Reference,
    params: Params,
    sample: str,
) -> Path:
    """Write the analysis TSV and return its path."""
    meta = reference.metadata
    shared = [
        sample,
        meta.bio_project,
        meta.bio_sample,
        meta.sra,
        meta.assembly_method,
        meta.sequencing_technology,
        meta.organism,
    ]
    settings = [
        str(params.min_read_length),
        str(params.min_match_length),
        str(params.division_cut),
        str(params.cov_max),
    ]

    rows: list[list[str]] = [
        [
            "run",
            *shared,
            "*",
            str(stats.reference_length),
            f"{stats.mean_coverage:.4f}",
            str(stats.total_reads),
            str(stats.counts[ReadClass.UNDIVIDED]),
            str(stats.counts[ReadClass.SHORT_DIVIDED]),
            str(stats.counts[ReadClass.LONG_DIVIDED]),
            str(stats.inverted_reads),
            str(len(stats.events())),
            *settings,
        ]
    ]

    for contig in reference.contigs:
        cs = stats.per_contig.get(contig.name)
        if cs is None:
            continue
        rows.append(
            [
                "contig",
                *shared,
                contig.name,
                str(contig.length),
                f"{cs.mean_coverage:.4f}",
                str(cs.reads),
                str(cs.counts[ReadClass.UNDIVIDED]),
                str(cs.counts[ReadClass.SHORT_DIVIDED]),
                str(cs.counts[ReadClass.LONG_DIVIDED]),
                str(cs.inverted),
                str(len(cs.events)),
                *settings,
            ]
        )

    target = Path(path)
    with open(target, "w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            handle.write("\t".join(_clean(cell) for cell in row) + "\n")

    return target


def _clean(value: str) -> str:
    """Keep the file parseable: no embedded tabs or newlines."""
    return value.replace("\t", " ").replace("\r", " ").replace("\n", " ").strip()
