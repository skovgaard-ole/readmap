"""Synthetic fixtures.

Everything is generated rather than checked in, so the expected classification
of every read can be derived from ``code_structure.md`` §5 by hand and stated
in one place -- :data:`EXPECTED` below.

The reads deliberately include the cases the Perl got wrong:

* ``0a1b2c3d-...``  an ONT-style UUID name, which the legacy ``R_<digits>``
  pattern never matches (finding B01);
* ``read(x).9``     a name containing a parenthesis, which broke the
  PostScript (finding B09);
* ``read.11``       mapped to an accession the reference does not have
  (finding B11);
* ``read.12``       three HSPs where two compete for the same stretch of the
  read, exercising the repeat-collapse and overlap rules.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from readrift.models import Junction, ReadClass
from readrift.params import Params

CONTIG_LENGTHS = {"ctgA": 40_000, "ctgB": 12_000}

# qseqid, sframe, qstart, qend, sstart, send, qlen, sseqid, btop
BTOP_ROWS: list[tuple] = [
    # undivided, forward
    ("read.1", 1, 10, 5010, 1000, 6000, 5200, "ctgA", "5000"),
    # undivided, reverse (sstart > send, as BLAST writes it)
    ("read.2", -1, 20, 4020, 16000, 12000, 4200, "ctgA", "4000"),
    # short-divided: a ~300 bp deletion
    ("read.3", 1, 10, 2010, 20000, 22000, 4200, "ctgA", "2000"),
    ("read.3", 1, 2100, 4100, 22300, 24300, 4200, "ctgA", "2000"),
    # short-divided with an inverted second piece
    ("read.4", 1, 10, 1510, 26000, 27500, 3200, "ctgA", "1500"),
    ("read.4", -1, 1600, 3100, 29200, 27700, 3200, "ctgA", "1500"),
    # long-divided: the two halves sit 26 kb apart
    ("read.5", 1, 10, 3010, 1000, 4000, 6200, "ctgA", "3000"),
    ("read.5", 1, 3100, 6100, 30000, 33000, 6200, "ctgA", "3000"),
    # contig join
    ("read.6", 1, 10, 2010, 36000, 38000, 4200, "ctgA", "2000"),
    ("read.6", 1, 2100, 4100, 1000, 3000, 4200, "ctgB", "2000"),
    # circular join: spans the whole of ctgB
    ("read.7", 1, 10, 2010, 10000, 12000, 4200, "ctgB", "2000"),
    ("read.7", 1, 2100, 4100, 10, 2010, 4200, "ctgB", "2000"),
    # below --min-read-length, dropped
    ("read.8", 1, 10, 500, 5000, 5490, 600, "ctgA", "490"),
    # ONT-style UUID name (finding B01)
    (
        "0a1b2c3d-e4f5-6789-abcd-ef0123456789",
        1, 10, 3010, 8000, 11000, 3200, "ctgA", "3000",
    ),
    # a parenthesis in the name (finding B09)
    ("read(x).9", 1, 10, 2010, 14000, 16000, 2200, "ctgA", "2000"),
    # accession the reference does not have (finding B11)
    ("read.11", 1, 10, 2010, 1000, 3000, 2200, "ctgZZZ", "2000"),
    # three HSPs, the last competing with the second for the same read bases
    ("read.12", 1, 10, 1010, 1000, 2000, 3200, "ctgA", "1000"),
    ("read.12", 1, 1100, 2100, 20000, 21000, 3200, "ctgA", "1000"),
    ("read.12", 1, 1100, 2100, 30000, 31000, 3200, "ctgA", "1000"),
]

#: read -> (class, junction, inverted).  Derived by hand from the algorithm.
EXPECTED: dict[str, tuple[ReadClass, Junction, bool]] = {
    "read.1": (ReadClass.UNDIVIDED, Junction.PLAIN, False),
    "read.2": (ReadClass.UNDIVIDED, Junction.PLAIN, False),
    "read.3": (ReadClass.SHORT_DIVIDED, Junction.PLAIN, False),
    "read.4": (ReadClass.SHORT_DIVIDED, Junction.PLAIN, True),
    "read.5": (ReadClass.LONG_DIVIDED, Junction.PLAIN, False),
    "read.6": (ReadClass.LONG_DIVIDED, Junction.CONTIG_JOIN, False),
    "read.7": (ReadClass.LONG_DIVIDED, Junction.CIRCLE, False),
    "0a1b2c3d-e4f5-6789-abcd-ef0123456789": (
        ReadClass.UNDIVIDED, Junction.PLAIN, False,
    ),
    "read(x).9": (ReadClass.UNDIVIDED, Junction.PLAIN, False),
    "read.12": (ReadClass.LONG_DIVIDED, Junction.PLAIN, False),
}

#: read.8 is below the length cut, read.11 is on an unknown accession.
DROPPED = {"read.8", "read.11"}

GENBANK = """\
LOCUS       ctgA                   40000 bp    DNA     circular BCT 01-JAN-2026
DEFINITION  Synthetic test contig A.
ACCESSION   ctgA
VERSION     ctgA.1
DBLINK      BioProject: PRJNA000001
            BioSample: SAMN00000001
            Sequence Read Archive: SRR0000001
SOURCE      Testus exemplaris
  ORGANISM  Testus exemplaris
            Bacteria; Proteobacteria.
COMMENT     ##Genome-Assembly-Data-START##
            Assembly Method       :: Flye v. 2.9
            Sequencing Technology :: Oxford Nanopore (MinION)
            ##Genome-Assembly-Data-END##
FEATURES             Location/Qualifiers
     source          1..40000
                     /organism="Testus exemplaris"
     gene            1000..2500
                     /gene="thrA"
                     /locus_tag="TST_0001"
     CDS             1000..2500
                     /gene="thrA"
                     /product="aspartokinase"
     gene            complement(3000..4200)
                     /locus_tag="TST_0002"
     CDS             complement(3000..4200)
                     /product="hypothetical protein"
     gene            5000..5075
                     /locus_tag="TST_0003"
     tRNA            5000..5075
                     /product="tRNA-Ala"
     gene            8000..9500
                     /locus_tag="TST_0004"
     CDS             8000..9500
                     /product="IS3 family transposase with a
                     wrapped product name"
     gene            12000..13500
                     /locus_tag="TST_0005"
                     /pseudo
     gene            20000..22000
                     /locus_tag="TST_0006"
     CDS             20000..22000
                     /product="phage tail protein"
     repeat_region   25000..25500
                     /rpt_type=dispersed
//
LOCUS       ctgB                   12000 bp    DNA     circular BCT 01-JAN-2026
DEFINITION  Synthetic test contig B.
ACCESSION   ctgB
VERSION     ctgB.1
FEATURES             Location/Qualifiers
     source          1..12000
     gene            500..1400
                     /gene="repB"
     CDS             500..1400
                     /product="replication protein"
//
"""


def _sequence(length: int) -> str:
    return ("ACGTTGCANRYACGT" * (length // 15 + 1))[:length]


def _wrap(text: str, width: int = 60) -> str:
    return "\n".join(text[i : i + width] for i in range(0, len(text), width))


@pytest.fixture
def fixture_dir(tmp_path: Path) -> Path:
    """A directory holding ref.fa, ref.gb, reads.btop, reads.fastq(.gz)."""
    fasta = "\n".join(
        f">{name} synthetic contig\n{_wrap(_sequence(length))}"
        for name, length in CONTIG_LENGTHS.items()
    )
    (tmp_path / "ref.fa").write_text(fasta + "\n", encoding="utf-8")
    (tmp_path / "ref.gb").write_text(GENBANK, encoding="utf-8")

    btop = "\n".join("\t".join(str(cell) for cell in row) for row in BTOP_ROWS)
    (tmp_path / "reads.btop").write_text(btop + "\n", encoding="utf-8")

    names = sorted({row[0] for row in BTOP_ROWS})
    fastq_parts = []
    for index, name in enumerate(names):
        seq = _sequence(200 + index)
        fastq_parts.append(f"@{name} runid=test\n{seq}\n+\n{'I' * len(seq)}\n")
    fastq = "".join(fastq_parts)
    (tmp_path / "reads.fastq").write_text(fastq, encoding="utf-8")
    with gzip.open(tmp_path / "reads.fastq.gz", "wt", encoding="utf-8") as handle:
        handle.write(fastq)

    return tmp_path


@pytest.fixture
def params(fixture_dir: Path) -> Params:
    """Parameters tuned to the fixture's coordinate scale."""
    return Params(
        seq_file=str(fixture_dir / "ref.fa"),
        btop_file=str(fixture_dir / "reads.btop"),
        out=str(fixture_dir / "out"),
        min_read_length=1_000,
        min_match_length=500,
        min_hsp_length=100,
        micro_match_length=100,
        map_scale=10.0,
        space=20,
    )
