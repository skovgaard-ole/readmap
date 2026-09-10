# ReadRift

[![CI](https://github.com/skovgaard-ole/readrift/actions/workflows/ci.yml/badge.svg)](https://github.com/skovgaard-ole/readrift/actions/workflows/ci.yml)

Visualise long sequence reads (Oxford Nanopore, PacBio) mapped onto a reference
by BLAST, classified by how their alignment divides.

Reads whose pieces map far apart, in opposite orientation, or to a different
contig are evidence of structural variation, inversions, transposition, phage
excision or circularity. `readrift` draws them all on one wide map, plus a
summary figures section.

Python port of `legacy/read_print_23.pl`, which is kept in this repository as
the reference the port was audited against. See `code_structure.md` for the
structure of this package, and `CHANGES.md` for every way the output differs
from that original.

## Install

Requires **Python 3.10 or newer**. The only dependencies are matplotlib and
numpy, and `pip` fetches both. There is **no Ghostscript, no BLAST and no
compiler** in the picture — the PDF is written directly, and nothing shells out.

Get the code, then pick the line for your platform:

```bash
git clone https://github.com/skovgaard-ole/readrift.git
cd readrift
```

**Linux / macOS**

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e .
```

**Windows (PowerShell)**

```powershell
py -m venv .venv; .\.venv\Scripts\Activate.ps1
pip install -e .
```

**conda**, on any platform — if that is how you already work:

```bash
conda create -n readrift python=3.11 && conda activate readrift
pip install -e .
```

Any of them gives you a `readrift` command on your PATH. `python -m readrift`
works identically and needs no install beyond the dependencies, so both spellings
appear below.

To run the checks as well, install the development extras and use them:

```bash
pip install -e ".[dev]"
python -m pytest        # 44 tests, about ten seconds
python -m ruff check .
```

Those two commands are exactly what CI runs on Linux, macOS and Windows across
Python 3.10–3.13.

> If PowerShell refuses to run the activation script, Windows is blocking local
> scripts. `Set-ExecutionPolicy -Scope CurrentUser RemoteSigned` once, in a
> normal (non-admin) PowerShell, fixes it for good.

## Use

```bash
python -m readrift AP027148.gb DRR325755.btop -x 10
```

Produces `DRR325755.pdf`.

The BTOP file comes from BLAST with this exact `-outfmt`:

```bash
blastn -db blastdb -query all_reads.fa -out reads.btop -outfmt '6 delim=	 qseqid sframe qstart qend sstart send qlen sseqid btop'
```

The reference may be GenBank or FASTA; either may be gzipped. With GenBank,
gene annotations are drawn under the axis and coloured by their likelihood of
causing genomic instability.

### Sample data

**A clone contains no data**, only source. The alignment used in the examples
above is ~985 MB, an order of magnitude past what a git remote accepts, and the
reference is a public accession — neither belongs in this history. Both are
fetched rather than cloned:

| | |
|---|---|
| Reference | GenBank accession [`AP027148`](https://www.ncbi.nlm.nih.gov/nuccore/AP027148) → `AP027148.gb` |
| Reads | Sequence Read Archive run `DRR325755`, e.g. `fasterq-dump DRR325755` from the SRA Toolkit |

Then align the reads against the reference with the `blastn -outfmt` above to
produce the `.btop`. That alignment is the slow step, and it is BLAST's, not
`readrift`'s.

Point it at your own reference and reads instead and nothing changes — the two
file names in the examples carry no special meaning.

### Common runs

Include shorter reads, require longer matches, cap the coverage drawn:

```bash
python -m readrift ref.fa reads.btop -r 4000 -m 1000 -x 20
```

Add the alignment-identity figure and the batch summary table:

```bash
python -m readrift ref.gb reads.btop --identity -a
```

Pull the reads covering a region back out of the original reads file:

```bash
python -m readrift ref.fa reads.btop -e reads.fastq.gz,chr1,50000,100000
```

Maps only, no figures:

```bash
python -m readrift ref.fa reads.btop --no-plots
```

## Browse

The PDF is a poster — on a real dataset, several pages 149 inches wide. To look
around it interactively:

```bash
python -m readrift browse DRR325755.readriftdb.npz
```

That opens a genome browser in your web browser: pan and zoom the reference,
with a coverage track, the gene annotations, and every read coloured by its
class. Click a read for its real name and every place it aligned; press `n` to
jump to the next structural event; export the region on screen to PDF or PNG.

Every normal run writes that `.readriftdb.npz` cache alongside the PDF, so
browsing never re-reads the BTOP file. Reading and classifying a gigabyte-scale
BTOP takes minutes; opening the cache takes about a second. Use `--no-cache` to
skip writing it, or `--no-pdf` to write only the cache.

If no cache exists yet, `browse` builds one from the same arguments a normal run
takes, then serves it:

```bash
python -m readrift browse ref.gb reads.btop -x 10
```

The cache is rebuilt automatically when an input file or a filtering option
changes; `--rebuild` forces it.

The server binds `127.0.0.1` on a free port and serves a read-only view. It has
no authentication, so leave `--host` alone unless you mean it. `--port N` pins
the port, `--no-open` prints the URL instead of launching a browser window.

| Key | |
|---|---|
| drag / wheel | pan / zoom at the cursor |
| `←` `→` `+` `−` | pan, zoom |
| `f` | whole contig |
| `n` `p` | next / previous structural event |
| `/` | focus the search box |
| `esc` | clear the selection |

The search box takes a coordinate (`ctgA:120000-160000`), a contig name, a gene
name, or a read name.

Reads carry the same arrow vocabulary as the printed map:

| Symbol | |
|---|---|
| single chevron | the read runs on past here — pointing the way it runs |
| double chevron | it continues onto another contig (`Contig_Join`) |
| chevron + ring | it bridges the origin of a circular replicon (`Circle`) |
| chevron against the run | that piece is inverted |
| short bar | a junction inside a divided read |

Arrowheads are hidden when the lanes are too thin to read them; zoom in or
narrow the class filter to bring them back.

## Output

| File | When |
|---|---|
| `<prefix>.pdf` | always, unless `--no-pdf` |
| `<prefix>.readriftdb.npz` | always, unless `--no-cache` — the browser cache |
| `<prefix>_Analysis.tsv` | with `-a` |
| `extract-reads-list_<contig>_<start>_<end>.txt` | with `-e` — one line per read |
| `extract-reads-list_<contig>_<start>_<end>.fastq` / `.fasta` | with `-e` — the sequences |

The PDF is: front page (parameters, legend, project information, statistics,
provenance) → figures → the read maps, one contig at a time.

## Options

`python -m readrift --help` lists everything. Every short flag from the Perl
version keeps its letter and meaning.

Three options were split apart because the Perl overloaded one:

| Option | Controls |
|---|---|
| `-m --min-match-length` | minimum aligned span of a **drawn** segment |
| `--division-threshold` | the short/long divided cut-off (defaults to `-m`) |
| `--index-window` | granularity of the `-e` extraction index |

and `--min-hsp-length` is now the minimum query span of a single BLAST HSP,
which is what `-y --micro-match-length` claimed to be but was not.
`--micro-match-length` still exists and still sets the junction tolerance
(`value / 4`), which is all it ever did.

### `--keep-fold-back`

Some long reads are sequenced through their template and then straight back
along it — adapter or end ligation joining the molecule to its own reverse
complement. They map as one long alignment followed by a run of alignments
re-reading the *same* stretch of reference in the opposite direction, and left
alone they masquerade as a long-distance junction with an inversion at every
fold, while counting their bases towards coverage twice.

These are dropped by default, and the run reports how many reads were affected.
Pass `--keep-fold-back` to keep them — for comparing against an older run, or
if your library really can produce genuine hairpins. See `CHANGES.md` §1.8 for
the exact rule and what it is careful *not* to remove.

## Layout

```
readrift/
├── cli.py          argparse, built from the same table the front page reads
├── params.py       option table + resolved Params
├── models.py       Hit, ReadGroup, Contig, Annotation, Segment, ...
├── labels.py       short display names for reads
├── inputs/         btop, fasta, genbank, reads (FASTA/FASTQ/gz)
├── classify.py     HSP selection and read classification  <- the science
├── layout.py       lane packing and page splitting
├── stats.py        counts, coverage, plot inputs
├── btop_trace.py   BTOP alignment-string parser
├── extract.py      -e region extraction
├── report.py       the analysis TSV
├── render/         theme, mapfig, frontpage, plots, PDF assembly
├── browser/        the interactive view
│   ├── store.py      the .readriftdb.npz cache: build and query by region
│   ├── region.py     one window, packed into lanes  <- shared by API and export
│   ├── server.py     stdlib HTTP server, loopback only
│   ├── export.py     the region on screen -> PDF/PNG
│   └── static/       index.html, app.js, style.css  (canvas front end)
└── pipeline.py     orchestration
```

The classification algorithm and the palette have exactly one implementation
each. `browser/region.py` packs lanes with `layout.LanePacker`, the same class
the printed map uses, and the browser's colours are served from
`render/theme.py` rather than restated in CSS — so a region on screen matches
the corresponding slice of the poster.

The classification algorithm is documented step by step in
`code_structure.md` §5 and implemented in `classify.py`. It is unchanged from
the Perl — only its implementation is different.

## Licence

MIT — see [`LICENSE`](LICENSE). © 2026 Ole Skovgaard.
