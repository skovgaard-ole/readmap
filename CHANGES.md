# Changes from `read_print_23.pl`

Every way the Python port's output differs from the Perl's, and why. Finding
IDs (`B01` …) refer to the audit in `port_plan.md` §4.

The Perl original is in this repository at **`legacy/read_print_23.pl`**, so the
line citations below (`read_print_23.pl:1141`) can be followed.

If you have Perl-generated results on file, **read §1 first** — those are the
numbers that move.

---

## 1. Numbers that change

### 1.1 Read counts — `B01`

The Perl rewrote every read name with

```perl
$btop_read[0] =~ /\w+[\._](\d+)/;
$btop_read[0] = "R_" . $1;
```

and never checked whether the match succeeded. On failure `$1` still held the
**previous** read's capture, so the read silently inherited the previous read's
identity and its HSPs were merged into that read's group. Reads were grouped by
"did this rewritten name change?", so a run of non-matching names collapsed into
a single read.

Names like `DRR325755.12345` match and were fine. Oxford Nanopore UUIDs
(`0a1b2c3d-e4f5-…`) never match, so on ONT data the whole file collapsed.

**Now:** read names are never rewritten. Grouping uses the real `qseqid`.
`R_<n>` survives only as a display label on the map, and colliding labels get a
`~2` suffix rather than silently merging two reads.

**Effect:** total read counts, and the split between classes, change on any
dataset with non-matching names. The run banner and the PDF front page report
**how many names failed the legacy pattern** — if that number is zero, this
change did not affect you.

### 1.2 Inversion count — `B17`

Three things were wrong. The Perl counted inversions where it *drew* them, so a
segment straddling a page boundary was neither marked red nor counted; it
counted drawn **segments**, not reads, while labelling the result "reads
indicating inversions"; and it only looked at short-divided reads.

**Now:** the inversion flag is computed during classification. It counts
**reads** with at least one segment opposite to their first, across all divided
classes, regardless of where the page folds fall.

**Effect:** the number changes in both directions — down where one read had
several inverted segments, up where inversions crossed a page or sat in a
long-divided read.

### 1.3 Coverage — `B10`, decision D3

The Perl summed matched span per HSP for undivided and long-divided reads, but
the whole first-to-last span *including the skipped reference* for
short-divided reads. The reported coverage therefore depended on how rearranged
the sample was.

**Now:** coverage is **matched bases only**, uniformly — the sum of aligned
reference spans. Reference that a divided read jumps over is not covered by it.

**Effect:** coverage is lower than before, by an amount proportional to the
short-divided fraction.

A second, smaller change in the same direction: long-divided segments shorter
than `--min-match-length` are still hidden from the map (as before) but now
**count towards coverage**, because they are matched bases.

### 1.4 `Circle` / `Contig_Join` decoration on short-divided reads — `B04`

In the short-divided drawing block, `@map` was left over from the preceding
grouping loop and held the **last record of the whole contig**. `$map[8]` was
then used to pick the arrow for every read on that contig, so the decoration
was applied all-or-nothing per contig.

**Now:** the junction tag belongs to the read and is decided during
classification.

**Effect:** circular-join and contig-join arrows appear on the reads that
actually earn them. On a contig whose last short-divided read happened to be a
`Circle`, this removes a lot of spurious ring markers.

### 1.5 Reads on accessions missing from the reference — `B11`

The Perl collected the BTOP accessions, wrote them to the debug log with the
note "They must match the Accession numbers from the GB/Fasta file", and never
compared them. A mismatch produced blank pages and a zero-read report with no
error.

**Now:** accessions are resolved through an alias table that also knows each
record's `ACCESSION` and `VERSION` (so `NC_000913` matches `NC_000913.3`).
Reads on genuinely unknown accessions are counted and reported, and a *total*
mismatch is a hard error with both lists printed, rather than an empty PDF.

**Effect:** runs that silently produced nothing now say why. Runs where the
BLAST database used versioned accessions and the GenBank `LOCUS` line did not
now work instead of producing nothing.

### 1.6 Contig-joining reads are drawn on both contigs

Related to 1.4. A read whose pieces map to two contigs is now placed on each,
showing the piece that belongs there. Previously the drawing code took the
whole read's coordinate range from its first HSP's contig.

### 1.7 Sequence lengths containing IUPAC ambiguity codes — `B22`

`s/[^acgtACGTnN]//g` deleted R, Y, S, W, K, M, B, D, H and V from the
reference before measuring it, so a contig containing them was recorded shorter
than it is. Lengths now count every letter.

**Effect:** contig lengths, page counts and the circularity test shift slightly
on references with ambiguity codes. Most bacterial assemblies have none.

### 1.8 Fold-back reads are dropped — new, not a Perl finding

A long read is sometimes sequenced through its template and then straight back
along it: adapter or end ligation joins the molecule to its own reverse
complement. It maps as one long HSP followed by a run of HSPs that re-read the
*same stretch of reference* in the opposite direction:

```
ctgA  1,248,192-1,269,037  -   read      1- 21,000     <- the molecule
ctgA  1,248,422-1,249,382  +   read 21,145- 22,076     <- reading back
ctgA  1,249,517-1,250,021  +   read 22,129- 22,578
...
```

The Perl's `@taken` rule cannot see this, and neither could this port before
now: every one of those HSPs claims **fresh read bases**, so all of them are
legitimately kept. It is the reference that is covered twice.

**Now:** `classify.drop_fold_back` (§5 step 3b) drops an HSP that runs opposite
to the read's first HSP *and* whose reference span is already at least
`FOLD_BACK_OVERLAP` (50%) covered by HSPs kept for the same read on the same
contig. Both halves of the test are needed — the strand test alone would delete
genuine inversions, the overlap test alone would delete genuine tandem
duplications, and a 50% majority rather than a small tolerance is what lets a
read turn a corner inside an inverted repeat without being mistaken for an
artefact.

**Effect:** the read above was `LONG_DIVIDED`, flagged inverted, and produced a
junction event at every gap and an inversion event at every fold; it is now
`UNDIVIDED`, `PLAIN`, not inverted, with no events. Its 21 kb was also being
counted towards coverage roughly twice, and now is not. So on data carrying
these artefacts the long-divided count, the inversion count, the event count
and the mean coverage all fall. **How many reads were affected is reported** as
a front-page note and in the console summary — if that note is absent, this
change did not affect your run.

`--keep-fold-back` restores the previous behaviour exactly.

---

## 2. Behaviour that changes

| Was | Now | Finding |
|---|---|---|
| `*_Analysis.txt` written on **every** run — the default `[]` is true in Perl | written only with `-a`, and named `*_Analysis.tsv` | `B02` |
| Summary row 14 columns, contig rows put the accession in column 16, no header | one fixed 17-column layout with a header row and a `record_type` column; loads with `pandas.read_csv(sep="\t")` | `B19` |
| `-y --micro-match-length` documented as the HSP filter; the filter was hardcoded `100` | `--min-hsp-length` is the HSP filter; `--micro-match-length` is documented as what it actually does (the junction tolerance) | `B05` |
| `--space2reads` documented, `--space-to-reads` implemented | both work; help and front page are generated from the same table as the parser | `B18` |
| `--min-match-length` controlled three unrelated things | `--min-match-length` (segment filter), `--division-threshold` (short/long cut, defaults to the same value), `--index-window` (extraction granularity) | `B15` |
| Crash on high coverage: the lane search ran past the lane array and died on `substr(undef,…)` | overflowing reads are drawn on the outermost lane and counted, with a note telling you to raise `--page-scale` | `B03` |
| BTOP filename had to contain `btop`, else an interactive STDIN prompt | any path; `.gz` accepted | `B14` |
| `ps2pdf` / Ghostscript required; `.ps` intermediates left in the working directory | PDF written directly, nothing left behind | D2, `B29` |
| A `(` in any read name, product, path or header produced broken PostScript | not applicable — no PostScript | `B09` |
| `-e` list file contained only a two-line header | contains one line per extracted read with its coordinates | `B06` |
| `-e` on a `.gz` file lost the first record and desynchronised the parse | decompressed in-process, no `seek` on a pipe | `B07` |
| `-e` shelled out to `gzip -dc $file` unquoted; unavailable on Windows | Python's `gzip` module | `B08` |
| Extraction index built on every run whether or not `-e` was given | built only when `-e` is given | `B25` |
| Whole reference sequence held in memory although only its length was used | lengths only | `B26` |
| A file starting with neither `>` nor `LOCUS` died with an uninitialized-value error | clear message naming the file and showing its first line | `B12` |
| GenBank without `BioProject:` etc. died | missing metadata is simply absent from the front page | `B13` |
| Reads longer than one page were drawn wrong | any number of page crossings handled | — |
| Annotations crossing a page boundary ran off the edge | drawn on both pages | `B23` |
| `misc_feature` / `repeat_region` / orphan `CDS` folded into the preceding gene | recognised as their own features | `B24` |
| `ncRNA`, `tmRNA`, `misc_RNA` not typed as RNA | whole RNA family recognised | `B24` |
| `>gi\|12345\|ref\|NC_000913.3\|` parsed as accession `12345\|ref\|NC_000913.3\|` | parsed as `NC_000913.3` | `B21` |

---

## 3. Cosmetic changes

- **Lane assignment.** Reads are sorted by start position before packing, so
  the map reads as a pileup rather than file-order scatter. Class order is
  unchanged — long-divided reads still take the lanes nearest the axis.
- **Lanes are packed by interval, not by rightmost end.** Reads are placed in
  class order, so the sweep restarts at the left edge twice. Knowing only where
  a lane last ended, one long-divided read near the end of a contig closed that
  lane to every later class; a read bridging the origin of a circular replicon
  closed one from base 1 to the last base, because it has a piece at each end
  by definition. A handful of such reads reserved a band of lanes across the
  whole reference. Short-divided and undivided reads now fill those holes.
  **Effect:** the map uses fewer lanes and `overflow` fires less often. Within a
  class nothing moves — while requests ascend the two rules agree exactly.
- **Front page and plots are A3 landscape**; map pages keep the original
  10750 × 7600 pt poster format. The old front page was a 149 × 105 inch sheet
  with text in one corner.
- **Read-class colours changed.** The original used `0.5 0.0 0.0` for
  long-divided and `0.8 0.0 0.0` for inversions — two dark reds that are
  near-indistinguishable under colour-vision deficiency and both too dark
  against white. Long-divided reads are now orange and inversions a lighter
  red; undivided stays green and short-divided stays blue. Annotation colours
  are unchanged (they are already redundant with their `*`, `#`, `+`, `^`, `?`
  symbols).
- **Inversion line width** is `4 × --line-width`, capped so it cannot swamp the
  neighbouring lane. The Perl used `6 ×` in-page and `2 ×` across a page break.
- **The "read starts at the contig start" arrow** is drawn only on page 1. The
  Perl tested a page-local coordinate, so it drew that arrow on every page
  whose first read began near the fold.
- **Duplicate junction tick marks** removed (the first segment's coordinates
  were pushed into the mark list twice).
- **`Begin mapping`** is logged once, not twice.
- **Gene names on the map are no longer cut to their last five characters.**
  The Perl did that with an uncommented `substr($name,-5)` in the middle of its
  drawing loop (line 1141). It is right for a locus tag only by accident — a
  tag ends in its number — and on anything else it eats the name from the
  front, turning `fadE16_1` into `E16_1`. Map pages now draw
  `labels.gene_label(name)`, which is exactly what the browser draws, so a gene
  is named the same way on every surface. The full name is kept in the data
  model and used by the plots and the browser's detail panel.
- **Features with no gene name are labelled by their number.** Most bacterial
  features have no `/gene`, so their name is a locus tag — one prefix repeated
  across the whole assembly, with only a number to tell them apart, and on a map
  that number is pushed to the right into the neighbouring label.
  `SS37A_41660` is now drawn as `41660`; real names (`recA`) are untouched.
  The map and the browser use the identical rule. The full name is kept
  everywhere and search still matches it.

  A numeric suffix on a *real* name is left alone: `ftsH_5` is the fifth
  `ftsH`, not tag number 5, and stays `ftsH_5`. The two shapes are identical,
  so the test is the width of the number — a locus tag is zero-padded to four
  or five places, a paralog suffix is a small integer
  (`labels.LOCUS_TAG_MIN_DIGITS`).

---

## 4. New

- **A figures section** (`--no-plots` to skip): read classification, coverage
  profile per contig, read-length distribution, division-distance distribution
  with the cut-off drawn on it, a one-page structural-event map, and a
  per-contig table.
- **`--identity`** parses the BTOP trace column — which the Perl read out of
  every line and never used — into per-read alignment identity, and adds an
  identity histogram.
- **Warnings surface on the front page**, not only in debug files: unparsable
  read names, unknown accessions, lane overflow, `--cov-max` truncation,
  malformed BTOP lines, out-of-order records.
- **`--compat`** restores the Perl's hardcoded `100` HSP filter for
  side-by-side comparison. It does *not* restore the bugs in §1 — those are
  fixed unconditionally.
- **`--keep-fold-back`** keeps the fold-back tail of an end-ligation artefact,
  which is dropped by default (§1.8). The only switch in this program that
  turns a classification rule *off*; it exists because that rule is the one
  judgement call here that is not in the Perl, and a run that needs to be
  compared against an older one has to be able to make it.
- **An interactive genome browser**: `readmap browse <cache>` opens the run in a
  web browser — pan and zoom the reference, with a coverage track, the gene
  annotations, the read pileup coloured by class and carrying the map's own
  arrow vocabulary (direction chevron, double chevron for a contig join, ring
  for a circular join, reversed chevron for an inverted piece, bar for a
  junction), click-a-read detail, search, jump-to-next-structural-event, and
  export of the region on screen to PDF or PNG. The server binds `127.0.0.1`
  only.
- **GenBank feature details in the browser.** Hovering a feature bar under the
  coordinate axis names it, gives its key, strand, position and product;
  clicking it opens the full entry — gene, locus tag, product, note — in the
  side panel. `/product`, `/locus_tag` and `/note` are now retained by the
  parser (`/product` was already being read to decide the annotation group and
  then discarded) and stored in the cache. **This is cache format 2**: a cache
  written by an earlier version is rebuilt rather than misread.
- **A new output file, written by default.** Every run now also writes
  `<prefix>.readmapdb.npz`, the cache that browser reads. It appears beside the
  PDF and is typically tens of MB. Reading and classifying a gigabyte-scale
  BTOP file takes minutes; the cache makes every later browsing session start in
  about a second. `--no-cache` skips it; `--no-pdf` writes only it. The numbers
  in the PDF are unaffected either way.

---

## 5. Deliberately unchanged

The classification algorithm itself, step for step: the HSP length filter, the
repeat-collapse comparison (including its dimensionally odd
`|sstart - origin + qlen|` term — flagged as `B20`, preserved because changing
it would change published numbers), the greedy non-overlapping selection with
its `micro_match_length / 4` tolerance, the `|reference span − read span|`
discriminator, the 20 bp circularity window, and the 50 bp repeat window.

Only the implementation changed: interval arithmetic instead of per-base loops.

The one addition to the algorithm is the fold-back rule (§1.8), which is a new
step rather than a changed one: every step above runs exactly as it did, and the
new step then discards HSPs that all of them accepted. `--keep-fold-back` takes
it back out.
