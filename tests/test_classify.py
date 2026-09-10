"""Golden tests for the classification algorithm (``code_structure.md`` §5).

If any of these change, the scientific output of the program has changed.
That is allowed only with a corresponding entry in ``CHANGES.md``.
"""

from __future__ import annotations

import random

from readmap.classify import (
    ClassifyStats,
    classify_hits,
    classify_stream,
    drop_fold_back,
    select_hits,
    select_non_overlapping,
)
from readmap.inputs.btop import BtopStats, iter_read_groups
from readmap.inputs.reference import load_reference
from readmap.labels import LabelAssigner
from readmap.models import Hit, ReadClass
from readmap.params import Params
from tests.conftest import DROPPED, EXPECTED


def _classify(params: Params):
    reference = load_reference(params.seq_file)
    labeller = LabelAssigner()
    stats = ClassifyStats()
    groups = iter_read_groups(params.btop_file, params.min_read_length, BtopStats())
    result = list(classify_stream(groups, reference, params, labeller, stats))
    return result, labeller, stats


def test_every_read_classified_as_expected(params: Params) -> None:
    groups, _labeller, _stats = _classify(params)
    seen = {g.read: (g.cls, g.junction, g.inverted) for g in groups}

    assert seen == EXPECTED


def test_dropped_reads_are_dropped(params: Params) -> None:
    groups, _labeller, stats = _classify(params)
    assert DROPPED.isdisjoint({g.read for g in groups})
    assert stats.unknown_contigs == {"ctgZZZ"}
    assert stats.unknown_contig_reads == 1


def test_class_totals(params: Params) -> None:
    groups, _labeller, _stats = _classify(params)
    counts = {cls: sum(1 for g in groups if g.cls is cls) for cls in ReadClass}

    assert counts[ReadClass.UNDIVIDED] == 4
    assert counts[ReadClass.SHORT_DIVIDED] == 2
    assert counts[ReadClass.LONG_DIVIDED] == 4
    assert sum(1 for g in groups if g.inverted) == 1


def test_read_names_are_never_rewritten(params: Params) -> None:
    """Finding B01: the Perl merged any read whose name failed one regex."""
    groups, labeller, _stats = _classify(params)
    names = {g.read for g in groups}

    assert "0a1b2c3d-e4f5-6789-abcd-ef0123456789" in names
    assert "read(x).9" in names
    # Both of those fail the legacy pattern and are reported as such.
    assert labeller.unmatched == 2
    assert labeller.warnings()


def test_labels_are_unique(params: Params) -> None:
    groups, _labeller, _stats = _classify(params)
    labels = [g.label for g in groups]
    assert len(labels) == len(set(labels))


def test_competing_hsps_are_dropped(params: Params) -> None:
    """read.12's third HSP claims read bases the second already took."""
    groups, _labeller, _stats = _classify(params)
    read12 = next(g for g in groups if g.read == "read.12")
    assert len(read12.hits) == 2
    assert [h.sstart for h in read12.hits] == [1000, 20000]


def test_interval_selection_matches_a_naive_per_base_implementation() -> None:
    """The interval sweep must agree with the Perl's per-base ``@taken`` array."""

    def naive(hits: list[Hit], tolerance: float) -> list[Hit]:
        taken: set[int] = set()
        kept: list[Hit] = []
        for hit in hits:
            lo, hi = hit.qlow, hit.qhigh
            if any(k in taken for k in range(lo, hi + 1)):
                continue
            import math

            start = math.floor(lo + tolerance)
            end = math.floor(hi - tolerance)
            if start <= end:
                taken.update(range(start, end + 1))
            kept.append(hit)
        return kept

    rng = random.Random(20260729)
    for _ in range(2_000):
        tolerance = rng.choice([0.0, 2.5, 25.0, 60.0])
        hits = []
        for _ in range(rng.randint(1, 6)):
            start = rng.randint(1, 400)
            length = rng.randint(1, 200)
            hits.append(
                Hit(
                    read="r", strand=1, qstart=start, qend=start + length,
                    sstart=start, send=start + length, qlen=1000,
                    contig="c", btop="",
                )
            )
        fast, _dropped = select_non_overlapping(hits, tolerance)
        assert fast == naive(hits, tolerance)


# --------------------------------------------------------------------------
# Fold-back (end-ligation) artefacts
# --------------------------------------------------------------------------

# The real read this rule was written for, from AP027142.  One 21 kb HSP on the
# minus strand, then eight HSPs re-reading that same 21 kb forwards.  Written as
# BLAST writes it: sstart > send on the minus strand.
_FOLD_BACK = [
    (-1, 1, 21_000, 1_269_037, 1_248_192),
    (1, 21_145, 22_076, 1_248_422, 1_249_382),
    (1, 22_129, 22_578, 1_249_517, 1_250_021),
    (1, 22_818, 22_936, 1_250_364, 1_250_482),
    (1, 23_310, 26_379, 1_250_989, 1_254_540),
    (1, 26_566, 30_037, 1_254_769, 1_258_475),
    (1, 31_817, 32_302, 1_260_730, 1_261_261),
    (1, 32_365, 37_550, 1_261_356, 1_266_879),
    (1, 37_760, 39_181, 1_267_529, 1_268_998),
]


def _hits(rows, contig="AP027142", read="fold.1", qlen=39_500):
    return [
        Hit(read=read, strand=s, qstart=qs, qend=qe, sstart=ss, send=se,
            qlen=qlen, contig=contig, btop="")
        for s, qs, qe, ss, se in rows
    ]


def test_fold_back_tail_is_dropped(params: Params) -> None:
    """The tail re-reads the anchor's reference backwards -- an artefact."""
    kept, _short, _overlap, folded = select_hits(_hits(_FOLD_BACK), params)

    assert folded == 8
    assert len(kept) == 1
    assert kept[0].qend == 21_000

    group = classify_hits(kept, "R_1", {"AP027142": 3_000_000}, params)
    assert group.cls is ReadClass.UNDIVIDED
    assert group.inverted is False
    assert group.distance == 0


def test_fold_back_is_kept_with_the_flag(params: Params) -> None:
    """--keep-fold-back restores every version before this rule existed."""
    from dataclasses import replace as replace_dc

    kept, _short, _overlap, folded = select_hits(
        _hits(_FOLD_BACK), replace_dc(params, keep_fold_back=True)
    )

    assert folded == 0
    assert len(kept) == 9

    group = classify_hits(kept, "R_1", {"AP027142": 3_000_000}, params)
    assert group.cls is ReadClass.LONG_DIVIDED
    assert group.inverted is True


def test_fold_back_rule_spares_a_genuine_inversion() -> None:
    """An inverted piece over *fresh* reference is evidence, not an artefact.

    This is fixture read.4's shape: the second HSP runs the other way but lands
    beyond the first, so nothing is re-read and both HSPs survive.  The golden
    map in EXPECTED asserts the same thing end to end -- read.4 stays
    SHORT_DIVIDED and inverted.
    """
    hits = _hits([
        (1, 10, 1_510, 26_000, 27_500),
        (-1, 1_600, 3_100, 29_200, 27_700),
    ], contig="ctgA", read="inv.1", qlen=3_200)

    kept, dropped = drop_fold_back(hits)
    assert dropped == 0
    assert kept == hits


def test_fold_back_rule_spares_a_tandem_duplication() -> None:
    """Two passes over one stretch on the *same* strand are a duplication.

    The read carries two copies of a region the reference has once.  The
    reference overlap is total, so only the strand test tells this apart from a
    fold-back -- which is why the rule needs both halves.
    """
    hits = _hits([
        (1, 10, 2_010, 5_000, 7_000),
        (1, 2_100, 4_100, 5_000, 7_000),
    ], contig="ctgA", read="dup.1", qlen=4_200)

    kept, dropped = drop_fold_back(hits)
    assert dropped == 0
    assert kept == hits


def test_fold_back_rule_needs_a_majority_overlap() -> None:
    """A breakpoint inside an inverted repeat shares reference but is not a fold.

    The inverted piece re-reads 400 bp of a 2 000 bp HSP -- the repeat the
    inversion occurred between -- and then carries on into reference no part of
    the read has seen. Below FOLD_BACK_OVERLAP, so it is kept.
    """
    hits = _hits([
        (1, 10, 2_010, 10_000, 12_000),
        (-1, 2_100, 4_100, 14_000, 11_600),
    ], contig="ctgA", read="rep.1", qlen=4_200)

    kept, dropped = drop_fold_back(hits)
    assert dropped == 0
    assert kept == hits


def test_genbank_and_fasta_agree(params: Params, fixture_dir) -> None:
    """The same BTOP file must classify identically against either reference."""
    from dataclasses import replace as replace_dc

    fasta_groups, _l, _s = _classify(params)
    gb_params = replace_dc(params, seq_file=str(fixture_dir / "ref.gb"))
    gb_groups, _l2, _s2 = _classify(gb_params)

    assert {g.read: g.cls for g in fasta_groups} == {
        g.read: g.cls for g in gb_groups
    }
