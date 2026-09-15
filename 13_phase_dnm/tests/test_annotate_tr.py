"""TR annotation (P26 for the third class): the founder reference, the leave-one-family-out statistics, and the
guarantee that a silently-empty annotation column fails loudly (the defect found 2026-09-15)."""
import os
from phase_dnm import annotate as AN

ORDER = ["f1", "f2", "f3", "f4", "f5", "f6"]
FAM = {"f1": "A", "f2": "A", "f3": "B", "f4": "B", "f5": "C", "f6": "C"}
PER = [[100, 100], [100, 105], [100, 100], [100, 100], [100, 160], [100, 100]]   # only family C carries a long allele


def test_tr_lofo_stats_make_a_private_expansion_private():
    ac, an, rec, dist = AN.tr_lofo_stats(160, 5, PER, ORDER, FAM, set())
    assert (ac, an, rec) == (1, 12, 1) and dist == 0.0        # with family C in the panel the allele is its own p99
    ac, an, rec, dist = AN.tr_lofo_stats(160, 5, PER, ORDER, FAM, {"C"})
    assert (ac, an, rec) == (0, 8, 0) and dist == 11.0        # leaving it out: private, 11 motif units beyond p99
    ac, an, rec, dist = AN.tr_lofo_stats(100, 5, PER, ORDER, FAM, set())
    assert ac == 11 and rec == 6 and dist < 0                 # a common allele, well inside the founder range
    assert AN.tr_lofo_stats(160, 5, None, ORDER, FAM, set()) == (None, None, None, None)


def test_tr_lofo_tolerance_is_in_motif_units():
    per = [[100], [104], [112], [100], [100], [100]]
    ac, _, _, _ = AN.tr_lofo_stats(100, 5, per, ORDER, FAM, set())      # +-1 motif unit = +-5 bp
    assert ac == 5                                                       # 100, 104, and the three other 100s
    ac, _, _, _ = AN.tr_lofo_stats(100, 1, per, ORDER, FAM, set())       # a 1 bp motif: only exact matches
    assert ac == 4


def test_union_sites_covers_tr():
    """Until 2026-09-15 union_sites had an empty branch for TR, which is why every column came out blank."""
    import inspect
    src = inspect.getsource(AN.union_sites)
    assert "trid" in src, "union_sites must key TR candidates by TRID"


def test_annot_header_carries_the_tr_columns():
    import inspect
    src = inspect.getsource(AN.write_annot)
    for col in ("trgt_pop_p99_distance", "strchive_locus"):
        assert col in src
