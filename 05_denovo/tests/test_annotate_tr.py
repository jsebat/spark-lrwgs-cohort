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


def test_union_sites_keys_tr_by_trid(tmp_path):
    """Until 2026-09-15 union_sites had an empty branch for TR, which is why every column came out blank. Behaviour, not
    source text: two children sharing a TRID collapse to one site keyed by the TRID; a different TRID is a second site."""
    from phase_dnm.records import CandidateRecord, write_candidates
    recs = [CandidateRecord("fam", "kid", "tr1", "chr1", 1000, 1060, ".", "CAG", "TR", "trgt", class_payload={"trid": "chr1_1000_1060_CAG"}),
            CandidateRecord("fam", "kid2", "tr1b", "chr1", 1000, 1060, ".", "CAG", "TR", "trgt", class_payload={"trid": "chr1_1000_1060_CAG"}),
            CandidateRecord("fam", "kid", "tr2", "chr2", 500, 530, ".", "A", "TR", "trgt", class_payload={"trid": "chr2_500_530_A"})]
    p1 = str(tmp_path / "a.tsv"); p2 = str(tmp_path / "b.tsv")
    write_candidates(recs[:1] + recs[2:], p1); write_candidates(recs[1:2], p2)
    sites, _ = AN.union_sites([p1, p2], "tr")
    assert sites == [("chr1", 1000, "chr1_1000_1060_CAG", "1060"), ("chr2", 500, "chr2_500_530_A", "530")]


def test_write_annot_emits_the_tr_columns_and_a_row_per_candidate(tmp_path):
    """The annotation table for a TR child carries the TR-specific columns (P26) and one row per candidate even when
    no founder table is available (columns empty, never absent)."""
    from phase_dnm.records import CandidateRecord, write_candidates
    recs = [CandidateRecord("fam", "kid", "tr1", "chr1", 1000, 1060, ".", "CAG", "TR", "trgt", class_payload={"trid": "chr1_1000_1060_CAG"})]
    cp = str(tmp_path / "c.tsv"); write_candidates(recs, cp)
    out = str(tmp_path / "annot.tsv")
    AN.write_annot(cp, out, "tr", {}, {}, [], {}, set(), set())
    lines = open(out).read().splitlines()
    hdr = lines[0].split("	")
    assert hdr[:2] == ["variant_id", "gnomad_af"] and "trgt_pop_p99_distance" in hdr and "strchive_locus" in hdr
    assert len(lines) == 2 and lines[1].split("	")[0] == "tr1"
