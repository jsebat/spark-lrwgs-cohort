"""Rule-layer tests on synthetic read observations: every P8 class, for SNV-like, SV-like and TR-like inputs.
The matrix and rules are class-agnostic by design; the TR variant of each scenario carries per-read allele
lengths to check the TR summaries flow through."""
from phase_dnm.evidence import hapmatrix as H

CHROM, POS, PS = "chr1", 1_000_000, 900_000


def labels(child_hap1="P", f_trans=1, m_trans=1):
    t = H.LabelTables()
    t.orient[(CHROM, PS)].append((1, 2_000_000, "HAP1_PAT" if child_hap1 == "P" else "HAP1_MAT"))
    t.trans[("F", CHROM, PS)].append((1, 2_000_000, "HAP1" if f_trans == 1 else "HAP2"))
    t.trans[("M", CHROM, PS)].append((1, 2_000_000, "HAP1" if m_trans == 1 else "HAP2"))
    return t


def reads(role, hp, n_alt, n_ref, al_alt=None, al_ref=None, mapq=60):
    out = [H.ReadObs(role, hp, PS, "ALT", mapq=mapq, al=al_alt, nm_rate=0.01) for _ in range(n_alt)]
    out += [H.ReadObs(role, hp, PS, "REF", mapq=mapq, al=al_ref, nm_rate=0.01) for _ in range(n_ref)]
    return out


def scenario(child_A=(10, 0), child_O=(0, 10), F1=(0, 10), F2=(0, 10), M1=(0, 10), M2=(0, 10), tr=False, hap1_is="P", f_trans=1):
    """child hap1 = A (alt-carrying), hap2 = O; parent rows (alt, ref). TR: give the alt allele length 120, ref 30."""
    aa, ar = (120, 30) if tr else (None, None)
    rs = reads("C", 1, *child_A, al_alt=aa, al_ref=ar) + reads("C", 2, *child_O, al_alt=aa, al_ref=ar)
    rs += reads("F", 1, *F1, al_alt=aa, al_ref=ar) + reads("F", 2, *F2, al_alt=aa, al_ref=ar)
    rs += reads("M", 1, *M1, al_alt=aa, al_ref=ar) + reads("M", 2, *M2, al_alt=aa, al_ref=ar)
    L = labels(child_hap1=hap1_is, f_trans=f_trans)
    hp, cp = H.HapParams(), H.ClassParams()
    m = H.build_matrix(rs, L, CHROM, POS, hp)
    f = H.features(m, hp)
    t = H.transmission_features(m, hp)
    c = H.classify(m, f, t, hp, cp, L, CHROM, POS)
    return m, f, t, c


def test_germline_phased_all_classes():
    for tr in (False, True):
        m, f, t, c = scenario(tr=tr)
        assert c["phase_class"] == "germline_DNM_phased", (tr, c)
        assert f["hap_obs_k5"] == 6 and f["c_alt_confined"] == 1 and f["c_alt_hap_frac"] == 1.0
        assert t["parent_of_origin"] == "paternal" and t["poo_confidence"] == 1.0 and t["t_alt_reads"] == 0 and t["t_hap_resolved"] == 1
        assert c["rule_score"] == 6 and c["flags"] == ""
        if tr:
            assert f["c_tr_al_hapA_mean"] == 120 and f["c_tr_al_hapO_mean"] == 30


def test_parent_of_origin_follows_orientation_and_transmission():
    # child alt on hap1, block says hap1 is MATERNAL; mother's transmitted haplotype is HAP2 -> t_* read from M2
    rs = reads("C", 1, 10, 0) + reads("C", 2, 0, 10) + reads("F", 1, 0, 10) + reads("F", 2, 0, 10) + reads("M", 1, 0, 10) + reads("M", 2, 0, 10)
    L = labels(child_hap1="M", m_trans=2)
    hp, cp = H.HapParams(), H.ClassParams()
    m = H.build_matrix(rs, L, CHROM, POS, hp)
    t = H.transmission_features(m, hp)
    assert t["parent_of_origin"] == "maternal" and t["t_dp"] == 10 and t["u_dp"] == 10 and t["nt_parent_alt_reads"] == 0


def test_inherited_missed_in_parent():
    # father's transmitted haplotype (F1) carries the alt in 6 of 10 reads: the parent is really het there
    m, f, t, c = scenario(F1=(6, 4))
    assert c["phase_class"] == "inherited_missed_in_parent"
    assert t["t_alt_reads"] == 6 and f["p_max_alt_any_hap"] == 6


def test_parental_mosaic_transmitted_and_its_floor():
    # 3 alt reads of 20 on the transmitted haplotype (15%): mosaic, floors met
    m, f, t, c = scenario(F1=(3, 17))
    assert c["phase_class"] == "parental_mosaic_transmitted", c
    # 2 alt of 8 reads on T: above the error allowance (max(1, 0.05*8) = 1) but depth below the floor -> inconclusive + flag (P10)
    m, f, t, c = scenario(F1=(2, 6))
    assert c["phase_class"] == "inconclusive" and "PARENTAL_ALT_LOW_DEPTH" in c["flags"], c
    # 1 alt of 8 is WITHIN the allowance: not evidence of anything, the germline call stands
    m, f, t, c = scenario(F1=(1, 7))
    assert c["phase_class"] == "germline_DNM_phased"


def test_child_postzygotic_mosaic_and_its_floor():
    # alt on 10 of 20 reads of the child's alt haplotype, nothing anywhere else
    m, f, t, c = scenario(child_A=(10, 10), child_O=(0, 20), F1=(0, 12), F2=(0, 12), M1=(0, 12), M2=(0, 12))
    assert c["phase_class"] == "child_postzygotic_mosaic", c
    assert f["c_alt_hap_frac"] == 0.5
    # 4 of 8 reads: fraction alike, depth below the floor
    m, f, t, c = scenario(child_A=(4, 4), child_O=(0, 8))
    assert c["phase_class"] == "inconclusive" and "MOSAIC_UNDERPOWERED" in c["flags"]


def test_phase_conflict_artifact():
    # alt on both child haplotypes
    m, f, t, c = scenario(child_A=(6, 4), child_O=(5, 5))
    assert c["phase_class"] == "phase_conflict_artifact" and t["poo_reason"] == "ALT_ON_BOTH_CHILD_HAPLOTYPES"
    # alt on two parental haplotypes with hom-ref parental GTs
    m, f, t, c = scenario(F2=(4, 6), M1=(4, 6))
    assert c["phase_class"] == "phase_conflict_artifact"


def test_germline_unphased_when_haplotypes_missing():
    # father's haplotype 2 unobserved (0 reads) -> hap_obs 5 -> cannot be 'phased', falls to unphased
    m, f, t, c = scenario(F2=(0, 0))
    assert c["phase_class"] == "germline_DNM_unphased" and "LOW_HAP_DEPTH" in c["flags"] and f["hap_obs_k5"] == 5
    # child block unoriented -> parent of origin undetermined with the reason, class still unphased
    rs = reads("C", 1, 10, 0) + reads("C", 2, 0, 10) + reads("F", 1, 0, 10) + reads("F", 2, 0, 10) + reads("M", 1, 0, 10) + reads("M", 2, 0, 10)
    L = H.LabelTables()          # no labels at all
    hp, cp = H.HapParams(), H.ClassParams()
    m = H.build_matrix(rs, L, CHROM, POS, hp)
    f = H.features(m, hp); t = H.transmission_features(m, hp); c = H.classify(m, f, t, hp, cp, L, CHROM, POS)
    assert t["parent_of_origin"] == "undetermined" and t["poo_reason"] == "CHILD_BLOCK_UNORIENTED"
    assert c["phase_class"] == "germline_DNM_unphased" and "CHILD_BLOCK_UNORIENTED" in c["flags"]


def test_untagged_reads_count_for_child_support_but_not_haplotypes():
    rs = reads("C", 1, 5, 0) + reads("C", 2, 0, 10) + [H.ReadObs("C", None, None, "ALT") for _ in range(4)]
    rs += reads("F", 1, 0, 10) + reads("F", 2, 0, 10) + reads("M", 1, 0, 10) + reads("M", 2, 0, 10)
    hp = H.HapParams(); L = labels()
    m = H.build_matrix(rs, L, CHROM, POS, hp); f = H.features(m, hp)
    assert f["c_untagged_alt"] == 4 and f["c_alt_hapA"] == 5 and abs(f["c_alt_tagged_frac"] - 5 / 9) < 1e-3


def test_low_mapq_reads_are_dropped():
    rs = reads("C", 1, 10, 0) + reads("C", 2, 0, 10, mapq=5) + reads("F", 1, 0, 10) + reads("F", 2, 0, 10) + reads("M", 1, 0, 10) + reads("M", 2, 0, 10)
    hp = H.HapParams(); L = labels()
    m = H.build_matrix(rs, L, CHROM, POS, hp)
    assert m.row("C", 2).dp == 0 and m.n_reads == 50


def test_single_alt_read_is_not_a_germline_call():
    # 1 alt read on A (frac 1.0), everything else clean: an evidence floor, not a DNM
    m, f, t, c = scenario(child_A=(1, 9))
    assert c["phase_class"] == "inconclusive" and "TOO_FEW_ALT_READS" in c["flags"]
    # 3 of 10 clears the floor but is a mosaic-like fraction at depth 10: under-powered mosaic, not germline
    m, f, t, c = scenario(child_A=(3, 7))
    assert c["phase_class"] == "inconclusive" and "MOSAIC_UNDERPOWERED" in c["flags"]
    # 8 of 10 on A, clean elsewhere: germline
    m, f, t, c = scenario(child_A=(8, 2))
    assert c["phase_class"] == "germline_DNM_phased"


def test_germline_phased_requires_all_parental_haplotypes_clean():
    # non-transmitting parent (mother) carries alt on one haplotype: not phased-germline
    m, f, t, c = scenario(M2=(4, 6))
    assert c["phase_class"] != "germline_DNM_phased"
    # non-transmitting parent's haplotype unobserved: cannot be 'phased' -> unphased
    m, f, t, c = scenario(M2=(0, 0))
    assert c["phase_class"] == "germline_DNM_unphased" and f["hap_obs_k5"] == 5


def test_label_lookup_uses_nearest_segment_within_margin():
    t = H.LabelTables()
    t.orient[(CHROM, PS)].append((1_000_000, 1_100_000, "HAP1_MAT"))
    assert t.child_hap1_is(CHROM, PS, 1_050_000) == "M"          # inside
    assert t.child_hap1_is(CHROM, PS, 1_120_000) == "M"          # 20 kb past the last het, same block: still the block
    assert t.child_hap1_is(CHROM, PS, 1_200_000) is None         # 100 kb away: not this block
    assert t.child_hap1_is(CHROM, PS + 1, 1_050_000) is None     # different phase set


def test_change_point_flags():
    L = labels()
    L.changes[("F", CHROM)] = [POS + 20_000]
    L.child_switches[CHROM] = [POS - 10_000]
    rs = reads("C", 1, 10, 0) + reads("C", 2, 0, 10) + reads("F", 1, 0, 10) + reads("F", 2, 0, 10) + reads("M", 1, 0, 10) + reads("M", 2, 0, 10)
    hp, cp = H.HapParams(), H.ClassParams()
    m = H.build_matrix(rs, L, CHROM, POS, hp); f = H.features(m, hp); t = H.transmission_features(m, hp)
    c = H.classify(m, f, t, hp, cp, L, CHROM, POS)
    assert "NEAR_CHANGE_POINT_F" in c["flags"] and "NEAR_CHILD_SWITCH" in c["flags"]
    assert c["phase_class"] == "germline_DNM_phased"     # flags inform, they do not reclassify


def reads_amb(role, hp, n_amb, mapq=60):
    return [H.ReadObs(role, hp, PS, "AMB", mapq=mapq) for _ in range(n_amb)]


def test_unreadable_reads_do_not_observe_a_haplotype():
    # the cohort pattern (2026-09-13): father's transmitted haplotype has 7 reads, 1 ALT + 6 unreadable (deleted at the
    # site). dp 7 >= k but only ONE readable read: T is unobserved, this cannot be 'germline_DNM_phased'
    rs = reads("C", 1, 10, 0) + reads("C", 2, 0, 10) + reads("F", 1, 1, 0) + reads_amb("F", 1, 6) + reads("F", 2, 0, 10)
    rs += reads("M", 1, 0, 10) + reads("M", 2, 0, 10)
    L = labels(); hp, cp = H.HapParams(), H.ClassParams()
    m = H.build_matrix(rs, L, CHROM, POS, hp); f = H.features(m, hp); t = H.transmission_features(m, hp)
    c = H.classify(m, f, t, hp, cp, L, CHROM, POS)
    assert m.row("F", 1).dp == 7 and m.row("F", 1).n == 1 and abs(m.row("F", 1).amb_frac - 6 / 7) < 1e-6
    assert f["hap_obs_k5"] == 5 and f["p_n_haps_obs_k5"] == 3 and f["p_amb_frac_max"] == 0.857 and f["c_amb_frac_hapA"] == 0.0
    assert c["phase_class"] != "germline_DNM_phased"
    assert "AMBIGUOUS_READS" in c["flags"] and "LOW_HAP_DEPTH" in c["flags"]
    assert t["t_dp"] == 7 and t["t_alt_reads"] == 1
    # the same haplotype with 7 readable ref reads and one alt: within the allowance, observed, germline
    rs2 = reads("C", 1, 10, 0) + reads("C", 2, 0, 10) + reads("F", 1, 1, 7) + reads("F", 2, 0, 10) + reads("M", 1, 0, 10) + reads("M", 2, 0, 10)
    m2 = H.build_matrix(rs2, L, CHROM, POS, hp); f2 = H.features(m2, hp); t2 = H.transmission_features(m2, hp)
    c2 = H.classify(m2, f2, t2, hp, cp, L, CHROM, POS)
    assert c2["phase_class"] == "germline_DNM_phased" and "AMBIGUOUS_READS" not in c2["flags"]


def test_error_allowance_is_on_readable_reads():
    # 20 reads on T of which 18 unreadable: 2 alt of 2 readable is NOT within max(1, 0.05*20) = 1 - it is alt on the
    # transmitted haplotype with no ref read to contradict it; T is unobserved (n=2 < k) and the class cannot be phased
    rs = reads("C", 1, 10, 0) + reads("C", 2, 0, 10) + reads("F", 1, 2, 0) + reads_amb("F", 1, 18) + reads("F", 2, 0, 10)
    rs += reads("M", 1, 0, 10) + reads("M", 2, 0, 10)
    L = labels(); hp, cp = H.HapParams(), H.ClassParams()
    m = H.build_matrix(rs, L, CHROM, POS, hp); f = H.features(m, hp); t = H.transmission_features(m, hp)
    c = H.classify(m, f, t, hp, cp, L, CHROM, POS)
    assert c["phase_class"] not in ("germline_DNM_phased", "germline_DNM_unphased"), c
    assert f["p_n_haps_with_alt"] == 1


def _review_style_row(m, f, t, c):
    """Mimic review.py's row construction for the columns reclassify_row reads."""
    row = {}
    row.update(c); row.update(t); row.update(f)
    row.update(child_hap1_is=m.child_hap1_is or ".", F_transmitted_hap=m.transmitted.get("F") or ".",
               M_transmitted_hap=m.transmitted.get("M") or ".", n_reads_used=m.n_reads)
    for (role, hp), pre in ((("C", 1), "C1"), (("C", 2), "C2"), (("F", 1), "F1"), (("F", 2), "F2"), (("M", 1), "M1"), (("M", 2), "M2")):
        s = m.row(role, hp).summary()
        for fld in ("dp", "alt", "ref", "amb", "mapq_mean"):
            row["%s_%s" % (pre, fld)] = s.get(fld)
    for role in ("C", "F", "M"):
        row["%s_untagged_dp" % role] = m.untagged[role].dp
        row["%s_untagged_alt" % role] = m.untagged[role].alt
    return {k: ("" if v is None else str(v)) for k, v in row.items()}


def test_reclassify_from_evidence_row_reproduces_the_rule_layer():
    hp, cp = H.HapParams(), H.ClassParams()
    cases = [dict(), dict(F1=(6, 4)), dict(F1=(3, 17)), dict(child_A=(10, 10), child_O=(0, 20), F1=(0, 12), F2=(0, 12), M1=(0, 12), M2=(0, 12)),
             dict(child_A=(6, 4), child_O=(5, 5)), dict(F2=(0, 0)), dict(child_A=(1, 9)), dict(M2=(4, 6)), dict(hap1_is="M", f_trans=2)]
    for kw in cases:
        m, f, t, c = scenario(**kw)
        # 1. without a positional flag the rule layer and the parent of origin round-trip exactly
        row = _review_style_row(m, f, t, c)
        out = H.reclassify_row(row, hp, cp, thresholds_version="test")
        assert out["phase_class"] == c["phase_class"], (kw, out["phase_class"], c["phase_class"])
        assert out["rule_score"] == c["rule_score"] and out["parent_of_origin"] == t["parent_of_origin"]
        assert out["t_alt_reads"] == t["t_alt_reads"] and out["hap_obs_k5"] == f["hap_obs_k5"]
        assert out["thresholds_version"] == "test"
        # 2. a positional flag must survive, and inside a phase-switch window the parent of origin is deliberately
        #    cleared (P34, 2026-09-16): the labels may be the wrong way round there, so the assignment is a coin toss
        row = _review_style_row(m, f, t, c)
        row["flags"] = (row["flags"] + ";NEAR_CHILD_SWITCH").strip(";")
        out = H.reclassify_row(row, hp, cp, thresholds_version="test")
        assert "NEAR_CHILD_SWITCH" in out["flags"]
        assert out["phase_class"] == c["phase_class"] and out["rule_score"] == c["rule_score"]
        if t["parent_of_origin"] in ("paternal", "maternal"):
            assert out["parent_of_origin"] == "undetermined", (kw, out["parent_of_origin"])
            assert str(out["poo_reason"]).startswith("PHASE_SWITCH_RISK"), out["poo_reason"]
        else:
            assert out["parent_of_origin"] == t["parent_of_origin"]
        assert out["c_alt_mapq_mean"] == row["c_alt_mapq_mean"]      # read-quality summary carried over, not recomputed
    # and a row from a table written before the ambiguity columns existed (no *_amb) still reclassifies
    m, f, t, c = scenario()
    row = _review_style_row(m, f, t, c)
    for k in [k for k in row if k.endswith("_amb")]:
        del row[k]
    assert H.reclassify_row(row, hp, cp)["phase_class"] == "germline_DNM_phased"
