"""The quality rubric, pinned on three synthetic rows shaped like the cases that define it.

These rows carry the evidence values of the cohort's two known pathogenic deletions (rounded) and of a typical
inherited deletion, with no identifiers. They pin BEHAVIOUR, not truth: P31 records that the deletion-path constants
were chosen with these events in view, so passing here is a regression guard, not validation.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from clinical import quality as Q  # noqa: E402


def mecp2_like():
    # 35 kb heterozygous deletion: depth 0.59, LOH inside, junction reads both ends, parents intact
    return dict(svtype="DEL", svlen="35454", start="1000000", end="1035454", phase_class="germline_DNM_unphased",
                rule_score="6", rf_prob="0.64", cohort_AC_loo="0", pon_founder_recurrence_loo="0", mask_overlap="0",
                sv_depth_ratio_inside_flank="0.5872", sv_het_snv_persistence="0.0", C_sv_junc_both_ends="1",
                C_sv_junc_hap1="2", C_sv_junc_hap2="3", C_sv_junc_hapu="5", F_sv_ratio_all="1.02", M_sv_ratio_all="0.98",
                F_sv_junc_hap1="0", F_sv_junc_hap2="0", F_sv_junc_hapu="0", M_sv_junc_hap1="0", M_sv_junc_hap2="0",
                M_sv_junc_hapu="0", child_sv_support="10", max_parent_sv_support="0")


def dnmt3a_like():
    # 302 bp coding deletion: depth 0.38 and phased, but no junction read at both ends and no het measurement (< 5 kb)
    return dict(svtype="DEL", svlen="302", start="2000000", end="2000302", phase_class="germline_DNM_phased",
                rule_score="5", rf_prob="0.72", cohort_AC_loo="0", pon_founder_recurrence_loo="0", mask_overlap="0",
                sv_depth_ratio_inside_flank="0.3814", sv_het_snv_persistence="", C_sv_junc_both_ends="0",
                C_sv_junc_hap1="4", C_sv_junc_hap2="0", C_sv_junc_hapu="1", F_sv_ratio_all="0.97", M_sv_ratio_all="1.01",
                F_sv_junc_hap1="0", F_sv_junc_hap2="0", F_sv_junc_hapu="0", M_sv_junc_hap1="0", M_sv_junc_hap2="0",
                M_sv_junc_hapu="0", child_sv_support="5", max_parent_sv_support="0")


def inherited_like():
    r = mecp2_like()
    r.update(phase_class="inherited_missed_in_parent", F_sv_ratio_all="0.52")
    return r


def test_mecp2_like_passes():
    v = Q.sv_quality(mecp2_like())
    assert v.verdict == Q.PASS, v.reasons
    assert "+CHILD_DEPTH_HALVED" in v.reasons and "+LOSS_OF_HETEROZYGOSITY" in v.reasons


def test_dnmt3a_like_is_review_not_fail():
    v = Q.sv_quality(dnmt3a_like())
    assert v.verdict == Q.REVIEW, v.reasons
    assert "-JUNCTION_BOTH_ENDS" in v.reasons and "+CHILD_DEPTH_HALVED" in v.reasons


def test_inherited_fails_on_parent_depth_and_phase_class():
    v = Q.sv_quality(inherited_like())
    assert v.verdict == Q.FAIL
    assert "INHERITED_MISSED_IN_PARENT" in v.reasons
    assert any(x.startswith("PARENT_DEPTH_DEPLETED") for x in v.reasons)


def test_seen_in_another_family_fails_regardless_of_evidence():
    r = mecp2_like(); r["cohort_AC_loo"] = "2"
    assert Q.sv_quality(r).verdict == Q.FAIL


def test_hets_persisting_with_depth_loss_is_review_mosaic():
    r = mecp2_like(); r["sv_het_snv_persistence"] = "0.9"
    v = Q.sv_quality(r)
    assert v.verdict == Q.REVIEW and any("HETS_PERSIST" in x for x in v.reasons)


def test_smallvar_pass_and_parent_alt_fail():
    ok = dict(phase_class="germline_DNM_phased", rule_score="6", rf_prob="0.99", c_alt_hapA="9", c_untagged_alt="1",
              c_alt_confined="1", p_max_alt_any_hap="0", t_alt_reads="0", child_DP="22", cohort_AC_loo="0",
              pon_founder_recurrence_loo="0", gnomad_af="", mask_overlap="0")
    assert Q.smallvar_quality(ok).verdict == Q.PASS
    bad = dict(ok, p_max_alt_any_hap="4")
    v = Q.smallvar_quality(bad)
    assert v.verdict == Q.FAIL and any(x.startswith("PARENT_ALT_READS") for x in v.reasons)
    common = dict(ok, gnomad_af="0.01")
    assert Q.smallvar_quality(common).verdict == Q.FAIL


def test_mosaic_is_review():
    r = dict(phase_class="child_postzygotic_mosaic", rule_score="3", rf_prob="0.5", c_alt_hapA="5", c_untagged_alt="0",
             c_alt_confined="1", p_max_alt_any_hap="0", child_DP="30", cohort_AC_loo="0", pon_founder_recurrence_loo="0")
    v = Q.smallvar_quality(r)
    assert v.verdict == Q.REVIEW and "POSTZYGOTIC_MOSAIC" in v.reasons


if __name__ == "__main__":
    import inspect
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and inspect.isfunction(fn):
            try:
                fn(); print("ok  ", name)
            except AssertionError as e:
                fails += 1; print("FAIL", name, e)
    sys.exit(1 if fails else 0)
