"""SV interval evidence (P17): depth across the interval per haplotype, which alt-read confinement cannot supply."""
from phase_dnm.evidence import hapmatrix as H


def test_sv_depth_features_order_haplotypes_by_evidence_not_by_label():
    # heterozygous deletion: one child haplotype loses depth, the other keeps it, both parents normal
    raw = dict(C_sv_ratio_hap1=0.05, C_sv_ratio_hap2=0.98, C_sv_ratio_all=0.52,
               F_sv_ratio_hap1=1.02, F_sv_ratio_hap2=0.97, M_sv_ratio_hap1=0.99, M_sv_ratio_hap2=1.03,
               C_sv_junc_hap1=4, C_sv_junc_hap2=0, C_sv_junc_hapu=3, C_sv_junc_both_ends=1)
    f = H.sv_depth_features(raw)
    assert f["c_sv_hap_depth_change_A"] == 0.05 and f["c_sv_hap_depth_change_O"] == 0.98
    assert f["sv_depth_ratio_inside_flank"] == 0.52
    assert f["p_sv_max_hap_depth_change"] == 0.03          # no parent loses depth -> de novo consistent
    assert f["c_sv_junction_hap_concentration"] == 1.0     # the tagged junction reads sit on one haplotype
    assert f["c_sv_junc_untagged_frac"] == round(3 / 7, 4)
    # the same event inherited: a parent loses depth too
    raw2 = dict(raw, F_sv_ratio_hap1=0.04)
    assert H.sv_depth_features(raw2)["p_sv_max_hap_depth_change"] == 0.96
    # ordering is by deviation from 1, whichever haplotype it is
    raw3 = dict(raw, C_sv_ratio_hap1=0.98, C_sv_ratio_hap2=0.05)
    assert H.sv_depth_features(raw3)["c_sv_hap_depth_change_A"] == 0.05
    assert H.sv_depth_features({}) == {}


def test_sv_depth_features_cope_with_missing_haplotypes():
    # chrX in a female child: the father contributes one X, so a parental haplotype ratio is absent
    raw = dict(C_sv_ratio_hap1=0.06, C_sv_ratio_hap2=None, C_sv_ratio_all=0.55,
               F_sv_ratio_hap1=1.0, F_sv_ratio_hap2=None, M_sv_ratio_hap1=0.98, M_sv_ratio_hap2=1.01,
               C_sv_junc_hap1=0, C_sv_junc_hap2=0, C_sv_junc_hapu=5)
    f = H.sv_depth_features(raw)
    assert f["c_sv_hap_depth_change_A"] == 0.06 and "c_sv_hap_depth_change_O" not in f
    assert f["p_sv_max_hap_depth_change"] == 0.02
    assert "c_sv_junction_hap_concentration" not in f      # no junction read carries a haplotype tag
    assert f["c_sv_junc_untagged_frac"] == 1.0


def test_deletion_is_classified_on_depth_not_on_alt_read_confinement():
    """P8 amendment: the cohort's four prioritised de novo SVs, with their MEASURED interval evidence."""
    from phase_dnm.evidence import hapmatrix as H
    cp = H.ClassParams()
    # MECP2, 35 kb, chrX: depth 0.59 of flank, complete loss of heterozygosity over 47 sites, neither parent depleted.
    # Under alt-read confinement this was `inconclusive` (5 junction reads, 3 untagged, 2 on opposite haplotypes).
    mecp2 = dict(sv_depth_ratio_inside_flank=0.5872, sv_het_snv_persistence=0.0, c_sv_junc_both_ends=1,
                 F_sv_ratio_all=1.2, M_sv_ratio_all=1.0191, C_sv_tagged_loss=0.3089)
    cls, flags = H.classify_sv_interval(mecp2, {"parent_of_origin": "undetermined"}, cp)
    assert cls == "germline_DNM_unphased" and "LOSS_OF_HETEROZYGOSITY" in flags
    assert H._rule_score_sv(mecp2, cls, cp) == 6
    assert H.classify_sv_interval(mecp2, {"parent_of_origin": "paternal"}, cp)[0] == "germline_DNM_phased"
    # CELSR1, 38 kb: the CHILD is not depleted and both PARENTS are -> inherited, not de novo
    celsr = dict(sv_depth_ratio_inside_flank=0.9291, sv_het_snv_persistence=0.8231, c_sv_junc_both_ends=1,
                 F_sv_ratio_all=0.6958, M_sv_ratio_all=0.5935)
    assert H.classify_sv_interval(celsr, {}, cp)[0] == "inherited_missed_in_parent"
    # FBRSL1: no depth loss, junction at one end only -> the depth evidence does not decide; alt-read rules apply
    assert H.classify_sv_interval(dict(sv_depth_ratio_inside_flank=1.1003, sv_het_snv_persistence=0.3491,
                                       c_sv_junc_both_ends=0, F_sv_ratio_all=1.1231, M_sv_ratio_all=0.9956), {}, cp) is None
    # DNMT3A, 302 bp: far too small for depth to move; the junction reads carry it through the alt-read rules
    assert H.classify_sv_interval({}, {}, cp) is None


def test_depth_loss_with_persisting_hets_is_not_constitutional():
    """Depth falls but heterozygous sites persist inside: both haplotypes are present, so mosaic or chimeric."""
    from phase_dnm.evidence import hapmatrix as H
    cp = H.ClassParams()
    sv = dict(sv_depth_ratio_inside_flank=0.55, sv_het_snv_persistence=0.9, c_sv_junc_both_ends=1,
              F_sv_ratio_all=1.0, M_sv_ratio_all=1.0)
    cls, flags = H.classify_sv_interval(sv, {}, cp)
    assert cls == "child_postzygotic_mosaic" and "HET_SNVS_PERSIST_INSIDE" in flags
    # and a deletion without junction reads at both ends is not called on depth alone
    sv2 = dict(sv, sv_het_snv_persistence=0.0, c_sv_junc_both_ends=0)
    assert H.classify_sv_interval(sv2, {}, cp)[0] == "inconclusive"


def test_rule_led_tier1_for_deletions_the_classifier_cannot_score():
    """P15 amendment 2026-09-15: the SV classifier has 43 positives between 10 and 50 kb, so its score is uninformative
    for large deletions — it ranked the pathogenic MECP2 deletion at 0.66 while the deterministic evidence scored 6 of 6."""
    from phase_dnm import integrate as I
    P = I.FinalParams(tau={"SV": 0.99}, tau_tier2={"SV": 0.97})
    mecp2 = dict(phase_class="germline_DNM_unphased", variant_class="SV", rule_score="6",
                 flags="LOW_HAP_DEPTH;AMBIGUOUS_READS;CHILD_BLOCK_UNORIENTED;SV_DEPTH_EVIDENCE",
                 gnomad_af="-1", cohort_AC_loo="0", pon_founder_recurrence_loo="0", segdup_overlap="0")
    d = I.decide(mecp2, P, 0.66)
    assert d["dnm_call"] == "YES" and d["dnm_tier"] == 1 and d["decision_reason"] == "TIER1_SV_DEPTH"
    # the path is gated: the rule layer still applies, the evidence must be complete, and it must BE depth evidence
    assert I.decide(dict(mecp2, cohort_AC_loo="9"), P, 0.66)["decision_reason"] == "RULES:cohort"
    assert I.decide(dict(mecp2, rule_score="3"), P, 0.66)["decision_reason"] == "BELOW_TAU"
    assert I.decide(dict(mecp2, flags="LOW_HAP_DEPTH"), P, 0.66)["decision_reason"] == "BELOW_TAU"
    assert I.decide(dict(mecp2, phase_class="inconclusive"), P, 0.66)["dnm_call"] == "NO"
    # it does not apply to the other classes, and a high-scoring SV is still called on its score
    assert I.decide(dict(mecp2, variant_class="SNV"), I.FinalParams(tau={"SNV": 0.99}), 0.66)["decision_reason"] == "BELOW_TAU"
    dnmt3a = dict(phase_class="germline_DNM_phased", variant_class="SV", rule_score="6", flags="NEAR_CHILD_SWITCH",
                  gnomad_af="-1", cohort_AC_loo="0", pon_founder_recurrence_loo="0", segdup_overlap="0")
    assert I.decide(dnmt3a, P, 0.998)["decision_reason"] == "TIER1"
    # and it can be switched off
    assert I.decide(mecp2, I.FinalParams(tau={"SV": 0.99}, sv_depth_rule_tier1=False), 0.66)["decision_reason"] == "BELOW_TAU"
