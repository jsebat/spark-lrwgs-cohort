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
