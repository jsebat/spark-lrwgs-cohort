"""The heuristic arms must reproduce the pipeline's decisions and order rows by the knob with failures below passes."""
from phase_dnm.eval import heuristics as Hx


def rec(**kw):
    base = dict(caller_gt="0/1", father_gt="0/0", mother_gt="0/0", caller_gq=40, father_gq=45, mother_gq=50, caller_dp=30, father_dp=28,
                mother_dp=31, child_ad="15,15", father_ad="28,0", mother_ad="31,0", class_payload={"allele_index": 1})
    base.update(kw)
    return base


def test_slivar_h1_verbatim():
    p, s = Hx.slivar_h1(rec())
    assert p and s == 40                                            # knob = min GQ of the trio
    assert not Hx.slivar_h1(rec(caller_gq=19))[0]                    # GQ below 20
    assert Hx.slivar_h1(rec(caller_gq=19))[1] == 19                  # still ordered by the knob among fixed-criteria passes
    assert not Hx.slivar_h1(rec(father_ad="28,1"))[0]                # one parental alt read fails
    assert Hx.slivar_h1(rec(father_ad="28,1"))[1] < Hx.FAIL / 2      # and scores below every passing row
    assert not Hx.slivar_h1(rec(child_ad="26,4"))[0]                 # AB 0.13 < 0.25
    assert not Hx.slivar_h1(rec(father_dp=9))[0]                     # DP < 10
    assert not Hx.slivar_h1(rec(father_gt="0/1"))[0]
    assert not Hx.slivar_h1(rec(caller_gt="1/1"))[0]
    # DP falls back to sum(AD) when the caller gives no DP
    assert Hx.slivar_h1(rec(caller_dp=None, father_dp=None, mother_dp=None))[0]
    # multi-allelic: the allele index selects the AD column
    assert Hx.slivar_h1(rec(caller_gt="0/2", child_ad="15,0,15", father_ad="28,0,0", mother_ad="31,0,0", class_payload={"allele_index": 2}))[0]


def test_snv_cohort_arms():
    r = rec()
    assert Hx.snv_h2(r, {"segdup_overlap": "0", "gnomad_af": ""})[0]              # absent from gnomAD counts as rare
    assert not Hx.snv_h2(r, {"segdup_overlap": "1", "gnomad_af": ""})[0]
    assert not Hx.snv_h2(r, {"segdup_overlap": "0", "gnomad_af": "0.01"})[0]
    assert Hx.snv_h3(r, {"segdup_overlap": "0", "gnomad_af": "", "pon_founder_recurrence_loo": "2", "sib_shared": "0"})[0]
    assert not Hx.snv_h3(r, {"segdup_overlap": "0", "gnomad_af": "", "pon_founder_recurrence_loo": "3", "sib_shared": "0"})[0]
    assert not Hx.snv_h3(r, {"segdup_overlap": "0", "gnomad_af": "", "pon_founder_recurrence_loo": "0", "sib_shared": "1"})[0]
    # ordering: a masked row scores below every unmasked passing row whatever its GQ
    assert Hx.snv_h2(rec(caller_gq=99), {"segdup_overlap": "1"})[1] < Hx.snv_h2(rec(caller_gq=21), {"segdup_overlap": "0"})[1]


def test_sv_arms():
    p, s = Hx.sv_h1(dict(caller_gt="0/1", father_gt="0/0", mother_gt="0/0", caller_gq=12))
    assert p and s == 12
    assert Hx.sv_h1(dict(caller_gt="1/1", father_gt="0/0", mother_gt="0/0", caller_gq=None))[0]        # no quality filter
    assert not Hx.sv_h1(dict(caller_gt="0/1", father_gt="0/1", mother_gt="0/0", caller_gq=99))[0]
    r = dict(caller_gt="0/1", father_gt="0/0", mother_gt="0/0", caller_gq=30)
    assert Hx.sv_h2(r, {"mask_frac": "0.2", "pon_founder_recurrence_loo": "1"})[0]
    assert not Hx.sv_h2(r, {"mask_frac": "0.6", "pon_founder_recurrence_loo": "1"})[0]
    assert not Hx.sv_h2(r, {"mask_frac": "0.0", "pon_founder_recurrence_loo": "6"})[0]


def test_tr_arms():
    pl = {"child_AL": [30, 120], "father_AL": [30, 33], "mother_AL": [28, 30], "motif_unit_bp": 3, "outlier_allele_idx": 1,
          "child_SD": "12,9", "father_SD": "10,11", "mother_SD": "8,9"}
    r = dict(class_payload=pl)
    p, s = Hx.tr_h1(r)
    assert p and s == 87                                                  # 120 - 33
    assert not Hx.tr_h1(dict(class_payload=dict(pl, child_AL=[30, 60])))[0]    # +27 bp < 50
    assert not Hx.tr_h1(dict(class_payload=dict(pl, child_AL=[30, 83], father_AL=[30, 60])))[0]   # +23 and < 1.5x
    p, s = Hx.tr_h2(r, {"trgt_founder_max_loo": "100"})
    assert p and s == 87
    assert not Hx.tr_h2(r, {"trgt_founder_max_loo": "125"})[0]           # not above the founder maximum
    assert not Hx.tr_h2(dict(class_payload=dict(pl, father_SD="4,11")), {})[0]     # a parental allele with SD < 5
    small = dict(class_payload=dict(pl, child_AL=[30, 40]))               # +7 bp < max(6, 3*3) = 9
    assert not Hx.tr_h2(small, {})[0] and Hx.tr_h2(small, {})[1] == 7
    assert Hx.tr_h2(dict(class_payload=dict(pl, child_AL=[30, 43])), {})[0]        # +10 >= 9


def test_score_row_names():
    out = Hx.score_row("snv_indel", rec(), {"segdup_overlap": "0"})
    assert set(out) == {"H1_slivar_pass", "H1_slivar_score", "H2_hiconf_pass", "H2_hiconf_score", "H3_cohort_pass", "H3_cohort_score"}
    assert out["H1_slivar_pass"] == 1 and out["H3_cohort_pass"] == 1
