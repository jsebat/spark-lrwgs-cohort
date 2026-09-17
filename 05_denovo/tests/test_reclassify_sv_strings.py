"""reclassify_row must survive an SV evidence row exactly as csv.DictReader delivers it: every value a string.

Regression for the 2026-09-16 review finding: classify_sv_interval and _rule_score_sv compare the sv_ columns with
floats, and the strings raised TypeError. The existing reclassify test built rows without sv_ columns and could not see
it; the one cohort run that reclassified SV tables happened while those columns were blank.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from phase_dnm.evidence import hapmatrix as H  # noqa: E402


def _row():
    r = {"chrom": "chr9", "start": "1000", "child_hap1_is": ".", "F_transmitted_hap": ".", "M_transmitted_hap": ".",
         "n_reads_used": "40", "flags": "", "C_untagged_dp": "6", "C_untagged_alt": "2", "F_untagged_dp": "2",
         "F_untagged_alt": "0", "M_untagged_dp": "2", "M_untagged_alt": "0"}
    for pre, alt, dp in (("C1", "3", "8"), ("C2", "0", "9"), ("F1", "0", "10"), ("F2", "0", "9"), ("M1", "0", "11"), ("M2", "0", "10")):
        r.update({pre + "_dp": dp, pre + "_alt": alt, pre + "_ref": str(int(dp) - int(alt)), pre + "_amb": "0"})
    # the deletion-path evidence, as strings
    r.update({"sv_interval_len": "35454", "sv_depth_ratio_inside_flank": "0.5872", "sv_het_snv_persistence": "0.0",
              "C_sv_junc_both_ends": "1", "c_sv_junc_both_ends": "1", "F_sv_ratio_all": "1.02", "M_sv_ratio_all": "0.98",
              "C_sv_tagged_loss": "0.3", "C_sv_ratio_hap1": "0.1", "C_sv_ratio_hap2": "1.0"})
    return r


def test_reclassify_sv_row_with_string_values_uses_the_depth_path():
    out = H.reclassify_row(_row(), H.HapParams(), H.ClassParams())
    assert out["phase_class"].startswith("germline_DNM"), out["phase_class"]
    assert "SV_DEPTH_EVIDENCE" in out["flags"]
    assert int(out["rule_score"]) == 6


def test_reclassify_sv_row_inherited_when_a_parent_is_depleted():
    r = _row(); r["F_sv_ratio_all"] = "0.48"
    out = H.reclassify_row(r, H.HapParams(), H.ClassParams())
    assert out["phase_class"] == "inherited_missed_in_parent"


if __name__ == "__main__":
    test_reclassify_sv_row_with_string_values_uses_the_depth_path()
    test_reclassify_sv_row_inherited_when_a_parent_is_depleted()
    print("ok")
