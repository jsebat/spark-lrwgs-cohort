"""Pure parts of the M4 harness: metrics, the phase re-ranking, tau selection, and a tiny end-to-end nested CV on
synthetic matrices (xgboost / sklearn are optional locally: the CV test is skipped when they are absent)."""
import numpy as np
import pandas as pd
import pytest

from phase_dnm.eval import harness as HZ
from phase_dnm.train import nested_cv as CV


def test_roc_and_pr_auc_basic():
    y = np.array([0, 0, 1, 1]); s = np.array([0.1, 0.4, 0.35, 0.8])
    assert abs(CV.roc_auc(y, s) - 0.75) < 1e-9
    assert CV.roc_auc(np.array([0, 0]), s[:2]) is None
    assert abs(CV.pr_auc(y, s) - (1.0 + 2 / 3) / 2) < 1e-9        # AP: precision at each positive (1.0 at rank 1, 2/3 at rank 3)
    assert abs(CV.roc_auc(np.array([0, 1, 1, 0]), np.array([1, 1, 1, 1])) - 0.5) < 1e-9     # ties -> 0.5
    assert abs(CV.brier(y, np.array([0, 0, 1, 1.0])) - 0.0) < 1e-9


def test_phase_rerank_demotes_and_rescues():
    s = np.array([0.9, 0.8, 0.7, 0.6])
    pc = ["germline_DNM_phased", "phase_conflict_artifact", "inconclusive", "germline_DNM_phased"]
    ho = ["6", "6", "6", "5"]; ps = ["0.95", "0.99", "0.2", "0.99"]
    r = HZ.phase_rerank(s, pc, ho, ps, allow_rescue=False)
    assert r[1] < min(r[0], r[2], r[3]) and r[0] > r[2] > r[3]                 # demotion only: order otherwise intact
    r2 = HZ.phase_rerank(s, pc, ho, ps, allow_rescue=True)
    assert r2[0] > r2[2] > r2[3] > r2[1] and r2[0] > 1e5                        # row 0 rescued above everything; row 3 (5 haplotypes) not


def test_operating_point_and_tau():
    y = np.array([1, 1, 0, 0, 0]); p = np.array([1, 0, 1, 0, 0])
    op = HZ.operating_point(y, p)
    assert op["tpr"] == 0.5 and abs(op["fpr"] - 1 / 3) < 1e-9 and op["precision"] == 0.5 and op["n_pass"] == 2
    y = np.array([1, 1, 1, 0, 1, 0, 0, 0]); prob = np.array([0.95, 0.9, 0.85, 0.8, 0.7, 0.3, 0.2, 0.1])
    assert HZ.choose_tau(y, prob, 1.0) == 0.85          # smallest threshold with precision 1.0
    assert HZ.choose_tau(y, prob, 0.8) == 0.7           # top 5: 4 positives -> 0.8
    assert HZ.choose_tau(np.zeros(3), prob[:3], 0.5) is None


def test_nested_cv_end_to_end_small():
    try:
        import xgboost  # noqa: F401
        import sklearn.isotonic  # noqa: F401
        import sklearn.ensemble  # noqa: F401
    except Exception as e:  # a blocked DLL raises ImportError after a partial import; importorskip alone does not catch it
        pytest.skip("ML backend unavailable here: %s" % e)
    rng = np.random.default_rng(0)
    n = 600
    fam = np.array(["f%02d" % (i % 12) for i in range(n)])
    y = (np.arange(n) % 3 == 0).astype(int)
    x1 = y * 2 + rng.normal(size=n); x2 = rng.normal(size=n); x3 = y + rng.normal(size=n) * 2
    X = pd.DataFrame({"c_alt_hap_frac": x1, "child_GQ": x2, "p_max_alt_any_hap": x3})
    ids = pd.DataFrame({"family_id": fam, "sample_id": fam, "variant_id": [str(i) for i in range(n)], "variant_class": "SNV", "origin": "real"})
    d = CV.Data(X=X, y=y, family=fam, ids=ids, features=list(X.columns), n_real=n, n_synth=0)
    assign = {"f%02d" % i: i % 4 for i in range(12)}
    res, prob, raw, fm = CV.nested_cv(d, assign, "snv_indel", seed=0, kind="xgb", inner_folds=2)
    assert res.auc is not None and res.auc > 0.8 and res.n_rows == n and len(fm) == 4
    assert np.isnan(prob).sum() == 0
    res2, _, raw2, _ = CV.nested_cv(d, assign, "snv_indel", seed=0, kind="xgb", ablation="no_phase", inner_folds=2, calibrate=False)
    assert res2.auc < res.auc                                                     # the signal lives in the phase columns here
    p_cal, p_raw = fm[0].predict(X)
    assert len(p_cal) == n and (p_cal >= 0).all() and (p_cal <= 1).all()


def test_choose_tau_fpr_on_real_rows():
    y = np.array([1] * 6 + [0] * 10); p = np.concatenate([np.full(6, 0.99), np.linspace(0.0, 0.9, 10)])
    assert HZ.choose_tau_fpr(y, p, 0.0) == 0.9           # no real row above threshold
    assert HZ.choose_tau_fpr(y, p, 0.1) == 0.8           # one of ten real rows may pass
    assert HZ.choose_tau_fpr(np.ones(3), p[:3], 0.1) is None


def test_presence_leak_guard_drops_label_encoding_columns(tmp_path):
    import csv
    # real matrix: column q empty; synthetic matrix: q filled -> q must be dropped, x kept
    cols = ["family_id", "sample_id", "variant_id", "variant_class", "x", "q"]
    def write(path, sid, fam, fill_q):
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t"); w.writerow(cols)
            for i in range(50):
                w.writerow([fam, sid, "v%d" % i, "SNV", str(i % 7), str(0.5) if fill_q else ""])
    (tmp_path / "real").mkdir(); (tmp_path / "syn").mkdir()
    write(str(tmp_path / "real" / "a.rf.tsv"), "kidA", "famA", False)
    write(str(tmp_path / "syn" / "b.rf.tsv"), "kidB", "famB", True)
    d = CV.load_matrices(str(tmp_path / "real" / "*.rf.tsv"), str(tmp_path / "syn" / "*.rf.tsv"), "snv_indel", {"kidA": "famA", "kidB": "famB"})
    assert d.features == ["x"] and "q" in CV.load_matrices.last_dropped
    assert CV.load_matrices.last_dropped["q"] == (0.0, 1.0)


def test_fold_quantile_score_is_a_pass_rate():
    from phase_dnm.train import rescore as RS
    ref = RS.ecdf_ref(np.array([0.1, 0.2, 0.3, 0.4, np.nan, 0.9]))
    assert len(ref) == 5
    q = RS.quantile(ref, np.array([0.95, 0.4, 0.05, 0.25]))
    assert list(q) == [1.0, 0.6, 0.0, 0.4]                      # 0.95 above all 5 -> pass rate 0 -> q 1.0; 0.4 -> 3 of 5 below
    assert np.isnan(RS.quantile(np.array([]), np.array([0.5]))).all()


def test_attribution_feature_families():
    from phase_dnm.eval import attribution as AT
    assert AT.family_of("c_alt_hap_frac") == "D_phase" and AT.family_of("hap_obs_k5") == "D_phase" and AT.family_of("p_min_hap_dp") == "D_phase"
    assert AT.family_of("c_alt_mapq_mean") == "C_reads" and AT.family_of("c_alt_rq_mean") == "C_reads"
    assert AT.family_of("child_GQ") == "A_caller" and AT.family_of("min_PL0") == "A_caller"
    assert AT.family_of("gc_200bp") == "B_context" and AT.family_of("segdup_overlap") == "B_context" and AT.family_of("delta_units") == "B_context"


def test_positive_rarity_gate_keeps_only_private_synthetic_rows(tmp_path):
    """P24 correction: SynthDNM keeps a synthetic positive only at AC = 2 (child + the one transmitting real parent).
    Ported as cohort_AC_loo == 0 with a population-frequency ceiling, applied to the POSITIVES only."""
    from phase_dnm.train import nested_cv as CV
    cols = ["family_id", "sample_id", "variant_id", "variant_class", "x"]
    def write(d, name, rows):
        d.mkdir(parents=True, exist_ok=True)
        with open(d / name, "w", newline="") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                fh.write("\t".join(str(v) for v in r) + "\n")
    real = tmp_path / "ev" / "famA" / "features"
    write(real, "kidA.snv_indel.features.rf.tsv", [("famA", "kidA", "v%d" % i, "SNV", i) for i in range(6)])
    syn = tmp_path / "tr" / "seed0" / "syn1" / "features"
    write(syn, "kidA.snv_indel.features.rf.tsv", [("famA", "kidA", "s%d" % i, "SNV", i) for i in range(4)])
    ann = tmp_path / "tr" / "seed0" / "syn1" / "annot"
    ann.mkdir(parents=True, exist_ok=True)
    with open(ann / "kidA.snv_indel.annot.tsv", "w", newline="") as fh:
        fh.write("variant_id\tgnomad_af\tcohort_AC_loo\n")
        fh.write("s0\t-1.0\t0\n")        # private, absent from gnomAD  -> KEEP
        fh.write("s1\t0.0005\t0\n")      # private, rare in gnomAD      -> KEEP
        fh.write("s2\t-1.0\t7\n")        # common in the cohort         -> drop
        fh.write("s3\t0.08\t0\n")        # common in gnomAD             -> drop
    d = CV.load_matrices(str(real / "*.rf.tsv"), str(syn / "*.rf.tsv"), "snv_indel", {"kidA": "famA"})
    assert d.n_synth == 2 and d.n_real == 6
    g = CV.load_matrices.last_positive_gate
    assert g["before"] == 4 and g["after"] == 2 and g["unannotated"] == 0
    d2 = CV.load_matrices(str(real / "*.rf.tsv"), str(syn / "*.rf.tsv"), "snv_indel", {"kidA": "famA"}, rare_positives=False)
    assert d2.n_synth == 4


def test_site_qual_is_not_a_registry_feature():
    """It is a joint-callset site statistic that rises with carrier count; SynthDNM uses QD / AQ instead (P24 correction)."""
    import os, yaml
    reg = yaml.safe_load(open(os.path.join(os.path.dirname(__file__), "..", "config", "features.yaml")))
    items = [it for sec in reg.values() if isinstance(sec, list) for it in sec if isinstance(it, dict)]
    sq = [it for it in items if it.get("name") == "site_qual"]
    assert sq and sq[0]["status"] == "drop" and sq[0]["classes"] == []
