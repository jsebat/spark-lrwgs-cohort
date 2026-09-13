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
