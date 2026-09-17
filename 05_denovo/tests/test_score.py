"""P21 transfer path: the frozen model scores a cohort and rf_q is the per-variant-class pass rate within that cohort."""
import hashlib
import json
import os

import numpy as np
import pandas as pd
import pytest


def _ml():
    try:
        import xgboost  # noqa: F401
    except Exception as e:  # a blocked DLL raises after a partial import
        pytest.skip("ML backend unavailable here: %s" % e)


def test_frozen_manifests_in_repo_match_models():
    root = os.path.join(os.path.dirname(__file__), "..", "models")
    for cls in ("snv_indel", "sv", "tr"):
        man = json.load(open(os.path.join(root, "%s.training_manifest.json" % cls)))
        got = hashlib.sha256(open(os.path.join(root, "%s.xgb.json" % cls), "rb").read()).hexdigest()
        assert got == man["model_sha256"], cls
        assert man["class_group"] == cls and man["n_features"] == len(man["features"]) and man["notes"]["seeds"] == [0, 1, 2, 3, 4]
        assert not any(f.startswith(("t_", "u_", "nt_parent", "parent_of_origin", "gnomad", "cohort_AC", "pon_founder", "lr_sv_catalog")) for f in man["features"])


def test_score_class_writes_rf_q_per_variant_class(tmp_path):
    _ml()
    import xgboost as xgb
    from phase_dnm.train import score as SC
    rng = np.random.default_rng(0)
    cols = ["a", "b"]
    X = rng.normal(size=(400, 2)); y = (X[:, 0] + 0.3 * rng.normal(size=400) > 0).astype(int)
    m = xgb.XGBClassifier(n_estimators=20, max_depth=2, tree_method="hist"); m.fit(X, y)
    mp = str(tmp_path / "snv_indel.xgb.json"); m.save_model(mp)
    json.dump({"class_group": "snv_indel", "features": cols, "model_sha256": hashlib.sha256(open(mp, "rb").read()).hexdigest()},
              open(str(tmp_path / "snv_indel.training_manifest.json"), "w"))
    ev = tmp_path / "ev"
    for fam, sid in (("famA", "kidA"), ("famB", "kidB")):
        d = ev / fam / "features"; d.mkdir(parents=True)
        n = 60
        df = pd.DataFrame({"variant_id": ["%s:%d" % (sid, i) for i in range(n)], "variant_class": ["SNV"] * 40 + ["INDEL"] * 20,
                           "a": rng.normal(size=n), "b": rng.normal(size=n)})
        df.to_csv(d / ("%s.snv_indel.features.rf.tsv" % sid), sep="\t", index=False)
    rep = SC.score_class("snv_indel", mp, str(ev), str(tmp_path / "out"), [("famA", "kidA"), ("famB", "kidB")])
    assert rep["children"] == 2 and rep["rows"] == 120 and rep["reference"] == {"SNV": 80, "INDEL": 40}
    q = pd.read_csv(tmp_path / "out" / "rf_probs" / "kidA.snv_indel.rf_probs.tsv", sep="\t")
    assert list(q.columns) == ["variant_id", "rf_prob", "rf_q", "n_seeds"] and q["rf_q"].between(0, 1).all()
    tau = json.load(open(tmp_path / "out" / "tau.snv_indel.json"))
    assert tau["score_column"] == "rf_q" and tau["model_sha256"] and tau["reference"] == {"SNV": 80, "INDEL": 40}
    # a tampered model is refused
    with open(mp, "ab") as fh:
        fh.write(b" ")
    with pytest.raises(ValueError):
        SC.load_frozen(mp)
