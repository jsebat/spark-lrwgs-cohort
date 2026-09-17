"""Apply a FROZEN classifier (P21) to a cohort: rf_prob from the model, rf_q against the scored cohort's own candidates.

This is the transfer path for a cohort that did not train the model: the three frozen models (models/<class>.xgb.json +
<class>.training_manifest.json) score every child's features.rf.tsv, and rf_q is the fraction of the scored cohort's real
candidates OF THE SAME VARIANT CLASS that score below the row - the same pass-rate semantics the thresholds were set on
(thresholds.yaml final.tau_q / tau_q_tier2), so `phase-dnm integrate` runs unchanged on the output. Calibration is not
needed for the tiers (rf_q is a ranking); rf_prob is the raw model probability. Caveat (P21): a frozen model transfers only to
the same callers and versions at similar coverage; otherwise retrain (swap -> synthetic -> train -> rescore).
"""
from __future__ import annotations

import glob
import hashlib
import json
import os
from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from . import nested_cv as CV
from .rescore import ecdf_ref, quantile


def load_frozen(model_path: str, manifest_path: Optional[str] = None):
    """(XGBClassifier, feature columns, manifest). Verifies the manifest's sha256 when a manifest is given."""
    import xgboost as xgb
    manifest_path = manifest_path or model_path.replace(".xgb.json", ".training_manifest.json")
    man = json.load(open(manifest_path)) if os.path.exists(manifest_path) else {}
    if man.get("model_sha256"):
        got = hashlib.sha256(open(model_path, "rb").read()).hexdigest()
        if got != man["model_sha256"]:
            raise ValueError("frozen model %s: sha256 %s != manifest %s" % (model_path, got[:16], man["model_sha256"][:16]))
    m = xgb.XGBClassifier()
    m.load_model(model_path)
    cols = list(man.get("features") or [])
    if not cols:
        cols = list(m.get_booster().feature_names or [])
    return m, cols, man


def score_class(class_group: str, model_path: str, evidence_dir: str, out_dir: str, children: Sequence[Tuple[str, str]],
                manifest_path: Optional[str] = None, log=None) -> Dict[str, object]:
    """children: (family_id, sample_id) pairs to score. Writes out_dir/rf_probs/<child>.<class>.rf_probs.tsv (rf_prob, rf_q,
    n_seeds=1) and out_dir/tau.<class>.json (score_column rf_q)."""
    log = log or (lambda m: None)
    m, cols, man = load_frozen(model_path, manifest_path)
    per_child: Dict[str, Tuple[pd.DataFrame, np.ndarray, np.ndarray]] = {}
    parts: Dict[str, List[np.ndarray]] = defaultdict(list)
    for fam, sid in children:
        paths = glob.glob(os.path.join(evidence_dir, fam, "features", "%s.%s.features.rf.tsv" % (sid, class_group)))
        if not paths:
            log("score %s: no rf matrix for %s" % (class_group, sid)); continue
        # _read_matrix, not read_csv: derived features (qd_child) are computed AT LOAD from columns the matrix
        # already holds, so a raw read leaves the column absent and the reindex below silently makes it NaN.
        # The model was trained with it, so every production score was computed with that feature missing
        # (train/score skew, found 2026-09-15).
        df = CV._read_matrix(paths[0])
        if df.empty:
            continue
        X = df.reindex(columns=cols).apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        p = m.predict_proba(X)[:, 1]
        vc = df["variant_class"].astype(str).values if "variant_class" in df.columns else np.array(["ALL"] * len(p))
        per_child[sid] = (df, p, vc)
        for c in np.unique(vc):
            parts[c].append(p[vc == c])
    refs = {c: ecdf_ref(np.concatenate(v)) for c, v in parts.items()}
    os.makedirs(os.path.join(out_dir, "rf_probs"), exist_ok=True)
    n_rows = 0
    for sid, (df, p, vc) in per_child.items():
        q = np.full(len(p), np.nan)
        for c, ref in refs.items():
            mk = vc == c
            if mk.any():
                q[mk] = quantile(ref, p[mk])
        with open(os.path.join(out_dir, "rf_probs", "%s.%s.rf_probs.tsv" % (sid, class_group)), "w", newline="") as fh:
            fh.write("variant_id\trf_prob\trf_q\tn_seeds\n")
            for vid, pv, qv in zip(df["variant_id"], p, q):
                fh.write("%s\t%.6f\t%.6f\t1\n" % (vid, pv, qv))
        n_rows += len(p)
    tau = dict(class_group=class_group, score_column="rf_q", tau_q=0.999, tau_q_rescue=0.99, frozen_model=os.path.abspath(model_path),
               model_sha256=man.get("model_sha256"), n_features=len(cols), reference={c: int(len(r)) for c, r in refs.items()},
               basis_q="rf_q = 1 - pass rate among the SCORED cohort's candidates of the same variant class under the frozen model; "
                       "operating points from thresholds.yaml final.tau_q / tau_q_tier2 (P21 transfer path)")
    with open(os.path.join(out_dir, "tau.%s.json" % class_group), "w") as fh:
        json.dump(tau, fh, indent=1)
    log("score %s: %d children, %d rows, reference %s" % (class_group, len(per_child), n_rows, tau["reference"]))
    return {"class_group": class_group, "children": len(per_child), "rows": n_rows, "reference": tau["reference"]}
