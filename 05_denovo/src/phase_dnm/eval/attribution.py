"""P14 attribution — how much of each decision is the phase block? TreeSHAP (XGBoost `pred_contribs`, exact, no extra
dependency) on the saved fold models, on each fold's HELD-OUT real rows and synthetic rows, |SHAP| summed by feature
family: A caller/genotype, B context, C reads (child read quality), D phase (the rf_safe part of the six-haplotype block).
Reported cohort-wide, per class, and for the rows the decision layer called YES via the classifier path. Two caveats
travel with every number (DESIGN P14): the classifier's phase share UNDERCOUNTS phase (the transmitted-haplotype test is
rf_safe: false and lives in the rule layer), and a high phase share on synthetic labels can be leakage - the external-truth
arm (P27) is the check.
"""
from __future__ import annotations

import csv
import glob
import json
import os
from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from ..train import nested_cv as CV

PHASE = ("c_", "p_", "hap_obs")
CONTEXT = ("gc_", "homopolymer", "seq_entropy", "segdup", "repeat_class", "mappab", "in_str", "indel_len", "dist_nearest", "cpg", "tr_", "sv", "delta_", "motif", "n_motif")
READS = ("child_alt_", "child_ref_", "c_alt_mapq", "c_alt_nm", "c_ref_nm", "c_alt_clip", "c_alt_supp", "c_alt_readlen", "c_alt_rq", "c_alt_mapq0")


def family_of(col: str) -> str:
    if col.startswith(READS):
        return "C_reads"
    if col.startswith(PHASE):
        return "D_phase"
    if col.startswith(CONTEXT):
        return "B_context"
    return "A_caller"


def shap_by_family(fm: CV.FoldModel, X: pd.DataFrame) -> Tuple[Dict[str, float], Dict[str, float]]:
    """(mean |SHAP| share per family, mean |SHAP| per feature) for the rows of X under the fold model."""
    import xgboost as xgb
    Xn = X.reindex(columns=fm.cols).apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    booster = fm.model.get_booster()
    contrib = booster.predict(xgb.DMatrix(Xn, feature_names=None), pred_contribs=True)   # (n, n_features + 1 bias)
    a = np.abs(contrib[:, :-1]).mean(axis=0)
    tot = a.sum() or 1.0
    per_feat = {c: float(v) for c, v in zip(fm.cols, a)}
    fam: Dict[str, float] = defaultdict(float)
    for c, v in per_feat.items():
        fam[family_of(c)] += float(v) / tot
    return dict(fam), per_feat


def run(class_group: str, harness_dir: str, evidence_dir: str, train_dir: str, folds_dir: str, seed: int, family_of_child: Dict[str, str],
        rf_q_threshold: float = 0.999, max_rows_per_child: int = 5000, log=None) -> Dict[str, object]:
    log = log or (lambda m: None)
    fold_of = {r["family_id"]: int(r["outer_fold"]) for r in csv.DictReader(open(os.path.join(folds_dir, "folds.seed%d.tsv" % seed)), delimiter="\t")}
    syn_dir_of: Dict[str, str] = {}
    t = os.path.join(folds_dir, "synthetic_trios.seed%d.tsv" % seed)
    for r in csv.DictReader(open(t), delimiter="\t"):
        syn_dir_of[r["child"]] = os.path.join(train_dir, "seed%d" % seed, r["synthetic_id"])
    rng = np.random.default_rng(seed)
    acc = {"real": defaultdict(list), "synthetic": defaultdict(list), "real_called": defaultdict(list)}
    per_feat_acc: Dict[str, List[float]] = defaultdict(list)
    n = {"real": 0, "synthetic": 0, "real_called": 0}
    for fam, k in sorted(fold_of.items()):
        prefix = os.path.join(harness_dir, "fold_models", "%s.seed%d.fold%d" % (class_group, seed, k))
        if not os.path.exists(prefix + ".meta.json"):
            continue
        fm = CV.FoldModel.load(prefix)
        for sid, f in family_of_child.items():
            if f != fam:
                continue
            real = glob.glob(os.path.join(evidence_dir, fam, "features", "%s.%s.features.rf.tsv" % (sid, class_group)))
            if real:
                X = pd.read_csv(real[0], sep="\t", dtype=str, keep_default_na=False)
                if len(X) > max_rows_per_child:
                    X = X.iloc[np.sort(rng.choice(len(X), max_rows_per_child, replace=False))]
                if len(X):
                    fam_share, pf = shap_by_family(fm, X)
                    for kk, v in fam_share.items():
                        acc["real"][kk].append(v)
                    for kk, v in pf.items():
                        per_feat_acc[kk].append(v)
                    n["real"] += len(X)
                # rows the classifier path called: rf_q >= threshold from rf_probs (if present)
                rp = os.path.join(harness_dir, "rf_probs", "%s.%s.rf_probs.tsv" % (sid, class_group))
                if os.path.exists(rp):
                    q = pd.read_csv(rp, sep="\t", dtype=str, keep_default_na=False)
                    if "rf_q" in q.columns:
                        hi = set(q.loc[pd.to_numeric(q["rf_q"], errors="coerce") >= rf_q_threshold, "variant_id"])
                        Xa = pd.read_csv(real[0], sep="\t", dtype=str, keep_default_na=False)
                        Xa = Xa[Xa["variant_id"].isin(hi)]
                        if len(Xa):
                            fs, _ = shap_by_family(fm, Xa)
                            for kk, v in fs.items():
                                acc["real_called"][kk].append(v)
                            n["real_called"] += len(Xa)
            sd = syn_dir_of.get(sid)
            if sd:
                sp = glob.glob(os.path.join(sd, "features", "%s.%s.features.rf.tsv" % (sid, class_group)))
                if sp:
                    Xs = pd.read_csv(sp[0], sep="\t", dtype=str, keep_default_na=False)
                    if len(Xs):
                        fs, _ = shap_by_family(fm, Xs)
                        for kk, v in fs.items():
                            acc["synthetic"][kk].append(v)
                        n["synthetic"] += len(Xs)
        log("attribution %s: fold %d done" % (class_group, k))
    out: Dict[str, object] = {"class_group": class_group, "seed": seed, "n_rows": n, "share_by_family": {}, "top_features": {}}
    for kind, d in acc.items():
        out["share_by_family"][kind] = {k: round(float(np.mean(v)), 4) for k, v in sorted(d.items())}
    top = sorted(((k, float(np.mean(v))) for k, v in per_feat_acc.items()), key=lambda kv: -kv[1])[:15]
    tot = sum(v for _, v in top) or 1.0
    out["top_features"] = {k: round(v, 5) for k, v in top}
    out["caveats"] = ["classifier phase share undercounts phase: the transmitted-haplotype test is rf_safe:false (rule layer)",
                      "a high phase share on synthetic labels can be leakage; the external-truth arm (P27) is the check"]
    return out
