"""Fold-quantile score (M4, 2026-09-13). A calibrated probability from one fold model is not on the same scale as the
next fold's: with one pooled tau, fold 0 passed 0.46 % of its held-out real candidates while the other folds passed
<= 0.04 % (per-fold AUCs were all ~0.997, so the RANKING was fine, the SCALE was not). The decision layer therefore uses

    rf_q = 1 - (fraction of the fold's held-out REAL candidates with a probability >= this row's probability)

computed per fold model from its own held-out real rows (label-free, so no leak), averaged over seeds. tau_q is then a
pass rate on real candidates by definition (P15: 0.1 % -> tau_q 0.999; rescue 1 % -> 0.99), identical across folds.
This module recomputes rf_q post hoc from the saved fold models, so a finished training run need not be repeated.
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

from . import nested_cv as CV


def ecdf_ref(scores: np.ndarray) -> np.ndarray:
    """Sorted reference scores of a fold's held-out real rows."""
    s = np.asarray(scores, dtype=float)
    return np.sort(s[~np.isnan(s)])


def quantile(ref: np.ndarray, p: np.ndarray) -> np.ndarray:
    """rf_q = fraction of reference scores strictly below p (1 - pass rate at threshold p)."""
    if len(ref) == 0:
        return np.full(len(p), np.nan)
    return np.searchsorted(ref, np.asarray(p, dtype=float), side="left") / len(ref)


def rescore_class(class_group: str, harness_dir: str, evidence_dir: str, folds_dir: str, seeds: Sequence[int],
                  family_of_child: Dict[str, str], max_ref_per_child: int = 20000, log=None) -> Dict[str, object]:
    """For every seed and fold: load the saved fold model, build the reference from that fold's held-out families' real
    rows (thinned per child, seeded), then score EVERY real row of those families -> rf_q. Average rf_q (and rf_prob)
    over seeds and write rf_probs/<child>.<class>.rf_probs.tsv with both columns; write tau.<class>.json with tau_q."""
    log = log or (lambda m: None)
    acc_q: Dict[Tuple[str, str], List[float]] = defaultdict(list)
    acc_p: Dict[Tuple[str, str], List[float]] = defaultdict(list)
    children_by_fam: Dict[str, List[str]] = defaultdict(list)
    for sid, fam in family_of_child.items():
        children_by_fam[fam].append(sid)
    n_models = 0
    for seed in seeds:
        fold_of = {r["family_id"]: int(r["outer_fold"]) for r in csv.DictReader(open(os.path.join(folds_dir, "folds.seed%d.tsv" % seed)), delimiter="\t")}
        rng = np.random.default_rng(seed)
        for k in sorted(set(fold_of.values())):
            prefix = os.path.join(harness_dir, "fold_models", "%s.seed%d.fold%d" % (class_group, seed, k))
            if not os.path.exists(prefix + ".meta.json"):
                log("rescore %s: no fold model for seed %d fold %d" % (class_group, seed, k)); continue
            fm = CV.FoldModel.load(prefix)
            n_models += 1
            fams = [f for f, kk in fold_of.items() if kk == k]
            per_child: Dict[str, Tuple[pd.DataFrame, np.ndarray, np.ndarray]] = {}
            ref_parts: Dict[str, List[np.ndarray]] = defaultdict(list)     # reference per VARIANT CLASS (SNV / INDEL separately; JS 2026-09-14)
            for fam in fams:
                for sid in children_by_fam.get(fam, []):
                    paths = glob.glob(os.path.join(evidence_dir, fam, "features", "%s.%s.features.rf.tsv" % (sid, class_group)))
                    if not paths:
                        continue
                    # see score.py: derived features must be added at load or the fold model scores them as NaN
                    df = CV._read_matrix(paths[0])
                    if df.empty:
                        continue
                    p, _ = fm.predict(df)
                    vc = df["variant_class"].astype(str).values if "variant_class" in df.columns else np.array(["ALL"] * len(p))
                    per_child[sid] = (df, p, vc)
                    for c in np.unique(vc):
                        pc = p[vc == c]
                        idx = np.arange(len(pc)) if len(pc) <= max_ref_per_child else np.sort(rng.choice(len(pc), max_ref_per_child, replace=False))
                        ref_parts[c].append(pc[idx])
            refs = {c: ecdf_ref(np.concatenate(parts)) for c, parts in ref_parts.items()}
            for sid, (df, p, vc) in per_child.items():
                q = np.full(len(p), np.nan)
                for c, ref in refs.items():
                    m = vc == c
                    if m.any():
                        q[m] = quantile(ref, p[m])
                for vid, pv, qv in zip(df["variant_id"], p, q):
                    acc_p[(sid, vid)].append(float(pv)); acc_q[(sid, vid)].append(float(qv))
            log("rescore %s seed %d fold %d: %d children, reference real rows %s" % (class_group, seed, k, len(per_child), {c: len(v) for c, v in refs.items()}))
    out_dir = os.path.join(harness_dir, "rf_probs")
    os.makedirs(out_dir, exist_ok=True)
    by_child: Dict[str, List[Tuple[str, float, float, int]]] = defaultdict(list)
    for key in acc_q:
        by_child[key[0]].append((key[1], float(np.mean(acc_p[key])), float(np.mean(acc_q[key])), len(acc_q[key])))
    for sid, rows in by_child.items():
        with open(os.path.join(out_dir, "%s.%s.rf_probs.tsv" % (sid, class_group)), "w", newline="") as fh:
            fh.write("variant_id\trf_prob\trf_q\tn_seeds\n")
            for vid, pv, qv, n in rows:
                fh.write("%s\t%.6f\t%.6f\t%d\n" % (vid, pv, qv, n))
    tau_path = os.path.join(harness_dir, "tau.%s.json" % class_group)
    tau = json.load(open(tau_path)) if os.path.exists(tau_path) else {"class_group": class_group}
    tau.update(score_column="rf_q", tau_q=0.999, tau_q_rescue=0.99,
               basis_q="rf_q = 1 - pass rate among the fold's held-out real candidates OF THE SAME VARIANT CLASS (per fold model, seed-averaged); "
                       "thresholds.yaml final.tau_q / tau_q_tier2 set the operating points (0.3.0: tier 1 / tier 2 + rule layer)")
    with open(tau_path, "w") as fh:
        json.dump(tau, fh, indent=1)
    return {"class_group": class_group, "fold_models": n_models, "children": len(by_child), "rows": len(acc_q)}
