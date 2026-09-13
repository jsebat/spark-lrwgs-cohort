"""M4 — swap-closed, family-grouped nested cross-validation per class (DESIGN P12, P14, P21, P24).

Rows: the rf_safe matrices (`*.features.rf.tsv`) of real trios (label 0, family = the child's family) and of synthetic
trios (label 1, family = the child's family - the same outer fold as the surrogate parents by construction). The
registry gate has already refused any transmission / parent-of-origin column (P11); this module refuses them again
by name before fitting.

Per outer fold and seed: inner folds over the training families choose the XGBoost hyperparameters (small grid,
PR-AUC), the model is refit on all training families, isotonic calibration is fitted on the inner out-of-sample
predictions, and the held-out families are scored once. Pooled held-out scores -> ROC-AUC, PR-AUC, Brier, reliability
bins; per-row `rf_prob` for every real candidate (the M3 input). Baselines on identical folds: RF (sklearn), logistic
regression, no-phase (D-block columns dropped) and phase-only ablations. The three shipped models are refit on all
families with the inner-selected hyperparameters and frozen with a manifest.

Everything here is plain numpy + xgboost/sklearn; matrices are read with pandas. Identifiers are values, never code.
"""
from __future__ import annotations

import glob
import hashlib
import json
import math
import os
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

ID_COLS = ["family_id", "sample_id", "variant_id", "variant_class"]
UNSAFE_PREFIXES = ("parent_of_origin", "poo_", "t_alt", "t_dp", "u_alt", "u_dp", "nt_parent", "t_hap", "t_tr", "u_tr", "child_hap1_is", "F_transmitted", "M_transmitted")
PHASE_PREFIXES = ("c_", "p_", "hap_obs")            # the D block: what "phase information" means in the ablation (P24)

GRID = [dict(max_depth=4, learning_rate=0.1, n_estimators=300, min_child_weight=5, subsample=0.8, colsample_bytree=0.8),
        dict(max_depth=6, learning_rate=0.05, n_estimators=600, min_child_weight=10, subsample=0.8, colsample_bytree=0.8),
        dict(max_depth=3, learning_rate=0.2, n_estimators=200, min_child_weight=20, subsample=0.9, colsample_bytree=0.9)]


@dataclass
class Data:
    X: pd.DataFrame
    y: np.ndarray
    family: np.ndarray
    ids: pd.DataFrame
    features: List[str]
    n_real: int
    n_synth: int


def _read_matrix(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    return df


def load_matrices(real_glob: str, synth_glob: str, class_group: str, family_of_child: Dict[str, str],
                  max_real_per_child: Optional[int] = None, seed: int = 0, allowed: Optional[Iterable[str]] = None,
                  max_presence_gap: float = 0.5, min_real_presence: float = 0.01) -> Data:
    """Concatenate real (label 0) and synthetic (label 1) rf matrices for one class group. Real rows may be thinned per
    child (seeded) for speed; the fold assignment uses the CHILD's family for both kinds of rows."""
    rng = np.random.default_rng(seed)
    frames = []
    n_real = n_synth = 0
    for path in sorted(glob.glob(real_glob)):
        df = _read_matrix(path)
        if df.empty:
            continue
        if max_real_per_child and len(df) > max_real_per_child:
            df = df.iloc[np.sort(rng.choice(len(df), max_real_per_child, replace=False))]
        df["label"] = 0; df["origin"] = "real"
        n_real += len(df); frames.append(df)
    for path in sorted(glob.glob(synth_glob)):
        df = _read_matrix(path)
        if df.empty:
            continue
        df["label"] = 1; df["origin"] = "synthetic"
        n_synth += len(df); frames.append(df)
    if not frames:
        raise ValueError("no matrices for %s" % class_group)
    all_ = pd.concat(frames, ignore_index=True)
    # TR candidate ids were the TRID until 2026-09-13 (two outlier alleles of one locus share an id): keep one row per key
    all_ = all_.drop_duplicates(subset=["sample_id", "variant_id", "origin"], keep="first").reset_index(drop=True)
    feats = [c for c in all_.columns if c not in ID_COLS + ["label", "origin"]]
    if allowed is not None:
        allowed = set(allowed)
        feats = [c for c in feats if c in allowed]        # the registry's CURRENT rf_safe set (a matrix may predate a registry change)
    bad = [c for c in feats if c.startswith(UNSAFE_PREFIXES)]
    if bad:
        raise ValueError("rf_safe violation in matrices: %s" % bad)
    X = all_[feats].apply(pd.to_numeric, errors="coerce")
    # presence leak guard (2026-09-13): a column that is filled for one label class and empty for the other encodes the
    # label through its missingness (e.g. read-quality summaries that exist only in reviews run after a code change).
    # Drop any feature whose non-missing fraction differs between classes by more than `max_presence_gap`, and any
    # feature that is essentially never present in the real rows; report both.
    lab = all_["label"].to_numpy().astype(int)
    present = X.notna()
    p_real = present[lab == 0].mean() if (lab == 0).any() else present.mean()
    p_syn = present[lab == 1].mean() if (lab == 1).any() else present.mean()
    dropped = {c: (round(float(p_real[c]), 3), round(float(p_syn[c]), 3)) for c in feats
               if abs(float(p_real[c]) - float(p_syn[c])) > max_presence_gap or float(p_real[c]) < min_real_presence}
    if dropped:
        feats = [c for c in feats if c not in dropped]
        X = X[feats]
    load_matrices.last_dropped = dropped  # type: ignore[attr-defined]
    fam = all_["sample_id"].map(lambda s: family_of_child.get(s, "?")).to_numpy()
    if (fam == "?").any():
        raise ValueError("%d rows whose child is not in the manifest" % int((fam == "?").sum()))
    return Data(X=X, y=all_["label"].to_numpy().astype(int), family=fam, ids=all_[ID_COLS + ["origin"]].reset_index(drop=True),
                features=feats, n_real=n_real, n_synth=n_synth)


# ----------------------------------------------------------------------------------------------
# metrics (no sklearn dependency for the basics, so a report can be recomputed anywhere)
# ----------------------------------------------------------------------------------------------
def roc_auc(y: np.ndarray, s: np.ndarray) -> Optional[float]:
    y = np.asarray(y); s = np.asarray(s, dtype=float)
    ok = ~np.isnan(s)
    y, s = y[ok], s[ok]
    n1, n0 = int(y.sum()), int((1 - y).sum())
    if n1 == 0 or n0 == 0:
        return None
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(len(s), dtype=float)
    ss = s[order]
    i = 0
    while i < len(ss):
        j = i
        while j + 1 < len(ss) and ss[j + 1] == ss[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2 + 1
        i = j + 1
    return float((ranks[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def pr_auc(y: np.ndarray, s: np.ndarray) -> Optional[float]:
    y = np.asarray(y); s = np.asarray(s, dtype=float)
    ok = ~np.isnan(s)
    y, s = y[ok], s[ok]
    if y.sum() == 0:
        return None
    order = np.argsort(-s, kind="mergesort")
    y = y[order]
    tp = np.cumsum(y); fp = np.cumsum(1 - y)
    prec = tp / (tp + fp); rec = tp / y.sum()
    # average precision: sum over positives of precision at each positive
    return float(np.sum(prec[y == 1]) / y.sum())


def brier(y: np.ndarray, p: np.ndarray) -> Optional[float]:
    p = np.asarray(p, dtype=float); ok = ~np.isnan(p)
    return float(np.mean((np.asarray(y)[ok] - p[ok]) ** 2)) if ok.any() else None


def reliability(y: np.ndarray, p: np.ndarray, bins: int = 10) -> List[dict]:
    p = np.asarray(p, dtype=float); y = np.asarray(y)
    out = []
    edges = np.linspace(0, 1, bins + 1)
    for i in range(bins):
        m = (p >= edges[i]) & (p < edges[i + 1] if i < bins - 1 else p <= edges[i + 1])
        if m.sum():
            out.append(dict(lo=float(edges[i]), hi=float(edges[i + 1]), n=int(m.sum()), mean_pred=float(p[m].mean()), frac_pos=float(y[m].mean())))
    return out


# ----------------------------------------------------------------------------------------------
# models
# ----------------------------------------------------------------------------------------------
def make_model(kind: str, params: Optional[dict] = None, seed: int = 0, spw: float = 1.0):
    if kind == "xgb":
        import xgboost as xgb
        p = dict(GRID[0]); p.update(params or {})
        return xgb.XGBClassifier(objective="binary:logistic", eval_metric="aucpr", random_state=seed, n_jobs=4, tree_method="hist",
                                 scale_pos_weight=spw, **p)
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.impute import SimpleImputer
    from sklearn.preprocessing import StandardScaler
    if kind == "rf":
        return make_pipeline(SimpleImputer(strategy="median"), RandomForestClassifier(n_estimators=400, min_samples_leaf=5, n_jobs=4, random_state=seed, class_weight="balanced_subsample"))
    if kind == "lr":
        return make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), LogisticRegression(max_iter=2000, class_weight="balanced", C=0.5))
    raise ValueError(kind)


def _fit_predict(kind, params, Xtr, ytr, Xte, seed):
    spw = float((ytr == 0).sum() / max(1, (ytr == 1).sum()))
    m = make_model(kind, params, seed, spw)
    m.fit(Xtr, ytr)
    return m, m.predict_proba(Xte)[:, 1]


def _split(fams: np.ndarray, assign: Dict[str, int], k: int) -> Tuple[np.ndarray, np.ndarray]:
    f = np.array([assign[x] for x in fams])
    return f != k, f == k


def _folds_from_families(families: Sequence[str], n: int, seed: int) -> Dict[str, int]:
    from .folds import Family, outer_folds
    return outer_folds([Family(f, ("x",)) for f in sorted(set(families))], n_folds=n, seed=seed, min_fold_size=1)


@dataclass
class CVResult:
    class_group: str
    seed: int
    kind: str
    ablation: str
    auc: Optional[float]
    pr: Optional[float]
    brier: Optional[float]
    n_rows: int
    n_pos: int
    chosen: List[dict] = field(default_factory=list)
    reliability: List[dict] = field(default_factory=list)
    per_fold_auc: List[Optional[float]] = field(default_factory=list)


def ablation_columns(features: Sequence[str], ablation: str) -> List[str]:
    if ablation == "no_phase":
        return [c for c in features if not c.startswith(PHASE_PREFIXES)]
    if ablation == "phase_only":
        return [c for c in features if c.startswith(PHASE_PREFIXES)]
    return list(features)


class FoldModel:
    """A fitted fold model with its calibrator and column list: scores any matrix with the same registry columns."""

    def __init__(self, model, iso, cols: List[str]):
        self.model, self.iso, self.cols = model, iso, cols

    def predict(self, X: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
        Xn = X.reindex(columns=self.cols)          # a per-child matrix may lack (or add) registry columns: missing -> NaN
        raw = self.model.predict_proba(Xn.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float))[:, 1]
        return (self.iso.predict(raw) if self.iso is not None else raw), raw


def nested_cv(d: Data, assign: Dict[str, int], class_group: str, seed: int, kind: str = "xgb", ablation: str = "full",
              inner_folds: int = 3, calibrate: bool = True, log=None) -> Tuple[CVResult, np.ndarray, np.ndarray, Dict[int, FoldModel]]:
    """Returns (result, held-out calibrated probabilities aligned with d rows, raw scores, fold -> FoldModel)."""
    cols = ablation_columns(d.features, ablation)
    X = d.X[cols].to_numpy(dtype=float)
    prob = np.full(len(d.y), np.nan); raw = np.full(len(d.y), np.nan)
    chosen = []; per_fold = []
    fold_models: Dict[int, FoldModel] = {}
    for k in sorted(set(assign.values())):
        tr, te = _split(d.family, assign, k)
        if te.sum() == 0 or d.y[tr].sum() == 0 or d.y[te].sum() == 0:
            per_fold.append(None); continue
        best, best_pr, inner_oos = None, -1.0, None
        if kind == "xgb":
            inner = _folds_from_families(d.family[tr], inner_folds, seed + 17 * k)
            for params in GRID:
                oos = np.full(int(tr.sum()), np.nan)
                Xtr_all, ytr_all, ftr_all = X[tr], d.y[tr], d.family[tr]
                for kk in sorted(set(inner.values())):
                    itr, ite = _split(ftr_all, inner, kk)
                    if ytr_all[itr].sum() == 0 or ite.sum() == 0:
                        continue
                    _, p = _fit_predict(kind, params, Xtr_all[itr], ytr_all[itr], Xtr_all[ite], seed)
                    oos[ite] = p
                score = pr_auc(ytr_all[~np.isnan(oos)], oos[~np.isnan(oos)]) or -1.0
                if score > best_pr:
                    best, best_pr, inner_oos = params, score, oos
        m, p_te = _fit_predict(kind, best, X[tr], d.y[tr], X[te], seed)
        raw[te] = p_te
        iso = None
        if calibrate and inner_oos is not None and (~np.isnan(inner_oos)).sum() > 50:
            from sklearn.isotonic import IsotonicRegression
            ok = ~np.isnan(inner_oos)
            iso = IsotonicRegression(out_of_bounds="clip").fit(inner_oos[ok], d.y[tr][ok])
            prob[te] = iso.predict(p_te)
        else:
            prob[te] = p_te
        fold_models[k] = FoldModel(m, iso, cols)
        per_fold.append(roc_auc(d.y[te], p_te))
        chosen.append(dict(fold=k, params=best, inner_pr_auc=round(best_pr, 4) if best_pr >= 0 else None, n_train=int(tr.sum()), n_test=int(te.sum())))
        if log:
            log("cv %s seed %d %s/%s fold %d: n_train=%d n_test=%d auc=%s" % (class_group, seed, kind, ablation, k, tr.sum(), te.sum(), per_fold[-1]))
    ok = ~np.isnan(raw)
    res = CVResult(class_group=class_group, seed=seed, kind=kind, ablation=ablation, auc=roc_auc(d.y[ok], raw[ok]), pr=pr_auc(d.y[ok], raw[ok]),
                   brier=brier(d.y[ok], prob[ok]), n_rows=int(ok.sum()), n_pos=int(d.y[ok].sum()), chosen=chosen,
                   reliability=reliability(d.y[ok], prob[ok]), per_fold_auc=per_fold)
    return res, prob, raw, fold_models


def fit_final(d: Data, class_group: str, params: dict, seed: int, out_dir: str, registry_manifest: dict, notes: dict) -> str:
    """Refit on ALL families with the chosen hyperparameters; freeze model + calibration + manifest (P21)."""
    import xgboost as xgb
    os.makedirs(out_dir, exist_ok=True)
    X = d.X[d.features].to_numpy(dtype=float)
    spw = float((d.y == 0).sum() / max(1, (d.y == 1).sum()))
    m = make_model("xgb", params, seed, spw)
    m.fit(X, d.y)
    path = os.path.join(out_dir, "%s.xgb.json" % class_group)
    m.save_model(path)
    man = dict(class_group=class_group, features=d.features, n_features=len(d.features), n_rows=int(len(d.y)), n_pos=int(d.y.sum()),
               params=params, seed=seed, xgboost_version=xgb.__version__, registry=registry_manifest, notes=notes,
               model_sha256=hashlib.sha256(open(path, "rb").read()).hexdigest())
    with open(os.path.join(out_dir, "%s.training_manifest.json" % class_group), "w") as fh:
        json.dump(man, fh, indent=1, sort_keys=True)
    return path


def attribution(d: Data, prob_model, seed: int = 0, n_perm: int = 3) -> Dict[str, object]:
    """Grouped permutation importance by feature family (A caller / B context / C reads / D phase) on the given rows, and
    TreeSHAP share by family when the model is XGBoost (P14 attribution plan)."""
    X = d.X[d.features].to_numpy(dtype=float)
    base = roc_auc(d.y, prob_model(X)) or 0.0
    groups = {"D_phase": [i for i, c in enumerate(d.features) if c.startswith(PHASE_PREFIXES)],
              "A_caller": [i for i, c in enumerate(d.features) if c.startswith(("child_", "min_", "max_", "AR_", "AB_", "indel_", "hap")) and not c.startswith(PHASE_PREFIXES)],
              "B_context": [i for i, c in enumerate(d.features) if c.startswith(("gc_", "homopolymer", "seq_", "segdup", "repeat", "mappab", "in_str", "tr_", "sv_", "svlen", "svtype", "delta_"))]}
    groups["other"] = [i for i in range(len(d.features)) if not any(i in g for g in groups.values())]
    rng = np.random.default_rng(seed)
    out: Dict[str, object] = {"base_auc": base, "grouped_permutation_auc_drop": {}}
    for g, idx in groups.items():
        if not idx:
            continue
        drops = []
        for _ in range(n_perm):
            Xp = X.copy()
            perm = rng.permutation(len(Xp))
            Xp[:, idx] = Xp[perm][:, idx]
            drops.append(base - (roc_auc(d.y, prob_model(Xp)) or 0.0))
        out["grouped_permutation_auc_drop"][g] = dict(n_features=len(idx), mean_drop=float(np.mean(drops)), range=[float(min(drops)), float(max(drops))])
    return out
