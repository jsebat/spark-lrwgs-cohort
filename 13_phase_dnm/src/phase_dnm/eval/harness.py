"""ONE harness for every arm (DESIGN P22, P25): per class group and fold seed, the same rows, the same swap-closed
family folds, the same metrics for

  RF          the phase-aware XGBoost classifier (nested CV, calibrated)            + ablations no_phase / phase_only
  RF+P        RF re-ranked by the P15 phase layer
  H*          the original pipeline's heuristic arms as sweep scores (eval/heuristics.py), with their operating point
  H*+P        the heuristic sweep re-ranked by the phase layer
  baselines   sklearn RandomForest and logistic regression on identical folds

On SYNTHETIC labels the phase layer can only DEMOTE (a synthetic trio has no transmitted haplotype, so
`germline_DNM_phased` — the rescue condition — never occurs for a positive); the rescue branch is evaluated on the
external truth sets (spike-ins, WES-confirmed exonic calls), never on synthetic labels. Every table says which.

Outputs per class group: harness.<cls>.tsv (arm x seed x metrics), cv_report.<cls>.json, rf_probs/<child>.<cls>.rf_probs.tsv
(every real candidate scored by the model of the fold that held its family out; seed-averaged), tau.<cls>.json
(tau / tau_rescue at fixed held-out precision, provisional until the external-truth FDR is available), and the frozen
per-class model (refit on all families) with its training manifest (P21).
"""
from __future__ import annotations

import csv
import glob
import json
import os
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from . import heuristics as Hx
from ..train import nested_cv as CV

CAND_COLS = ["variant_id", "caller_gt", "father_gt", "mother_gt", "caller_gq", "father_gq", "mother_gq", "caller_dp", "father_dp", "mother_dp",
             "child_ad", "father_ad", "mother_ad", "class_payload"]
EV_COLS = ["variant_id", "phase_class", "hap_obs_k5", "phase_score"]
ANNOT_COLS = ["variant_id", "segdup_overlap", "gnomad_af", "pon_founder_recurrence_loo", "sib_shared", "mask_frac", "trgt_founder_max_loo"]
DEMOTE_CLASSES = ("phase_conflict_artifact", "inherited_missed_in_parent", "child_postzygotic_mosaic", "parental_mosaic_transmitted")


# ----------------------------------------------------------------------------------------------
# the P15 layer as a monotone re-ranking (pure)
# ----------------------------------------------------------------------------------------------
def phase_rerank(score: np.ndarray, phase_class: Sequence[str], hap_obs: Sequence, phase_score: Sequence, allow_rescue: bool,
                 min_score: float = 0.9) -> np.ndarray:
    """Demoted rows (phase-conflict / inherited-missed / mosaic) drop below every other row, keeping their order; with
    allow_rescue, phased-germline rows with six haplotypes observed and posterior >= min_score rise above every other row."""
    s = np.asarray(score, dtype=float).copy()
    pc = np.asarray(phase_class, dtype=object)
    demote = np.isin(pc, list(DEMOTE_CLASSES))
    s[demote] = s[demote] - 1e6
    if allow_rescue:
        ho = pd.to_numeric(pd.Series(hap_obs), errors="coerce").to_numpy()
        ps = pd.to_numeric(pd.Series(phase_score), errors="coerce").to_numpy()
        resc = (pc == "germline_DNM_phased") & (ho >= 6) & (ps >= min_score) & ~demote
        s[resc] = s[resc] + 1e6
    return s


def operating_point(y: np.ndarray, pass_: np.ndarray) -> Dict[str, float]:
    y = np.asarray(y).astype(bool); p = np.asarray(pass_).astype(bool)
    tp = int((y & p).sum()); fp = int((~y & p).sum()); fn = int((y & ~p).sum()); tn = int((~y & ~p).sum())
    return dict(tpr=tp / max(1, tp + fn), fpr=fp / max(1, fp + tn), precision=tp / max(1, tp + fp), n_pass=int(p.sum()))


def choose_tau(y: np.ndarray, p: np.ndarray, precision: float) -> Optional[float]:
    """Smallest probability threshold whose held-out precision on these labels is >= `precision`."""
    y = np.asarray(y); p = np.asarray(p, dtype=float)
    ok = ~np.isnan(p); y, p = y[ok], p[ok]
    if y.sum() == 0:
        return None
    order = np.argsort(-p, kind="mergesort")
    ys, ps = y[order], p[order]
    tp = np.cumsum(ys); prec = tp / np.arange(1, len(ys) + 1)
    good = np.where(prec >= precision)[0]
    if len(good) == 0:
        return None
    return float(ps[good.max()])


# ----------------------------------------------------------------------------------------------
# loaders
# ----------------------------------------------------------------------------------------------
def _read_cols(path: str, cols: List[str]) -> pd.DataFrame:
    head = pd.read_csv(path, sep="\t", nrows=0).columns
    use = [c for c in cols if c in head]
    df = pd.read_csv(path, sep="\t", usecols=use, dtype=str, keep_default_na=False)
    for c in cols:
        if c not in df.columns:
            df[c] = ""
    return df


def side_tables(ids: pd.DataFrame, evidence_dir: str, train_dir: str, class_group: str, seed: int, family_of_child: Dict[str, str],
                synth_dir_of: Dict[str, str]) -> pd.DataFrame:
    """Candidate fields, phase columns and annotation columns for every (sample_id, variant_id, origin) row."""
    parts = []
    for (sid, origin), sub in ids.groupby(["sample_id", "origin"], sort=False):
        if origin == "real":
            base_c = os.path.join(evidence_dir, family_of_child[sid], "candidates", "%s.%s.candidates.tsv" % (sid, class_group))
            base_e = os.path.join(evidence_dir, family_of_child[sid], "evidence", "%s.%s.evidence.lik.tsv" % (sid, class_group))
            base_f = os.path.join(evidence_dir, family_of_child[sid], "features", "%s.%s.features.tsv" % (sid, class_group))
        else:
            d = synth_dir_of[(sid, seed)]
            base_c = os.path.join(d, "candidates", "%s.%s.candidates.tsv" % (sid, class_group))
            base_e = os.path.join(d, "evidence", "%s.%s.evidence.tsv" % (sid, class_group))
            base_f = os.path.join(d, "features", "%s.%s.features.tsv" % (sid, class_group))
        c = _read_cols(base_c, CAND_COLS) if os.path.exists(base_c) else pd.DataFrame(columns=CAND_COLS)
        e = _read_cols(base_e, EV_COLS) if os.path.exists(base_e) else pd.DataFrame(columns=EV_COLS)
        f = _read_cols(base_f, ANNOT_COLS) if os.path.exists(base_f) else pd.DataFrame(columns=ANNOT_COLS)
        m = sub[["sample_id", "variant_id", "origin"]].merge(c, on="variant_id", how="left").merge(e, on="variant_id", how="left").merge(f, on="variant_id", how="left")
        parts.append(m)
    out = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    return out.fillna("")


def heuristic_scores(side: pd.DataFrame, class_group: str) -> pd.DataFrame:
    rows = []
    for r in side.to_dict("records"):
        pl = r.get("class_payload") or "{}"
        try:
            r["class_payload"] = json.loads(pl) if isinstance(pl, str) else pl
        except ValueError:
            r["class_payload"] = {}
        rows.append(Hx.score_row(class_group, r, r))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------------------------
# the run
# ----------------------------------------------------------------------------------------------
def run_class(class_group: str, evidence_dir: str, train_dir: str, folds_dir: str, seeds: Sequence[int], out_dir: str,
              family_of_child: Dict[str, str], max_real_per_child: Optional[int], registry_manifest: dict, log=None,
              baselines: bool = True, freeze: bool = True) -> Dict[str, object]:
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(os.path.join(out_dir, "rf_probs"), exist_ok=True)
    log = log or (lambda m: None)
    table: List[dict] = []
    report: Dict[str, object] = {"class_group": class_group, "seeds": list(seeds), "arms": {}}
    prob_acc: Dict[Tuple[str, str], List[float]] = defaultdict(list)
    chosen_params: List[dict] = []
    # synthetic trio directories per (child, seed)
    synth_dir_of: Dict[Tuple[str, int], str] = {}
    for seed in seeds:
        t = os.path.join(folds_dir, "synthetic_trios.seed%d.tsv" % seed)
        if os.path.exists(t):
            for r in csv.DictReader(open(t), delimiter="\t"):
                synth_dir_of[(r["child"], seed)] = os.path.join(train_dir, "seed%d" % seed, r["synthetic_id"])
    for seed in seeds:
        assign = {r["family_id"]: int(r["outer_fold"]) for r in csv.DictReader(open(os.path.join(folds_dir, "folds.seed%d.tsv" % seed)), delimiter="\t")}
        real_glob = os.path.join(evidence_dir, "*", "features", "*.%s.features.rf.tsv" % class_group)
        synth_glob = os.path.join(train_dir, "seed%d" % seed, "*", "features", "*.%s.features.rf.tsv" % class_group)
        d = CV.load_matrices(real_glob, synth_glob, class_group, family_of_child, max_real_per_child=max_real_per_child, seed=seed)
        log("harness %s seed %d: %d real + %d synthetic rows, %d features" % (class_group, seed, d.n_real, d.n_synth, len(d.features)))
        side = side_tables(d.ids, evidence_dir, train_dir, class_group, seed, family_of_child, synth_dir_of)
        side = d.ids[["sample_id", "variant_id", "origin"]].merge(side, on=["sample_id", "variant_id", "origin"], how="left").fillna("")
        pc, ho, ps = side["phase_class"].to_numpy(), side["hap_obs_k5"].to_numpy(), side["phase_score"].to_numpy()

        def add(arm: str, score: np.ndarray, prob: Optional[np.ndarray] = None, pass_: Optional[np.ndarray] = None, note: str = ""):
            ok = ~np.isnan(np.asarray(score, dtype=float))
            row = dict(class_group=class_group, arm=arm, seed=seed, n=int(ok.sum()), n_pos=int(d.y[ok].sum()),
                       roc_auc=CV.roc_auc(d.y[ok], score[ok]), pr_auc=CV.pr_auc(d.y[ok], score[ok]),
                       brier=CV.brier(d.y[ok], prob[ok]) if prob is not None else None, note=note)
            if pass_ is not None:
                row.update({"op_" + k: v for k, v in operating_point(d.y[ok], pass_[ok]).items()})
            table.append(row)
            log("  %-22s auc=%s pr=%s%s" % (arm, None if row["roc_auc"] is None else round(row["roc_auc"], 4), None if row["pr_auc"] is None else round(row["pr_auc"], 4),
                                            "" if pass_ is None else " op(tpr=%.3f fpr=%.4f)" % (row["op_tpr"], row["op_fpr"])))

        # classifier arms
        res, prob, raw, fold_models = CV.nested_cv(d, assign, class_group, seed, kind="xgb", ablation="full", log=log)
        add("RF", raw, prob)
        add("RF+P(demote)", phase_rerank(raw, pc, ho, ps, allow_rescue=False), note="synthetic labels: demotion only")
        chosen_params += res.chosen
        report["arms"].setdefault("RF", []).append(res.__dict__)
        for abl in ("no_phase", "phase_only"):
            r2, p2, raw2, _ = CV.nested_cv(d, assign, class_group, seed, kind="xgb", ablation=abl, calibrate=False, log=log)
            add("RF_" + abl, raw2)
            report["arms"].setdefault("RF_" + abl, []).append(dict(auc=r2.auc, pr=r2.pr, n_rows=r2.n_rows))
        if baselines:
            for kind in ("rf", "lr"):
                r3, p3, raw3, _ = CV.nested_cv(d, assign, class_group, seed, kind=kind, calibrate=False, log=log)
                add("baseline_" + kind, raw3)
        # heuristic arms on the same rows
        hs = heuristic_scores(side, class_group)
        for arm in sorted({c[:-6] for c in hs.columns if c.endswith("_score")}):
            sc = hs[arm + "_score"].to_numpy(dtype=float); pa = hs[arm + "_pass"].to_numpy()
            add(arm, sc, pass_=pa)
            add(arm + "+P(demote)", phase_rerank(sc, pc, ho, ps, allow_rescue=False), pass_=pa.astype(bool) & ~np.isin(pc, list(DEMOTE_CLASSES)), note="synthetic labels: demotion only")
        # rf_prob for EVERY real candidate of the held-out families (full per-child matrices), seed-averaged later
        for path in sorted(glob.glob(real_glob)):
            df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
            if df.empty:
                continue
            fam = family_of_child.get(df["sample_id"].iloc[0])
            fm = fold_models.get(assign.get(fam, -1))
            if fm is None:
                continue
            p_cal, _ = fm.predict(df)
            for vid, sid, pv in zip(df["variant_id"], df["sample_id"], p_cal):
                prob_acc[(sid, vid)].append(float(pv))
        # tau on held-out synthetic-vs-real labels (provisional, P15)
        okp = ~np.isnan(prob)
        report.setdefault("tau_by_seed", []).append(dict(seed=seed, tau_p95=choose_tau(d.y[okp], prob[okp], 0.95), tau_p80=choose_tau(d.y[okp], prob[okp], 0.80)))
    # write outputs
    with open(os.path.join(out_dir, "harness.%s.tsv" % class_group), "w", newline="") as fh:
        cols = sorted({k for r in table for k in r}, key=lambda c: (c not in ("class_group", "arm", "seed"), c))
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in table:
            w.writerow({c: ("" if r.get(c) is None else r.get(c)) for c in cols})
    by_child: Dict[str, List[Tuple[str, float, int]]] = defaultdict(list)
    for (sid, vid), ps_ in prob_acc.items():
        by_child[sid].append((vid, float(np.mean(ps_)), len(ps_)))
    for sid, rows in by_child.items():
        with open(os.path.join(out_dir, "rf_probs", "%s.%s.rf_probs.tsv" % (sid, class_group)), "w", newline="") as fh:
            fh.write("variant_id\trf_prob\tn_seeds\n")
            for vid, pv, n in rows:
                fh.write("%s\t%.6f\t%d\n" % (vid, pv, n))
    taus = [t for t in report.get("tau_by_seed", []) if t["tau_p95"] is not None]
    tau = dict(class_group=class_group, tau=float(np.median([t["tau_p95"] for t in taus])) if taus else None,
               tau_rescue=float(np.median([t["tau_p80"] for t in taus if t["tau_p80"] is not None])) if taus else None,
               basis="held-out synthetic-vs-real precision 0.95 / 0.80, median over seeds; PROVISIONAL until the external-truth FDR (P15)")
    with open(os.path.join(out_dir, "tau.%s.json" % class_group), "w") as fh:
        json.dump(tau, fh, indent=1)
    report["tau"] = tau
    report["summary"] = _summarise(table)
    if freeze and chosen_params:
        # the most often chosen hyperparameters across folds/seeds
        keyf = lambda c: json.dumps(c["params"], sort_keys=True)
        best = max({keyf(c) for c in chosen_params}, key=lambda k: sum(1 for c in chosen_params if keyf(c) == k))
        d_all = CV.load_matrices(os.path.join(evidence_dir, "*", "features", "*.%s.features.rf.tsv" % class_group),
                                 os.path.join(train_dir, "seed*", "*", "features", "*.%s.features.rf.tsv" % class_group),
                                 class_group, family_of_child, max_real_per_child=max_real_per_child, seed=0)
        path = CV.fit_final(d_all, class_group, json.loads(best), 0, os.path.join(out_dir, "models"), registry_manifest,
                            notes=dict(seeds=list(seeds), max_real_per_child=max_real_per_child, tau=tau))
        report["frozen_model"] = path
        log("frozen model -> %s" % path)
    with open(os.path.join(out_dir, "cv_report.%s.json" % class_group), "w") as fh:
        json.dump(report, fh, indent=1, sort_keys=True, default=str)
    return report


def _summarise(table: List[dict]) -> Dict[str, dict]:
    out: Dict[str, dict] = {}
    by_arm: Dict[str, List[dict]] = defaultdict(list)
    for r in table:
        by_arm[r["arm"]].append(r)
    for arm, rows in by_arm.items():
        aucs = [r["roc_auc"] for r in rows if r["roc_auc"] is not None]
        prs = [r["pr_auc"] for r in rows if r["pr_auc"] is not None]
        out[arm] = dict(n_seeds=len(rows), roc_auc_mean=float(np.mean(aucs)) if aucs else None, roc_auc_min=min(aucs) if aucs else None,
                        roc_auc_max=max(aucs) if aucs else None, pr_auc_mean=float(np.mean(prs)) if prs else None,
                        op_tpr_mean=float(np.mean([r["op_tpr"] for r in rows])) if "op_tpr" in rows[0] else None,
                        op_fpr_mean=float(np.mean([r["op_fpr"] for r in rows])) if "op_fpr" in rows[0] else None)
    return out
