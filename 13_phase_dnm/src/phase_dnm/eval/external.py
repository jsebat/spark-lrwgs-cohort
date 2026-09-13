"""P27 — external-truth arm of the harness: rows whose truth is known by construction or by an orthogonal assay, scored
by the fold model that held the child's family out, with the same arms as the synthetic harness.

Spike-ins (the only arm implemented so far; WES-confirmed exonic calls and GIAB follow the same interface):
  positives  = planted germline (G) rows, all classes, from sim/spike.py's per-child tables
  negatives  = the child's real raw candidates (P13: ~99.7 % non-DNM) plus the planted inherited-missed (IM) rows
  sensitivity = planted mosaics (CM / PM) reported as the fraction above the class tau, never as positives (P10)
Arms: RF (fold model), RF+P with demotion AND rescue (legitimate here: real trios, transmitted haplotypes known), the
heuristic sweeps incl. H2/H3 (population filters are legitimate on these labels), and recall at the class tau.
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

from . import harness as HZ
from . import heuristics as Hx
from ..train import nested_cv as CV
from ..train import rescore as RS

GROUP_OF = {"SNV": "snv_indel", "INDEL": "snv_indel", "SV": "sv", "TR": "tr"}


def load_fold_models(harness_dir: str, class_group: str, seed: int) -> Dict[int, CV.FoldModel]:
    out = {}
    for p in glob.glob(os.path.join(harness_dir, "fold_models", "%s.seed%d.fold*.meta.json" % (class_group, seed))):
        k = int(p.rsplit(".fold", 1)[1].split(".")[0])
        out[k] = CV.FoldModel.load(p[:-len(".meta.json")])
    return out


def spike_rows(spike_dir: str, child: str, class_group: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """(rf matrix rows of the spiked candidates, plan+evidence side table) for one child and class group."""
    rf = os.path.join(spike_dir, "%s.features.rf.tsv" % class_group)
    plan = os.path.join(spike_dir, "plan.tsv")
    ev = os.path.join(spike_dir, "%s.evidence.lik.tsv" % class_group)
    if not (os.path.exists(rf) and os.path.exists(plan) and os.path.exists(ev)):
        return pd.DataFrame(), pd.DataFrame()
    X = pd.read_csv(rf, sep="\t", dtype=str, keep_default_na=False)
    pl = pd.read_csv(plan, sep="\t", dtype=str, keep_default_na=False)[["variant_id", "scenario", "variant_class", "child_frac", "parent_frac"]]
    e = HZ._read_cols(ev, HZ.EV_COLS + HZ.CAND_COLS)
    side = X[["variant_id"]].merge(pl, on="variant_id", how="left").merge(e, on="variant_id", how="left").fillna("")
    return X, side


def evaluate_spikes(evidence_dir: str, harness_dir: str, class_group: str, seed: int, family_of_child: Dict[str, str],
                    fold_of_family: Dict[str, int], taus: Dict[str, float], max_real_per_child: int = 20000, log=None) -> Dict[str, object]:
    log = log or (lambda m: None)
    fms = load_fold_models(harness_dir, class_group, seed)
    if not fms:
        raise ValueError("no fold models for %s seed %d in %s" % (class_group, seed, harness_dir))
    rows: List[dict] = []
    rng = np.random.default_rng(seed)
    for sd in sorted(glob.glob(os.path.join(evidence_dir, "*", "spike", "*"))):
        child = os.path.basename(sd)
        fam = family_of_child.get(child)
        if fam is None or fam not in fold_of_family:
            continue
        fm = fms.get(fold_of_family[fam])
        if fm is None:
            continue
        Xs, side = spike_rows(sd, child, class_group)
        if Xs.empty:
            continue
        p_s, raw_s = fm.predict(Xs)
        for i, r in side.reset_index(drop=True).iterrows():
            rows.append(dict(child=child, fold=fold_of_family[fam], origin="spike", scenario=r["scenario"], variant_class=r["variant_class"], variant_id=r["variant_id"],
                             prob=float(p_s[i]), raw=float(raw_s[i]), phase_class=r["phase_class"], hap_obs_k5=r["hap_obs_k5"], phase_score=r["phase_score"],
                             cand={k: r.get(k, "") for k in HZ.CAND_COLS}))
        # the child's real raw candidates as the negative pool (thinned)
        real_rf = glob.glob(os.path.join(evidence_dir, fam, "features", "%s.%s.features.rf.tsv" % (child, class_group)))
        if real_rf:
            Xr = pd.read_csv(real_rf[0], sep="\t", dtype=str, keep_default_na=False)
            if len(Xr) > max_real_per_child:
                Xr = Xr.iloc[np.sort(rng.choice(len(Xr), max_real_per_child, replace=False))].reset_index(drop=True)
            ev = HZ._read_cols(os.path.join(evidence_dir, fam, "evidence", "%s.%s.evidence.lik.tsv" % (child, class_group)), HZ.EV_COLS + HZ.CAND_COLS)
            sr = Xr[["variant_id"]].merge(ev.drop_duplicates("variant_id"), on="variant_id", how="left").fillna("")
            p_r, raw_r = fm.predict(Xr)
            for i, r in sr.iterrows():
                rows.append(dict(child=child, fold=fold_of_family[fam], origin="real", scenario="", variant_class="", variant_id=r["variant_id"], prob=float(p_r[i]), raw=float(raw_r[i]),
                                 phase_class=r["phase_class"], hap_obs_k5=r["hap_obs_k5"], phase_score=r["phase_score"], cand={k: r.get(k, "") for k in HZ.CAND_COLS}))
        log("external %s: child scored (%d spike rows)" % (class_group, len(Xs)))
    df = pd.DataFrame(rows)
    if df.empty:
        return {"class_group": class_group, "n_rows": 0}
    # fold-quantile score against the fold's REAL rows (the same construction as train/rescore.py)
    df["rf_q"] = np.nan
    for k, g in df.groupby("fold"):
        ref = RS.ecdf_ref(g.loc[g["origin"] == "real", "prob"].to_numpy())
        df.loc[g.index, "rf_q"] = RS.quantile(ref, g["prob"].to_numpy())
    lab = np.where((df["origin"] == "spike") & (df["scenario"] == "G"), 1, np.where((df["origin"] == "real") | (df["scenario"] == "IM"), 0, -1))
    keep = lab >= 0
    y = lab[keep]
    out: Dict[str, object] = {"class_group": class_group, "seed": seed, "n_rows": int(keep.sum()), "n_pos": int(y.sum()),
                              "n_real_neg": int((df["origin"][keep] == "real").sum()), "n_IM_neg": int((df["scenario"][keep] == "IM").sum()), "arms": {}}
    prob = df["prob"].to_numpy()[keep]; raw = df["raw"].to_numpy()[keep]; rfq = df["rf_q"].to_numpy()[keep]
    pc = df["phase_class"].to_numpy()[keep]; ho = df["hap_obs_k5"].to_numpy()[keep]; ps = df["phase_score"].to_numpy()[keep]
    use_q = taus.get("score_column") == "rf_q"
    tau = taus.get("tau_q", 0.999) if use_q else taus.get(class_group)
    score_for_tau = rfq if use_q else prob
    out["decision_score"] = "rf_q" if use_q else "rf_prob"; out["tau"] = tau
    def arm(name, score, pass_=None):
        rec = dict(roc_auc=CV.roc_auc(y, score), pr_auc=CV.pr_auc(y, score))
        if pass_ is not None:
            rec.update(HZ.operating_point(y, pass_))
        out["arms"][name] = rec
        log("  %-22s auc=%s pr=%s%s" % (name, None if rec["roc_auc"] is None else round(rec["roc_auc"], 4), None if rec["pr_auc"] is None else round(rec["pr_auc"], 4),
                                        "" if pass_ is None else " op(tpr=%.3f fpr=%.4f)" % (rec["tpr"], rec["fpr"])))
    arm("RF", raw, (score_for_tau >= tau) if tau is not None else None)
    arm("RF+P(demote)", HZ.phase_rerank(raw, pc, ho, ps, allow_rescue=False))
    arm("RF+P(demote+rescue)", HZ.phase_rerank(raw, pc, ho, ps, allow_rescue=True))
    cand = pd.DataFrame(list(df["cand"][keep]))
    hs = []
    for r in cand.to_dict("records"):
        try:
            r["class_payload"] = json.loads(r.get("class_payload") or "{}")
        except ValueError:
            r["class_payload"] = {}
        hs.append(Hx.score_row(class_group, r, r))
    hs = pd.DataFrame(hs)
    # heuristic arms need caller fields; planted candidates carry none (caller = spike), so an arm whose positives are all
    # "fixed-criterion failures" is NOT evaluable on this truth set and is reported as such (TR H1 works from allele lengths)
    for a in sorted({c[:-6] for c in hs.columns if c.endswith("_score")}):
        sc = hs[a + "_score"].to_numpy(dtype=float)
        if (sc[y == 1] < Hx.FAIL / 2).all():
            out["arms"][a] = dict(roc_auc=None, pr_auc=None, note="not evaluable on spike-ins: planted candidates carry no caller genotype/GQ")
            log("  %-22s not evaluable (no caller fields on planted candidates)" % a)
            continue
        arm(a, sc, hs[a + "_pass"].to_numpy().astype(bool))
    # mosaic sensitivity at tau
    col = "rf_q" if use_q else "prob"
    if tau is not None:
        for sc in ("CM", "PM"):
            m = (df["origin"] == "spike") & (df["scenario"] == sc)
            if m.any():
                out["mosaic_sensitivity_at_tau_" + sc] = dict(n=int(m.sum()), frac_above_tau=float((df[col][m] >= tau).mean()))
    # recall of G at tau per variant class (spike classes are finer than the group); real pass rate alongside
    g = (df["origin"] == "spike") & (df["scenario"] == "G")
    if tau is not None and g.any():
        out["recall_at_tau_by_class"] = {vc: float((df[col][g & (df["variant_class"] == vc)] >= tau).mean()) for vc in sorted(set(df["variant_class"][g]))}
        out["real_pass_rate_at_tau"] = float((df[col][df["origin"] == "real"] >= tau).mean())
    return out
