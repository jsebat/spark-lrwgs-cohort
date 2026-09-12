#!/usr/bin/env python
"""Step 08 - panel similarity and ranking; proband report card against the pre-registered predictions.

Inputs : results/06_layer1/scores.tsv (long), results/05_qc/sample_qc.tsv, results/07_layer2/features_resid.tsv (optional)
Outputs: results/08_panel/shape_r.<null>.<tier>.tsv, mean_z.<null>.<tier>.tsv   (samples x signatures)
         results/08_panel/ranks.tsv        long: sample, signature, tier, null, metric, value, rank_children, n_children,
                                                 rank_all, n_all, pct_children, pct_all, is_proband, flag_gene, qc_outlier
         results/08_panel/proband_report_card.tsv   one row per prediction: predicted vs observed, PASS/FAIL
         results/08_panel/global_outliers.tsv       samples top-5 (children) on >= 4 signatures (shape_r)
         results/08_panel/flagged_and_family.tsv    ranks of flagged samples and of the proband's family
Ranks are 1 = highest shape_r (or highest mean_z magnitude in the signature's expected direction).
Nothing is softened: PASS/FAIL is mechanical from config.predictions.
"""
import argparse
import pathlib

import numpy as np
import pandas as pd
import yaml


def wide(scores, metric, null, tier):
    s = scores[(scores["metric"] == metric) & (scores["null"] == null) & (scores["tier"] == tier)]
    return s.pivot(index="sample_id", columns="signature", values="value")


def rank_desc(v):
    """1 = largest; NaN stays NaN; ties -> min rank."""
    return v.rank(ascending=False, method="min")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", required=True)
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config))
    out = pathlib.Path(cfg["paths"]["output_root"]); od = out / "08_panel"; od.mkdir(parents=True, exist_ok=True)
    scores = pd.read_csv(out / "06_layer1" / "scores.tsv", sep="\t", dtype={"sample_id": str})
    qc = pd.read_csv(out / "05_qc" / "sample_qc.tsv", sep="\t", dtype={"sample_id": str}).set_index("sample_id")
    children = qc.index[qc["role"] == "child"]
    proband = qc.index[qc["is_proband"].astype(str).str.lower() == "true"].tolist()
    fam = qc.index[qc["proband_family"].astype(str).str.lower() == "true"].tolist()
    print(f"samples {len(qc)}, children {len(children)}, proband rows {len(proband)}, proband-family rows {len(fam)}")

    rows = []
    for null in scores["null"].unique():
        for tier in scores["tier"].unique():
            for metric in ("shape_r", "mean_z"):
                W = wide(scores, metric, null, tier)
                if W.empty:
                    continue
                W.to_csv(od / f"{metric}.{null}.{tier}.tsv", sep="\t")
                for sig in W.columns:
                    v = W[sig]
                    if metric == "mean_z":  # rank magnitude in the expected direction: hypo signatures -> most negative first
                        sign = -1.0 if (scores[(scores.signature == sig)].shape[0] and cfg["predictions"].get(sig, {}).get("shape_r_sign") != "-") else 1.0
                        v_rank = v * sign  # TBRS-like (hypo) => more negative mean_z ranks first
                    else:
                        v_rank = v
                    rc = rank_desc(v_rank.loc[v_rank.index.intersection(children)]); ra = rank_desc(v_rank)
                    for s in W.index:
                        rows.append(dict(sample_id=s, signature=sig, tier=tier, null=null, metric=metric, value=v[s],
                                         rank_children=rc.get(s, np.nan), n_children=int(rc.notna().sum()),
                                         rank_all=ra.get(s, np.nan), n_all=int(ra.notna().sum()),
                                         is_proband=s in proband, proband_family=s in fam,
                                         flag_gene=qc.loc[s, "flag_gene"] if pd.notna(qc.loc[s, "flag_gene"]) else "",
                                         qc_outlier=str(qc.loc[s, "qc_outlier"]).lower() == "true", role=qc.loc[s, "role"]))
    R = pd.DataFrame(rows)
    R["pct_children"] = (R["rank_children"] - 1) / (R["n_children"] - 1).clip(lower=1)
    R["pct_all"] = (R["rank_all"] - 1) / (R["n_all"] - 1).clip(lower=1)
    R.to_csv(od / "ranks.tsv", sep="\t", index=False)

    # ---- global outliers: top-5 among children on >= 4 signatures (shape_r, children null, tier2)
    top = R[(R.metric == "shape_r") & (R.null == "children") & (R.tier == "tier2") & (R.rank_children <= 5)]
    g = top.groupby("sample_id")["signature"].agg(["count", lambda x: ",".join(sorted(x))]); g.columns = ["n_top5", "signatures"]
    g = g[g["n_top5"] >= 4].reset_index()
    for c in ("epi_frac", "mean_depth", "qc_outlier", "flag_gene", "role"):
        g[c] = g["sample_id"].map(qc[c])
    g["label"] = "global outlier, likely composition/QC"
    g.to_csv(od / "global_outliers.tsv", sep="\t", index=False)

    # ---- flagged samples and proband family
    ff = R[(R.flag_gene != "") | R.proband_family].copy()
    ff.to_csv(od / "flagged_and_family.tsv", sep="\t", index=False)

    # ---- report card
    card = []
    if not proband:
        print("WARNING: no proband row in the manifest/QC table; report card not produced")
    for p in proband:
        for sig, pred in cfg["predictions"].items():
            for tier in sorted(R[R.signature == sig]["tier"].unique(), reverse=True):  # tier2 first, then tier1 if present
                sub = R[(R.sample_id == p) & (R.signature == sig) & (R.tier == tier) & (R.metric == "shape_r")]
                for null in ("children", "all"):
                    r = sub[sub.null == null]
                    if r.empty:
                        continue
                    r = r.iloc[0]
                    # z of the proband's shape_r relative to the null set's shape_r distribution (excluding proband)
                    pool = R[(R.signature == sig) & (R.tier == tier) & (R.metric == "shape_r") & (R.null == null) & (R.sample_id != p)]
                    pool = pool[pool.role == "child"] if null == "children" else pool
                    pool = pool[(pool.flag_gene == "") & (~pool.qc_outlier)]
                    zr = (r.value - pool.value.mean()) / pool.value.std() if pool.value.std() > 0 else np.nan
                    perm = scores[(scores.sample_id == p) & (scores.signature == sig) & (scores.tier == tier) & (scores.null == null) & (scores.metric == "perm_p_shape_r")]
                    perm_p = float(perm.value.iloc[0]) if len(perm) else np.nan
                    obs_sign = "+" if r.value > 0 else "-" if r.value < 0 else "0"
                    if np.isnan(zr):
                        sign_ok = rank_ok = None
                    else:
                        sign_ok = {"+": bool(r.value > 0 and zr > 2), "-": bool(r.value < 0 and zr < -2), "0": bool(abs(zr) < 2)}.get(str(pred["shape_r_sign"]))
                        rk = str(pred["rank"])
                        rank_ok = {"1": bool(r.rank_children == 1), "top_quartile": bool(r.pct_children <= 0.25), "null": bool(abs(zr) < 2)}.get(rk)
                    rk = str(pred["rank"])
                    verdict = "PASS" if (sign_ok is True and rank_ok is True) else ("FAIL" if (sign_ok is False or rank_ok is False) else "not_evaluable")
                    card.append(dict(sample="proband", signature=sig, tier=tier, null=null, predicted_sign=pred["shape_r_sign"], predicted_rank=rk,
                                     observed_shape_r=round(r.value, 4), observed_sign=obs_sign, shape_r_z_vs_null=round(zr, 2) if not np.isnan(zr) else np.nan,
                                     rank_children=r.rank_children, n_children=r.n_children, rank_all=r.rank_all, n_all=r.n_all,
                                     perm_p_shape_r=perm_p, sign_ok=sign_ok, rank_ok=rank_ok, verdict=verdict, note=pred.get("note", "")))
        # Layer 2 predictions
        l2 = out / "07_layer2" / "features_resid.tsv"
        if l2.exists():
            fr = pd.read_csv(l2, sep="\t", dtype={"sample_id": str}).set_index("sample_id")
            for feat, pred in cfg.get("layer2_predictions", {}).items():
                if feat not in fr.columns or p not in fr.index:
                    continue
                v = fr.loc[p, feat]
                if "resid_within" in pred:
                    ok = bool(abs(v) < pred["resid_within"]); predicted = f"|resid| < {pred['resid_within']}"
                else:
                    ok = bool((v > 0) if pred["sign"] == "+" else (v < 0)); predicted = f"sign {pred['sign']}"
                card.append(dict(sample="proband", signature=f"L2:{feat}", tier="-", null="all", predicted_sign=predicted, predicted_rank="-",
                                 observed_shape_r=round(float(v), 3), observed_sign="+" if v > 0 else "-", shape_r_z_vs_null=np.nan,
                                 rank_children=np.nan, n_children=np.nan, rank_all=np.nan, n_all=np.nan, perm_p_shape_r=np.nan,
                                 sign_ok=ok, rank_ok=ok, verdict="PASS" if ok else "FAIL", note=pred.get("note", "")))
    C = pd.DataFrame(card)
    C.to_csv(od / "proband_report_card.tsv", sep="\t", index=False)
    if len(C):
        show = C[(C.null == "children") | C.signature.str.startswith("L2:")][["signature", "tier", "predicted_sign", "predicted_rank", "observed_shape_r", "shape_r_z_vs_null", "rank_children", "n_children", "perm_p_shape_r", "verdict"]]
        print(show.to_string(index=False))
        print(f"\nglobal outliers: {len(g)}; flagged/family rows: {ff.sample_id.nunique()} samples")


if __name__ == "__main__":
    main()
