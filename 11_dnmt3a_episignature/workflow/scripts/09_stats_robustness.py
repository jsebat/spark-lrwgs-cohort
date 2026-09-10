#!/usr/bin/env python
"""Step 09 - statistics and robustness for the proband's TBRS_LOF result.

1. Empirical rank p: proband shape_r rank among children (primary) and among all samples (secondary), p = rank / n.
2. Sensitivity of the proband's TBRS_LOF shape_r rank (children null) to:
     - tier1 vs tier2 (from step 06 scores)
     - region_pad 0 / 250 / 500: TBRS_LOF regions rebuilt from results/02_signatures/TBRS_LOF.probes.tsv at each pad
       and re-scored from the per-sample parquets (same min_cpgs/min_depth rules)
     - composition-clean: drop regions whose beta correlates with epithelial fraction (|r| > 0.5) across clean controls
3. Bootstrap over control samples (n_boot resamples of the clean children with replacement): distribution of the
   proband's rank and shape_r; reports the fraction of bootstraps with rank 1.
Output: results/09_stats/sensitivity.tsv (one row per analysis), bootstrap.tsv, STEP_LOG.md.
"""
import argparse
import importlib.util
import pathlib

import numpy as np
import pandas as pd
import yaml

HERE = pathlib.Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("m06", HERE / "06_layer1_scoring.py"); m06 = importlib.util.module_from_spec(spec); spec.loader.exec_module(m06)
spec2 = importlib.util.spec_from_file_location("m02", HERE / "02_harmonize_signatures.py"); m02 = importlib.util.module_from_spec(spec2); spec2.loader.exec_module(m02)


def score_matrix(samples, out, reg, th):
    idx = m06.index_regions(reg.reset_index(drop=True))
    B = np.full((len(samples), len(reg)), np.nan)
    for i, s in enumerate(samples):
        b, _, _ = m06.score_regions(out / "04_methylation" / f"{s}.cpg.parquet", idx, len(reg), th["min_cpgs"], th["min_depth"])
        B[i] = b
    return B


def rank_result(B, target, samples, ctrl_mask, child_mask, proband):
    mu, sd, p5 = m06.loo_stats(B, ctrl_mask)
    r = np.array([m06.metrics(B[i], mu[i], sd[i], p5[i], target)["shape_r"] for i in range(len(samples))])
    z = np.array([m06.metrics(B[i], mu[i], sd[i], p5[i], target)["mean_z"] for i in range(len(samples))])
    pi = samples.index(proband)
    rc = pd.Series(r[child_mask], index=np.array(samples)[child_mask]).rank(ascending=False, method="min")
    ra = pd.Series(r, index=samples).rank(ascending=False, method="min")
    return dict(shape_r=r[pi], mean_z=z[pi], rank_children=rc.get(proband, np.nan), n_children=int(np.sum(~np.isnan(r[child_mask]))),
                rank_all=ra.get(proband, np.nan), n_all=int(np.sum(~np.isnan(r))), n_regions=B.shape[1]), r


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", required=True)
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config)); th = cfg["thresholds"]
    out = pathlib.Path(cfg["paths"]["output_root"]); od = out / "09_stats"; od.mkdir(parents=True, exist_ok=True)
    qc = pd.read_csv(out / "05_qc" / "sample_qc.tsv", sep="\t", dtype={"sample_id": str})
    samples = qc["sample_id"].tolist()
    proband = qc.loc[qc["is_proband"].astype(str).str.lower() == "true", "sample_id"].tolist()
    if len(proband) != 1:
        raise SystemExit(f"expected exactly one proband, found {len(proband)}")
    proband = proband[0]
    clean = ((qc["flag_gene"].fillna("") == "") & ~qc["qc_outlier"].astype(str).str.lower().eq("true")).to_numpy()
    child = (qc["role"] == "child").to_numpy()
    ctrl_children = clean & child; ctrl_all = clean
    rows = []

    # 1 + tier sensitivity from step 06 scores
    R = pd.read_csv(out / "08_panel" / "ranks.tsv", sep="\t", dtype={"sample_id": str})
    for tier in sorted(R.tier.unique(), reverse=True):
        for null in ("children", "all"):
            r = R[(R.sample_id == proband) & (R.signature == "TBRS_LOF") & (R.tier == tier) & (R.null == null) & (R.metric == "shape_r")]
            if r.empty:
                continue
            r = r.iloc[0]
            rows.append(dict(analysis=f"step06:{tier}", null=null, shape_r=r.value, rank_children=r.rank_children, n_children=r.n_children,
                             rank_all=r.rank_all, n_all=r.n_all, p_rank_children=r.rank_children / r.n_children, p_rank_all=r.rank_all / r.n_all, n_regions=np.nan))

    # 2 region_pad sensitivity (TBRS_LOF tier2 union rebuilt from the per-probe/DMR table)
    probes = pd.read_csv(out / "02_signatures" / "TBRS_LOF.probes.tsv", sep="\t")
    B_ref = None
    for pad in (0, 250, 500):
        reg = m02.merge_regions(probes, pad, "TBRS_LOF")
        reg["target"] = pd.to_numeric(reg["published_delta_beta"], errors="coerce").fillna(reg["direction"].map({"+": 1.0, "-": -1.0}).fillna(0.0))
        B = score_matrix(samples, out, reg, th)
        res, _ = rank_result(B, reg["target"].to_numpy(float), samples, ctrl_children, child, proband)
        rows.append(dict(analysis=f"region_pad={pad}", null="children", p_rank_children=res["rank_children"] / res["n_children"], **res))
        if pad == int(th["region_pad"]):
            B_ref, reg_ref = B, reg

    # 2a evidence-threshold sensitivity (pre-registered: min_cpgs=10, min_depth=15). NA-driven, so re-scored from parquets.
    for mc, md in ((10, 15), (10, 10), (5, 10), (5, 15)):
        th2 = dict(th, min_cpgs=mc, min_depth=md)
        B = B_ref if (mc, md) == (int(th["min_cpgs"]), int(th["min_depth"])) else score_matrix(samples, out, reg_ref, th2)
        na = float(np.isnan(B[ctrl_children]).mean())
        res, _ = rank_result(B, reg_ref["target"].to_numpy(float), samples, ctrl_children, child, proband)
        rows.append(dict(analysis=f"min_cpgs={mc},min_depth={md} (region_pad={int(th['region_pad'])}; control NA frac {na:.2f})", null="children",
                         p_rank_children=res["rank_children"] / res["n_children"], **res))

    # 2b composition-clean: drop regions correlated with epithelial fraction in clean controls
    epi = qc["epi_frac"].to_numpy(float)
    keep = np.ones(B_ref.shape[1], bool)
    for j in range(B_ref.shape[1]):
        x = B_ref[ctrl_children, j]; y = epi[ctrl_children]; ok = ~np.isnan(x) & ~np.isnan(y)
        if ok.sum() >= 8 and np.std(x[ok]) > 0 and np.std(y[ok]) > 0 and abs(np.corrcoef(x[ok], y[ok])[0, 1]) > 0.5:
            keep[j] = False
    res, _ = rank_result(B_ref[:, keep], reg_ref["target"].to_numpy(float)[keep], samples, ctrl_children, child, proband)
    rows.append(dict(analysis=f"composition_clean (dropped {int((~keep).sum())} of {len(keep)} regions with |r(epi_frac)|>0.5)", null="children",
                     p_rank_children=res["rank_children"] / res["n_children"], **res))

    # 3 bootstrap over clean children controls
    rng = np.random.default_rng(a.seed)
    ctrl_idx = np.where(ctrl_children)[0]
    target = reg_ref["target"].to_numpy(float)
    boot = []
    pi = samples.index(proband)
    for b in range(a.n_boot):
        draw = rng.choice(ctrl_idx, size=len(ctrl_idx), replace=True)
        idx = np.concatenate([[pi], draw]) if pi not in draw else draw
        Bb = B_ref[idx]; names = [f"{samples[i]}#{k}" for k, i in enumerate(idx)]  # duplicates need distinct names
        cm = np.array([i != pi for i in idx]); ch = np.ones(len(idx), bool)
        mu, sd, p5 = m06.loo_stats(Bb, cm)
        r = np.array([m06.metrics(Bb[k], mu[k], sd[k], p5[k], target)["shape_r"] for k in range(len(idx))])
        pk = list(idx).index(pi)
        rank = int(1 + np.sum(r[np.arange(len(idx)) != pk] > r[pk])) if not np.isnan(r[pk]) else np.nan
        boot.append(dict(boot=b, shape_r=r[pk], rank=rank, n=len(idx)))
    BT = pd.DataFrame(boot); BT.to_csv(od / "bootstrap.tsv", sep="\t", index=False)
    rows.append(dict(analysis=f"bootstrap_controls (n_boot={a.n_boot})", null="children", shape_r=BT.shape_r.mean(), rank_children=BT["rank"].median(),
                     n_children=int(BT["n"].median()), p_rank_children=float((BT["rank"] == 1).mean()),
                     rank_all=np.nan, n_all=np.nan, p_rank_all=np.nan, n_regions=B_ref.shape[1]))

    S = pd.DataFrame(rows)
    S.to_csv(od / "sensitivity.tsv", sep="\t", index=False)
    log = ["# Step 09 - statistics and robustness (proband, TBRS_LOF)", "", S.round(4).to_string(index=False), "",
           f"bootstrap: fraction of resamples with proband rank 1 = {(BT['rank'] == 1).mean():.3f}; rank quantiles (5/50/95%) = "
           f"{BT['rank'].quantile([.05, .5, .95]).tolist()}"]
    (od / "STEP_LOG.md").write_text("\n".join(log) + "\n"); print("\n".join(log))


if __name__ == "__main__":
    main()
