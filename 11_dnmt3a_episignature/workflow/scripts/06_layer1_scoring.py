#!/usr/bin/env python
"""Step 06 - Layer 1 region scoring: per-sample region betas, leave-one-out z, shape_r, permutation null.

Region betas: unweighted mean of per-CpG betas over CpGs in the region; NA if n_cpgs < min_cpgs or mean depth < min_depth.
Regions come from signatures/hg38/<ID>.bed (tier2) and TBRS_LOF.tier1.bed. A background pool of array-probe-centred
regions (HM450+EPIC union, +/- region_pad, merged) is scored the same way and used for CpG-density-matched permutations.

Per sample x signature (x tier) x null set:
  mean_z        mean over regions of (beta - LOO control mean) / LOO control SD
  shape_r       Pearson r between (beta - LOO control mean) and the published delta-beta (sign vector where NA)
  frac_below_p5 fraction of regions with beta below the LOO controls' 5th percentile
  na_frac       fraction of signature regions NA for this sample
  perm_p_mean_z, perm_p_shape_r   empirical p from n_permutations random region sets matched on CpG count and width
Null sets: children = non-flagged, non-QC-outlier children (primary); all = non-flagged, non-QC-outlier samples (secondary).
The scored sample itself is always excluded (leave-one-out). Flagged samples and QC outliers are scored but never used
as controls.

Outputs (results/06_layer1/): region_betas.parquet (samples x regions), region_ncpgs.parquet, region_depth.parquet,
pool_betas.parquet, region_info.tsv, scores.tsv (long), permutation_null.parquet (per sample x signature summary).
"""
import argparse
import pathlib
import sys

import numpy as np
import pandas as pd
import yaml

MAIN_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
BED_COLS = ["chrom", "start", "end", "signature", "source", "direction", "published_delta_beta", "n_array_probes"]
METRICS = ["mean_z", "shape_r", "frac_below_p5", "na_frac"]


# ----------------------------------------------------------------------------- regions
def load_signature_regions(sig_dir, sigs):
    frames = []
    for sig in sigs:
        for tier, fn in (("tier2", f"{sig}.bed"), ("tier1", f"{sig}.tier1.bed")):
            p = pathlib.Path(sig_dir, fn)
            if not p.exists() or p.stat().st_size == 0:
                continue
            df = pd.read_csv(p, sep="\t", header=None, names=BED_COLS, na_values=["NA"])
            df["tier"] = tier
            df["set"] = f"{sig}:{tier}"
            frames.append(df)
    r = pd.concat(frames, ignore_index=True)
    r["region_id"] = r["set"] + ":" + r["chrom"] + ":" + r["start"].astype(str) + "-" + r["end"].astype(str)
    r["dir_num"] = r["direction"].map({"+": 1.0, "-": -1.0}).fillna(0.0)
    r["target"] = r["published_delta_beta"].where(r["published_delta_beta"].notna(), r["dir_num"])
    return r


def build_pool(res_dir, pad):
    """Array-probe-centred background regions (probe +/- pad, merged), hg38."""
    parts = []
    for f in sorted(pathlib.Path(res_dir, "manifests").glob("*.hg38.manifest.tsv.gz")):
        if "EPICv2" in f.name:
            continue
        m = pd.read_csv(f, sep="\t", usecols=["Probe_ID", "CpG_chrm", "CpG_beg"], low_memory=False)
        m = m[m["CpG_chrm"].isin(MAIN_CHROMS) & m["CpG_beg"].notna() & m["Probe_ID"].str.startswith("cg")]
        parts.append(m[["CpG_chrm", "CpG_beg"]].rename(columns={"CpG_chrm": "chrom", "CpG_beg": "pos"}))
    p = pd.concat(parts).drop_duplicates().astype({"pos": np.int64}).sort_values(["chrom", "pos"])
    out = []
    for chrom, g in p.groupby("chrom", sort=False):
        s = (g["pos"].to_numpy() - pad).clip(0); e = g["pos"].to_numpy() + 2 + pad
        cs, ce = s[0], e[0]
        for a, b in zip(s[1:], e[1:]):
            if a <= ce:
                ce = max(ce, b)
            else:
                out.append((chrom, cs, ce)); cs, ce = a, b
        out.append((chrom, cs, ce))
    pool = pd.DataFrame(out, columns=["chrom", "start", "end"])
    pool["region_id"] = "pool:" + pool["chrom"] + ":" + pool["start"].astype(str) + "-" + pool["end"].astype(str)
    return pool


def score_regions(parquet, regions_by_chrom, n_regions, min_cpgs, min_depth):
    """Return (beta, ncpg, depth) arrays aligned to region index for one sample. `parquet` may be a path or a loaded DataFrame."""
    d = parquet if isinstance(parquet, pd.DataFrame) else pd.read_parquet(parquet)
    beta = np.full(n_regions, np.nan, dtype=np.float32); ncpg = np.zeros(n_regions, dtype=np.int32); dep = np.zeros(n_regions, dtype=np.float32)
    for chrom, g in d.groupby("chrom", observed=True, sort=False):
        if chrom not in regions_by_chrom:
            continue
        starts, ends, idxs = regions_by_chrom[chrom]
        pos = g["pos"].to_numpy(); b = g["beta"].to_numpy(np.float64); dp = g["depth"].to_numpy(np.float64)
        k = np.searchsorted(starts, pos, side="right") - 1
        ok = (k >= 0) & (pos < ends[np.clip(k, 0, len(ends) - 1)])
        k = k[ok]; ridx = idxs[k]
        n = np.bincount(k, minlength=len(starts)); sb = np.bincount(k, weights=b[ok], minlength=len(starts)); sd = np.bincount(k, weights=dp[ok], minlength=len(starts))
        has = n > 0
        ncpg[idxs[has]] = n[has]; dep[idxs[has]] = sd[has] / n[has]; beta[idxs[has]] = sb[has] / n[has]
    bad = (ncpg < min_cpgs) | (dep < min_depth)
    beta[bad] = np.nan
    return beta, ncpg, dep


def index_regions(reg):
    """Per-chromosome sorted arrays (start, end, global index); regions within a frame must not overlap."""
    by = {}
    for chrom, g in reg.groupby("chrom", sort=False):
        g = g.sort_values("start")
        by[chrom] = (g["start"].to_numpy(), g["end"].to_numpy(), g.index.to_numpy())
    return by


# ----------------------------------------------------------------------------- statistics
def loo_stats(B, ctrl_mask):
    """B: samples x regions. For each sample i: mean/SD/p5 over controls excluding i. Returns arrays (samples x regions)."""
    n_s, n_r = B.shape
    C = B[ctrl_mask]
    cnt = np.sum(~np.isnan(C), axis=0).astype(float); s1 = np.nansum(C, axis=0); s2 = np.nansum(C ** 2, axis=0)
    mu = np.full((n_s, n_r), np.nan); sd = np.full((n_s, n_r), np.nan); p5 = np.full((n_s, n_r), np.nan)
    ctrl_idx = np.where(ctrl_mask)[0]
    for i in range(n_s):
        if ctrl_mask[i]:
            own = B[i]; ok = ~np.isnan(own)
            c = cnt - ok; a = s1 - np.where(ok, own, 0); q = s2 - np.where(ok, own ** 2, 0)
        else:
            c, a, q = cnt, s1, s2
        with np.errstate(invalid="ignore", divide="ignore"):
            m = a / c; v = (q - c * m ** 2) / (c - 1)
        mu[i] = np.where(c >= 3, m, np.nan); sd[i] = np.where(c >= 3, np.sqrt(np.clip(v, 0, None)), np.nan)
        sub = C if not ctrl_mask[i] else B[ctrl_idx[ctrl_idx != i]]
        p5[i] = np.nanpercentile(sub, 5, axis=0) if len(sub) >= 3 else np.nan
    return mu, sd, p5


def metrics(b, mu, sd, p5, target):
    """b, mu, sd, p5, target: region vectors for one sample/signature."""
    ok = ~np.isnan(b) & ~np.isnan(mu) & ~np.isnan(sd) & (sd > 0)
    out = dict(na_frac=float(np.mean(np.isnan(b))), n_used=int(ok.sum()))
    if ok.sum() < 5:
        out.update(mean_z=np.nan, shape_r=np.nan, frac_below_p5=np.nan); return out
    z = (b[ok] - mu[ok]) / sd[ok]; delta = b[ok] - mu[ok]
    out["mean_z"] = float(np.mean(z))
    t = target[ok]
    out["shape_r"] = float(np.corrcoef(delta, t)[0, 1]) if np.std(t) > 0 and np.std(delta) > 0 else np.nan
    okp = ok & ~np.isnan(p5)
    out["frac_below_p5"] = float(np.mean(b[okp] < p5[okp])) if okp.sum() else np.nan
    return out


def match_bins(reg_n, reg_w, pool_n, pool_w):
    """Joint bins on CpG count (deciles of pool) x width (3 bins). Returns bin id per region and per pool region."""
    qn = np.unique(np.nanquantile(pool_n, np.linspace(0, 1, 11)[1:-1])); qw = np.unique(np.nanquantile(pool_w, [1 / 3, 2 / 3]))
    def b(n, w): return np.searchsorted(qn, n) * 10 + np.searchsorted(qw, w)
    return b(reg_n, reg_w), b(pool_n, pool_w)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", required=True)
    ap.add_argument("--seed", type=int, default=20260909)
    ap.add_argument("--skip-perm", action="store_true")
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config)); th = cfg["thresholds"]
    out = pathlib.Path(cfg["paths"]["output_root"]); od = out / "06_layer1"; od.mkdir(parents=True, exist_ok=True)
    qc = pd.read_csv(out / "05_qc" / "sample_qc.tsv", sep="\t", dtype={"sample_id": str})
    samples = qc["sample_id"].tolist()
    reg = load_signature_regions(cfg["paths"]["signatures_dir"], cfg["signatures"]["layer1"])
    reg.to_csv(od / "region_info.tsv", sep="\t", index=False)
    pool = build_pool(cfg["paths"]["resources_dir"], int(th["region_pad"]))
    print(f"{len(reg)} signature regions ({reg['set'].nunique()} sets), {len(pool):,} pool regions, {len(samples)} samples")

    # ---- region betas (signature sets scored one set at a time so overlapping regions across sets are fine)
    B = np.full((len(samples), len(reg)), np.nan, np.float32); N = np.zeros_like(B, dtype=np.int32); D = np.zeros_like(B)
    P = np.full((len(samples), len(pool)), np.nan, np.float32); PN = np.zeros_like(P, dtype=np.int32)
    set_index = {s: index_regions(reg[reg["set"] == s]) for s in reg["set"].unique()}
    pool_index = index_regions(pool)
    for i, s in enumerate(samples):
        pq = pd.read_parquet(out / "04_methylation" / f"{s}.cpg.parquet")   # read once, score every set + the pool
        for sname, idx in set_index.items():
            b, n, d = score_regions(pq, idx, len(reg), th["min_cpgs"], th["min_depth"])
            sel = reg.index[reg["set"] == sname].to_numpy()
            B[i, sel], N[i, sel], D[i, sel] = b[sel], n[sel], d[sel]
        pb, pn, _ = score_regions(pq, pool_index, len(pool), th["min_cpgs"], th["min_depth"])
        P[i], PN[i] = pb, pn
        if i % 10 == 0:
            print(f"  scored {i + 1}/{len(samples)}")
    pd.DataFrame(B, index=samples, columns=reg["region_id"]).to_parquet(od / "region_betas.parquet")
    pd.DataFrame(N, index=samples, columns=reg["region_id"]).to_parquet(od / "region_ncpgs.parquet")
    pd.DataFrame(D, index=samples, columns=reg["region_id"]).to_parquet(od / "region_depth.parquet")
    pd.DataFrame(P, index=samples, columns=pool["region_id"]).to_parquet(od / "pool_betas.parquet")

    # ---- null sets
    clean = (qc["flag_gene"].fillna("") == "") & (~qc["qc_outlier"].astype(str).str.lower().eq("true"))
    nulls = {"children": (clean & (qc["role"] == "child")).to_numpy(), "all": clean.to_numpy()}
    print("null set sizes:", {k: int(v.sum()) for k, v in nulls.items()})
    rng = np.random.default_rng(a.seed)
    reg_n = np.nanmedian(np.where(N > 0, N, np.nan), axis=0); reg_w = (reg["end"] - reg["start"]).to_numpy()
    pool_n = np.nanmedian(np.where(PN > 0, PN, np.nan), axis=0); pool_w = (pool["end"] - pool["start"]).to_numpy()
    pool_ok = ~np.isnan(pool_n)
    rb, pb_ = match_bins(reg_n, reg_w, pool_n[pool_ok], pool_w[pool_ok]); pool_ids = np.where(pool_ok)[0]
    pool_by_bin = {bb: pool_ids[pb_ == bb] for bb in np.unique(pb_)}

    rows, perm_rows = [], []
    for null_name, cmask in nulls.items():
        mu, sd, p5 = loo_stats(B.astype(np.float64), cmask)
        pmu, psd, _ = loo_stats(P.astype(np.float64), cmask)
        for sname in reg["set"].unique():
            sel = reg.index[reg["set"] == sname].to_numpy(); target = reg.loc[sel, "target"].to_numpy(float)
            sig, tier = sname.split(":")
            # permutation sets: same size, bin-matched draw from the pool
            n_perm = 0 if a.skip_perm else int(th["n_permutations"])
            perm_sets = []
            for _ in range(n_perm):
                draw = np.empty(len(sel), dtype=int)
                for j, bb in enumerate(rb[sel]):
                    cand = pool_by_bin.get(bb)
                    if cand is None or len(cand) == 0:
                        cand = pool_ids
                    draw[j] = rng.choice(cand)
                perm_sets.append(draw)
            for i, s in enumerate(samples):
                m = metrics(B[i, sel].astype(float), mu[i, sel], sd[i, sel], p5[i, sel], target)
                if n_perm and not np.isnan(m["mean_z"]):
                    pz, pr = [], []
                    for draw in perm_sets:
                        mm = metrics(P[i, draw].astype(float), pmu[i, draw], psd[i, draw], np.full(len(draw), np.nan), target)
                        pz.append(mm["mean_z"]); pr.append(mm["shape_r"])
                    pz, pr = np.array(pz), np.array(pr)
                    # two-sided for mean_z (magnitude in the expected direction is judged in step 08), one-sided for shape_r > 0
                    m["perm_p_mean_z"] = float((np.sum(np.abs(pz[~np.isnan(pz)]) >= abs(m["mean_z"])) + 1) / (np.sum(~np.isnan(pz)) + 1))
                    m["perm_p_shape_r"] = float((np.sum(pr[~np.isnan(pr)] >= m["shape_r"]) + 1) / (np.sum(~np.isnan(pr)) + 1)) if not np.isnan(m["shape_r"]) else np.nan
                    perm_rows.append(dict(sample_id=s, signature=sig, tier=tier, null=null_name, perm_mean_z_mean=float(np.nanmean(pz)),
                                          perm_mean_z_sd=float(np.nanstd(pz)), perm_shape_r_mean=float(np.nanmean(pr)), perm_shape_r_sd=float(np.nanstd(pr))))
                for k, v in m.items():
                    rows.append(dict(sample_id=s, signature=sig, tier=tier, null=null_name, metric=k, value=v))
            print(f"  {null_name:8s} {sname:22s} done")
    scores = pd.DataFrame(rows)
    scores.to_csv(od / "scores.tsv", sep="\t", index=False)
    if perm_rows:
        pd.DataFrame(perm_rows).to_parquet(od / "permutation_null.parquet")
    # acceptance: no signature with >30% NA regions in the median sample
    na = scores[(scores["metric"] == "na_frac") & (scores["null"] == "children")].groupby(["signature", "tier"])["value"].median()
    nu = scores[(scores["metric"] == "n_used") & (scores["null"] == "children")].groupby(["signature", "tier"])["value"].median()
    print("median NA fraction per set:\n" + na.round(3).to_string())
    bad = na[na > 0.30]
    if len(bad):
        print(f"ACCEPTANCE WARNING: {len(bad)} set(s) exceed 30% NA in the median sample: {list(bad.index)}")
    ncpg_med = np.nanmedian(np.where(N > 0, N, np.nan), axis=0)
    structural = pd.Series(ncpg_med < th["min_cpgs"], index=reg["set"]).groupby(level=0).mean()
    log = ["# Step 06 Layer 1 scoring", "", f"samples {len(samples)}; null sizes " + str({k: int(v.sum()) for k, v in nulls.items()}) +
           f"; min_cpgs={th['min_cpgs']} min_depth={th['min_depth']} region_pad={th['region_pad']} n_perm={th['n_permutations']}", "",
           "| set | n regions | median NA frac (children null) | median regions used | frac regions with < min_cpgs CpGs (structural NA) |", "|---|---|---|---|---|"]
    for (sig, tier), v in na.items():
        key = f"{sig}:{tier}"
        log.append(f"| {key} | {int((reg['set'] == key).sum())} | {v:.3f} | {int(nu[(sig, tier)])} | {structural.get(key, float('nan')):.3f} |")
    log += ["", f"Acceptance (median NA <= 0.30 for every set): {'PASS' if not len(bad) else 'FAIL - ' + ', '.join(f'{a}:{b}' for a, b in bad.index)}"]
    (od / "STEP_LOG.md").write_text("\n".join(log) + "\n")


if __name__ == "__main__":
    main()
