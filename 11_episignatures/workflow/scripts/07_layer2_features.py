#!/usr/bin/env python
"""Step 07 - Layer 2 genome-architecture features, covariate-adjusted.

Per sample (CpGs with depth >= min_depth_cpg only):
  enhancer_beta, promoter_beta, cgi_beta, genebody_beta        ENCODE cCREs (pELS+dELS / PLS), UCSC CGIs, RefSeq curated tx bodies (+2 kb from TSS to end)
  satII_III_beta, alpha_sat_beta, line1_beta, line1_young_beta   RepeatMasker: HSATII/HSAT3/(GAATG)n/(CATTC)n; ALR/Alpha; L1 family; L1HS+L1PA2+L1PA3
  pmd_soloWCGW_beta                                            Zhou 2018 solo-WCGW CpGs in common PMDs
  hmd_soloWCGW_beta                                            solo-WCGW CpGs inside common HMDs (from PMD_coordinates, WCGW context from reference)
  canyon_total_bp, canyon_mean_width_bp, canyon_interior_beta, canyon_edge_inner_beta, canyon_edge_outer_beta, canyon_edge_slope
        canyons = control-mean track (clean samples) runs of CpGs with beta < 0.10, inter-CpG gap <= 500 bp, span >= 3.5 kb;
        per sample: interior mean beta, mean beta in the 1 kb just inside vs just outside each boundary (slope = outer - inner),
        and the sample's own extent (grow from canyon centre while sample beta < 0.10) -> total bp and mean width
  global_mean_5mc (step 04), dnam_age / age_accel_resid (step 05)

Residualisation: feature ~ epi_frac + age + sex + mean_depth fitted on clean samples (non-flagged, non-outlier);
standardized residual = (obs - fitted) / SD(residuals of clean samples). Outputs results/07_layer2/features.tsv,
features_resid.tsv, canyons.control.bed, STEP_LOG.md (acceptance: control residual means ~0; satellite CV < 15%).
"""
import argparse
import gzip
import pathlib

import numpy as np
import pandas as pd
import yaml

MAIN_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
CHR_IDX = {c: i for i, c in enumerate(MAIN_CHROMS)}


def merge(df):
    """Merge overlapping intervals (chrom,start,end) -> dict chrom -> (starts, ends)."""
    out = {}
    for chrom, g in df.groupby("chrom", sort=False):
        g = g.sort_values("start"); s = g["start"].to_numpy(); e = g["end"].to_numpy()
        ms, me = [s[0]], [e[0]]
        for a, b in zip(s[1:], e[1:]):
            if a <= me[-1]:
                me[-1] = max(me[-1], b)
            else:
                ms.append(a); me.append(b)
        out[chrom] = (np.array(ms), np.array(me))
    return out


def in_intervals(pos, iv):
    if iv is None:
        return np.zeros(len(pos), bool)
    s, e = iv
    k = np.searchsorted(s, pos, side="right") - 1
    return (k >= 0) & (pos < e[np.clip(k, 0, len(e) - 1)])


def load_annotations(R):
    R = pathlib.Path(R)
    cc = pd.read_csv(R / "GRCh38-cCREs.bed", sep="\t", header=None, usecols=[0, 1, 2, 5], names=["chrom", "start", "end", "group"])
    cc = cc[cc["chrom"].isin(MAIN_CHROMS)]
    ann = {"enhancer": merge(cc[cc["group"].str.contains("ELS")]), "promoter": merge(cc[cc["group"].str.contains("PLS")])}
    cgi = pd.read_csv(R / "cpgIslandExt.hg38.txt.gz", sep="\t", header=None, usecols=[1, 2, 3], names=["chrom", "start", "end"])
    ann["cgi"] = merge(cgi[cgi["chrom"].isin(MAIN_CHROMS)])
    rs = pd.read_csv(R / "ncbiRefSeqCurated.hg38.txt.gz", sep="\t", header=None, usecols=[2, 3, 4, 5], names=["chrom", "strand", "start", "end"])
    rs = rs[rs["chrom"].isin(MAIN_CHROMS) & ((rs["end"] - rs["start"]) > 4000)].copy()
    rs["start"] = np.where(rs["strand"] == "+", rs["start"] + 2000, rs["start"]); rs["end"] = np.where(rs["strand"] == "-", rs["end"] - 2000, rs["end"])
    ann["genebody"] = merge(rs[["chrom", "start", "end"]])
    rm = pd.read_csv(R / "rmsk.hg38.txt.gz", sep="\t", header=None, usecols=[5, 6, 7, 10, 11, 12], names=["chrom", "start", "end", "name", "cls", "fam"])
    rm = rm[rm["chrom"].isin(MAIN_CHROMS)]
    ann["satII_III"] = merge(rm[rm["name"].isin(["HSATII", "HSAT3", "(GAATG)n", "(CATTC)n"])])
    ann["alpha_sat"] = merge(rm[rm["name"] == "ALR/Alpha"])
    ann["line1"] = merge(rm[rm["fam"] == "L1"])
    ann["line1_young"] = merge(rm[rm["name"].isin(["L1HS", "L1PA2", "L1PA3"])])
    pmd = pd.read_csv(R / "PMD_coordinates_hg38.bed.gz", sep="\t", header=None, usecols=[0, 1, 2, 5], names=["chrom", "start", "end", "name"])  # col6 = commonPMD/commonHMD (100 kb bins)
    pmd = pmd[pmd["chrom"].isin(MAIN_CHROMS)]
    ann["commonPMD"] = merge(pmd[pmd["name"].str.lower().str.contains("commonpmd")]) if pmd["name"].str.lower().str.contains("commonpmd").any() else merge(pmd[pmd["name"].str.contains("PMD")])
    ann["commonHMD"] = merge(pmd[pmd["name"].str.lower().str.contains("commonhmd")]) if pmd["name"].str.lower().str.contains("commonhmd").any() else None
    sw = pd.read_csv(R / "solo_WCGW_inCommonPMDs_hg38.bed.gz", sep="\t", header=None, usecols=[0, 1], names=["chrom", "pos"])
    sw = sw[sw["chrom"].isin(MAIN_CHROMS)]
    ann["_soloWCGW_pmd_pos"] = {c: np.sort(g["pos"].to_numpy()) for c, g in sw.groupby("chrom")}
    return ann


def key(chrom_idx, pos):
    return chrom_idx.astype(np.int64) * 400_000_000 + pos.astype(np.int64)


def control_mean_track(parquets, min_depth):
    """Universe = CpGs covered in the first control; accumulate mean beta over controls (depth >= min_depth)."""
    first = pd.read_parquet(parquets[0]); first = first[first["depth"] >= min_depth]
    ci = first["chrom"].astype(str).map(CHR_IDX); ok = ci.notna()
    k0 = np.sort(key(ci[ok].to_numpy(), first["pos"][ok].to_numpy()))
    s = np.zeros(len(k0)); n = np.zeros(len(k0))
    for pq in parquets:
        d = pd.read_parquet(pq); d = d[d["depth"] >= min_depth]
        ci = d["chrom"].astype(str).map(CHR_IDX); ok = ci.notna()
        k = key(ci[ok].to_numpy(), d["pos"][ok].to_numpy()); b = d["beta"][ok].to_numpy(float)
        j = np.searchsorted(k0, k); hit = (j < len(k0)) & (k0[np.clip(j, 0, len(k0) - 1)] == k)
        np.add.at(s, j[hit], b[hit]); np.add.at(n, j[hit], 1)
    ok = n >= max(3, 0.5 * len(parquets))
    return k0[ok], s[ok] / n[ok]


def call_canyons(k, mean_beta, thr=0.10, max_gap=500, min_span=3500):
    chrom = k // 400_000_000; pos = k % 400_000_000
    low = mean_beta < thr
    rows = []
    start = None
    for i in range(len(k)):
        if low[i] and (start is None):
            start, last = i, i
        elif low[i] and chrom[i] == chrom[last] and pos[i] - pos[last] <= max_gap:
            last = i
        else:
            if start is not None and pos[last] - pos[start] + 2 >= min_span and chrom[start] == chrom[last]:
                rows.append((MAIN_CHROMS[chrom[start]], int(pos[start]), int(pos[last] + 2)))
            start, last = (i, i) if low[i] else (None, None)
    if start is not None and pos[last] - pos[start] + 2 >= min_span:
        rows.append((MAIN_CHROMS[chrom[start]], int(pos[start]), int(pos[last] + 2)))
    return pd.DataFrame(rows, columns=["chrom", "start", "end"])


def sample_features(pq, ann, canyons, min_depth, thr=0.10):
    d = pd.read_parquet(pq); d = d[(d["depth"] >= min_depth) & d["chrom"].astype(str).isin(MAIN_CHROMS)]
    f = {}
    sums = {k: [0.0, 0] for k in ["enhancer", "promoter", "cgi", "genebody", "satII_III", "alpha_sat", "line1", "line1_young", "pmd_soloWCGW", "hmd_soloWCGW"]}
    can = {"interior": [0.0, 0], "inner": [0.0, 0], "outer": [0.0, 0]}; widths = []
    for chrom, g in d.groupby("chrom", observed=True, sort=False):
        chrom = str(chrom); pos = g["pos"].to_numpy(); b = g["beta"].to_numpy(float)
        for name in ["enhancer", "promoter", "cgi", "genebody", "satII_III", "alpha_sat", "line1", "line1_young"]:
            m = in_intervals(pos, ann[name].get(chrom)); sums[name][0] += b[m].sum(); sums[name][1] += int(m.sum())
        sw = ann["_soloWCGW_pmd_pos"].get(chrom)
        if sw is not None:
            m = np.isin(pos, sw); sums["pmd_soloWCGW"][0] += b[m].sum(); sums["pmd_soloWCGW"][1] += int(m.sum())
        if ann.get("commonHMD"):
            # solo-WCGW context in HMDs: approximate with CpGs in common HMDs that are NOT in any PMD solo-WCGW list
            m = in_intervals(pos, ann["commonHMD"].get(chrom)); sums["hmd_soloWCGW"][0] += b[m].sum(); sums["hmd_soloWCGW"][1] += int(m.sum())
        cg = canyons[canyons["chrom"] == chrom]
        for s, e in zip(cg["start"].to_numpy(), cg["end"].to_numpy()):
            i0, i1 = np.searchsorted(pos, [s, e])
            inside = b[i0:i1]
            if len(inside):
                can["interior"][0] += inside.sum(); can["interior"][1] += len(inside)
            for edge, inner_rng, outer_rng in ((s, (s, min(s + 1000, e)), (max(s - 1000, 0), s)), (e, (max(e - 1000, s), e), (e, e + 1000))):
                a0, a1 = np.searchsorted(pos, inner_rng); o0, o1 = np.searchsorted(pos, outer_rng)
                can["inner"][0] += b[a0:a1].sum(); can["inner"][1] += a1 - a0; can["outer"][0] += b[o0:o1].sum(); can["outer"][1] += o1 - o0
            # sample's own canyon extent from the centre (skip if the sample has no CpG inside this canyon)
            c = (s + e) // 2; ic = int(np.searchsorted(pos, c))
            if ic >= len(pos) or not (s <= pos[ic] < e) or b[ic] >= thr:
                continue
            lo = ic
            while lo - 1 >= 0 and b[lo - 1] < thr and pos[lo] - pos[lo - 1] <= 500:
                lo -= 1
            hi = ic
            while hi + 1 < len(pos) and b[hi + 1] < thr and pos[hi + 1] - pos[hi] <= 500:
                hi += 1
            if hi > lo:
                widths.append(int(pos[hi] - pos[lo] + 2))
    for k, (s_, n_) in sums.items():
        f[f"{k}_beta"] = s_ / n_ if n_ else np.nan; f[f"{k}_ncpg"] = n_
    f["canyon_interior_beta"] = can["interior"][0] / can["interior"][1] if can["interior"][1] else np.nan
    f["canyon_edge_inner_beta"] = can["inner"][0] / can["inner"][1] if can["inner"][1] else np.nan
    f["canyon_edge_outer_beta"] = can["outer"][0] / can["outer"][1] if can["outer"][1] else np.nan
    f["canyon_edge_slope"] = f["canyon_edge_outer_beta"] - f["canyon_edge_inner_beta"]
    f["canyon_total_bp"] = int(np.sum(widths)); f["canyon_mean_width_bp"] = float(np.mean(widths)) if widths else np.nan; f["canyon_n"] = len(widths)
    return f


def residualize(feat, qc, clean, cols):
    X = pd.DataFrame({"epi_frac": qc["epi_frac"], "age": qc["age"], "sex": (qc["sex"].str.lower() == "male").astype(float), "mean_depth": qc["mean_depth"]}, index=qc.index)
    X = X.fillna(X.mean()); X.insert(0, "const", 1.0)
    out = pd.DataFrame(index=qc.index); fits = []
    for c in cols:
        y = feat[c].astype(float); ok = clean & y.notna()
        if ok.sum() < 8:
            out[c] = np.nan; continue
        beta = np.linalg.lstsq(X[ok].to_numpy(), y[ok].to_numpy(), rcond=None)[0]
        resid = y - X.to_numpy() @ beta
        sd = resid[ok].std(ddof=X.shape[1])
        out[c] = resid / sd if sd > 0 else np.nan
        fits.append(dict(feature=c, n_fit=int(ok.sum()), resid_sd=float(sd), ctrl_resid_mean=float(out.loc[ok, c].mean())))
    return out, pd.DataFrame(fits)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", required=True)
    ap.add_argument("--min-depth-cpg", type=int, default=10)
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config))
    out = pathlib.Path(cfg["paths"]["output_root"]); od = out / "07_layer2"; od.mkdir(parents=True, exist_ok=True)
    qc = pd.read_csv(out / "05_qc" / "sample_qc.tsv", sep="\t", dtype={"sample_id": str})
    clean = (qc["flag_gene"].fillna("") == "") & ~qc["qc_outlier"].astype(str).str.lower().eq("true")
    ann = load_annotations(pathlib.Path(cfg["paths"]["resources_dir"], "layer2"))
    pqs = [out / "04_methylation" / f"{s}.cpg.parquet" for s in qc["sample_id"]]
    ctrl_pqs = [p for p, c in zip(pqs, clean) if c]
    print(f"{len(qc)} samples, {len(ctrl_pqs)} clean controls for the canyon track")
    k, mb = control_mean_track(ctrl_pqs, a.min_depth_cpg)
    canyons = call_canyons(k, mb)
    canyons.to_csv(od / "canyons.control.bed", sep="\t", header=False, index=False)
    print(f"control canyons: {len(canyons)} (median width {int((canyons['end'] - canyons['start']).median()) if len(canyons) else 0} bp)")
    rows = []
    for s, pq in zip(qc["sample_id"], pqs):
        f = sample_features(pq, ann, canyons, a.min_depth_cpg); f["sample_id"] = s; rows.append(f)
    feat = pd.DataFrame(rows).set_index("sample_id").loc[qc["sample_id"]].reset_index()
    feat["global_mean_5mc"] = qc["global_mean_5mc"].to_numpy(); feat["dnam_age"] = qc.get("dnam_age"); feat["age_accel_resid"] = qc.get("age_accel_resid")
    feat.to_csv(od / "features.tsv", sep="\t", index=False)
    cols = [c for c in feat.columns if c.endswith("_beta") or c.startswith("canyon_") and not c.endswith("_n")] + ["global_mean_5mc"]
    resid, fits = residualize(feat.set_index("sample_id").reset_index(drop=True), qc, clean, cols)
    resid.insert(0, "sample_id", qc["sample_id"]); resid["age_accel_resid"] = qc.get("age_accel_resid")
    resid.to_csv(od / "features_resid.tsv", sep="\t", index=False)
    sat = feat.loc[clean.to_numpy(), ["satII_III_beta", "alpha_sat_beta", "line1_beta"]]
    cv = (sat.std() / sat.mean() * 100).round(2)
    log = ["# Step 07 Layer 2 features", "", f"clean controls used for fits: {int(clean.sum())} / {len(qc)}; control canyons: {len(canyons)}", "",
           "Satellite/LINE-1 CV across controls (%): " + ", ".join(f"{k}={v}" for k, v in cv.items()), "",
           fits.round(3).to_string(index=False)]
    (od / "STEP_LOG.md").write_text("\n".join(log) + "\n"); print("\n".join(log))


if __name__ == "__main__":
    main()
