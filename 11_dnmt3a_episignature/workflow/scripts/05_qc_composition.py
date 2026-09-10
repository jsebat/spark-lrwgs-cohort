#!/usr/bin/env python
"""Step 05 - cohort QC: array-probe beta matrix, HEpiDISH composition, DNAm age, global 5mC, outlier flags.

1. For every sample parquet from step 04, take the beta at each Illumina probe CpG (HM450 + EPIC union, hg38 via
   the Zhou manifests; pb-CpG-tools `pos` and manifest `CpG_beg` are both the 0-based C of the CpG). Depth < min_depth_probe
   -> NA. Writes results/05_qc/probe_betas.parquet (probes x samples) and probe_missing.parquet (bool mask).
2. Cohort-mean imputation of missing probe betas (probes missing in >50% of samples are dropped instead);
   per-sample % imputed is reported overall and, by the R step, for the clock CpGs specifically.
3. Calls workflow/scripts/05b_qc_r.R (EpiDISH HEpiDISH; methylclock DNAmAge) on the imputed matrix.
4. Merges with step 04 qc.json (depth, global mean 5mC) and the manifest -> results/05_qc/sample_qc.tsv with
   epi_frac, fib_frac, immune_frac, dnam_age_<clock>, age_accel_resid_<clock>, dnam_age, age_accel_resid,
   global_mean_5mc, pct_clock_cpgs_imputed, qc_outlier (+ reason). Outliers: epi_frac > median + 3*MAD, or mean depth
   < min_depth_sample, or |global 5mC z| > 4. Outliers are kept and marked, never dropped here.
"""
import argparse
import json
import pathlib
import subprocess
import sys

import numpy as np
import pandas as pd
import yaml

MAIN_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def load_probe_positions(res_dir):
    frames = []
    for f in sorted(pathlib.Path(res_dir, "manifests").glob("*.hg38.manifest.tsv.gz")):
        if "EPICv2" in f.name:
            continue
        df = pd.read_csv(f, sep="\t", usecols=["Probe_ID", "CpG_chrm", "CpG_beg"], low_memory=False)
        df = df[df["CpG_chrm"].isin(MAIN_CHROMS) & df["CpG_beg"].notna() & df["Probe_ID"].str.startswith("cg")]
        frames.append(df.rename(columns={"Probe_ID": "probe", "CpG_chrm": "chrom", "CpG_beg": "pos"}))
    p = pd.concat(frames).drop_duplicates("probe")
    p["pos"] = p["pos"].astype(np.int64)
    return p.set_index(["chrom", "pos"])["probe"]


def probe_betas(parquet, probes, min_depth):
    d = pd.read_parquet(parquet)
    d = d[d["depth"] >= min_depth]
    d = d.set_index([d["chrom"].astype(str), "pos"])["beta"]
    hit = d.reindex(probes.index)
    return pd.Series(hit.values, index=probes.values, dtype="float32")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", required=True)
    ap.add_argument("--min-depth-probe", type=int, default=10)
    ap.add_argument("--min-depth-sample", type=float, default=12.0)
    ap.add_argument("--rscript", default="Rscript")
    ap.add_argument("--r-only", action="store_true", help="reuse results/05_qc/probe_betas.parquet + probe_missing.parquet; rerun only the R step and merge")
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config))
    out = pathlib.Path(cfg["paths"]["output_root"]); qc_dir = out / "05_qc"; qc_dir.mkdir(parents=True, exist_ok=True)
    man = pd.read_csv(out / "03_manifest" / "samples.frozen.tsv", sep="\t", dtype=str)
    man["age"] = pd.to_numeric(man["age"], errors="coerce")
    probes = load_probe_positions(cfg["paths"]["resources_dir"])
    print(f"{len(probes):,} array probe CpGs on hg38")

    cols, qcs, missing_samples, failed = {}, [], [], []
    qcj = {s: json.load(open(str(out / "04_methylation" / f"{s}.cpg.parquet.qc.json"))) for s in man["sample_id"]
           if (out / "04_methylation" / f"{s}.cpg.parquet.qc.json").exists()}
    med_cpgs = float(np.median([q["n_cpgs"] for q in qcj.values()])) if qcj else 0
    for s in man["sample_id"]:
        pq = out / "04_methylation" / f"{s}.cpg.parquet"
        if not pq.exists():
            missing_samples.append(s); continue
        q = qcj[s]
        qcs.append(dict(sample_id=s, **{k: v for k, v in q.items() if k != "sample"}))
        if q["n_cpgs"] < 0.5 * med_cpgs:      # degenerate extraction: keep the row, exclude from every matrix
            failed.append(s); continue
        if not a.r_only:
            cols[s] = probe_betas(pq, probes, a.min_depth_probe)
    if missing_samples:
        print(f"WARNING: {len(missing_samples)} manifest samples have no step-04 parquet; excluded from QC")
    if failed:
        print(f"WARNING: {len(failed)} sample(s) with < 50% of the median CpG count (extraction failure); kept in sample_qc.tsv as qc_outlier, excluded from all matrices")
    if a.r_only:
        B = pd.read_parquet(qc_dir / "probe_betas.parquet"); miss = pd.read_parquet(qc_dir / "probe_missing.parquet")
        print(f"--r-only: reusing saved probe matrix {B.shape}")
    else:
        B = pd.DataFrame(cols)                   # probes x samples
        miss = B.isna()
        keep = miss.mean(axis=1) <= 0.5
        B, miss = B[keep], miss[keep]
        B.to_parquet(qc_dir / "probe_betas.parquet"); miss.to_parquet(qc_dir / "probe_missing.parquet")
    imputed = B.apply(lambda r: r.fillna(r.mean()), axis=1)
    imputed.index.name = "ProbeID"
    imputed.reset_index().to_csv(qc_dir / "probe_betas.imputed.csv", index=False)
    miss.astype(int).reset_index().rename(columns={"index": "ProbeID"}).to_csv(qc_dir / "probe_missing.csv", index=False)
    man[["sample_id", "age", "sex", "role"]].to_csv(qc_dir / "pheno.csv", index=False)
    print(f"probe matrix {imputed.shape[0]:,} probes x {imputed.shape[1]} samples; dropped {int((~keep).sum()):,} probes missing in >50% of samples")

    r = subprocess.run([a.rscript, "workflow/scripts/05b_qc_r.R", str(qc_dir)], capture_output=True, text=True)
    (qc_dir / "r_step.log").write_text(r.stdout + "\n" + r.stderr)
    if r.returncode != 0:
        sys.exit(f"R step failed; see {qc_dir / 'r_step.log'}\n{r.stderr[-3000:]}")
    comp = pd.read_csv(qc_dir / "r_composition.csv").rename(columns={"sample": "sample_id"})
    clocks = pd.read_csv(qc_dir / "r_clocks.csv").rename(columns={"sample": "sample_id"})

    q = pd.DataFrame(qcs).merge(man[["sample_id", "family_id", "role", "age", "sex", "flag_gene", "is_proband", "proband_family"]], on="sample_id")
    q = q.merge(comp, on="sample_id", how="left").merge(clocks, on="sample_id", how="left")
    q["pct_probes_imputed"] = q["sample_id"].map((miss.mean(axis=0) * 100).round(2))
    q = q.rename(columns={"global_mean_beta_depth_ge10": "global_mean_5mc"})
    # age acceleration = residual of DNAm age on chronological age, fit on all samples with both values
    clock_cols = [c for c in q.columns if c.startswith("dnam_age_")]
    for c in clock_cols:
        ok = q[c].notna() & q["age"].notna()
        if ok.sum() >= 10:
            X = np.c_[np.ones(ok.sum()), q.loc[ok, "age"]]
            beta = np.linalg.lstsq(X, q.loc[ok, c], rcond=None)[0]
            q.loc[ok, "age_accel_resid_" + c[9:]] = q.loc[ok, c] - X @ beta
            print(f"{c}: r(age) = {np.corrcoef(q.loc[ok, 'age'], q.loc[ok, c])[0, 1]:.3f} (n={int(ok.sum())})")
    primary = cfg.get("qc", {}).get("primary_clock", "Horvath")
    if f"dnam_age_{primary}" in q:
        q["dnam_age"] = q[f"dnam_age_{primary}"]; q["age_accel_resid"] = q.get(f"age_accel_resid_{primary}")
    # outliers
    med, mad = q["epi_frac"].median(), (q["epi_frac"] - q["epi_frac"].median()).abs().median() * 1.4826
    z5 = (q["global_mean_5mc"] - q["global_mean_5mc"].mean()) / q["global_mean_5mc"].std()
    reasons = []
    for i, r_ in q.iterrows():
        why = []
        if r_["sample_id"] in failed: why.append("extraction_failed(n_cpgs<50%median)")
        if pd.notna(r_["epi_frac"]) and r_["epi_frac"] > med + 3 * mad: why.append("epi_frac>median+3MAD")
        if r_["mean_depth"] < a.min_depth_sample: why.append(f"mean_depth<{a.min_depth_sample}")
        if r_.get("frac_cpgs_depth_ge10", 1) < cfg["thresholds"].get("min_frac_cpgs_10x", 0.70): why.append("frac_cpgs_10x<0.70(step04 acceptance)")
        if pd.notna(z5.iloc[i]) and abs(z5.iloc[i]) > 4: why.append("global_5mc_z>4")
        reasons.append(";".join(why))
    q["qc_outlier"] = [bool(x) for x in reasons]; q["qc_outlier_reason"] = reasons
    q.to_csv(qc_dir / "sample_qc.tsv", sep="\t", index=False)
    summ = q[["epi_frac", "fib_frac", "immune_frac", "mean_depth", "global_mean_5mc", "pct_probes_imputed"]].describe().round(3)
    print(summ.to_string()); print(f"qc_outliers: {int(q['qc_outlier'].sum())} / {len(q)}")
    (qc_dir / "STEP_LOG.md").write_text("# Step 05 QC summary (counts and distributions only)\n\n```\n" + summ.to_string() +
                                        f"\n```\n\nqc_outliers: {int(q['qc_outlier'].sum())} / {len(q)}\n" +
                                        "".join(f"\n- {c}: r(age) computed in run log" for c in clock_cols) + "\n")


if __name__ == "__main__":
    main()
