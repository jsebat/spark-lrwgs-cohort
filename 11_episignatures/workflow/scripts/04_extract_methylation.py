#!/usr/bin/env python
"""Step 04 - per-sample CpG methylation from HiFi BAM MM/ML tags via pb-CpG-tools.

For one sample: run `aligned_bam_to_cpg_scores` on each of its BAMs (one merged BAM, or several per-movie BAMs),
then combine across BAMs by summing estimated modified / unmodified counts per CpG, and write a parquet:
    chrom, pos (0-based CpG start, hg38), beta, depth
where beta = est_mod / (est_mod + est_unmod) and depth = total coverage. Sites with depth 0 are dropped.

Combined (haplotype-agnostic) scores only; haplotype-specific outputs are kept alongside if the tool produced them.
Also writes <out>.qc.json: n CpGs, mean depth, fraction of CpGs with depth >= 10 (step 04 acceptance), tool version.
Fallback (config.tools.fallback = modkit pileup) is used only if pb-CpG-tools is absent.
"""
import argparse
import glob
import json
import os
import pathlib
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
import yaml

# pb-CpG-tools combined.bed columns (v2.x / v3.x): chrom start end mod_score hap coverage est_mod est_unmod [discretized]
PB_COLS = ["chrom", "start", "end", "mod_score", "hap", "coverage", "est_mod", "est_unmod"]
MAIN_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def sh(cmd, log):
    with open(log, "a") as fh:
        fh.write(f"\n$ {cmd}\n")
        r = subprocess.run(cmd, shell=True, stdout=fh, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        sys.exit(f"command failed (see {log}): {cmd}")


def tool_version(exe):
    try:
        return subprocess.run([exe, "--version"], capture_output=True, text=True).stdout.strip()
    except FileNotFoundError:
        return None


def run_pbcpg(bam, prefix, ref, threads, cfg, log):
    exe = cfg["tools"]["cpg_scores"]
    ver = tool_version(exe)
    if ver is None:
        return None
    args = cfg["tools"].get("cpg_scores_args", "")
    # v3 takes --ref; v2 does not (model-based pileup needs the bundled model path in v2 -> --model)
    extra = f"--ref {ref}" if "3." in ver else ""
    model = cfg["tools"].get("cpg_scores_model")
    if model and "2." in ver:
        extra += f" --model {model}"
    sh(f"{exe} --bam {bam} --output-prefix {prefix} --threads {threads} {args} {extra}", log)
    beds = sorted(b for b in glob.glob(f"{prefix}*combined*.bed*") if not b.endswith(".bw"))
    if not beds:
        sys.exit(f"pb-CpG-tools produced no combined bed for prefix {prefix}")
    return ver, beds[0]


def read_pbcpg_bed(path):
    comp = "gzip" if str(path).endswith(".gz") else None   # explicit: extension inference failed on Expanse
    df = pd.read_csv(path, sep="\t", header=None, comment="#", low_memory=False, compression=comp)
    if df.shape[1] < 8:
        sys.exit(f"{path}: expected >=8 columns, got {df.shape[1]}")
    df = df.iloc[:, :8]; df.columns = PB_COLS
    df = df[df["chrom"].isin(MAIN_CHROMS)]
    return df[["chrom", "start", "coverage", "est_mod", "est_unmod"]].astype({"start": np.int64, "coverage": np.float64, "est_mod": np.float64, "est_unmod": np.float64})


def run_modkit(bam, prefix, ref, threads, cfg, log):
    exe = "modkit"
    ver = tool_version(exe)
    if ver is None:
        sys.exit("neither pb-CpG-tools nor modkit is available")
    out = f"{prefix}.modkit.bed"
    sh(f"{exe} pileup --cpg --combine-strands --ref {ref} --threads {threads} {bam} {out}", log)
    # bedMethyl: chrom start end code score strand ... Nvalid_cov percent_mod Nmod Ncanonical ...
    df = pd.read_csv(out, sep="\t", header=None, low_memory=False)
    df = df[df[3] == "m"]
    d = pd.DataFrame({"chrom": df[0], "start": df[1].astype(np.int64), "coverage": df[9].astype(float),
                      "est_mod": df[11].astype(float), "est_unmod": df[12].astype(float)})
    return ver, d[d["chrom"].isin(MAIN_CHROMS)]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--bams", required=True, help="semicolon-separated BAM paths")
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--workdir", default=None, help="scratch dir for tool outputs (default: config.paths.scratch_dir/04/<sample>)")
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config))
    ref = cfg["paths"]["reference_fasta"]
    bams = [b for b in a.bams.split(";") if b]
    if not bams:
        sys.exit(f"{a.sample}: no BAMs")
    work = pathlib.Path(a.workdir or os.path.join(cfg["paths"]["scratch_dir"], "04", a.sample)); work.mkdir(parents=True, exist_ok=True)
    log = work / "tool.log"
    out = pathlib.Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)

    # some SPARK staging movie BAMs are unaligned (no @SQ in header) despite having a .bai next to them: skip those
    aligned, skipped = [], []
    for bam in bams:
        hdr = subprocess.run(["samtools", "view", "-H", bam], capture_output=True, text=True).stdout
        (aligned if "\n@SQ\t" in hdr or hdr.startswith("@SQ\t") else skipped).append(bam)
    if skipped:
        print(f"{a.sample}: skipping {len(skipped)} unaligned BAM(s) of {len(bams)}", file=sys.stderr)
    if not aligned:
        sys.exit(f"{a.sample}: no aligned BAM among {len(bams)} inputs")
    bams = aligned
    parts, versions, mode = [], set(), None
    for i, bam in enumerate(bams):
        prefix = str(work / f"part{i}")
        res = run_pbcpg(bam, prefix, ref, a.threads, cfg, log)
        if res is None:
            ver, d = run_modkit(bam, prefix, ref, a.threads, cfg, log); mode = "modkit"
        else:
            ver, bed = res; d = read_pbcpg_bed(bed); mode = "pb-CpG-tools"
        versions.add(ver); parts.append(d)
    d = pd.concat(parts, ignore_index=True)
    if len(parts) > 1:
        d = d.groupby(["chrom", "start"], as_index=False, sort=False)[["coverage", "est_mod", "est_unmod"]].sum()
    tot = d["est_mod"] + d["est_unmod"]
    d = d[tot > 0].copy()
    d["beta"] = (d["est_mod"] / (d["est_mod"] + d["est_unmod"])).clip(0, 1).astype(np.float32)
    d["depth"] = d["coverage"].round().astype(np.int32)
    d = d.rename(columns={"start": "pos"})[["chrom", "pos", "beta", "depth"]]
    d["chrom"] = pd.Categorical(d["chrom"], categories=MAIN_CHROMS)
    d = d.sort_values(["chrom", "pos"]).reset_index(drop=True)
    d.to_parquet(out, index=False, compression="zstd")
    qc = dict(sample=a.sample, n_bams=len(bams), n_bams_skipped_unaligned=len(skipped), tool=mode, tool_version=sorted(v for v in versions if v),
              n_cpgs=int(len(d)), mean_depth=float(d["depth"].mean()), median_depth=float(d["depth"].median()),
              frac_cpgs_depth_ge10=float((d["depth"] >= 10).mean()), frac_cpgs_depth_ge15=float((d["depth"] >= 15).mean()),
              global_mean_beta=float(d["beta"].mean()), global_mean_beta_depth_ge10=float(d.loc[d["depth"] >= 10, "beta"].mean()))
    pathlib.Path(str(out) + ".qc.json").write_text(json.dumps(qc, indent=1))
    print(json.dumps({k: v for k, v in qc.items() if k != "sample"}))
    if not cfg.get("keep_tool_output", False):
        for f in work.glob("part*"):
            f.unlink()


if __name__ == "__main__":
    main()
