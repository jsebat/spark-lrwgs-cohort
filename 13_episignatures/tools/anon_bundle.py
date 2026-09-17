#!/usr/bin/env python
"""Build an anonymised bundle of summary tables for report writing off-cluster.

Replaces every sample_id with the label from results/report/sample_key.tsv and family_id with F01.. codes; drops the
raw identifiers. Only derived, aggregate-level tables are included (no per-CpG data, no genotypes, no BAM paths).
Run on Expanse from the repo root:  python tools/anon_bundle.py  -> results/report/anon/ and anon_bundle.tgz
"""
import pathlib
import shutil
import tarfile

import numpy as np
import pandas as pd

out = pathlib.Path("results/report/anon")
out.mkdir(exist_ok=True)
key = pd.read_csv("results/report/sample_key.tsv", sep="\t", dtype=str).set_index("sample_id")["label"]

qc = pd.read_csv("results/05_qc/sample_qc.tsv", sep="\t", dtype={"sample_id": str})
qc["label"] = qc.sample_id.map(key)
fam_codes = {f: f"F{i + 1:02d}" for i, f in enumerate(sorted(qc.family_id.unique()))}
qc["family_label"] = qc.family_id.map(fam_codes)
qc.drop(columns=["sample_id", "family_id"]).to_csv(out / "sample_qc.anon.tsv", sep="\t", index=False)

for src, name in [("results/08_panel/ranks.tsv", "ranks.anon.tsv"), ("results/06_layer1/scores.tsv", "scores.anon.tsv"),
                  ("results/07_layer2/features_resid.tsv", "features_resid.anon.tsv"), ("results/07_layer2/features.tsv", "features.anon.tsv")]:
    df = pd.read_csv(src, sep="\t", dtype={"sample_id": str})
    df["label"] = df.sample_id.map(key)
    df.drop(columns=["sample_id"]).to_csv(out / name, sep="\t", index=False)

for f in ["results/08_panel/proband_report_card.tsv", "results/09_stats/sensitivity.tsv", "results/09_stats/bootstrap.tsv",
          "results/03_manifest/STEP_LOG.md", "results/05_qc/STEP_LOG.md", "results/07_layer2/STEP_LOG.md", "results/09_stats/STEP_LOG.md",
          "results/06_layer1/region_info.tsv", "results/01_variant/annotation.json", "results/07_layer2/canyons.control.bed"]:
    shutil.copy(f, out / pathlib.Path(f).name)

# TBRS region-level: proband beta vs clean-children distribution (aggregates only)
rb = pd.read_parquet("results/06_layer1/region_betas.parquet")
info = pd.read_csv("results/06_layer1/region_info.tsv", sep="\t")
is_out = qc.qc_outlier.astype(str).str.lower().eq("true")
pro = qc.loc[qc.is_proband.astype(str).str.lower() == "true", "sample_id"].iloc[0]
clean_ch = qc[(qc.role == "child") & (qc.flag_gene.fillna("") == "") & ~is_out].sample_id
for setname in ["TBRS_LOF:tier1", "TBRS_LOF:tier2"]:
    ids = info[info.set == setname].region_id.tolist()
    C = rb.loc[clean_ch, ids]
    d = pd.DataFrame({"region_id": ids, "target": info.set_index("region_id").loc[ids, "target"].to_numpy(),
                      "ctrl_mean": C.mean().to_numpy(), "ctrl_q05": C.quantile(.05).to_numpy(), "ctrl_q95": C.quantile(.95).to_numpy(),
                      "ctrl_n": C.notna().sum().to_numpy(), "proband_beta": rb.loc[pro, ids].to_numpy()})
    d.to_csv(out / (setname.replace(":", "_") + "_regions.tsv"), sep="\t", index=False)

# NA fraction by threshold (step 06 acceptance diagnostics)
N = pd.read_parquet("results/06_layer1/region_ncpgs.parquet")
D = pd.read_parquet("results/06_layer1/region_depth.parquet")
ok = qc.loc[~is_out, "sample_id"]
N, D = N.loc[ok], D.loc[ok]
rows = []
for st in sorted(info.set.unique()):
    m = (info.set == st).to_numpy()
    ncpg = N.loc[:, m].replace(0, np.nan).median()
    for mc, md in ((10, 15), (10, 10), (5, 10)):
        na = ((N.loc[:, m] < mc) | (D.loc[:, m] < md)).mean(axis=1).median()
        rows.append(dict(set=st, n_regions=int(m.sum()), frac_lt10cpg=float((ncpg < 10).mean()), min_cpgs=mc, min_depth=md, median_na=float(na)))
pd.DataFrame(rows).to_csv(out / "na_by_threshold.tsv", sep="\t", index=False)

with tarfile.open("results/report/anon_bundle.tgz", "w:gz") as tf:
    tf.add(out, arcname="anon")
print("bundle files:", len(list(out.iterdir())))
