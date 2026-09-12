#!/usr/bin/env python
"""Step 10 - assemble results/report/report.html (self-contained; matplotlib figures embedded as PNG).

Every number and figure is read from pipeline outputs; nothing is typed in. Samples are shown under anonymised labels
(S001.., with role) - the label->sample_id key is written to results/report/sample_key.tsv (Expanse only, gitignored).
Sections follow PLAN.md step 10.
"""
import argparse
import base64
import html
import io
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml


def fig_to_html(fig, width=900):
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=110, bbox_inches="tight"); plt.close(fig)
    return f'<img style="max-width:{width}px;width:100%" src="data:image/png;base64,{base64.b64encode(buf.getvalue()).decode()}">'


def table(df, max_rows=200):
    return df.head(max_rows).to_html(index=False, border=0, classes="tbl", na_rep="", float_format=lambda x: f"{x:.3f}")


def md_to_html(text):
    out = []
    for line in text.splitlines():
        if line.startswith("|"):
            cells = [c.strip() for c in line.strip("|").split("|")]
            if set("".join(cells)) <= set("-: "):
                continue
            out.append("<tr>" + "".join(f"<td>{html.escape(c)}</td>" for c in cells) + "</tr>")
        elif line.startswith("#"):
            lvl = len(line) - len(line.lstrip("#")); out.append(f"<h{lvl + 2}>{html.escape(line.lstrip('# '))}</h{lvl + 2}>")
        elif line.strip():
            out.append(f"<p>{html.escape(line)}</p>")
    s = "\n".join(out)
    return s.replace("</tr>\n<tr>", "</tr><tr>").replace("<tr>", '<table class="tbl"><tr>', 1).replace("</tr>", "</tr></table>", s.count("<tr>")) if "<tr>" in s else s


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", required=True)
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config))
    out = pathlib.Path(cfg["paths"]["output_root"]); rd = out / "report"; rd.mkdir(parents=True, exist_ok=True)
    qc = pd.read_csv(out / "05_qc" / "sample_qc.tsv", sep="\t", dtype={"sample_id": str})
    qc = qc.sort_values(["role", "sample_id"]).reset_index(drop=True)
    qc["label"] = [f"S{i + 1:03d}-{r[:1].upper()}" for i, r in enumerate(qc["role"])]
    qc[["label", "sample_id"]].to_csv(rd / "sample_key.tsv", sep="\t", index=False)
    lab = qc.set_index("sample_id")["label"]
    proband = qc.loc[qc["is_proband"].astype(str).str.lower() == "true", "sample_id"].tolist()
    flagged = qc.loc[qc["flag_gene"].fillna("") != "", "sample_id"].tolist()
    ranks = pd.read_csv(out / "08_panel" / "ranks.tsv", sep="\t", dtype={"sample_id": str})
    card = pd.read_csv(out / "08_panel" / "proband_report_card.tsv", sep="\t")
    sens = pd.read_csv(out / "09_stats" / "sensitivity.tsv", sep="\t")
    sigs = cfg["signatures"]["layer1"]
    S = ["<style>body{font-family:system-ui,sans-serif;max-width:1100px;margin:auto;padding:1em;color:#222} .tbl{border-collapse:collapse;font-size:12px} .tbl td,.tbl th{border:1px solid #ccc;padding:2px 6px} h2{border-bottom:2px solid #444;margin-top:2em} .pass{color:#1a7f37;font-weight:bold} .fail{color:#b91c1c;font-weight:bold} .warn{background:#fff7ed;padding:.5em;border-left:4px solid #f97316}</style>",
         "<h1>DNMT3A episignature panel - SPARK lrWGS pilot</h1>",
         f"<p>Generated from pipeline outputs under <code>{out}</code>. Samples are anonymised (key: results/report/sample_key.tsv). "
         f"n = {len(qc)} samples ({int((qc.role == 'child').sum())} children, {int(qc.role.isin(['mother', 'father']).sum())} parents); "
         f"{len(flagged)} flagged; {int(qc.qc_outlier.astype(str).str.lower().eq('true').sum())} QC outliers.</p>"]

    # 1 variant
    S.append("<h2>1. Variant annotation</h2>")
    for f in ("CHECKPOINT.md", "annotation.md"):
        p = out / "01_variant" / f
        if p.exists():
            S.append(f"<details {'open' if f == 'CHECKPOINT.md' else ''}><summary>{f}</summary>{md_to_html(p.read_text(encoding='utf-8'))}</details>")
    # 2 inventory
    S.append("<h2>2. Signature inventory</h2>" + md_to_html(pathlib.Path("signatures/INVENTORY.md").read_text(encoding="utf-8")))
    # 3 QC
    S.append("<h2>3. Cohort QC</h2>")
    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    for c, col in zip(ax, ["epi_frac", "mean_depth", "global_mean_5mc"]):
        for role, g in qc.groupby("role"):
            c.scatter(np.random.default_rng(0).normal(["child", "father", "mother", "other"].index(role) if role in ["child", "father", "mother", "other"] else 3, 0.05, len(g)), g[col], s=14, label=role, alpha=.7)
        pb = qc[qc.sample_id.isin(proband)]
        if len(pb):
            c.scatter([0] * len(pb), pb[col], s=90, facecolors="none", edgecolors="red", linewidths=2, label="proband")
        c.set_title(col); c.set_xticks([])
    ax[0].legend(fontsize=8); S.append(fig_to_html(fig, 1000))
    clocks = [c for c in qc.columns if c.startswith("dnam_age_")]
    if clocks:
        fig, ax = plt.subplots(1, len(clocks), figsize=(4 * len(clocks), 4), squeeze=False)
        for c, col in zip(ax[0], clocks):
            ok = qc[col].notna() & qc.age.notna()
            r = np.corrcoef(qc.age[ok], qc[col][ok])[0, 1] if ok.sum() > 2 else np.nan
            c.scatter(qc.age[ok], qc[col][ok], s=14); c.plot([0, 80], [0, 80], "k--", lw=.7)
            pb = qc[qc.sample_id.isin(proband) & ok]
            c.scatter(pb.age, pb[col], s=90, facecolors="none", edgecolors="red", linewidths=2)
            c.set_title(f"{col[9:]}  r={r:.2f}"); c.set_xlabel("chronological age"); c.set_ylabel("DNAm age")
        S.append(fig_to_html(fig, 1100))
    S.append(table(qc[["label", "role", "age", "sex", "mean_depth", "epi_frac", "immune_frac", "global_mean_5mc", "dnam_age", "age_accel_resid", "pct_probes_imputed", "flag_gene", "qc_outlier"]].round(3)))

    # 4 heatmap shape_r
    S.append("<h2>4. Samples x signatures: shape_r (children null, tier2)</h2>")
    W = ranks[(ranks.metric == "shape_r") & (ranks.null == "children") & (ranks.tier == "tier2")].pivot(index="sample_id", columns="signature", values="value").reindex(columns=[s for s in sigs if s in ranks.signature.unique()])
    W = W.loc[qc.sample_id[qc.sample_id.isin(W.index)]]
    fig, ax = plt.subplots(figsize=(8, max(6, 0.16 * len(W))))
    v = np.nanmax(np.abs(W.values)) if W.notna().values.any() else 1
    im = ax.imshow(W.values, aspect="auto", cmap="RdBu_r", vmin=-v, vmax=v)
    ax.set_xticks(range(W.shape[1])); ax.set_xticklabels(W.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(W))); ax.set_yticklabels([lab[s] + (" *PROBAND*" if s in proband else " [flag]" if s in flagged else "") for s in W.index], fontsize=6)
    for i, s in enumerate(W.index):
        if s in proband:
            ax.add_patch(plt.Rectangle((-0.5, i - 0.5), W.shape[1], 1, fill=False, edgecolor="red", lw=2))
    plt.colorbar(im, ax=ax, label="shape_r"); S.append(fig_to_html(fig, 800))

    # 5 radar
    S.append("<h2>5. Proband radar across the panel (shape_r, children null, tier2)</h2>")
    if proband and len(W.columns):
        ang = np.linspace(0, 2 * np.pi, len(W.columns), endpoint=False).tolist(); ang += ang[:1]
        fig = plt.figure(figsize=(6, 6)); ax = plt.subplot(111, polar=True)
        ctrl = W.loc[[s for s in W.index if s not in proband and s not in flagged and qc.set_index("sample_id").loc[s, "role"] == "child"]]
        for q_, st in ((0.05, ":"), (0.5, "--"), (0.95, ":")):
            vals = ctrl.quantile(q_).tolist(); ax.plot(ang, vals + vals[:1], "k" + st, lw=.8, label=f"children q{int(q_ * 100)}")
        vals = W.loc[proband[0]].tolist(); ax.plot(ang, vals + vals[:1], "r-", lw=2, label="proband"); ax.fill(ang, vals + vals[:1], "r", alpha=.15)
        ax.set_xticks(ang[:-1]); ax.set_xticklabels(W.columns, fontsize=8); ax.legend(fontsize=7, loc="lower right", bbox_to_anchor=(1.3, -0.1))
        S.append(fig_to_html(fig, 600))

    # 6 rank plots
    S.append("<h2>6. Rank plots per signature (shape_r among children; proband in red; empirical p = rank/n)</h2>")
    sub = ranks[(ranks.metric == "shape_r") & (ranks.null == "children") & (ranks.tier == "tier2") & (ranks.role == "child")]
    fig, axes = plt.subplots(2, 4, figsize=(16, 7)); axes = axes.ravel()
    for axx, sig in zip(axes, [s for s in sigs if s in sub.signature.unique()]):
        g = sub[sub.signature == sig].sort_values("value", ascending=False).reset_index(drop=True)
        cols = ["red" if s in proband else "orange" if s in flagged else "grey" for s in g.sample_id]
        axx.bar(range(len(g)), g.value, color=cols); axx.set_title(sig, fontsize=10); axx.set_xlabel("children, ranked"); axx.set_ylabel("shape_r")
        pr = g[g.sample_id.isin(proband)]
        if len(pr):
            axx.text(0.98, 0.95, f"proband rank {int(pr.rank_children.iloc[0])}/{int(pr.n_children.iloc[0])}\np={pr.rank_children.iloc[0] / pr.n_children.iloc[0]:.3f}", transform=axx.transAxes, ha="right", va="top", fontsize=8)
    S.append(fig_to_html(fig, 1100))

    # 7 TBRS tier1 delta-beta heatmap
    S.append("<h2>7. TBRS_LOF tier1 regions: proband delta-beta vs control distribution</h2>")
    rb = pd.read_parquet(out / "06_layer1" / "region_betas.parquet"); info = pd.read_csv(out / "06_layer1" / "region_info.tsv", sep="\t")
    t1 = info[info.set == "TBRS_LOF:tier1"].region_id.tolist()
    if t1 and proband:
        ctrl_ids = [s for s in rb.index if s not in proband and s not in flagged and qc.set_index("sample_id").loc[s, "role"] == "child"]
        C = rb.loc[ctrl_ids, t1]; mu = C.mean(); d_pro = rb.loc[proband[0], t1] - mu; d_ctrl = C - mu
        order = np.argsort(info.set_index("region_id").loc[t1, "target"].to_numpy())
        fig, ax = plt.subplots(figsize=(12, 4))
        x = np.arange(len(t1))
        ax.fill_between(x, d_ctrl.quantile(.05).to_numpy()[order], d_ctrl.quantile(.95).to_numpy()[order], color="lightgrey", label="children 5-95%")
        ax.plot(x, d_pro.to_numpy()[order], "r.", ms=4, label="proband delta-beta"); ax.axhline(0, color="k", lw=.5)
        ax.set_xlabel("tier1 regions (sorted by published delta)"); ax.set_ylabel("beta - control mean"); ax.legend(fontsize=8)
        ax.set_title(f"proband: {int((d_pro < d_ctrl.quantile(.05)).sum())} of {len(t1)} tier1 regions below the children 5th percentile")
        S.append(fig_to_html(fig, 1100))

    # 8 Layer 2
    S.append("<h2>8. Layer 2 standardized residuals (covariate-adjusted)</h2>")
    l2p = out / "07_layer2" / "features_resid.tsv"
    if l2p.exists():
        fr = pd.read_csv(l2p, sep="\t", dtype={"sample_id": str}).set_index("sample_id")
        cols = [c for c in fr.columns if c != "sample_id"]
        ctrl = fr.loc[[s for s in fr.index if s not in proband and s not in flagged]]
        fig, ax = plt.subplots(figsize=(13, 4))
        ax.boxplot([ctrl[c].dropna() for c in cols], positions=range(len(cols)), widths=.5, showfliers=False)
        if proband:
            ax.plot(range(len(cols)), fr.loc[proband[0], cols], "ro", label="proband")
        ax.axhline(0, color="k", lw=.5); ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=60, ha="right", fontsize=8); ax.set_ylabel("std. residual"); ax.legend()
        S.append(fig_to_html(fig, 1100))
        S.append("<details><summary>Step 07 log</summary><pre>" + html.escape((out / "07_layer2" / "STEP_LOG.md").read_text()) + "</pre></details>")

    # 9 predictions
    S.append("<h2>9. Pre-registered predictions: observed vs predicted</h2>")
    c2 = card.copy(); c2["verdict"] = c2["verdict"].map(lambda v: f'<span class="{v.lower() if v in ("PASS", "FAIL") else ""}">{v}</span>')
    S.append(c2.to_html(index=False, border=0, classes="tbl", escape=False, na_rep=""))
    S.append("<h3>Sensitivity (step 09)</h3>" + table(sens.round(4)))
    go = out / "08_panel" / "global_outliers.tsv"
    if go.exists() and go.stat().st_size > 0:
        g = pd.read_csv(go, sep="\t", dtype={"sample_id": str})
        if len(g):
            g["sample"] = g.sample_id.map(lab); S.append('<div class="warn">Global outliers (top-5 on >= 4 signatures):</div>' + table(g.drop(columns=["sample_id"])))

    # 10 limitations
    S.append("<h2>10. Limitations</h2><ul>"
             "<li>Episignatures were derived from blood (array) or fibroblast (RRBS) data and are applied here to saliva, whose epithelial fraction varies by sample; composition is a covariate in Layer 2 and a sensitivity analysis in Layer 1, not a correction of the signatures themselves.</li>"
             "<li>~20x HiFi coverage gives noisy single-CpG calls; all statistics are region-level with >=10 CpGs and >=15x mean depth, and regions failing that are NA (counted, never imputed).</li>"
             "<li>Small cohort: the primary null has only the non-flagged children; empirical p-values are bounded below by 1/n.</li>"
             "<li>The EpiSign classifiers themselves are not reproduced; we test similarity of direction and magnitude to published probe sets.</li>"
             "<li>Rare-variant flags use cohort-internal allele counts (no population frequencies offline) and canonical-transcript consequences only.</li>"
             "<li>Genome build and coordinate conventions of supplementary tables were verified where possible (see signatures/INVENTORY.md); Smith 2021 build was inferred and checked against CpG positions.</li></ul>")
    (rd / "report.html").write_text("\n".join(S), encoding="utf-8")
    print(f"wrote {rd / 'report.html'}")


if __name__ == "__main__":
    main()
