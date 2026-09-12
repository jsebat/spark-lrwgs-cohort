#!/usr/bin/env python
"""Publication figures from the anonymised bundle (results/report/anon) -> results/report/figures/*.png.
Every value is read from pipeline outputs. Run from results/report:  python ../../tools/report_figures.py
"""
import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

A = pathlib.Path("anon"); F = pathlib.Path("figures"); F.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 200, "savefig.dpi": 300})
SIGS = ["TBRS_LOF", "DNMT3A_R882", "DNMT3A_GOF", "DNMT1", "DNMT3B", "NSD1_SOTOS", "NEG_KMT2D", "NEG_CHD7"]
NICE = {"TBRS_LOF": "TBRS (DNMT3A LOF)", "DNMT3A_R882": "DNMT3A R882", "DNMT3A_GOF": "DNMT3A GOF (HESJAS)", "DNMT1": "DNMT1 (ADCA-DN)",
        "DNMT3B": "DNMT3B (ICF1)", "NSD1_SOTOS": "NSD1 (Sotos)", "NEG_KMT2D": "KMT2D (Kabuki)", "NEG_CHD7": "CHD7 (CHARGE)"}
RED, GREY, ORANGE, BLUE = "#c0392b", "#8c8c8c", "#e67e22", "#2c6fbb"

q = pd.read_csv(A / "sample_qc.anon.tsv", sep="\t")
q["qc_outlier"] = q.qc_outlier.astype(str).str.lower() == "true"; q["flag"] = q.flag_gene.fillna("") != ""
q["is_proband"] = q.is_proband.astype(str).str.lower() == "true"
pro = q[q.is_proband].label.iloc[0]
clean = q[(~q.flag) & (~q.qc_outlier)]
r = pd.read_csv(A / "ranks.anon.tsv", sep="\t")

# ---------------- Figure 1: gene model + deletion + protein domains
ann = json.load(open(A / "annotation.json"))
ex = ann["exons"]; m = ann["meta"]; a = [x for x in ann["analyses"] if x["convention"] == "bed"][0]
fig, axes = plt.subplots(2, 1, figsize=(7.2, 2.6), gridspec_kw={"height_ratios": [1, 1]})
ax = axes[0]
for e in ex:
    ax.add_patch(Rectangle((e["start"], 0.3), e["end"] - e["start"], 0.4, color=BLUE if e["coding"] else "#b0c4de", lw=0))
gs, ge = min(e["start"] for e in ex), max(e["end"] for e in ex)
ax.plot([gs, ge], [0.5, 0.5], color="k", lw=0.6, zorder=0)
d0, d1 = [int(v) for v in a["deleted_interval_1based"].split(":")[1].split("-")]
ax.add_patch(Rectangle((d0, 0.15), d1 - d0, 0.7, color=RED, alpha=0.9, lw=0))
ax.annotate(f"302-bp deletion, exons 14-15 (chr2:{d0:,}-{d1:,})", xy=(d0, 0.85), xytext=(d0 + 25000, 1.05), fontsize=7, color=RED, va="center",
            arrowprops=dict(arrowstyle="-", color=RED, lw=0.6))
ax.set_xlim(gs - 2000, ge + 2000); ax.set_ylim(0, 1.7); ax.set_yticks([]); ax.spines["left"].set_visible(False)
ax.set_xlabel("chr2 (GRCh38), minus strand; transcript " + m["transcript"], fontsize=7)
ax.text(gs, 1.5, "A  DNMT3A gene model (23 exons; coding exons dark)", fontsize=8, fontweight="bold")
ax.ticklabel_format(style="plain", axis="x"); ax.tick_params(labelsize=6)
ax = axes[1]
L = m["protein_length_aa"]
ax.add_patch(Rectangle((1, 0.35), L, 0.3, color="#dddddd", lw=0))
doms = {"PWWP domain": ("PWWP", "#7fb3d5", 0.5), "GATA1-like zinc finger": ("ADD (GATA)", "#f5b041", 0.24), "PHD zinc finger": ("ADD (PHD)", "#eb984e", 0.09), "C-5 cytosine methyltransferase": ("MTase", "#58d68d", 0.5)}
for dm in m["domains_all"]:
    for k, (lab, col, ty) in doms.items():
        if dm["db"] == "Pfam" and k in dm["description"]:
            ax.add_patch(Rectangle((dm["start"], 0.35), dm["end"] - dm["start"], 0.3, color=col, lw=0)); ax.text((dm["start"] + dm["end"]) / 2, ty, lab, ha="center", va="center", fontsize=6)
c1, c2 = [int(v) for v in a["codons_affected"].replace("p.", "").split("-")]
ax.add_patch(Rectangle((c1, 0.25), c2 - c1, 0.5, fill=False, edgecolor=RED, lw=1.2))
sc = a["scenario_fused_partial_exons"]
ax.annotate(f"frameshift at codon {c1}; stop after {sc['aberrant_residues']} aberrant residues\npredicted {sc['protein_length']}-aa product, NMD predicted", xy=(c1, 0.75), xytext=(c1 + 40, 1.2), fontsize=7, color=RED, arrowprops=dict(arrowstyle="-", color=RED, lw=0.6))
ax.axvline(882, color="k", lw=0.6, ls=":"); ax.text(882, 0.1, "R882", ha="center", fontsize=6)
ax.set_xlim(0, L + 10); ax.set_ylim(0, 1.6); ax.set_yticks([]); ax.spines["left"].set_visible(False); ax.set_xlabel("amino acid", fontsize=7); ax.tick_params(labelsize=6)
ax.text(0, 1.45, f"B  DNMT3A protein ({L} aa) with Pfam domains", fontsize=8, fontweight="bold")
plt.tight_layout(); fig.savefig(F / "fig1_variant.png"); plt.close(fig)

# ---------------- Figure 2: cohort QC
fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.4))
roles = ["child", "mother", "father"]
ax = axes[0]
for i, ro in enumerate(roles):
    g = q[(q.role == ro) & ~q.qc_outlier]; go = q[(q.role == ro) & q.qc_outlier]
    ax.scatter(np.random.default_rng(i).normal(i, 0.07, len(g)), g.epi_frac, s=9, color=GREY, alpha=.7)
    ax.scatter(np.random.default_rng(i + 9).normal(i, 0.07, len(go)), go.epi_frac, s=9, facecolors="none", edgecolors=ORANGE, lw=.8)
ax.scatter([0], q[q.is_proband].epi_frac, s=60, facecolors="none", edgecolors=RED, lw=1.5)
ax.set_xticks(range(3)); ax.set_xticklabels(["children", "mothers", "fathers"]); ax.set_ylabel("epithelial fraction (HEpiDISH)"); ax.set_title("A  Saliva composition", loc="left", fontsize=8, fontweight="bold")
ax = axes[1]
ax.scatter(q.mean_depth[~q.qc_outlier], q.frac_cpgs_depth_ge10[~q.qc_outlier], s=9, color=GREY, alpha=.7)
ax.scatter(q.mean_depth[q.qc_outlier], q.frac_cpgs_depth_ge10[q.qc_outlier], s=9, facecolors="none", edgecolors=ORANGE, lw=.8)
ax.scatter(q[q.is_proband].mean_depth, q[q.is_proband].frac_cpgs_depth_ge10, s=60, facecolors="none", edgecolors=RED, lw=1.5)
ax.axhline(0.70, color="k", ls=":", lw=.6); ax.set_xlim(0, 60); ax.set_xlabel("mean CpG depth (x)"); ax.set_ylabel("fraction of CpGs at >=10x"); ax.set_title("B  Coverage", loc="left", fontsize=8, fontweight="bold")
ax = axes[2]
ok = ~q.qc_outlier
for col, lab, colr in (("dnam_age_Horvath", "Horvath", BLUE), ("dnam_age_skinHorvath", "skin & blood", "#16a085"), ("dnam_age_PedBE", "PedBE", ORANGE)):
    mm = ok & q[col].notna(); rr = np.corrcoef(q.age[mm], q[col][mm])[0, 1]
    ax.scatter(q.age[mm], q[col][mm], s=7, color=colr, alpha=.6, label=f"{lab} (r = {rr:.2f})")
ax.plot([0, 80], [0, 80], "k:", lw=.6); ax.set_xlabel("chronological age (y)"); ax.set_ylabel("DNAm age (y)"); ax.legend(fontsize=6, frameon=False); ax.set_title("C  Epigenetic clocks", loc="left", fontsize=8, fontweight="bold")
plt.tight_layout(); fig.savefig(F / "fig2_qc.png"); plt.close(fig)

# ---------------- Figure 3: heatmap + radar
W = r[(r.metric == "shape_r") & (r.null == "children") & (r.tier == "tier2")].pivot(index="label", columns="signature", values="value")[SIGS]
order = pd.concat([q[q.role == "child"].sort_values("label").label, q[q.role != "child"].sort_values("label").label]); W = W.loc[[l for l in order if l in W.index]]
fig = plt.figure(figsize=(7.2, 6.0)); gs_ = fig.add_gridspec(1, 2, width_ratios=[1.15, 1])
ax = fig.add_subplot(gs_[0, 0])
v = np.nanmax(np.abs(W.values)); im = ax.imshow(W.values, aspect="auto", cmap="RdBu_r", vmin=-v, vmax=v)
ax.set_xticks(range(len(SIGS))); ax.set_xticklabels([NICE[s] for s in SIGS], rotation=60, ha="right", fontsize=6)
ax.set_yticks([]); i_p = list(W.index).index(pro); ax.add_patch(Rectangle((-0.5, i_p - 0.5), len(SIGS), 1, fill=False, edgecolor=RED, lw=1.5))
ax.text(len(SIGS) - 0.4, i_p, "proband", color=RED, va="center", fontsize=7)
nch = (q.role == "child").sum(); ax.axhline(nch - 0.5, color="k", lw=.5); ax.text(-0.6, nch / 2, "children", rotation=90, va="center", ha="right", fontsize=7); ax.text(-0.6, nch + (len(W) - nch) / 2, "parents", rotation=90, va="center", ha="right", fontsize=7)
for i, l in enumerate(W.index):
    row = q.set_index("label").loc[l]
    if row.flag: ax.plot(-0.5, i, marker=">", color=ORANGE, ms=3, clip_on=False)
    if row.qc_outlier: ax.plot(len(SIGS) - 0.5, i, marker="<", color="k", ms=3, clip_on=False)
cb = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02); cb.set_label("shape r", fontsize=7); cb.ax.tick_params(labelsize=6)
ax.set_title("A  Shape correlation, all samples x panel", loc="left", fontsize=8, fontweight="bold")
ax = fig.add_subplot(gs_[0, 1], polar=True)
ang = np.linspace(0, 2 * np.pi, len(SIGS), endpoint=False).tolist(); ang += ang[:1]
ctrl = W.loc[[l for l in W.index if l in set(clean[clean.role == "child"].label)]]
for qq, st, lab in ((0.05, ":", "clean children 5th / 95th pct"), (0.5, "--", "clean children median"), (0.95, ":", None)):
    vals = ctrl.quantile(qq).tolist(); ax.plot(ang, vals + vals[:1], color="k", ls=st, lw=.8, label=lab)
vals = W.loc[pro].tolist(); ax.plot(ang, vals + vals[:1], color=RED, lw=1.8, label="proband"); ax.fill(ang, vals + vals[:1], color=RED, alpha=.12)
ax.set_xticks(ang[:-1]); ax.set_xticklabels([NICE[s].replace(" (", "\n(") for s in SIGS], fontsize=6); ax.tick_params(axis="y", labelsize=6)
ax.legend(fontsize=6, loc="upper right", bbox_to_anchor=(1.35, 1.15), frameon=False); ax.set_title("B  Proband vs clean children", loc="left", fontsize=8, fontweight="bold", pad=18)
plt.tight_layout(); fig.savefig(F / "fig3_panel.png"); plt.close(fig)

# ---------------- Figure 4: rank plots
sub = r[(r.metric == "shape_r") & (r.null == "children") & (r.tier == "tier2") & (r.role == "child")]
fig, axes = plt.subplots(2, 4, figsize=(7.2, 3.8), sharey=False); axes = axes.ravel()
for ax, sig in zip(axes, SIGS):
    g = sub[sub.signature == sig].sort_values("value", ascending=False).reset_index(drop=True)
    cols = [RED if p_ else (ORANGE if (f_ or o_) else GREY) for p_, f_, o_ in zip(g.is_proband, g.flag_gene.fillna("") != "", g.qc_outlier.astype(str).str.lower() == "true")]
    ax.bar(range(len(g)), g.value, color=cols, width=0.85); ax.axhline(0, color="k", lw=.5)
    pr_ = g[g.is_proband]
    ax.set_title(NICE[sig], fontsize=7); ax.set_xticks([]); ax.tick_params(labelsize=6)
    if len(pr_): ax.text(0.98, 0.95, f"rank {int(pr_.rank_children.iloc[0])}/{int(pr_.n_children.iloc[0])}", transform=ax.transAxes, ha="right", va="top", fontsize=6.5, color=RED)
axes[0].set_ylabel("shape r"); axes[4].set_ylabel("shape r")
fig.text(0.5, 0.01, "children ranked by shape r (red = proband; orange = variant-flagged or QC outlier; grey = clean)", ha="center", fontsize=7)
plt.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(F / "fig4_ranks.png"); plt.close(fig)

# ---------------- Figure 5: TBRS tier1 region deltas
t1 = pd.read_csv(A / "TBRS_LOF_tier1_regions.tsv", sep="\t").dropna(subset=["proband_beta", "ctrl_mean"]).sort_values("target")
x = np.arange(len(t1))
fig, ax = plt.subplots(figsize=(7.2, 2.6))
ax.fill_between(x, t1.ctrl_q05 - t1.ctrl_mean, t1.ctrl_q95 - t1.ctrl_mean, color="#d5d8dc", label="clean children, 5th-95th percentile")
ax.plot(x, t1.proband_beta - t1.ctrl_mean, ".", color=RED, ms=4, label="proband")
ax.axhline(0, color="k", lw=.6); ax.set_xlabel("TBRS tier-1 regions (sorted by published delta-beta)"); ax.set_ylabel("beta - control mean")
nb = int((t1.proband_beta < t1.ctrl_q05).sum()); ax.set_title(f"{nb} of {len(t1)} scored tier-1 regions below the children 5th percentile", fontsize=8, loc="left")
ax.legend(fontsize=6.5, frameon=False, loc="lower left"); plt.tight_layout(); fig.savefig(F / "fig5_tbrs_regions.png"); plt.close(fig)

# ---------------- Figure 6: Layer 2 residuals
fr = pd.read_csv(A / "features_resid.anon.tsv", sep="\t").set_index("label")
feats = ["enhancer_beta", "promoter_beta", "cgi_beta", "genebody_beta", "pmd_soloWCGW_beta", "hmd_soloWCGW_beta", "canyon_interior_beta", "canyon_edge_slope", "canyon_total_bp", "canyon_mean_width_bp", "satII_III_beta", "alpha_sat_beta", "line1_beta", "line1_young_beta", "global_mean_5mc", "age_accel_resid"]
labels = ["enhancer 5mC", "promoter 5mC", "CGI 5mC", "gene-body 5mC", "PMD solo-WCGW", "HMD 5mC", "canyon interior", "canyon edge slope", "canyon total bp", "canyon width", "SatII/III", "alpha satellite", "LINE-1", "young LINE-1", "global 5mC", "age accel. (Horvath)"]
ctrl_l = [l for l in fr.index if l in set(clean.label)]
# age acceleration is in years in the table; express it in control SD units like the other features
mu_a, sd_a = fr.loc[ctrl_l, "age_accel_resid"].mean(), fr.loc[ctrl_l, "age_accel_resid"].std()
fr["age_accel_resid"] = (fr["age_accel_resid"] - mu_a) / sd_a
fig, ax = plt.subplots(figsize=(7.2, 2.8))
ax.boxplot([fr.loc[ctrl_l, f].dropna() for f in feats], positions=range(len(feats)), widths=.5, showfliers=False, medianprops=dict(color="k"))
ax.plot(range(len(feats)), fr.loc[pro, feats], "o", color=RED, ms=5, label="proband")
ax.axhline(0, color="k", lw=.5); ax.axhline(2, color=GREY, lw=.5, ls=":"); ax.axhline(-2, color=GREY, lw=.5, ls=":")
ax.set_xticks(range(len(feats))); ax.set_xticklabels(labels, rotation=55, ha="right", fontsize=6.5); ax.set_ylabel("standardized residual\n(vs clean controls)")
ax.legend(fontsize=7, frameon=False); plt.tight_layout(); fig.savefig(F / "fig6_layer2.png"); plt.close(fig)

# ---------------- Figure 7: robustness
s = pd.read_csv(A / "sensitivity.tsv", sep="\t"); s = s[s.null == "children"]
lab = s.analysis.str.replace(r" \(.*", "", regex=True).str.replace("step06:", "step 06 ").str.replace("composition_clean", "composition-clean").str.replace("bootstrap_controls", "bootstrap (n=200)")
fig, ax = plt.subplots(figsize=(7.2, 2.6))
ax.barh(range(len(s)), s.shape_r, color=[RED if rk == 1 else GREY for rk in s.rank_children]); ax.set_yticks(range(len(s))); ax.set_yticklabels(lab, fontsize=6.5); ax.invert_yaxis()
for i, (v, rk, n) in enumerate(zip(s.shape_r, s.rank_children, s.n_children)):
    ax.text(v + 0.005, i, f"rank {int(rk)}/{int(n)}", va="center", fontsize=6.5)
ax.set_xlabel("proband TBRS shape r"); ax.set_xlim(0, max(s.shape_r) + 0.1); plt.tight_layout(); fig.savefig(F / "fig7_robustness.png"); plt.close(fig)
print("figures:", sorted(p.name for p in F.iterdir()))
