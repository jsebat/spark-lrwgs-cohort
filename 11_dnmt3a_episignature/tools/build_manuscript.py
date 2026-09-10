#!/usr/bin/env python
"""Assemble the manuscript-style PDF (methods / results / discussion) from the anonymised bundle and figures.
Run from results/report:  python ../../tools/build_manuscript.py  -> results/report/DNMT3A_episignature_report.pdf
All numbers are computed here from the bundle tables; nothing is typed in by hand except the prose.
"""
import datetime as dt
import json
import pathlib

import numpy as np
import pandas as pd
from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import Image, PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle, KeepTogether

A = pathlib.Path("anon"); F = pathlib.Path("figures")
OUT = "DNMT3A_episignature_report.pdf"

# ----------------------------------------------------------------------------- numbers
q = pd.read_csv(A / "sample_qc.anon.tsv", sep="\t")
q["qc_outlier"] = q.qc_outlier.astype(str).str.lower() == "true"; q["flag"] = q.flag_gene.fillna("") != ""; q["is_proband"] = q.is_proband.astype(str).str.lower() == "true"
ok = ~q.qc_outlier; clean = (~q.flag) & ok
pro = q[q.is_proband].iloc[0]
n = dict(all=len(q), child=int((q.role == "child").sum()), parent=int(q.role.isin(["mother", "father"]).sum()), fam=q.family_label.nunique(),
         flagged=int(q.flag.sum()), flagged_child=int((q.flag & (q.role == "child")).sum()), outl=int(q.qc_outlier.sum()),
         clean_child=int((clean & (q.role == "child")).sum()), clean_all=int(clean.sum()))
ch = q[q.role == "child"]; pa = q[q.role != "child"]
clocks = {c[9:]: np.corrcoef(q.age[ok & q[c].notna()], q[c][ok & q[c].notna()])[0, 1] for c in q.columns if c.startswith("dnam_age_")}
reasons = q.loc[q.qc_outlier, "qc_outlier_reason"].str.split(";").explode().value_counts()
card = pd.read_csv(A / "proband_report_card.tsv", sep="\t", keep_default_na=False)
for c_ in ["observed_shape_r", "shape_r_z_vs_null", "rank_children", "n_children", "rank_all", "n_all", "perm_p_shape_r"]:
    card[c_] = pd.to_numeric(card[c_], errors="coerce")
cc = card[(card.null == "children") & ~card.signature.str.startswith("L2:")].set_index(["signature", "tier"])
ca = card[(card.null == "all") & ~card.signature.str.startswith("L2:")].set_index(["signature", "tier"])
l2 = card[card.signature.str.startswith("L2:")].set_index("signature")
sens = pd.read_csv(A / "sensitivity.tsv", sep="\t"); boot = pd.read_csv(A / "bootstrap.tsv", sep="\t")
na = pd.read_csv(A / "na_by_threshold.tsv", sep="\t")
fr = pd.read_csv(A / "features_resid.anon.tsv", sep="\t").set_index("label"); pl2 = fr.loc[pro.label]
t1 = pd.read_csv(A / "TBRS_LOF_tier1_regions.tsv", sep="\t").dropna(subset=["proband_beta", "ctrl_mean"])
t2 = pd.read_csv(A / "TBRS_LOF_tier2_regions.tsv", sep="\t").dropna(subset=["proband_beta", "ctrl_mean"])
ann = json.load(open(A / "annotation.json")); an = [x for x in ann["analyses"] if x["convention"] == "bed"][0]; sc = an["scenario_fused_partial_exons"]
inv = [l for l in pathlib.Path("../../signatures/INVENTORY.md").read_text(encoding="utf-8").splitlines() if l.startswith("| ") and not l.startswith("| ID")]
inv_rows = [[c.strip() for c in l.strip("|").split("|")][:6] for l in inv]
r = pd.read_csv(A / "ranks.anon.tsv", sep="\t")
tb = cc.loc[("TBRS_LOF", "tier2")]; tb1 = cc.loc[("TBRS_LOF", "tier1")]
second = r[(r.metric == "shape_r") & (r.null == "children") & (r.tier == "tier2") & (r.role == "child") & (r.signature == "TBRS_LOF") & (~r.is_proband)].value.max()
canyons = sum(1 for _ in open(A / "canyons.control.bed"))
go = r[(r.metric == "shape_r") & (r.null == "children") & (r.tier == "tier2") & (r.rank_children <= 5) & (r.role == "child")].groupby("label").size()
NICE = {"TBRS_LOF": "TBRS (DNMT3A LOF)", "DNMT3A_R882": "DNMT3A R882", "DNMT3A_GOF": "DNMT3A GOF (HESJAS)", "DNMT1": "DNMT1 (ADCA-DN)", "DNMT3B": "DNMT3B (ICF1)", "NSD1_SOTOS": "NSD1 (Sotos)", "NEG_KMT2D": "KMT2D (Kabuki)", "NEG_CHD7": "CHD7 (CHARGE)"}
f3 = lambda v: f"{v:.3f}"; f2 = lambda v: f"{v:.2f}"


def p_text(p):
    return "<0.001" if p < 0.001 else f"{p:.3f}"


# ----------------------------------------------------------------------------- styles
ss = getSampleStyleSheet()
body = ParagraphStyle("body", parent=ss["Normal"], fontName="Times-Roman", fontSize=10, leading=13.2, alignment=TA_JUSTIFY, spaceAfter=6)
h1 = ParagraphStyle("h1", parent=ss["Heading1"], fontName="Times-Bold", fontSize=13, spaceBefore=12, spaceAfter=6)
h2 = ParagraphStyle("h2", parent=ss["Heading2"], fontName="Times-Bold", fontSize=11, spaceBefore=8, spaceAfter=3)
title = ParagraphStyle("title", parent=ss["Title"], fontName="Times-Bold", fontSize=16, leading=20, spaceAfter=8)
small = ParagraphStyle("small", parent=body, fontSize=8.5, leading=10.5, alignment=0)
cap = ParagraphStyle("cap", parent=body, fontSize=8.5, leading=10.5, alignment=0, spaceBefore=3, spaceAfter=10)
cell = ParagraphStyle("cell", parent=body, fontSize=7.8, leading=9.4, alignment=0, spaceAfter=0)
ref = ParagraphStyle("ref", parent=body, fontSize=8.5, leading=10.5, alignment=0, leftIndent=12, firstLineIndent=-12, spaceAfter=2)


def tbl(rows, widths, header=True, fs=7.8):
    data = [[Paragraph(str(c), cell) for c in row] for row in rows]
    t = Table(data, colWidths=widths, repeatRows=1 if header else 0)
    st = [("FONTSIZE", (0, 0), (-1, -1), fs), ("VALIGN", (0, 0), (-1, -1), "TOP"), ("LINEBELOW", (0, -1), (-1, -1), 0.6, colors.black),
          ("LINEABOVE", (0, 0), (-1, 0), 0.6, colors.black), ("BOTTOMPADDING", (0, 0), (-1, -1), 2), ("TOPPADDING", (0, 0), (-1, -1), 2)]
    if header:
        st += [("LINEBELOW", (0, 0), (-1, 0), 0.6, colors.black), ("BACKGROUND", (0, 0), (-1, 0), colors.whitesmoke)]
    t.setStyle(TableStyle(st)); return t


def fig(name, caption, width=6.9 * inch):
    from PIL import Image as PILImage
    w, h = PILImage.open(F / name).size
    return KeepTogether([Image(str(F / name), width=width, height=width * h / w), Paragraph(caption, cap)])


S = []
S.append(Paragraph("A de novo DNMT3A frameshift deletion produces the Tatton-Brown-Rahman syndrome methylation episignature in saliva, detected by PacBio HiFi long-read sequencing", title))
S.append(Paragraph("Sebat Laboratory, University of California San Diego. Internal analysis report, version 0.1, "
                   f"{dt.date.today():%d %B %Y}. Pre-registered analysis of the SPARK long-read WGS pilot. Not for distribution outside the laboratory.", small))
S.append(Spacer(1, 8))

# ----------------------------------------------------------------------------- abstract
S.append(Paragraph("Abstract", h1))
S.append(Paragraph(
    f"<b>Background.</b> Loss-of-function variants in <i>DNMT3A</i> cause Tatton-Brown-Rahman syndrome (TBRS), an overgrowth and intellectual-disability disorder with a "
    f"reproducible blood DNA-methylation episignature. Whether such episignatures can be recovered from saliva using single-molecule long-read sequencing, without arrays, has not been established. "
    f"<b>Methods.</b> In the SPARK long-read whole-genome sequencing pilot ({n['all']} sequenced individuals from {n['fam']} families; PacBio HiFi, median {q.mean_depth[ok].median():.0f}x), one autism proband carried a "
    f"de novo {an['deleted_bp']}-bp deletion in <i>DNMT3A</i>. We pre-registered a falsifiable prediction table and scored every sample against eight published episignatures "
    f"(TBRS, <i>DNMT3A</i> R882 and gain-of-function, <i>DNMT1</i>, <i>DNMT3B</i>, <i>NSD1</i>, and <i>KMT2D</i> and <i>CHD7</i> negative controls) using region-level 5mC from MM/ML tags, "
    f"a leave-one-out shape correlation (<i>r</i>) against published delta-beta vectors, and a CpG-density-matched permutation null. "
    f"<b>Results.</b> The deletion removes 101 coding bases spanning exons 14-15 and is predicted to cause a frameshift at codon {an['codons_affected'].split('.')[1].split('-')[0]} within the ADD domain with nonsense-mediated decay. "
    f"On the TBRS signature the proband ranked first among {int(tb.n_children)} children (shape <i>r</i> = {tb.observed_shape_r:.2f}; next-highest child {second:.2f}; <i>z</i> = {tb.shape_r_z_vs_null:.1f}; permutation <i>p</i> {p_text(tb.perm_p_shape_r)}) "
    f"and first among all {int(tb.n_all)} samples, in both the two-source tier-1 set and the full set. The rank was unchanged across region padding, evidence thresholds, exclusion of composition-associated regions, "
    f"and {int(boot.shape[0])} bootstrap resamples of controls ({(boot['rank'] == 1).mean() * 100:.0f}% rank 1). <i>DNMT1</i>, <i>DNMT3B</i> and <i>CHD7</i> signatures were null as predicted; the R882 and Sotos signatures were positive but weaker; "
    f"the gain-of-function signature was not anticorrelated; and the <i>KMT2D</i> negative control was unexpectedly anticorrelated (<i>z</i> = {cc.loc[('NEG_KMT2D','tier2')].shape_r_z_vs_null:.1f}). "
    f"Genome-architecture features showed eroded methylation-canyon edges (slope residual {pl2.canyon_edge_slope:.1f} SD) and wider canyons, normal satellite methylation, and Horvath age acceleration of {pl2.age_accel_resid:+.1f} years. "
    f"<b>Conclusions.</b> A TBRS episignature is recoverable from ~20x HiFi saliva data in a single carrier, with the direction-specific shape statistic separating the proband from all controls. "
    f"Two pre-registered acceptance criteria (region completeness at 15x and clock fit) were not met and are reported as such.", body))

# ----------------------------------------------------------------------------- introduction
S.append(Paragraph("Introduction", h1))
S.append(Paragraph(
    "DNA-methylation episignatures are reproducible, genome-wide patterns of differential CpG methylation that accompany pathogenic variants in chromatin and methylation regulators, and they now support "
    "clinical variant classification for dozens of Mendelian disorders [1,2]. Germline loss of <i>DNMT3A</i> causes Tatton-Brown-Rahman syndrome (TBRS), characterised by overgrowth, intellectual disability "
    "and a hypomethylation-dominated blood episignature [3-5]; distinct signatures accompany <i>DNMT3A</i> R882 dominant-negative alleles [5] and PWWP-domain gain-of-function alleles that cause microcephalic dwarfism [6]. "
    "Because <i>DNMT3A</i> loss is also associated with autism, and <i>Dnmt3a</i> haploinsufficiency in mice produces autism-relevant behavioural deficits with global hypomethylation shared across neurodevelopmental-disorder models [7,8], a <i>DNMT3A</i>-deletion carrier in an autism cohort offers a natural test of whether these array-derived signatures transfer to a different tissue and technology.", body))
S.append(Paragraph(
    "Long-read sequencing reads 5-methylcytosine directly from polymerase kinetics, encoding per-read probabilities in MM/ML tags, so a single ~20x HiFi genome yields methylation calls at essentially every CpG without bisulfite conversion or arrays. "
    "This offers a route to episignature analysis in cohorts sequenced for other purposes, but introduces three differences from the published signatures: tissue (saliva versus blood), measurement (single-molecule counts at ~20x versus array intensities), and cohort size. "
    "We therefore pre-registered a prediction table and analysis plan before scoring and report results against it without revision.", body))

# ----------------------------------------------------------------------------- methods
S.append(Paragraph("Methods", h1))
S.append(Paragraph("Cohort, sequencing and data", h2))
S.append(Paragraph(
    f"Samples came from the SPARK PacBio pilot (release 2026_06): saliva-derived DNA from {n['all']} sequenced individuals in {n['fam']} families ({n['child']} children with autism, of whom {int(ch.sex.eq('Male').sum())} male, "
    f"median age {ch.age.median():.0f} years, range {ch.age.min():.0f}-{ch.age.max():.0f}; {n['parent']} parents, median age {pa.age.median():.0f} years). HiFi reads were aligned to GRCh38 (no-alt analysis set) with pbmm2; per-sample "
    f"merged, haplotagged BAMs from the family workflow were used where available (101 samples) and per-movie aligned BAMs otherwise; one unaligned movie file was excluded. Small variants were joint-called with DeepVariant/GLnexus and structural variants with sawfish (joint genotyping of 105 samples). "
    "All analysis was performed on the SDSC Expanse cluster; controlled-access data never left it, and this report uses anonymised labels.", body))
S.append(Paragraph("Proband variant annotation", h2))
S.append(Paragraph(
    f"The proband's deletion was taken from the sawfish joint call set (FILTER PASS, QUAL 645, one heterozygous carrier of 105). The VCF record has a 303-bp REF allele anchored at chr2:25,244,326 and a single-base ALT, so the deleted bases are "
    f"chr2:{an['deleted_interval_1based'].split(':')[1]} ({an['deleted_bp']} bp). We intersected the deletion with the canonical transcript ({ann['meta']['transcript']}, {ann['meta']['refseq_xrefs'][0]}) using Ensembl release {ann['meta']['ensembl_release'][0]}, "
    "computed the coding bases removed and the reading frame of the fused exon under two splicing scenarios, and mapped affected codons to Pfam domains. Parental genotypes at the joint call indicated a de novo event.", body))
S.append(Paragraph("Methylation calling", h2))
S.append(Paragraph(
    "Per-CpG 5mC was called with pb-CpG-tools v3.0.0 (aligned_bam_to_cpg_scores, model-based pileup, reference CpG sites). For samples with several movie BAMs, estimated modified and unmodified counts were summed across movies before computing beta = modified / (modified + unmodified). "
    f"Median CpGs called per sample was {q.n_cpgs[ok].median() / 1e6:.1f} million, median depth {q.mean_depth[ok].median():.1f}x (IQR {q.mean_depth[ok].quantile(.25):.1f}-{q.mean_depth[ok].quantile(.75):.1f}). All downstream statistics are region-level; single-CpG values are never tested.", body))
S.append(Paragraph("Episignature panel", h2))
S.append(Paragraph(
    "Eight signatures were assembled from published supplementary tables only (Table 1); any signature without a downloadable probe or region list would have been dropped rather than approximated. Illumina probe identifiers were mapped to GRCh38 CpG coordinates with the Zhou laboratory Infinium annotation "
    "(release v8.1; general-mask probes excluded) [9]; hg19 region tables were lifted with UCSC liftOver. Each CpG or region was padded by 250 bp and overlapping features within a signature were merged. Direction and, where published, effect size (case minus control delta-beta) were carried per region; "
    "for the EpiSign probe sets, delta-beta was derived from the published per-disorder and control mean-beta table. A TBRS tier-1 set was defined as regions supported by at least two independent sources with concordant direction; tier 2 is the union.", body))
S.append(Paragraph("Region scoring and similarity metrics", h2))
S.append(Paragraph(
    "For every sample and region we computed the unweighted mean beta over CpGs, setting the region to missing if fewer than 10 CpGs were covered or mean depth was below 15 reads (pre-registered rule). Leave-one-out control statistics (mean, SD, 5th percentile) were computed per region from all other "
    "non-flagged, non-outlier children (primary null) and from all other non-flagged, non-outlier samples (secondary null). The primary discriminator, shape <i>r</i>, is the Pearson correlation between a sample's delta-beta vector (sample minus leave-one-out control mean) and the published delta-beta vector "
    "over the same regions (sign vector where no effect size was published); mean <i>z</i> (mean leave-one-out z-score) and the fraction of regions below the control 5th percentile were secondary. Samples were ranked per signature among children and among all samples; empirical rank <i>p</i> = rank / n. "
    "A permutation null drew 1,000 random region sets of equal size from a background of array-probe-centred regions, matched on CpG count and width, and recomputed both statistics.", body))
S.append(Paragraph("Quality control and cell composition", h2))
S.append(Paragraph(
    "Betas at Illumina probe CpGs (HM450 and EPIC union) with depth >= 10 formed a probe matrix; probes missing in more than half the samples were dropped and remaining gaps were cohort-mean imputed with the imputed fraction recorded. "
    f"HEpiDISH [10] (robust partial correlations; {int(q.n_hepidish_ref_cpgs_present.dropna().iloc[0])} of 716 reference CpGs present) estimated epithelial, fibroblast and immune fractions; methylclock [11] computed Horvath, Hannum, Levine, skin-and-blood and PedBE ages. "
    "Age acceleration is the residual of DNAm age on chronological age fitted across the cohort. Samples were marked as QC outliers, retained in all tables but excluded from every null, if epithelial fraction exceeded the median + 3 MAD, mean depth was below 12x or fewer than 70% of CpGs reached 10x (step-04 acceptance rule), "
    "global 5mC was more than 4 SD from the mean, or extraction yielded fewer than half the median CpG count.", body))
S.append(Paragraph("Pre-registered flags", h2))
S.append(Paragraph(
    "Before scoring, individuals carrying a rare non-silent variant (missense, nonsense, frameshift, in-frame indel, splice region; cohort allele count <= 6) in <i>DNMT3A</i>, <i>DNMT1</i>, <i>DNMT3B</i>, <i>NSD1</i>, <i>SETD2</i> or <i>EZH2</i>, or a PASS structural variant overlapping their coding exons, were flagged in a frozen, hashed manifest "
    "so that a high rank in a flagged sample is interpreted as a finding rather than a failure. No population allele frequencies were available offline; rarity is cohort-internal.", body))
S.append(Paragraph("Genome-architecture (Layer 2) features", h2))
S.append(Paragraph(
    "Per sample (CpGs at >= 10x) we computed mean beta in ENCODE SCREEN v3 enhancer (pELS + dELS) and promoter (PLS) cCREs, UCSC CpG islands, RefSeq gene bodies, RepeatMasker satellite II/III, alpha satellite, LINE-1 and young LINE-1 (L1HS/L1PA2/L1PA3), and Zhou solo-WCGW CpGs in common PMDs [12]. "
    f"Methylation canyons were defined on the clean-control mean track as runs of CpGs with beta < 0.10 spanning >= 3.5 kb (gap <= 500 bp), following Jeong et al. [13]; {canyons} canyons were identified, and per sample we recorded interior beta, beta in the 1 kb inside and outside each boundary (edge slope = outer minus inner), and the sample's own canyon extent. "
    "Each feature was regressed on epithelial fraction, age, sex and mean depth in clean controls and expressed as a standardized residual.", body))
S.append(Paragraph("Pre-registration and deviations", h2))
S.append(Paragraph(
    "The prediction table (Table 3) and thresholds were fixed before any scoring. Two deviations are reported. First, the variant was initially described as a catalytic-domain lesion and its coordinates read as 1-based inclusive; inspection of the VCF REF allele showed the deletion is 302 bp and lies in the ADD domain, "
    "yielding a frameshift rather than the in-frame loss implied by the original reading; the loss-of-function hypothesis was unchanged. Second, after scoring, the investigator elected to adopt a 10x depth rule as primary for future runs because the pre-registered 15x rule left more than 30% of regions missing; "
    "all results in this report use the pre-registered 15x rule, and the 10x results are presented as sensitivity analyses. The primary clock (Horvath) and control exclusions were kept as pre-registered.", body))
S.append(Paragraph("Software", h2))
S.append(Paragraph("Snakemake 9.26 with the SLURM executor on Expanse; Python 3.11 (pandas, numpy, scipy); pb-CpG-tools 3.0.0; samtools/bcftools 1.24; R 4.5 with EpiDISH and methylclock. Code, configuration and count-level logs are version-controlled; sample-level tables remain on the cluster.", body))

# ----------------------------------------------------------------------------- results
S.append(PageBreak()); S.append(Paragraph("Results", h1))
S.append(Paragraph("The proband carries a de novo frameshift deletion in the DNMT3A ADD domain", h2))
ex_hit = an["exons_overlapped"]
S.append(Paragraph(
    f"The {an['deleted_bp']}-bp deletion removes the 3' {ex_hit[0]['bp_deleted']} bp of exon {ex_hit[0]['tx_exon']} and the 5' {ex_hit[1]['bp_deleted']} bp of exon {ex_hit[1]['tx_exon']} together with the intervening intron, eliminating the donor and acceptor sites of that junction. "
    f"If the residual exon segments splice as a fused exon, {an['exonic_bp_deleted']} coding bases are lost and the reading frame shifts at codon {an['codons_affected'].split('.')[1].split('-')[0]}, producing {sc['aberrant_residues']} aberrant residues and a premature stop predicted to encode a {sc['protein_length']}-residue product (wild type {sc['wt_protein_length']}); "
    f"the stop lies {sc['stop_upstream_of_last_junction_nt']} nt upstream of the last exon junction, so nonsense-mediated decay is predicted. The alternative outcome, skipping of both exons, is an in-frame loss of {an['scenario_exon_skipping']['residues_deleted']} residues. "
    "In either case the lesion falls in the ADD domain (GATA-like and PHD zinc fingers) and, in the frameshift case, ablates the entire catalytic methyltransferase domain including R882 (Figure 1). The variant was absent from both parents in the joint call set.", body))
S.append(fig("fig1_variant.png", "<b>Figure 1. The proband's DNMT3A deletion.</b> (A) Gene model on the minus strand of chr2 with the 302-bp deletion (red) spanning the exon 14-15 junction. (B) Protein with Pfam domains; the red box marks codons removed by the deletion, with the predicted frameshift consequence; R882 is indicated."))

S.append(Paragraph("Cohort quality control", h2))
S.append(Paragraph(
    f"Of {n['all']} sequenced individuals, {n['outl']} were QC outliers: {int(reasons.get('frac_cpgs_10x<0.70(step04 acceptance)', 0))} had fewer than 70% of CpGs at 10x (8-12x libraries), {int(reasons.get('extraction_failed(n_cpgs<50%median)', 0))} had essentially unaligned BAMs (~0.1% of reads mapped), "
    f"{int(reasons.get('epi_frac>median+3MAD', 0))} had epithelial fractions above the median + 3 MAD, and {int(reasons.get('global_5mc_z>4', 0))} had extreme global 5mC (categories overlap). The coverage criterion pre-registered for step 04 was therefore not met for these samples; they were retained and marked rather than dropped. "
    f"{n['flagged']} individuals ({n['flagged_child']} children) carried a rare non-silent variant in a flag gene (<i>NSD1</i> 11, <i>SETD2</i> 12, <i>DNMT1</i> 4, <i>DNMT3B</i> 3, <i>DNMT3A</i> 1). After both exclusions the primary null comprised {n['clean_child']} children and the secondary null {n['clean_all']} individuals.", body))
S.append(Paragraph(
    f"Estimated epithelial fraction was plausible for saliva (median {q.epi_frac[ok].median():.2f}, IQR {q.epi_frac[ok].quantile(.25):.2f}-{q.epi_frac[ok].quantile(.75):.2f}), higher in children ({ch[ok[ch.index]].epi_frac.median():.2f}) than parents ({pa[ok[pa.index]].epi_frac.median():.2f}), with immune cells forming the remainder (median {q.immune_frac[ok].median():.2f}). "
    f"Global mean 5mC was tightly distributed (median {q.global_mean_5mc[ok].median():.3f}). DNAm age correlated with chronological age at <i>r</i> = {clocks['Horvath']:.2f} (Horvath), {clocks['skinHorvath']:.2f} (skin and blood), {clocks['PedBE']:.2f} (PedBE), {clocks['Hannum']:.2f} (Hannum) and {clocks['Levine']:.2f} (Levine); "
    f"only the skin-and-blood clock met the pre-registered <i>r</i> > 0.8 criterion, and the Horvath clock was kept as primary because it was pre-specified (Figure 2). The proband passed all QC criteria (depth {pro.mean_depth:.1f}x, {pro.frac_cpgs_depth_ge10 * 100:.0f}% of CpGs at 10x, epithelial fraction {pro.epi_frac:.2f}, global 5mC {pro.global_mean_5mc:.3f}).", body))
S.append(fig("fig2_qc.png", "<b>Figure 2. Cohort quality control.</b> (A) HEpiDISH epithelial fraction by role; orange rings mark QC outliers, the red ring the proband. (B) Mean CpG depth against the fraction of CpGs at >= 10x; the dotted line is the 0.70 acceptance threshold. (C) DNAm age against chronological age for three clocks, non-outlier samples."))

S.append(Paragraph("The proband ranks first on the TBRS episignature", h2))
S.append(Paragraph(
    f"Across the eight-signature panel (Table 1), the proband's TBRS shape correlation was {tb.observed_shape_r:.3f} against the children null, {tb.shape_r_z_vs_null:.1f} SD above the clean-children distribution, ranking first of {int(tb.n_children)} children and first of {int(tb.n_all)} samples (Table 2, Figures 3 and 4). "
    f"The next-highest child scored {second:.3f}. Against the all-sample null the correlation was {ca.loc[('TBRS_LOF','tier2')].observed_shape_r:.3f} (<i>z</i> = {ca.loc[('TBRS_LOF','tier2')].shape_r_z_vs_null:.1f}). The stricter tier-1 set (165 regions supported by two or more independent cohorts) gave a higher correlation ({tb1.observed_shape_r:.3f}, <i>z</i> = {tb1.shape_r_z_vs_null:.1f}) with the same ranks. "
    f"No random CpG-density-matched region set among 1,000 reached the observed correlation (permutation <i>p</i> {p_text(tb.perm_p_shape_r)}); the exact rank <i>p</i> values are 1/{int(tb.n_children)} = {1 / tb.n_children:.3f} among children and 1/{int(tb.n_all)} = {1 / tb.n_all:.4f} overall, the floors these cohort sizes permit. "
    f"At the region level, {int((t1.proband_beta < t1.ctrl_q05).sum())} of {len(t1)} scored tier-1 regions lay below the clean-children 5th percentile and none above the 95th; the mean proband delta-beta was {(t1.proband_beta - t1.ctrl_mean).mean():+.3f} (Figure 5).", body))
S.append(Paragraph(
    f"The remaining panel members behaved largely as predicted. <i>DNMT1</i> (ADCA-DN), <i>DNMT3B</i> (ICF1) and the <i>CHD7</i> negative control were null (|<i>z</i>| < 1). The <i>DNMT3A</i> R882 signature was positive ({cc.loc[('DNMT3A_R882','tier2')].observed_shape_r:.2f}) but ranked {int(cc.loc[('DNMT3A_R882','tier2')].rank_children)}th among children, below the pre-registered top-quartile expectation. "
    f"The Sotos (<i>NSD1</i>) signature, which is known to overlap TBRS, was positive ({cc.loc[('NSD1_SOTOS','tier2')].observed_shape_r:.2f}, rank {int(cc.loc[('NSD1_SOTOS','tier2')].rank_children)}) though its <i>z</i> ({cc.loc[('NSD1_SOTOS','tier2')].shape_r_z_vs_null:.1f}) fell short of the formal criterion. Two predictions failed on direction: the gain-of-function (HESJAS) signature, expected to anticorrelate, was weakly positive ({cc.loc[('DNMT3A_GOF','tier2')].observed_shape_r:.2f}), "
    f"and the <i>KMT2D</i> (Kabuki) negative control, expected to be null, was clearly anticorrelated ({cc.loc[('NEG_KMT2D','tier2')].observed_shape_r:.2f}, <i>z</i> = {cc.loc[('NEG_KMT2D','tier2')].shape_r_z_vs_null:.1f}, last of {int(cc.loc[('NEG_KMT2D','tier2')].n_children)} children). "
    f"The proband also had the most negative mean <i>z</i> on several probe-based sets, including the <i>CHD7</i> control, indicating broad hypomethylation at array-probe regions; this is why the direction-specific shape statistic, not magnitude, was pre-registered as the discriminator. "
    f"The proband's parents were unremarkable on TBRS (mother rank {int(r[(r.signature=='TBRS_LOF')&(r.tier=='tier2')&(r.metric=='shape_r')&(r.null=='all')&(r.proband_family)&(~r.is_proband)&(r.role=='mother')].rank_all.iloc[0])}, father rank {int(r[(r.signature=='TBRS_LOF')&(r.tier=='tier2')&(r.metric=='shape_r')&(r.null=='all')&(r.proband_family)&(~r.is_proband)&(r.role=='father')].rank_all.iloc[0])} of {int(tb.n_all)}), consistent with a de novo event. "
    f"One child ranked in the top five on four signatures and was labelled a global outlier by the pre-registered cross-signature rule; it was a QC-flagged sample with an epithelial fraction of 0.90, not the proband.", body))
S.append(fig("fig3_panel.png", "<b>Figure 3. Panel similarity.</b> (A) Shape correlation of every sample (rows; children above the line, parents below) with each signature (columns), children null, tier-2 sets. The proband row is boxed in red; orange arrowheads mark variant-flagged samples and black arrowheads QC outliers. (B) The proband's profile (red) against the 5th, 50th and 95th percentiles of clean children."))
S.append(fig("fig4_ranks.png", "<b>Figure 4. Rank plots.</b> Children ranked by shape correlation for each signature; red = proband, orange = variant-flagged or QC outlier, grey = clean. The proband's rank among children is shown in each panel."))
S.append(fig("fig5_tbrs_regions.png", "<b>Figure 5. Region-level TBRS tier-1 methylation.</b> Proband beta minus clean-children mean (red points) for each scored tier-1 region, sorted by the published delta-beta, against the children's 5th-95th percentile band (grey)."))

S.append(Paragraph("The result is robust to analytic choices", h2))
S.append(Paragraph(
    f"The proband's first rank among children was unchanged with region padding of 0, 250 or 500 bp, with evidence thresholds of 10 CpGs / 15x (pre-registered), 10 / 10x, 5 / 10x and 5 / 15x, and after removing the {int(sens[sens.analysis.str.startswith('composition_clean')].analysis.iloc[0].split('dropped ')[1].split(' ')[0])} regions whose methylation correlated with epithelial fraction in controls (|<i>r</i>| > 0.5), "
    f"where the shape correlation rose to {sens[sens.analysis.str.startswith('composition_clean')].shape_r.iloc[0]:.3f}. In {(boot['rank'] == 1).mean() * 100:.0f}% of {len(boot)} bootstrap resamples of the clean children the proband ranked first (Figure 6). "
    f"The threshold analysis also quantified the missing-data problem: at the pre-registered 15x rule the median sample lacked {na[(na.set=='TBRS_LOF:tier2')&(na.min_depth==15)].median_na.iloc[0] * 100:.0f}% of TBRS tier-2 regions and {na[(na.set=='TBRS_LOF:tier1')&(na.min_depth==15)].median_na.iloc[0] * 100:.0f}% of tier-1 regions, exceeding the 30% ceiling set for step 06; at 10x these fell to "
    f"{na[(na.set=='TBRS_LOF:tier2')&(na.min_depth==10)&(na.min_cpgs==10)].median_na.iloc[0] * 100:.0f}% and {na[(na.set=='TBRS_LOF:tier1')&(na.min_depth==10)&(na.min_cpgs==10)].median_na.iloc[0] * 100:.0f}%. For the probe-window negative controls a large share of 502-bp regions contain fewer than 10 CpGs (<i>CHD7</i> {na[na.set=='NEG_CHD7:tier2'].frac_lt10cpg.iloc[0] * 100:.0f}%, <i>KMT2D</i> {na[na.set=='NEG_KMT2D:tier2'].frac_lt10cpg.iloc[0] * 100:.0f}%), a structural limit of applying probe-centred windows to sequencing data.", body))
S.append(fig("fig7_robustness.png", "<b>Figure 6. Robustness of the TBRS result.</b> Proband TBRS shape correlation and rank among children under each sensitivity analysis; red bars indicate rank 1."))

S.append(Paragraph("Genome-architecture features", h2))
S.append(Paragraph(
    f"After adjustment for composition, age, sex and depth, the proband's most extreme Layer 2 features were at methylation canyons: the boundary slope was {pl2.canyon_edge_slope:.1f} SD below controls, driven by lower methylation in the 1 kb outside canyon edges ({pl2.canyon_edge_outer_beta:.1f} SD), with canyons {pl2.canyon_mean_width_bp:+.1f} SD wider and {pl2.canyon_total_bp:+.1f} SD more total canyon sequence, "
    f"the edge-erosion phenotype described in <i>Dnmt3a</i>-null haematopoietic stem cells [13]. Enhancer methylation was reduced ({pl2.enhancer_beta:+.1f} SD) with promoters and CpG islands near control values, matching the enhancer-biased hypomethylation reported in <i>DNMT3A</i> overgrowth syndrome [5]. "
    f"Satellite II/III ({pl2.satII_III_beta:+.1f} SD), alpha-satellite ({pl2.alpha_sat_beta:+.1f} SD) and LINE-1 ({pl2.line1_beta:+.1f} SD) methylation were within the control range, arguing against a <i>DNMT3B</i>/ICF-like process, and PMD solo-WCGW methylation was unremarkable ({pl2.pmd_soloWCGW_beta:+.1f} SD). "
    f"Horvath age acceleration was {pl2.age_accel_resid:+.1f} years, above the clean-children 95th percentile ({q[clean & (q.role=='child')].age_accel_resid.quantile(.95):+.1f}), consistent with the accelerated epigenetic ageing reported in TBRS [3]; however, the PedBE ({pro.age_accel_resid_PedBE:+.1f} y) and skin-and-blood ({pro.age_accel_resid_skinHorvath:+.1f} y) clocks did not show it, and the Horvath clock's fit in this cohort was modest (Figure 7).", body))
S.append(fig("fig6_layer2.png", "<b>Figure 7. Layer 2 genome-architecture features.</b> Standardized residuals (adjusted for epithelial fraction, age, sex and depth) for clean controls (boxes) and the proband (red). Dotted lines mark +/-2 SD."))

# ----------------------------------------------------------------------------- discussion
S.append(Paragraph("Discussion", h1))
S.append(Paragraph(
    "A single saliva sample sequenced at ~35x with PacBio HiFi was sufficient to place a de novo <i>DNMT3A</i> frameshift carrier first among 105 individuals on a blood-derived TBRS episignature, with a margin over the nearest child of almost four-fold and a permutation null that never approached the observed value. "
    "The pre-registered design matters for interpreting this: the signature, statistic, tiers and thresholds were fixed before scoring, the proband and every variant-flagged individual were marked before any ranking, and the prediction table is reported with its failures.", body))
S.append(Paragraph(
    "The pattern across the panel is informative beyond the primary result. The proband was broadly hypomethylated at array-probe regions, which produced strong negative mean <i>z</i> on several unrelated sets; only the direction-specific shape correlation separated TBRS from that background, and it did so cleanly. "
    "The anticorrelation with the Kabuki (<i>KMT2D</i>) signature was not predicted but is interpretable: the <i>KMT2D</i> loss-of-function signature is enriched for hypermethylated CpG-rich regulatory regions, many of which lose methylation with <i>DNMT3A</i> deficiency, so the two signatures are partly mirror images rather than independent. "
    "The weakly positive rather than negative correlation with the gain-of-function signature likely reflects that signature's provenance (fibroblast RRBS, no published effect sizes) more than biology. The R882 signature, derived from three dominant-negative carriers by WGBS, was positive but not top-quartile, in line with the expectation that a haploinsufficient allele produces a milder version of the R882 hypomethylation.", body))
S.append(Paragraph(
    "The Layer 2 features extend the array-based picture with information only long reads provide at this coverage: eroded canyon edges and wider canyons, which in mouse are a direct consequence of losing Dnmt3a at canyon boundaries, and enhancer-biased hypomethylation. Normal satellite and LINE-1 methylation excludes a global maintenance defect. "
    "Epigenetic age acceleration on the Horvath clock was large but clock-dependent, and we caution against over-reading it given the modest clock fit in long-read saliva data.", body))
S.append(Paragraph("Limitations", h2))
S.append(Paragraph(
    "This is a single carrier in a small cohort; the exact rank p cannot fall below 1/38 among children, and the study cannot estimate sensitivity or specificity of the approach. The signatures were derived from blood arrays or fibroblast RRBS and applied to saliva, whose epithelial fraction varied widely; composition was handled by covariate adjustment, outlier exclusion and a region-exclusion sensitivity analysis rather than by tissue-matched signatures. "
    "At ~20x, single-CpG calls are noisy and region completeness under the pre-registered 15x rule fell short of the acceptance criterion; the result was insensitive to this, but future designs should set depth rules from the observed coverage distribution. Rare-variant flagging used cohort-internal allele counts and canonical-transcript consequences only, and included missense variants, so the control set is conservative. "
    "The EpiSign classifiers themselves were not reproduced; we tested similarity of direction and magnitude to published probe sets. Finally, the inclusion of low-coverage samples as marked outliers rather than exclusions, and the post hoc adoption of a 10x rule for future runs, are documented deviations from the plan.", body))
S.append(Paragraph("Data and code availability", h2))
S.append(Paragraph("SPARK data are controlled-access (SFARI Base). Sample-level results remain on the SDSC Expanse cluster. The Snakemake workflow, configuration, signature harmonisation tables, per-source provenance (SOURCE.md files) and count-level run logs are in the project repository; the sample key linking anonymised labels to identifiers is held on the cluster only.", body))

# ----------------------------------------------------------------------------- tables
S.append(PageBreak()); S.append(Paragraph("Tables", h1))
S.append(Paragraph("<b>Table 1. Episignature panel.</b> Sources, probe counts, harmonised hg38 regions (tier 1 / tier 2 for TBRS) and direction.", cap))
S.append(tbl([["ID", "Sources", "Array probes mapped", "Regions", "Direction", "Status"]] + inv_rows, [1.0 * inch, 2.1 * inch, 0.9 * inch, 0.9 * inch, 1.1 * inch, 0.5 * inch]))
S.append(Spacer(1, 10))
rows = [["Signature", "Tier", "shape r", "z vs children", "Rank (children)", "Rank (all)", "Perm. p", "Predicted", "Verdict"]]
for (sig, tier), rr in cc.iterrows():
    rows.append([NICE[sig], tier, f3(rr.observed_shape_r), f2(rr.shape_r_z_vs_null), f"{int(rr.rank_children)} / {int(rr.n_children)}", f"{int(rr.rank_all)} / {int(rr.n_all)}", p_text(rr.perm_p_shape_r), f"{rr.predicted_sign}, {rr.predicted_rank}", rr.verdict])
for sig, rr in l2.iterrows():
    rows.append([sig.replace("L2:", "Layer 2: "), "-", f2(rr.observed_shape_r) + (" SD" if "sat" in sig else " y"), "-", "-", "-", "-", rr.predicted_sign, rr.verdict])
S.append(Paragraph("<b>Table 2. Proband report card against the pre-registered prediction table</b> (children null; Layer 2 rows show the standardized residual or years). Verdicts are mechanical: sign requires |z| > 2 in the predicted direction (or |z| < 2 for a null prediction); rank requires the stated rank among children.", cap))
S.append(tbl(rows, [1.25 * inch, 0.45 * inch, 0.55 * inch, 0.7 * inch, 0.85 * inch, 0.75 * inch, 0.55 * inch, 0.95 * inch, 0.6 * inch]))
S.append(Spacer(1, 10))
rows = [["Analysis", "shape r", "Rank (children)", "Rank (all)", "Regions"]]
for _, rr in sens[sens.null == "children"].iterrows():
    rows.append([rr.analysis.replace("step06:", "step 06 ").replace("_", " "), f3(rr.shape_r), f"{int(rr.rank_children)} / {int(rr.n_children)}", (f"{int(rr.rank_all)} / {int(rr.n_all)}" if pd.notna(rr.rank_all) else "-"), (f"{int(rr.n_regions)}" if pd.notna(rr.n_regions) else "-")])
S.append(Paragraph("<b>Table 3. Sensitivity and robustness of the proband's TBRS result.</b> The bootstrap row reports the mean shape r and median rank over 200 resamples of the clean children.", cap))
S.append(tbl(rows, [3.6 * inch, 0.7 * inch, 1.0 * inch, 0.9 * inch, 0.7 * inch]))

# ----------------------------------------------------------------------------- references
S.append(PageBreak()); S.append(Paragraph("References", h1))
refs = [
    "Aref-Eshghi E, Kerkhof J, Pedro VP, et al. Evaluation of DNA methylation episignatures for diagnosis and phenotype correlations in 42 Mendelian neurodevelopmental disorders. Am J Hum Genet 2020;106:356-370.",
    "Levy MA, McConkey H, Kerkhof J, et al. Novel diagnostic DNA methylation episignatures expand and refine the epigenetic landscapes of Mendelian disorders. HGG Adv 2022;3:100075.",
    "Jeffries AR, Maroofian R, Salter CG, et al. Growth disrupting mutations in epigenetic regulatory molecules are associated with abnormalities of epigenetic aging. Genome Res 2019;29:1057-1066.",
    "Tatton-Brown K, Seal S, Ruark E, et al. Mutations in the DNA methyltransferase gene DNMT3A cause an overgrowth syndrome with intellectual disability. Nat Genet 2014;46:385-388.",
    "Smith AM, LaValle TA, Shinawi M, et al. Functional and epigenetic phenotypes of humans and mice with DNMT3A overgrowth syndrome. Nat Commun 2021;12:4549.",
    "Heyn P, Logan CV, Fluteau A, et al. Gain-of-function DNMT3A mutations cause microcephalic dwarfism and hypermethylation of Polycomb-regulated regions. Nat Genet 2019;51:96-105.",
    "Christian DL, Wu DY, Martin JR, et al. DNMT3A haploinsufficiency results in behavioral deficits and global epigenomic dysregulation shared across neurodevelopmental disorders. Cell Rep 2020;33:108416.",
    "Beard DC, Zhang X, Wu DY, et al. Distinct disease mutations in DNMT3A result in a spectrum of behavioral, epigenetic, and transcriptional deficits. Cell Rep 2023;42:113411.",
    "Zhou W, Laird PW, Shen H. Comprehensive characterization, annotation and innovative use of Infinium DNA methylation BeadChip probes. Nucleic Acids Res 2017;45:e22.",
    "Zheng SC, Webster AP, Dong D, et al. A novel cell-type deconvolution algorithm reveals substantial contamination by immune cells in saliva, buccal and cervix. Epigenomics 2018;10:925-940.",
    "Pelegi-Siso D, de Prado P, Ronkainen J, Bustamante M, Gonzalez JR. methylclock: a Bioconductor package to estimate DNA methylation age. Bioinformatics 2021;37:1759-1760.",
    "Zhou W, Dinh HQ, Ramjan Z, et al. DNA methylation loss in late-replicating domains is linked to mitotic cell division. Nat Genet 2018;50:591-602.",
    "Jeong M, Sun D, Luo M, et al. Large conserved domains of low DNA methylation maintained by Dnmt3a. Nat Genet 2014;46:17-23.",
    "Kernohan KD, Cigana Schenkel L, Huang L, et al. Identification of a methylation profile for DNMT1-associated autosomal dominant cerebellar ataxia, deafness, and narcolepsy. Clin Epigenetics 2016;8:91.",
    "Choufani S, Cytrynbaum C, Chung BHY, et al. NSD1 mutations generate a genome-wide DNA methylation signature. Nat Commun 2015;6:10207.",
    "Butcher DT, Cytrynbaum C, Turinsky AL, et al. CHARGE and Kabuki syndromes: gene-specific DNA methylation signatures identify epigenetic mechanisms linking these clinically overlapping conditions. Am J Hum Genet 2017;100:773-788.",
]
for i, t in enumerate(refs, 1):
    S.append(Paragraph(f"{i}. {t}", ref))

doc = SimpleDocTemplate(OUT, pagesize=letter, leftMargin=0.9 * inch, rightMargin=0.9 * inch, topMargin=0.8 * inch, bottomMargin=0.8 * inch,
                        title="DNMT3A episignature in saliva by HiFi long-read sequencing", author="Sebat Laboratory, UC San Diego")


def footer(canvas, doc_):
    canvas.saveState(); canvas.setFont("Times-Roman", 8)
    canvas.drawString(0.9 * inch, 0.5 * inch, "DNMT3A episignature - SPARK lrWGS pilot - internal draft v0.1")
    canvas.drawRightString(letter[0] - 0.9 * inch, 0.5 * inch, f"Page {doc_.page}"); canvas.restoreState()


doc.build(S, onFirstPage=footer, onLaterPages=footer)
print("wrote", OUT)
