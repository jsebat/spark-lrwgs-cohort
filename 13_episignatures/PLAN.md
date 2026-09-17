# PLAN.md — implementation steps

Work through these in order. Each step has inputs, outputs, and an acceptance check.
Stop and ask Jonathan at any **[CHECKPOINT]**. Commit after every step.

---

## 00 — Environment & scaffold
- Create repo layout from CLAUDE.md; `envs/py.yaml`, `envs/r.yaml`; Snakefile skeleton with one rule per step.
- `config/config.yaml`: paths to BAM dir, reference FASTA, output root; thresholds (min_cpgs=10, min_depth=15, region_pad=250).
- `config/samples.tsv` template columns: `sample_id, family_id, role(proband|mother|father), age, sex, mean_depth, flag_gene, flag_note`. Gitignore it.
- **Accept:** `snakemake -n` runs clean; `pb-CpG-tools`, `bedtools`, `liftOver`, R packages all resolve.

## 01 — Annotate proband variant
- Intersect chr2:25,244,326-25,244,628 with DNMT3A exons of NM_022552 / ENST00000264709 (GENCODE GTF, hg38).
- Report: overlapped exon(s), reading frame, predicted consequence (frameshift/NMD vs in-frame vs intronic/splice).
- Also report distance to R882 codon and which protein domain is affected.
- Output: `results/01_variant/annotation.md`.
- **[CHECKPOINT]** If the deletion is intronic-only, the LOF hypothesis weakens — discuss before continuing.

## 02 — Signature inventory & harmonization
For each Layer 1 signature in CLAUDE.md:
1. Locate the supplementary table (probe IDs and/or DMR coordinates, with direction and effect size if given). Store in `signatures/raw/<ID>/` with a `SOURCE.md` (paper, DOI, table number, URL, download date).
2. Probe IDs → hg38 CpG coordinates via Illumina EPIC/450K manifest or Zhou lab InfiniumAnnotation (hg38). DMRs in hg19 → liftOver.
3. Expand each CpG/probe to a region: ±`region_pad` bp; merge overlapping regions within a signature.
4. Write `signatures/hg38/<ID>.bed` with columns: chrom, start, end, signature, source, direction(+/-), published_delta_beta (NA if absent), n_array_probes.
5. Build TBRS_LOF tiers: `tier1` = regions supported by >=2 sources with concordant direction; `tier2` = union.
- Output: `signatures/INVENTORY.md` — one row per signature: n regions, sources, direction summary, status (OK / DROPPED + reason).
- **Accept:** every kept signature has >=20 regions; no region with conflicting direction inside tier1.
- **[CHECKPOINT]** Review INVENTORY.md — which signatures survived, and does TBRS tier1 look sensible (expected: dozens to a few hundred regions, predominantly hypomethylated).

## 03 — Sample manifest & pre-registration of flags
- Fill `samples.tsv` from SPARK metadata. Confirm mean depth per BAM (`samtools coverage` or existing QC).
- From available variant calls, flag samples with rare coding/splice variants in DNMT3A, DNMT1, DNMT3B, NSD1, SETD2, EZH2. Record gene and variant class in `flag_gene`/`flag_note`. Mark the DNMT3A proband and its parents.
- Freeze the manifest (`results/03_manifest/samples.frozen.tsv`, hashed) before any scoring.
- **Accept:** n children, n parents, n flagged, and proband row all present. Report counts (not IDs) in the step log.

## 04 — Methylation extraction
- Per BAM: `aligned_bam_to_cpg_scores` (pb-CpG-tools, `--modsites-mode reference`, model-based) → per-CpG bed with beta and coverage. Fallback: `modkit pileup --cpg --combine-strands`.
- Store per-sample bedGraph-like parquet: chrom, pos, beta, depth.
- **Accept:** per-sample genome-wide CpG count and mean depth within 20% of expectation; no sample with <70% CpGs covered at >=10x.

## 05 — QC & cell composition
- Extract betas at HEpiDISH reference CpGs (need hg38 coordinates via manifest); run HEpiDISH → epithelial, fibroblast, immune fractions per sample. Also extract Horvath/PC-clock CpGs and compute DNAm age; impute missing clock CpGs with cohort mean and report % imputed.
- Global mean 5mC per sample.
- Output: `results/05_qc/sample_qc.tsv` (adds epi_frac, immune_frac, dnam_age, age_accel_resid, global_mean_5mc, pct_clock_cpgs_imputed).
- Flag QC outliers (epi_frac > cohort median + 3 MAD, or depth outliers) — keep them but mark `qc_outlier=TRUE`.
- **Accept:** epithelial fractions plausible for saliva (roughly 5–60%); DNAm age correlates with chronological age (r > 0.8 across cohort).

## 06 — Layer 1 region scoring
- For every sample x signature region: mean beta over CpGs in region, n_cpgs, mean depth. Apply min_cpgs/min_depth → NA.
- Matrix `results/06_layer1/region_betas.parquet` (samples x regions).
- Leave-one-out: for each sample, control mean/SD per region from all *other* non-flagged children (primary) and separately from all other non-flagged samples (secondary).
- Compute per signature: `mean_z`, `shape_r` (vs published delta_beta where available; where absent, use sign-only direction vector), fraction of regions below LOO 5th percentile, and NA fraction.
- Permutation null: 1000 random region sets matched on region count and CpG density → empirical p for `mean_z` and `shape_r`.
- Output: `results/06_layer1/scores.tsv` (sample x signature x metric), long format.
- **Accept:** score matrix has no signature with >30% NA regions in the median sample; proband row present.

## 07 — Layer 2 global features
Annotation resources (download to `resources/`, hg38): ENCODE cCREs (enhancer/promoter), UCSC CpG islands, RepeatMasker (satellite II/III, LINE-1), PMD/solo-WCGW set (Zhou 2018), canyon definitions (derive: merge hypomethylated (<10%) CpG runs >3.5 kb in the control-mean track).
Per sample compute: enhancer/promoter/CGI/gene-body mean beta; satellite II/III and L1 mean beta; PMD solo-WCGW mean beta; canyon metrics (total canyon bp, mean canyon width, edge slope over ±1 kb of each control-defined canyon boundary); age acceleration from step 05.
Fit `feature ~ epi_frac + age + sex + mean_depth` on non-flagged samples; report standardized residual per sample.
- Output: `results/07_layer2/features.tsv`, `features_resid.tsv`.
- **Accept:** residuals approximately centered on 0 for controls; satellite metrics stable across controls (CV < 15%).

## 08 — Panel similarity & ranking
- Assemble sample x signature `shape_r` and `mean_z` matrices; rank samples per signature (children-only and all).
- Proband report card: observed rank/percentile per signature vs. the pre-registered prediction table in CLAUDE.md. Do not soften mismatches.
- Cross-signature check: samples ranking top-5 on >=4 signatures → label "global outlier, likely composition/QC" and cross-reference step 05 flags.
- Flagged samples (step 03): report their ranks separately.
- Proband's parents: report TBRS rank; interpret against inheritance (de novo vs inherited) once genotype is known.
- Output: `results/08_panel/ranks.tsv`, `proband_report_card.tsv`.

## 09 — Statistics & robustness
- Empirical p for proband TBRS rank (children-only primary, all-sample secondary).
- Sensitivity: tier1 vs tier2; region_pad 0/250/500; with/without epithelial-fraction exclusion of composition-driven regions (regions where beta ~ epi_frac r > 0.5 in controls are dropped in the "composition-clean" run).
- Bootstrap over control samples for rank stability.
- Output: `results/09_stats/sensitivity.tsv`.

## 10 — Report
`results/report/report.html` (Jupyter → nbconvert or Quarto), containing:
1. Variant annotation summary.
2. Signature inventory table.
3. Cohort QC (composition, depth, DNAm age vs age).
4. Heatmap: samples x signatures (`shape_r`), proband and flagged samples marked.
5. Proband radar plot across the panel.
6. Rank plots per signature with proband highlighted; empirical p.
7. TBRS tier1 delta-beta heatmap (proband vs control distribution).
8. Layer 2 residual panel (canyons, enhancers, satellites, PMD, age acceleration).
9. Prediction table with observed vs predicted, pass/fail per row.
10. Limitations: blood-derived signatures applied to saliva, 20x single-CpG noise, n of cohort, EpiSign classifier not reproduced.
- **[CHECKPOINT]** Review with Jonathan before anything is shared beyond the lab.

---

## Open items to resolve with Jonathan before step 03
- Exact SPARK pilot sample count (children / parents) and whether variant calls are available for flagging.
- Compute environment and SLURM partition; where BAMs live.
- Whether the proband's deletion is de novo or inherited.
