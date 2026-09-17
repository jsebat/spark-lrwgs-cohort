# DNMT3A episignature panel — long-read methylation (SPARK lrWGS pilot)

## Purpose
Quantify similarity of every sample in the SPARK long-read WGS pilot to a panel of
published DNA-methylation episignatures, and test the pre-registered prediction that
the DNMT3A-deletion proband ranks first on the TBRS (DNMT3A loss-of-function) signature.

## Study design (fixed — do not change without asking Jonathan)
- **Proband variant:** heterozygous deletion, GRCh38 chr2:25,244,326-25,244,628 (DNMT3A, minus strand,
  3' half of gene / catalytic domain region). Working hypothesis: exon-disrupting → frameshift →
  haploinsufficiency (NOT R882 dominant-negative). Step 01 must confirm exon vs intron.
- **Tissue:** saliva (leukocyte-rich, variable buccal epithelial fraction).
- **Platform:** PacBio HiFi, ~20x. Methylation from MM/ML tags (5mC CpG).
- **Cohort:** SPARK lrWGS pilot trios. Probands (mostly non-DNMT3A ASD cases) + parents.
  Primary null = non-DNMT3A children. Secondary null = parents. Proband's own parents = within-family controls.
- **Reference genome:** GRCh38 only. All coordinates hg38, BED = 0-based half-open.
- **Methylation unit:** beta in [0,1]. Region-level, never single-CpG, for any statistic.
- **Minimum evidence per region:** >=10 CpGs AND mean depth >= `thresholds.min_depth` reads (config.yaml; 15 pre-registered, changed to 10 on 2026-09-10 after the step 06 acceptance check), else region = NA for that sample.

## Two-layer readout
**Layer 1 — probe/DMR-based episignature panel** (all from published supplementary tables):
| ID | Disorder / mechanism | Expected in proband |
|---|---|---|
| TBRS_LOF | DNMT3A haploinsufficiency (Jeffries 2019 Genome Res; EpiSign TBRS; Smith 2021 germline-LOF subgroup) | **Rank 1, hypomethylated, positive shape correlation** |
| DNMT3A_R882 | DNMT3A dominant-negative (Smith 2021 R882 subgroup) | positive correlation, lower magnitude than published |
| DNMT3A_GOF | Heyn-Sproul-Jackson (Heyn 2019 Nat Genet, PWWP missense, hypermethylation) | negative correlation |
| DNMT1 | ADCA-DN / HSAN1E (EpiSign ADCADN) | ~0 |
| DNMT3B | ICF1 (EpiSign ICF) | ~0 |
| NSD1_SOTOS | Sotos — known TBRS overlap | modest positive |
| NEG_KMT2D, NEG_CHD7 | Kabuki, CHARGE — negative controls | ~0 |

**Layer 2 — genome-architecture features** (long-read only, region-aggregated, covariate-adjusted):
canyon width/edge erosion (Jeong 2014), enhancer vs promoter vs CGI mean methylation (Smith 2021 pattern),
satellite II/III + LINE-1 methylation (DNMT3B/ICF readout — expected NORMAL in proband),
PMD / solo-WCGW methylation (Zhou 2018), epigenetic age acceleration (Horvath / PC-clock),
global mean 5mC.

## Similarity metrics (per sample x signature)
1. `mean_z`: mean leave-one-out z-score of region betas vs. all other samples (magnitude).
2. `shape_r`: Pearson r between the sample's delta-beta vector (sample minus LOO control mean)
   and the published case-vs-control delta-beta vector over the same regions (direction/shape).
   `shape_r` is the primary discriminator; `mean_z` is secondary.
3. Empirical p = rank of proband among all scored samples (report children-only and all-sample ranks).
4. Permutation null: 1000 random CpG-density-matched region sets of equal size.

## Pre-registered prediction (falsifiable; write results against this table, do not revise it)
TBRS_LOF shape_r rank 1 among children; R882 positive but weaker; GOF negative; DNMT1/DNMT3B/NEG ~0;
Sotos modest positive; satellites within control range; age acceleration > 0.

## Hard rules
- Controlled-access data: never copy SPARK BAMs/genotypes off the compute environment, never
  paste sample IDs or genotypes into chat, never commit data or sample manifests to git.
- Any signature without a downloadable probe/DMR list is **dropped**, not approximated. Record why in `signatures/INVENTORY.md`.
- No silent imputation. Missing regions stay NA and are counted in QC.
- Every figure is generated from real pipeline output. Never fabricate values, placeholder curves, or "example" numbers.
- Covariates in every Layer 2 model: epithelial fraction, age, sex, mean depth. Report residuals, not raw values.
- Flag before ranking: samples with variants in DNMT3A, DNMT1, DNMT3B, NSD1, SETD2, or EZH2 are marked in the
  manifest *before* scoring so a high rank there is interpreted as a finding, not a failure.
- Don't infer or hard-code coordinates for probes; always derive from Illumina manifest / Zhou lab hg38 annotation.
- If a published paper's supplementary content differs from what CLAUDE.md assumes, trust the paper and update INVENTORY.md.

## Tooling
Python 3.11: pandas, pyranges, numpy, scipy, statsmodels, matplotlib, pysam. CLI: pb-CpG-tools
(`aligned_bam_to_cpg_scores`), modkit (fallback), bedtools, UCSC liftOver (hg19→hg38 for old DMRs).
R (via rpy2 or separate scripts): EpiDISH (HEpiDISH epithelial/immune deconvolution), methylclock.
Workflow runner: Snakemake. Environment: conda `envs/*.yaml`. Compute: HPC (confirm — likely Expanse) via SLURM.

## Repo layout
```
config/            config.yaml (paths, thresholds), samples.tsv (NOT committed)
envs/              conda yamls
signatures/        raw/ (downloaded supp tables), hg38/ (harmonized BED), INVENTORY.md
resources/         manifests, annotations (ENCODE cCREs, CGIs, RepeatMasker satellites, PMDs, canyons)
workflow/          Snakefile + rules/ + scripts/ (numbered 01..10)
results/           per-step outputs; results/report/ is the deliverable
notebooks/         exploratory only, never load-bearing
```
