# spark-lrwgs-cohort

Cohort-level analysis workflow for long-read (PacBio HiFi) whole-genome sequencing of
autism families. Starts from joint-called cohort callsets (small variants, structural
variants, tandem repeats) and per-family pipeline outputs, and produces:

1. **Variant tiering** — LOFTEE high-confidence LoF tiered by GeneBayes s_het
   (`lof_t1/t2/t3`), missense tiered by four rankscores (`miss_t1..t4`).
2. **Panel of normals** — unaffected founders, leave-one-family-out, as an artifact filter.
3. **Burden** — per-sample tier counts; conditional logistic regression with genome-wide
   depth as the nuisance covariate; transmission test.
4. **De novo SVs** — coding-impact prioritisation in disease genes, short-read
   re-genotyping (GraphTyper) with a three-way outcome, cross-family recurrence filter, founder-panel recurrence matched by reciprocal overlap >= 0.5 (an empty panel aborts the run).
5. **Inheritance models** — recessive / X-linked / compound-het with DDG2P allelic
   requirement enforced.
6. **Ascertainment** — measured against an external de novo callset.
7. **Phenotype** — cohort severity vs the parent study; syndrome-match scoring of
   candidate carriers against published feature sets.

Methods for every step are in `docs/METHODS.md`. Configuration lives in
`config/cohort.env` (copy from `cohort.env.example`); scripts take family or sample
identifiers as arguments and never embed them. **Run `scripts/phi_scan.sh` before every
push.** This repository is distinct from the earlier per-family clinical-interpretation
workflow, which is not included.

Directory order is execution order: `01_qc` → `02_tiering` → `03_panel` → `04_burden`,
`05_denovo`, `06_inheritance`, `07_ascertainment`, `08_phenotype`.

## Relationship to the per-family workflow
The earlier **per-family (N-of-1) clinical interpretation** workflow lives in the separate
repository `longread-autism-workflow` and is not part of this one.

### Tandem repeats (02_tiering/tr_*)
`tr_cohort.sb` runs the STRchive known-locus screen on the joint TRGT callset (`known_repeats.py`, population-modal
gate) and the genome-wide scan `tr_outliers.py`: allele length = TRGT AL, QC on spanning reads only (allele purity is
not usable on a multi-motif catalog), leave-one-family-out founder reference, outlier > founder p99 + 3 motif units,
de novo expansion > both parents + the same margin, transmission test, then `tr_clogit.R` (depth-adjusted burden).

### Read-level SV review (05_denovo/sv_read_review.sh) — required before a candidate is validation-ready
Junction reads at both breakpoints, haplotype-resolved read counts inside vs flanks, and persistence of heterozygous
SNVs inside the interval. A joint-caller genotype on three reads passed every upstream filter once; this step is
what caught it.

### Status
Frozen for the first cohort report (September 2026). Analyses of the polygenic-score comparison are a separate
downstream project and are not included. No identifiers, names or data files are in this repository; run
`scripts/phi_scan.sh` before any commit.

### Methylation (09_methylation)
`meth_stage1.py`: imprinted-gene coordinates and promoter islands (GENCODE), cohort CpG-island matrix from methbat
profiles (combined / hap1 / hap2 / ASM p), haplotype-convention and parent-of-origin test at 18 imprinting control
regions, LoF sites in imprinted genes, marker means on the Loyfer U25 atlas. `deconv.R`: NNLS deconvolution and a
calibrated two-compartment epithelial fraction (atlas values are UNMETHYLATED fractions; sample methylation is
converted to 1 - m). `pofo_local.py`: parent of origin from the NEAREST informative phased SNVs with a phase-switch
flag (whole-block tallies are wrong about 10% of the time). `meth_stage2a.py`: per-child ASM islands, parent-of-origin
catalogue (imprinted-like / sequence-dependent / sporadic), loss-of-imprinting screen. `meth_stage2b.py`:
cohort-relative LOI screen, gene annotation, cis SNV / SV / repeat-length tests in founders with composition and DNA
source as covariates. Inputs: CPG_DIR and PHASED_VCF_DIR hold pb-CpG-tools beds and HiPhase VCFs linked BESIDE their
indexes (the pipeline writes indexes to sibling directories); METH_REFS holds the geneimprint list (included) and
Atlas.U25.l4.hg38.tsv from nloyfer/UXM_deconv.

### Transmission follow-up in short-read cohorts (10_transmission)
`gene_tdt.py`: parent-of-origin TDT of LoF alleles for one gene (GENE, REGION, RECURRENT env) over the rare-variant
pipeline's annotated family-genotype tables: four strata (mother / father to proband / sibling), GQ, DP and
allele-balance filters, WGS precedence for families on two platforms, per-cohort and pooled exact binomial tests.
`imprinted_tdt_all.py`: the same over all catalogued imprinted genes with expressed-allele test arms, silenced-allele
and sibling controls, Bonferroni ranking, collapsed burden, and clonal-hematopoiesis flags (gnomAD outlier_lof plus a
curated driver list) that exclude genes from the burden only. Under-transmission on every arm marks spurious parental
heterozygote calls; the control arms are not optional.
