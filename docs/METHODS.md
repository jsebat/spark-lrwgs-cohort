# METHODS (running) - long-read autism cohort
Maintained from 2026-09-06 so the final report (a cohort-level short publication, JS) can assemble
its Methods without reconstruction. Each entry cites the DECISIONS.md record (D#) holding the
numbers and the scripts.

## 1. Cohort, sequencing and variant calling
105 PacBio HiFi genomes from 35 families (33 trios; 2 quads: family-A affected+unaffected,
family-B affected+affected); 101 SPARK, 4 REACH. HiFi-human-WGS-WDL v3.3.1 via miniwdl on SDSC
Expanse; GRCh38 no-alt analysis set. Small variants DeepVariant with GLnexus joint calling
(Freeze 1: 21,780,907 sites); structural variants sawfish joint calling (347,631 records); tandem
repeats TRGT (937,195 loci). Composition from pipeline-written pedigrees: 68 founders (65
unaffected, 3 affected) and 37 offspring (36 affected, 1 unaffected). D54, D56.
*How it was run (module `00_upstream`):* HiFi-human-WGS-WDL v3.3.1 (commit 477ef39) `family.wdl`, one miniwdl run per
family on Slurm with the `slurm_singularity` backend (SingularityPro 4.1.2), PacBio resource bundle v3.1.0, backend
`HPC`, no preemption; tools pinned by container digest in the workflow's image manifest (pbmm2; DeepVariant 1.10.0 with
GLnexus; HiPhase 1.6.0; sawfish 2.2.1; TRGT 5.0.0; mosdepth, paraphase, mitorsaw, pbstarphase, MethBat, pb-cpg-tools,
svpack, slivar for tertiary annotation). Inputs are one five-key `inputs.json` per family generated from the cohort
manifest; a driver job runs miniwdl and a watchdog relaunches from the call cache on driver death. The per-task Slurm
time limit (6 h) is calibrated to ~22-24x and is raised per run for deeper inputs. The cohort-wide GLnexus joint call
(`cohort/glnexus_full.sb`) produces freeze 1 from the 105 per-sample gVCFs.

## 2. Quality control, identity and masking
somalier relatedness and sex concordance. Cross-platform identity: long-read and short-read
genotypes of the same individual related at ~1.0, because relatedness alone cannot detect a
same-sex sibling exchange (D57). Genotype QC for burden: GQ >= 20, DP >= 10, heterozygous allele
balance 0.25-0.75. Region mask: genomicSuperDups + simpleRepeat + RepeatMasker
Simple_repeat/Low_complexity (338 Mb, ~11% of the genome); point variants masked on any overlap,
SV intervals masked when >= 50% of the interval is repeat. D56, D61, D63.

## 3. Panel of normals
65 unaffected founders; leave-one-family-out (the family's own founders removed, N=63).
Recurrence (PON_AC) is an ARTIFACT filter; population rarity comes from gnomAD. Small variants
matched by exact allele after bcftools norm on both sides; SVs matched by the rule in section 5.
D53, D56, D60, D62.

## 4. Annotation and tiering
VEP 115 with LOFTEE (HC) and dbNSFP; gnomAD v4.1 joint GLOBAL allele frequency (JS decision).
Gene constraint: GeneBayes s_het joined by Ensembl gene id (GeneBayes carries no symbol). LoF
tiers (J. Guevara, rare-variant-pipeline tier_variants.py; Tier 3 added by JS): lof_t1 s_het >=
0.18; lof_t2 0.03 <= s_het < 0.18; lof_t3 the remainder INCLUDING genes without an s_het.
Missense tiers by n_flag = number of ClinPred >= 0.4298, AlphaMissense >= 0.9603, popEVE >=
0.9209, MPC >= 0.8947 (dbNSFP 5.3.1a; popEVE coverage 100%); miss_t1 = 4 flags, which is 1.12% of
all missense genome-wide. Rare = AF < 0.001. D59, D64.

## 5. Structural-variant validation and impact
De novo SV candidates re-genotyped in short reads with GraphTyper 2.7.2 genotype_sv across all
family members. Matching rule (JS): INS and DUP treated as one class; a match on the left
boundary, the right boundary, or the length (within 10%, minimum 500 bp), or reciprocal overlap
>= 0.5. Three-way outcome CONFIRMED / REFUTED / UNINFORMATIVE; short-read non-detection is never
scored as refutation. An identical "de novo" interval in more than one family is an artifact.
Coding impact by GENCODE v44 CDS intersection; breakends excluded. D61, D62, D63.

## 6. Repeat expansions
73 STRchive disease loci scored on the pathogenic motif with QC on spanning reads and purity. A
PATHOGENIC or INTERMEDIATE call whose allele equals the population-modal allele in the founder
panel is rendered NOT_CALLABLE (misconfigured threshold); benign calls are untouched. D55, D59.

## 7. Inheritance models
De novo: slivar trio model, gnomAD AF < 0.001, region mask, panel leave-one-out, exclusion of
variants shared by both siblings. Recessive, X-linked and compound heterozygous: slivar tags on
the tertiary VCF with the DDG2P (definitive|strong) allelic requirement enforced (a dominant call
must fall in a monoallelic gene, a recessive call in a biallelic gene); dominant AF < 1e-4;
recessive AF < 0.005 with nhomalt <= 2; X-linked AF < 1e-4; the unaffected sibling must not
carry the same qualifying genotype. D60, D63.

## 8. Association testing
Conditional logistic regression (survival::clogit, R 4.5.3): case ~ tier1 + tier2 + tier3 +
strata(family), 39 cases, 66 controls, 35 informative families, with genome-wide mean depth as
the nuisance covariate (cases 22.4x vs controls 24.4x; depth correlates 0.57 with detected
burden). Site-level low-quality count rejected as a covariate: collinear with depth (r = -0.73)
and circular. Transmission test: informative meioses (one heterozygous parent, other homozygous
reference) to affected offspring, binomial against 0.5. Both null. The transmission test is the
better-specified analysis because 65 of 66 controls are parents of the cases. D64.

## 9. Ascertainment
Measured against SPARK WES SynthDNM de novo calls: our probands carry 0.0000 miss_t1 de novo per
proband vs 0.0076 in 29,893 other SPARK probands (0/33 vs 0.7%), with no depletion of lower
tiers, indicating selection against likely-diagnostic variants specifically. Coding de novo
burden is therefore biased by construction; inherited burden is judged valid (JS). D64.

## 10. Clinical phenotype
SPARK phenotype release 2026-06-25: core_descriptive_variables and approximated_cognitive_
impairment. Predictors with coverage in lrWGS probands: cognitive_impairment_latest (parent-
reported ID diagnosis, 100%), language_level_latest (100%), derived_cog_impair (88%), SCQ (79%),
RBS-R (76%), Vineland ABC (61%); no IQ scores exist for any lrWGS proband. lrWGS vs SPARK: Fisher
exact and Mann-Whitney, with age-band matching for the age-accumulating ID label (D65).
Syndrome match: published feature sets mapped to SPARK fields, na_survey_logic coded as UNKNOWN (the dictionary defines it as not asked); absent only for an
explicit no or a blank checkbox under a positive gate, every proband scored on every syndrome, the carrier ranked on its own syndrome
with empirical p = fraction of probands scoring at least as high (D66). A like-for-like ranking (carrier vs probands assessable on the same features, weighted by published frequency) was evaluated and found unsupportable: comparison sets are conditioned on parental flagging and are small (n=11) or trivial (n=2); feature-level concordance is reported instead (D66 addenda).
