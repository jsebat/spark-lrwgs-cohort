# Phase-aware de novo mutation calling — Methods draft (13_phase_dnm)

Assembled from DESIGN.md decisions P1–P27 (one paragraph per decision, in pipeline order) for the de novo calling section
of the SPARK long-read WGS pipeline paper (DESIGN §5: pipeline paper, cohort as proof of concept). Cohort values are the
2026-09-14 state (thresholds 0.3.0, `harness_v2` models, P18 v6) and are marked *measured*. Citations marked **unverified**
have not been checked against the publisher's record and must be before submission. Upstream alignment and variant
calling are described once, in the repository-level methods (`docs/METHODS.md`, "Alignment and variant calling"), and
only referenced here.

## 1. Input data

Per-family outputs of the PacBio HiFi-human-WGS-WDL (v3.3.1, commit 477ef39) run per family on the cohort: pbmm2-aligned,
HiPhase 1.6.0 haplotagged BAMs (`HP`/`PS` tags), the family joint small-variant VCF (DeepVariant + GLnexus), sawfish 2.2.1
structural-variant VCFs with supporting-read lists, TRGT 5.0.0 tandem-repeat VCFs with spanning-read BAMs, and the
cohort-wide joint callsets (105 samples). Coverage is ~22–24× per genome (*measured*: per-haplotype depth median 10.0×;
88 % of depth haplotagged; three genomes at 8.3–10.0× total are reported separately in every cohort summary, R11). The
module operates on 33 complete-trio families (35 children, two quads); the two mother–child duos are excluded from
transmission, folds and cohort rates (P19). The module adds no aligner or caller (P17): every read-level quantity comes
from the three existing BAM types joined by read name.

## 2. Module 1 — Parent-of-origin phasing and transmission map

**Orientation (P2).** HiPhase read-based phasing is kept; no sample is re-phased. Each child phase block is oriented to a
parent of origin by a Mendelian vote at informative sites (one parent heterozygous, the other homozygous): the block is
`HAP1_PAT` or `HAP1_MAT` when the majority fraction is ≥ 0.95 over ≥ 20 informative sites, otherwise `AMBIGUOUS`; segments
of 10–19 sites are oriented only when unanimous. Phase-switch errors inside a block are located by an exact dynamic
programme over the position-ordered votes (at most three cuts, each raising the majority count by ≥ 2), so a lone
dissenting vote is treated as a genotype error and a run of two at a block edge as a switch; each segment is oriented on
its own, and a read's parent-of-origin label is keyed by (phase block, position). *Measured:* 9.9–16.4 k blocks and
1.9–2.2 M informative sites per child; Mendelian-inconsistent votes 0.06–0.11 % of informative sites; 313–367 located
switches per genome, ≈ 1 per 8 kb of phased heterozygotes; at the 18 imprinting control regions the orientation agreed
with the independent nearest-SNV parent-of-origin calls of the methylation module in 30 of 30 comparable cases.

**Transmission map (P3).** For each parent and parent phase block, the child's allele from that parent (known at child
homozygous sites, or at phased heterozygotes in an oriented segment) votes for the transmitted parental haplotype; votes
are segmented with the same change-point procedure (isolated runs ≥ 10 sites). A within-block change of transmitted
haplotype is either a meiotic crossover or a phase-switch error in the parent's read-based phasing, and the two are
indistinguishable at the VCF level; each change is therefore emitted as a candidate with its resolution interval and
classified from the parent's reads spanning the interval by allele concordance with the tagged haplotype (intact phase →
crossover; discordant reads → switch error). Crossovers between parent blocks are invisible, so counts are lower bounds.
*Measured (35 meioses per parent):* median 30 paternal and 43 maternal crossovers per meiosis excluding intervals with a
located child switch, ratio 1.43, against ~26 / ~43 in deCODE pedigrees (Halldorsson *et al.*, *Science* 2019 —
**unverified** here); ~80 % of change candidates remain ambiguous (gaps longer than a read, or disagreeing reads in
pericentromeric segmental duplications) and mark nearby candidates `NEAR_CHANGE_POINT` rather than being resolved. The
paper reports crossover candidates with the child-switch flag, not a recombination map.

**Tags (P4) and QC.** Module 1 rewrites no BAM: it emits block-level tables (child block → parent of origin; parent block →
transmitted / untransmitted per segment) that Module 2 resolves per read through the HiPhase read-name → block map.
Cohort QC gates: crossovers per meiosis within [15, 65], ≥ 90 % of parental heterozygotes inside resolved segments, total
depth ≥ 12×, tagged fraction ≥ 0.75, sex consistency. *Measured:* 31 of 35 children pass all gates; the four flags are all
`LOW_DEPTH` in the three low-coverage genomes.

## 3. Module 2 — Six-haplotype review of every candidate

**Candidates enter unfiltered (P5, P6).** Small variants: every family-VCF site where a child is heterozygous and both
parents homozygous reference, at any genotype quality (*measured*: ~25–35 k per child). SVs: every sawfish record where the
child carries an allele absent from both parents' genotypes. TRs: every TRGT locus where a child allele exceeds both
parental alleles by ≥ 1 motif unit. The lab's segmental-duplication / simple-repeat mask (10 % of the genome) and the
original pipeline's own call lists are carried as flags (`mask_overlap`, `source_tier`), so that rescue and demotion can be
reported against the existing pipeline (P18). Variants shared by both siblings of a quad are kept and flagged, since a
variant in both siblings and absent from both parents is the parental-mosaic signature.

**The six-haplotype matrix (P7, P17).** For each candidate the reads of the three haplotagged BAMs are partitioned into six
haplotypes — maternal 1/2, paternal 1/2, child maternal/paternal (the latter by the Module 1 orientation) — and counted as
REF, ALT or `AMB` (a read that is deleted, soft-clipped or carries a third base at the site). A haplotype is *observed*
when it carries ≥ k readable reads (k = 5 working, k = 3 fallback with a `LOW_HAP_DEPTH` flag); `AMB` reads count towards
depth but never observe an allele, and rows with > 30 % unreadable reads on any haplotype are flagged. SV alt-support comes
from sawfish's supporting-read lists plus an independent junction scan in all three BAMs (so a hom-ref parent's one or two
junction reads are counted); TR per-read repeat lengths come from the TRGT spanning-read BAMs (query length minus the two
flank lengths, within ±2 bp of the called allele in 98.4 % of reads, *measured*). Untagged reads count in child-level
support but in no haplotype row.

**Rule layer (P8).** A transparent, ordered rule set assigns each candidate a `phase_class`: `phase_conflict_artifact`
(alt on both child haplotypes, or on ≥ 2 parental haplotypes with no parental het genotype), `inherited_missed_in_parent`
(alt on the transmitted parental haplotype at ≥ max(3 reads, 30 %)), `parental_mosaic_transmitted`, `child_postzygotic_mosaic`,
`germline_DNM_phased` (alt confined to one child haplotype at ≥ 80 %, all four parental haplotypes observed and alt-free,
transmitted haplotype resolved), `germline_DNM_unphased` (same evidence with fewer than six haplotypes observed or an
unresolved transmission), or `inconclusive` with a reason. A per-read error allowance of max(1 read, 5 %) and a floor of
three child alt reads apply. All thresholds live in a versioned `thresholds.yaml` written into every output row, and the
rule layer can be re-run from the evidence tables without the BAMs. *Measured:* the rule layer enriches but does not
purify the raw set — the paternal fraction of rule-only "phased germline" SNVs is 0.66 (0.75 at GQ ≥ 20), the mixture
signature of real DNMs (~0.77) and artefacts (0.50) — which is why a rule-layer class alone is never a call (P15).

**Likelihood layer (P9) and mosaic floor (P10).** Per-haplotype counts are modelled as binomial under five hypotheses
(germline DNM on the child's alt haplotype; inherited from the transmitted haplotype; parental mosaic on it at fraction
m; child post-zygotic mosaic at fraction c; artefact at a per-class error rate ε estimated from cohort hom-ref sites),
with m and c integrated over Uniform(0.02, 0.5) and priors from expected rates; the germline posterior is `phase_score`
and the log10 likelihood ratio is reported against the best alternative. At ~10 reads per haplotype the nearest
alternative to a clean germline pattern is a low-fraction parental mosaic (posterior ~0.98), so the ratio is
depth-limited at cohort coverage. A mosaic class requires ≥ 15 reads on the relevant haplotype and ≥ 3 minority-allele
reads passing mapping and clipping checks; otherwise the row is `inconclusive: MOSAIC_UNDERPOWERED`. Mosaic sensitivity is
stated from spike-ins (below), never as a mosaic rate (R1).

**Spike-ins.** Germline, inherited-missed, child-mosaic and transmitted-parental-mosaic events of every class are planted
into the haplotagged reads of each child (and, where required, a parent) and run through the identical review, giving
class-exact recovery, parent-of-origin recovery and, for mosaics, a sensitivity curve. *Measured (35 children, 7,828
planted sites):* germline and inherited-missed events recover class and parent of origin at ≥ 0.91 exact / ≥ 0.98 lenient
with parent of origin 1.00; child mosaics at 15–30 % and transmitted parental mosaics at 10–25 % are recovered at 0.05–0.15,
the depth floor at ~10 reads per haplotype.

## 4. Module 4 — Per-class classifiers trained on synthetic de novo variants

**Construction (P11, P23, P24).** Following SynthDNM (Lian *et al.* 2021), positives are synthetic de novo variants: the
child of one family is paired with the parents of another, the candidate set of that synthetic trio is regenerated from the
cohort joint callsets, and each candidate is reviewed against the child's and the surrogate parents' haplotagged reads by
the same Module 2 code path as a real trio, without parent-of-origin or transmission labels (which do not exist for a
synthetic trio). Negatives are the raw, unfiltered candidates of the real trios (P13; contamination by true DNMs ≈ 0.3 %).
A feature is admitted to the classifier (`rf_safe`) only if it is a function of the child's evidence or of parental
evidence symmetric under swapping a parent's two haplotypes and swapping the two parents; child haplotypes are sorted
alt-first and never labelled paternal or maternal. Parent of origin, transmission, and population-frequency or
cohort-recurrence quantities are never features: the former are inferences made in real trios only, the latter are label
leaks under the synthetic construction (positives are inherited, hence common) and enter the call as rules instead (§5).
A presence-leak guard drops any feature whose missingness differs by > 0.5 between synthetic and real rows or that is
present in < 1 % of real rows; real and synthetic reviews are produced by the same code revision.

**Folds and models (P12, P14).** Families are split into five outer folds, each swap pair drawn inside a fold (so a
child's real-trio negatives and synthetic-trio positives never cross a boundary), the single blood-derived family always
with company; inner folds are built identically for hyperparameter selection. Three XGBoost classifiers (SNV/indel with
an indel flag, SV, TR; 55 / 38 / 41 features) are trained with nested cross-validation over five outer-fold seeds, with
isotonic calibration on inner out-of-sample predictions, and compared on identical folds with a no-phase ablation, a
phase-only ablation, a random forest and logistic regression. Attribution uses TreeSHAP summed by feature family (caller,
context, read-level, phase) on held-out rows (*measured*: phase share 0.17 / 0.29 / 0.68 for SNV/indel / SV / TR).

**Heuristic arms on the same harness (P22, P25).** The original pipeline's filters are implemented as scorers returning
the pipeline's actual pass/fail decision and a sweep of the rule's natural knob (GQ for the slivar trio rule and the SV
genotype rule, gain in bp/units for TR), so each heuristic has an operating point and a curve on the same rows and folds
as the classifier. Population and artefact filters are information both arms may use, so every comparison is reported with
and without that layer. Provenance: slivar (Pedersen *et al.*, *npj Genomic Medicine* 2021, **verified**); TRGT
(Dolzhenko *et al.*, *Nature Biotechnology* 2024, **verified**); sawfish, GraphTyper, STRetch and the GIAB
stratifications **unverified** here.

**Cohort annotation (P26).** Leave-one-family-out cohort allele counts and founder-panel recurrence are computed from the 65
unaffected founders excluding the child's and the surrogate parents' families, so real and synthetic rows see the same
kind of panel; gnomAD v4.1 allele frequency is a site property.

**Fold-quantile score (P15).** Because calibrated probability scales differ between fold models and seeds, the decision
score is `rf_q` = 1 − pass rate among the fold's held-out real candidates of the same variant class (SNV and INDEL
separately), averaged over the 25 fold models; a threshold on `rf_q` is a pass rate by construction.

**Frozen models and transfer (P21).** After cross-validation the three models are refit on all families with the selected
hyperparameters and released with training manifests, checksums and operating points (`models/FROZEN.md`). A frozen model
applied to a new cohort (`phase-dnm score`) computes `rf_q` against that cohort's own candidates, so the operating points
keep their pass-rate meaning; a frozen model transfers only to the same callers and versions at similar coverage, and the
swap-closed retraining workflow is part of the release.

## 5. Module 3 — Two-tier decision with a rule layer

Demotions and mosaics are decided first: `inherited_missed_in_parent` and `phase_conflict_artifact` rows are NO, mosaic
classes are NO and flagged. **Tier 1** (`dnm_call = YES`) requires `rf_q ≥ τ_q`, a germline review class
(`germline_DNM_phased` or `_unphased`; an `inconclusive` row cannot be called on the score alone), and the rule layer:
gnomAD AF < 0.001 (absent = rare; for SVs the long-read SV catalogue AF), leave-one-family-out founder-panel recurrence 0,
leave-one-family-out cohort allele count 0 (so sib-shared DNMs of a quad survive), and the region mask — a filter for
SNV/indel and TR, a flag for SV because 92 % of SV candidates lie in segmental duplications and simple repeats and the
original SV set never used the mask. **Tier 2** (`dnm_call = CANDIDATE`, reported for experimental validation with score,
phase class and parent of origin) covers `τ_q,tier2 ≤ rf_q < τ_q` under the same gate and rules, TR candidates additionally
requiring ≥ 3 motif units of expansion. Operating points (thresholds 0.3.0): SNV/indel 0.99 / 0.95, SV 0.99 / 0.97, TR
0.999 / 0.997, set as described in §6. Parent of origin, transmitted-haplotype support and crossover proximity are computed
for every row of the final tables after the decision, as inferences on real trios only. Final tables and sites VCFs carry
the score, tier, decision reason, phase class, six-haplotype counts, parent of origin and all flags per candidate.

**Large deletions are called without the classifier (P28).** For a structural deletion of at least 300 bp whose class
is decided by the interval evidence, the decision is deterministic and the classifier score is not used. The conditions,
all recorded per row, are: depth across the interval at or below 0.7 of flanking depth in the child; neither parent's
interval depth reduced (both at or above 0.85); junction reads at both breakpoints; heterozygous sites inside collapsed
towards zero relative to the flanks; and the population rule layer passed. Two measurements justify the split. First,
the pedigree swap cannot supply the training data: after the rarity gate the structural-variant positives contain 43
deletions between 10 and 50 kb and 13 above 50 kb per seed, about eight per held-out fold in the size class of the
cohort's pathogenic MECP2 deletion. Second, the classifier mis-ranks those events: adding the interval evidence to the
feature matrix moved that deletion's score down from 0.967 to 0.66 while leaving overall structural-variant ROC-AUC
unchanged at 0.922, because a depth ratio of 0.59 is unlike any positive in its training set. The same evidence scores
it 6 of 6 deterministically. Calls made this way carry `decision_reason = TIER1_SV_DEPTH` and are counted separately
from score-led calls throughout. Recall as a function of deletion size, for the classifier and for the deterministic
path, is measured on planted deletions from the extended spike-in planter, which removes the deleted haplotype's reads
inside the interval and clips crossing reads into junction reads, because no HiFi read spans a 35 kb event.

**Concordance with the original pipeline (P18).** The original de novo sets (slivar trio rule + gnomAD + mask + founder
recurrence for small variants; the unfiltered sawfish de novo list; the TRGT expansion table) are matched to the final
tables (small variants by position and alleles with representation-independent normalisation; SVs by type with breakpoints
within 500 bp or ≥ 50 % reciprocal overlap; TRs by locus) and reported as concordant, original-only (by phase class and
decision reason), recovered at tier 2, module-only, with per-proband rates before and after and the paternal fraction of
phased calls as the guard against buying sensitivity with artefacts (R9).

## 6. Setting the operating points and validation

**Exome-confirmed exonic truth (P27).** For the 33 children present in SPARK iWES v3 (DeepVariant pVCF), every exonic
long-read candidate is labelled from the exome trio genotypes: positive when the child is 0/1 with GQ ≥ 20, depth ≥ 10
and ≥ 2 alt reads and both parents are 0/0 with GQ ≥ 20, depth ≥ 10 and no alt read; negative when a parent carries the
allele or the child is 0/0 at depth ≥ 20; otherwise not evaluable (hom-alt children are not evaluable). This is a
confirmation of our own candidates on an orthogonal platform, not a curated truth list, and it is exonic only (~1 % of the
genome; R10). *Measured:* 56 positives and 2,926 negatives. The SNV/indel operating points were chosen on a sweep of τ_q
with the rule layer applied: tier 1 (0.99) recovers 41 of 56 with 0 false positives among 2,926 and yields 73 calls per
proband genome-wide (64 SNV + 7 indel) at a paternal fraction of 0.770; tier 2 (0.95) adds 7 positives against 10 false
candidates (cumulative recall 0.86, cumulative exonic FDR 0.17). The original pipeline's set recovers 34 of 56 at 2 false
positives with 38 calls per proband. The bare classifier without the rule layer does not beat the original set at any
threshold (recall 0.72 at precision 0.65 at τ_q 0.99), which is the empirical case for the rule layer.

**SV and TR operating points.** No orthogonal truth exists for these classes; thresholds were set from planted-germline
recall and planted-inherited pass rate on the spike-ins, together with calls per proband against expected rates, the
paternal fraction of phased calls, and overlap with the original sets. *Measured:* TR purity tracks the threshold smoothly
(paternal 0.80 at 7 per proband, 0.78 at 15, 0.69 at 37, 0.59 at 106) with planted recall 0.99 down to τ_q 0.998; SV tier 1
gives 10 calls cohort-wide (0.3 per genome, in line with published de novo SV rates of ~0.2–0.3 per genome, **unverified**
figure), tier 2 59. Planted-germline SV recall plateaus at 0.63 for the score alone; the six-haplotype demotion, not part
of the score-only spike arm, removes inherited SVs on real data (3,782 phase-conflict demotions in the original SV set).

**Genome-wide precision and GIAB (P16).** Exonic precision is an upper bound for the genome. The GIAB Ashkenazi trio
(HG002/HG003/HG004) is held out of every fold, calibration and threshold and is run through the identical WDL and module at
full depth (48 / 46 / 36×, *measured*) and, depth-matched, after subsampling the unaligned reads to ~23×; calls outside the
HG002 benchmark within confident regions give the genome-wide false-positive estimate that the exome cannot. Mask-stratified
recall and FDR (inside vs outside the lab mask; GIAB stratifications) are reported for every arm.

## 7. Reporting guards (DESIGN §4)

Mosaic events are reported as candidates with read counts, not rates (R1); sib-shared DNMs from two quads illustrate the
parental-mosaic class, they do not validate it (R2); SV and TR thresholds rest on ~33 trios and spike-ins and are labelled
provisional until GIAB (R3); the paternal fraction validates the cohort call set, not individual calls (R4); classifier
performance on synthetic held-out labels is a ceiling and the headline numbers are external-truth numbers (R7); the
per-proband rate is reported before and after against the expected ~60–80 DNMs per genome (**unverified** figure) with the
paternal fraction alongside (R9); exome non-detection never refutes a long-read call (R10); low observability in the three
low-coverage genomes is a coverage statement, and their families are reported separately (R11).

## 8. Software

`phase_dnm` (Python 3.11; pysam, numpy, pandas, XGBoost 3.2.0, scikit-learn 1.9.1), 108 unit tests, Snakemake workflow
over per-family Slurm scripts, versioned thresholds and feature registry, no identifiers in code or configuration.
Repository: https://github.com/jsebat/spark-lrwgs-cohort, directory `13_phase_dnm`.

## Results summary for the paper (cohort state 2026-09-14)

| item | value |
|---|---|
| Harness ROC-AUC, five seeds (RF / no-phase / phase-only vs best heuristic) | SNV/indel 0.9967 / 0.9958 / 0.953 vs 0.887; SV 0.9945 / 0.990 / 0.872 vs 0.61; TR 0.892 / 0.840 / 0.881 vs 0.59 |
| Phase attribution (TreeSHAP share) | 0.17 / 0.29 / 0.68 |
| Tier-1 calls per proband | 73 SNV/indel (64 + 7), 0.3 SV, 15 TR |
| Tier-2 candidates per proband | 25 SNV/indel (8 + 17), ~1.4 SV, 3 TR |
| Paternal fraction, tier 1 phased | 0.770 SNV/indel (original set 0.780); 0.703 TR |
| Exome truth (56 positives / 2,926 negatives) | tier 1: 41 recovered, 0 FP; tiers 1+2: 48 recovered, 10 FP (FDR 0.17); original set: 34 recovered, 2 FP |
| Concordance with original SNV/indel set (1,335) | 1,058 concordant at tier 1, 66 at tier 2, 211 lost (88 below τ, 96 demoted, 14 rules, 13 other) |
| Original SV set (8,129) | 3,782 phase-conflict, 3,829 below τ, 205 inherited-missed, 239 mosaic-flagged; 2 concordant |
| Long insertions in the original small-variant set | 118 (9 %) are inherited insertions with parental read support missed by the joint genotyper |
| Crossovers per meiosis (median) | 30 paternal / 43 maternal |
| Mosaic spike-in sensitivity at ~10 reads per haplotype | 0.05–0.15 |
