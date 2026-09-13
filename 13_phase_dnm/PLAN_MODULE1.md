# Module 1 build plan — phasing orientation, transmission map, QC

## Status
- **MODULE 1 GATE PASSED (2026-09-12).** `phase-dnm phase-qc` over all 33 complete-trio families / 35 children:
  **31 PASS**; the four flags are all `LOW_DEPTH` (children 9.4× and 11.5×, fathers 8.3× and 11.9×; R11) — every
  orientation, transmission, crossover and sex-consistency gate passes. Cohort table `cohort_phase_qc.tsv` +
  `cohort_phase_qc.summary.json` on the filer; per child `<child>.phase_qc.json`. Deliverables for M2 per child:
  `orientation.tsv` (read → parent of origin by block + position), `transmission.tsv` (parent haplotype →
  transmitted/untransmitted by block + position), `changepoints.resolved.tsv` (with `child_switch_in_interval`),
  and per sample `hapdepth.tsv.gz`. `haplotag --export-bam` (IGV) is deferred; the tables are the labels (P4).
- **Module 2 started (2026-09-12): shared `CandidateRecord` + unfiltered generators for SNV/indel, SV, TR
  (`phase-dnm candidates`), run cohort-wide (array 54273350, 33/33, ~100 s per family).** Per child: small variants
  31,757 raw (GQ ≥ 20: 3,008; GQ ≥ 30: 609) → P13 contamination 0.22 % / 2.3 % / 11 %; SV 627 (INS 374, DEL 219, BND 46);
  TR 38,921 at ≥ 1 motif unit (990 at ≥ 3 units). Details in DESIGN §0.2. Next: the six-haplotype extractor
  (`hapmatrix.py`) with the three alt-support adapters, on a read-level minitrio built with the cached WDL images.
- **M2 core built and smoke-tested on the first quad (2026-09-12; `hapmatrix.py`, `readers.py`, `review.py`,
  `phase-dnm review`, 38 tests).** Smoke #1 (1,500 per class per child) found: base lookup 0.45 s/candidate
  (aligned-pairs → CIGAR walk), SV junction test accepting any SA tag (448/627 conflicts → event-matched signature),
  TRGT `AL` tag misread as a length (it is the allele index; per-read length = `query_length − FL0 − FL1`, verified
  ±2 bp in 98.4 % of 4,385 reads), single-read "germline" calls (→ `min_alt_reads` 3, all four parental
  haplotypes observed + alt-free for `germline_DNM_phased`), 35 % unoriented from exact-cover label lookup
  (→ nearest segment within 30 kb). Smoke #2 after the fixes: **0.06 s per small variant, 0.1 s per SV, 3 ms per
  TR → ~35 min per child, ~70 min per quad**, 397 MB; ≥ 1 child ALT read in 98 % INDEL / 59 % SNV / 93 % SV / 99 %
  TR rows; 82 % of INDEL candidates with a resolved transmitted haplotype carry alt reads on it (systematic
  homopolymer errors, correctly rejected as inherited/conflict); raw SNV candidates sit in low-mappability
  regions (tagged depth median 7, `hap_obs_k5 = 0` in 42 %). Remaining TR defect — 56 % "alt on both child
  haplotypes" from a read-length tolerance as wide as the allele separation — fixed by capping the tolerance at
  half the distance to the nearest competing allele. Cohort review array follows.
- **Features + likelihood cohort arrays (2026-09-13, arrays 54278875 / 54278919, 33/33 + 33/33 COMPLETED, 1 core,
  14–30 s per child per class for features, 23–42 s for likelihood; 105 feature tables + 105 `rf_safe` matrices +
  105 `*.evidence.lik.tsv`).** Feature coverage report per class: SNV/indel 56/87 applicable features produced, SV 39/80,
  TR 42/78. **Never produced** (extractors do not exist yet; must be written before M4 trains): the read-level quality set
  (`child_alt_mapq_mean`, `child_alt_mapq0_frac`, `child_alt_rq_mean`, `child_alt_nm_rate`, `child_alt_clip_frac`,
  `child_alt_readlen_median`, `child_alt_readpos_frac_median`), the block-geometry set (`c_block_len_log10`,
  `c_dist_block_edge_log10`, `c_n_informative_1kb`, `c_local_rephase_agreement`, `c_hp_rephase_conflict`), SV
  `bp_precision` / `c_sv_hap_depth_change_*` / `c_sv_junction_hap_concentration`, TR `c_tr_expanded_one_hap_only`, and
  `child_DP`/`child_GQ` for SV and TR (sawfish/TRGT FORMAT differ from DeepVariant's).
- **Rules vs likelihood on the phased-germline SNV class (2026-09-13, job 54279056, 35 children, rows with a resolved
  parent of origin; paternal fraction as the purity read-out, P8):** rule `germline_DNM_phased` n = 6,400, paternal
  0.571 (GQ ≥ 20: 0.668); `phase_score ≥ 0.9` n = 9,521, 0.556 (0.659); **intersection n = 2,124 (median 59 per
  child), paternal 0.713 (GQ ≥ 20: 1,692 rows, 0.762)**; rule-only 4,276 at 0.501 and likelihood-only 7,397 at 0.511 —
  i.e. each layer alone admits a disjoint noise set with no parent-of-origin signal, and the intersection is where the
  DNMs are (a 0.71 paternal fraction against an expected ~0.8 for true SNV DNMs puts the intersection at roughly 70 %
  real, ~40 real phased SNV DNMs per child). Likelihood-only rows are mostly rule `germline_DNM_unphased` (7,442) and
  `inconclusive` (4,918): the posterior is a statement about hypotheses given whatever depth exists and does not
  encode observability (haplotypes below k reads); the rule layer does (P7/P8). Rule-phased rows have median
  `phase_score` 0.049; best alternative there is `parental_mosaic` in 65 % and `inherited` in 30 %. Diagnostic of
  what the likelihood rejects: see the next entry. Consequence for P15: the final call uses **both** layers as
  designed (rule observability + posterior), and ε/δ are to be re-estimated from cohort hom-ref sites rather than
  the 0.01/0.03 defaults before the layer is used as a feature.
- **Why the likelihood rejected the rule's phased-germline calls (2026-09-13, jobs 54279111 / 54279173):** among rule
  `germline_DNM_phased` SNVs with `phase_score` < 0.5 the winning hypothesis was *inherited* (3,143), and the rejected rows
  differ from the accepted ones in caller GQ (median 8 vs 41), DP (14 vs 23), alt reads on parental haplotypes (73 % vs
  17 % have any) and ambiguous reads (`C2_amb` 3 vs 0). Row-level view: the transmitted haplotype had `dp` 7–8 but
  `alt/ref` = 1/0 — six of seven reads UNREADABLE at the site (deleted / clipped: low-GQ SNVs next to indels). The rule
  counted depth as observation and one alt read as within the allowance → germline; the likelihood saw 1/1 alt on T →
  inherited (or, with 0 readable reads, the prior alone → inherited). Column semantics were consistent (`T alt mismatch
  rule vs lik: 0`). **Fix (thresholds 0.2.0, P7 clarified):** observability, allowance and mosaic floors on READABLE
  reads (`Row.n`); features `c_amb_frac_hapA`, `p_amb_frac_max`; flag `AMBIGUOUS_READS`; `phase-dnm reclassify` re-runs
  the rule layer from the evidence columns without BAMs (`*.evidence.review.tsv` immutable, `*.evidence.tsv` working);
  `workflow/m2_reclassify_family.sb` chains reclassify → likelihood → features. 3 new tests (54 passing).
- **Reclassify cohort-wide under thresholds 0.2.0 (2026-09-13, smoke 54279292 + array 54279293, 33/33 COMPLETED, ~20 s per
  child per step; `*.evidence.review.tsv` kept for all 105 tables) and the comparison re-run (job 54279309):** rule
  `germline_DNM_phased` SNVs with a resolved parent of origin **6,400 → 2,747 (median 72 per child), paternal fraction
  0.571 → 0.664 (GQ ≥ 20: 1,858 rows, 0.746)**; their median `phase_score` 0.049 → **0.976** (p10 0.74); 2,057 of them
  (75 %) are inside the likelihood ≥ 0.9 set (intersection paternal 0.722; 0.762 at GQ ≥ 20). The rule-only remainder is
  690 rows at 0.488 (was 4,276 at 0.501) — what is left of the disagreement is the depth-limited parental-mosaic
  alternative (best alternative for 2,596 of 2,747 rule-phased rows), not a definitional mismatch. Likelihood-only stays
  7,463 at 0.511 (unobservable haplotypes; the posterior does not encode observability). Feature coverage after the fix:
  58/89 applicable SNV/indel features.
- **Spike-in harness built and smoke-run (2026-09-13, `sim/edit.py` + `sim/spike.py`, `phase-dnm spike plan|apply|evaluate`,
  `workflow/m2_spike_family.sb`; smoke job 54279524 on family 1, N_PER = 4, ~5 min; cohort array 54279525, N_PER = 12).**
  Plants SNV / indel (1–12 bp) / SV (DEL 200–1,500, INS 100–400 bp) / TR (3–8 motif units) into the REAL haplotagged
  reads (CIGAR surgery on plain tuples, pysam only for I/O), on a child haplotype whose parent of origin and the parent's
  transmitted haplotype come from the M1 tables at the site; four scenarios G / CM (15, 30 %) / PM (10, 25 %) / IM; sites
  pass a read-level QC (>= 5 readable reads per haplotype per sample, unanimous base, indel-free window, outside the mask).
  Smoke: 82/96 planned sites placed (skips: haplotype QC 2,642, window 2,071, mask 1,659 of 6,751 tries); the review ran
  on the spiked slices unchanged. **Recovery among observable sites: G class 1.0 and parent of origin 1.0 for all four
  classes (posterior 0.98); IM 1.0 (posterior of germline 0.00); CM 0.0–0.25 — every miss is `TOO_FEW_ALT_READS` (15 % of
  ~10 reads = 2 alt reads) or `MOSAIC_UNDERPOWERED` (< 15 readable reads on A); PM 0.0 — m = 0.10 leaves 0–1 alt reads on
  T (within the allowance → germline, posterior 0.72–0.98), m = 0.25 gives 2–4 alt reads → `PARENTAL_ALT_LOW_DEPTH` or
  `inherited_missed_in_parent`.** This is P10 measured: at ~10 reads per haplotype a transmitted parental mosaic below
  ~25 % is indistinguishable from germline and above it from missed inheritance; the mosaic classes are reported with
  that caveat and their posteriors, not as calls (DESIGN P10). The array's cohort-wide summary follows.
- **Module 3 built (2026-09-13): `integrate.py` (final unfiltered table per child and class group: README §2.3 core →
  class columns → registered feature vector; P15 decision with the provisional phase-only mode until M4; demotion; mosaics
  flagged, never YES), `io/vcfinfo.py` (PDNM_* INFO block + sites VCF), `concordance.py` (P18 against the baselines copied
  from Lustre the same morning), CLI `integrate` / `concordance`, `workflow/m3_integrate_family.sb` + `m3_concordance.sb`;
  5 new tests (71 passing).** Cohort run next: integrate array over the 33 families, then the concordance job.
- **M3 cohort run (2026-09-13, smoke 54280056 + array 54280057, 33/33, ~30 s per family; concordance 54280058):** 105 final
  tables + sites VCFs in `$LRC/phase_dnm/final/`, provisional phase-only mode. Per child YES: SNV/indel median 65 (16–131),
  SV 3 (0–8), TR 88 (15–158; stutter-dominated raw set — the TR YES set is not a DNM count until M4/P8 class thresholds).
  **P18 first read-out, SNV/indel (35 probands):** concordant YES 666; original-only 546 — `germline_DNM_unphased` 484,
  `inconclusive` 35, demoted 13 (12 phase-conflict, 1 inherited), mosaic 4, low posterior 10; module-only 1,631 (all
  UNFILTERED tier, 0 in the mask); per-proband median 38 → 65; paternal fraction of module YES 0.72; **120 original calls
  never entered the unfiltered candidate set** (investigated below). SV against the pipeline's unfiltered 7,736-row list:
  concordant 46, original-only 8,083 (phase-conflict 3,782, inconclusive 3,119, unphased 719, inherited 205, mosaic 239),
  module-only 69, per-proband 247 → 3, paternal fraction 0.643. TR concordance OOM at 8 GB → rerun at 32 GB.
- **The 117 "unseen" originals resolved (2026-09-13):** all are insertions >= 30 bp (Alu-like sequence) whose ALT the tiered
  table truncates at 30 characters; the module has the same events at the same positions — classed mostly
  `phase_conflict_artifact` / `inconclusive`, i.e. not clean heterozygous insertions on one child haplotype. Concordance
  now matches a same-position indel whose baseline core is a >= 20 bp prefix of the module's. Paper-relevant: ~9 % of the
  original de novo SNV/indel set are long insertions that the read-level review does not support as germline DNMs
  (candidates for polymorphic mobile-element insertions missed in the parents' genotypes); to be examined against the
  parental read evidence (`t_alt_reads`, `p_max_alt_any_hap`) before any claim.
- **Spike-in cohort array (2026-09-13, 54279525, 32/32 COMPLETED + smoke; 34 children, 7,828 planted sites, N_PER = 12):**
  recovery among observable sites — **G: class-exact SNV 0.959 / INDEL 0.940 / SV 0.927 / TR 0.913 (lenient, i.e. incl.
  `germline_DNM_unphased`, 0.995 / 0.999 / 1.000 / 0.979), parent of origin 1.00 (TR 0.995); IM: 1.000 / 1.000 / 0.997 /
  0.984; CM (15–30 %): 0.138 / 0.137 / 0.145 / 0.099; PM (10–25 %): 0.105 / 0.128 / 0.067 / 0.054.** The mosaic numbers
  are the P10 depth floor measured cohort-wide and go in the paper as the module's mosaic sensitivity at ~10 reads per
  haplotype; the germline/inherited numbers are the review's recovery on real reads with known truth. Tables per child in
  `$LRC/phase_dnm/evidence/<FAM>/spike/<child>/`.
- **M4 design fixed before code (2026-09-13): P23 synthetic trios through the same chain (no labels, positives subsampled),
  P24 rf_safe-only training matrix with synthetic/real labels and XGBoost + RF + LR baselines, P25 heuristic arms as
  `pass` + `sweep_score` scorers.** Build order inside M4: heuristics (pure, testable) → synthetic-trio candidates from the
  cohort callsets → review/features for synthetic trios (array) → folds + nested CV → rf_probs → integrate (rf+phase) →
  concordance → harness table.
- **M4 build (2026-09-13, commits 7c56885 → b285375): swap-closed folds + within-fold pedigree swaps (`train/folds.py`,
  `train/swap.py`, `phase-dnm swap`: 33 families → 5 folds of 6–7, 35 synthetic trios per seed, 5 seeds, blood family kept
  with company); synthetic trios through the same chain (`m4_synthetic_trio.sb`: bcftools trio extraction from the cohort
  callsets 226 s, candidates thinned to 3,000 per variant class, review without label tables, features); nested CV +
  harness (`train/nested_cv.py`, `eval/harness.py`, `phase-dnm train`: XGBoost with inner-fold grid + isotonic calibration,
  no-phase / phase-only ablations, RF + LR baselines, the heuristic sweeps H1–H3 with operating points, +P demotion
  re-ranking on synthetic labels, rf_probs for every real candidate from the model that held its family out, τ at
  held-out precision 0.95/0.80, frozen per-class models with manifests; `integrate --tau-json`); P26 annotation
  (`annotate.py`: gnomAD v4.1 via the WDL's slivar gnotate zip, leave-one-FAMILY-out founder counts from the cohort
  BCF/SV VCF excluding both the child's and the parents' families, sib-shared for the quads; `m4_annotate.sb`,
  `m4_annotate_synthetic.sb`); registry clean-up (B-block read-quality names superseded by per-haplotype D-block
  summaries incl. new `c_alt_mapq0_frac`/`c_alt_supp_frac`/`c_alt_readlen_median`/`c_alt_rq_mean` — need a re-review to
  fill; phase-segment geometry `c_block_len_log10`/`c_dist_block_edge_log10` and crossover distance from the M1 tables).
  Smoke synthetic trio (54280890): SNV/indel stage 606 s for 6,000 rows. Seed-0 array (34 trios, %10) queued behind it;
  annotation jobs 54281886/7 running. 90 tests + 1 skipped locally (the ML backend's DLLs are blocked on this Windows
  host; the CV test runs on Expanse).
- **First harness table, SV, seed 0 (2026-09-13, job 54282295; 22,434 real + 105,000 synthetic rows, 64 features):** RF
  0.998 ROC-AUC (no-phase 0.994, phase-only 0.977; sklearn RF 0.999, LR 0.992); heuristic H1 (child carries, parents 0/0,
  GQ sweep) 0.61 with an operating point that passes every candidate (the original SV de novo rule IS the candidate
  definition); H2 (mask + founder-panel recurrence) 0.20. **Two things this table taught before it can be quoted:** (i)
  I had registered `gnomad_af` / `cohort_AC_loo` / `pon_founder_recurrence_loo` as classifier features that morning —
  under the SynthDNM construction a positive is an inherited, common variant, so population frequency separates the labels
  by construction (SynthDNM's universal set has none; JS's design was right, mine reintroduced the leak) → `rf_safe: false`,
  heuristic arms and final table only, training rerun; (ii) the same effect makes the population-filter heuristic arms
  (H2/H3) unfair on synthetic labels — the fair heuristic comparison on synthetic labels is H1, and H2/H3 vs RF+population
  belong on the external truth (P22's 2×2, now measured rather than anticipated). Also fixed: per-child matrices with
  different column sets (rf columns came from the first row's class), TR candidate ids shared by two alleles of one locus
  (now `TRID:a<idx>`; old tables deduplicated at load), and τ chosen at a fixed pass rate on REAL held-out candidates
  (0.1 % / 1 %) instead of precision on the positive-heavy mix (which gave τ = 0). SNV/indel and TR runs rerun after the fixes.
- **Seed-0 harness table, all classes (2026-09-13, jobs 54283205/6/7 after the fixes; ROC-AUC on identical swap-closed
  family folds, real raw candidates = 0, synthetic = 1; population-frequency features excluded from every classifier arm):**
  SNV/indel (698,885 real ≤ 20k/child + 210,000 synthetic, 68 features): RF 0.986, no-phase 0.997, phase-only 0.975,
  sklearn RF 0.999, LR 0.984; H1 slivar sweep 0.887 with operating point TPR 0.48 / FPR 0.0078; H2/H3 0.874 (synthetic-label
  caveat). SV (22,434 + 105,000, 61 features): RF 0.998, no-phase 0.990, phase-only 0.977, RF/LR 0.999/0.989; H1 genotype
  rule 0.61 (its operating point passes every candidate). TR (700,000 + 105,000, 61 features): **RF 0.977, no-phase 0.835,
  phase-only 0.975**, RF/LR 0.990/0.851; H1 family rule 0.59 (TPR 0.005), H2 cohort rule 0.59. Read-out: phase features are
  decisive for TR and additive for SV; for SNV/indel the caller features already saturate on synthetic labels (no-phase >
  full is inner-grid variance on a saturated task; the sklearn RF at 0.999 says the ceiling is the label construction, not
  the features) — the SNV/indel phase benefit has to be shown on the external truth (WES-confirmed exonic calls, spike-ins,
  P22). "+P(demote)" lowers every arm on synthetic labels, as expected from the construction (P22 caveat). τ (0.1 % pass
  rate on held-out real candidates): SNV/indel 0.348, SV 0.989, TR 0.258; τ_rescue (1 %): 0.054 / 0.764 / 0.059. Frozen
  seed-0 models + manifests in `$TRAIN_DIR/harness/models/`; rf_probs for every real candidate in `harness/rf_probs/`.
  Next: integrate in rf+phase mode → P18 re-issue; seeds 1–4 for the mean ± range.
- **Seed-0 RF / phase-only numbers above are INVALID — presence leak found 2026-09-13 (JS asked why no-phase beat the full
  model):** per-fold AUC showed the gap was one fold, then one child (full 0.843 / no-phase 0.979 / phase-only 0.685 for
  that child; the other six 0.9998–1.000). Not coverage (25.9×), not the surrogates (the deepest parents), not the feature
  marginals (identical). Cause: the four read-quality columns added to the D block that morning (`c_alt_rq_mean`,
  `c_alt_readlen_median`, `c_alt_mapq0_frac`, `c_alt_supp_frac`) are filled in 34/35 synthetic trios (reviewed after the
  code pull) and in 0 % of real rows (all real reviews predate it); the one synthetic trio reviewed before the pull (the
  smoke) is the collapsing child. The model learned presence ⇒ positive. The no-phase arm (immune) and the heuristic
  arms stand; RF and phase-only are rerun with a **presence-leak guard** in `load_matrices` (drop a feature whose
  non-missing fraction differs by > 0.5 between classes or is < 1 % present in real rows; dropped features logged in
  `cv_report`). Process rule: synthetic reviews and real reviews must come from the same code revision — the final cohort
  run re-reviews the real trios so the four columns are filled on both sides and the guard keeps them.
- **Corrected seed-0 harness table (2026-09-13, jobs 54283814/5/6, presence-leak guard active; ROC-AUC, identical
  swap-closed family folds):** SNV/indel (51 features after the guard): RF 0.997, no-phase 0.997, phase-only 0.951,
  sklearn RF 0.995, LR 0.984; H1 slivar 0.887 (TPR 0.48 @ FPR 0.0078), H2/H3 0.874. SV (34): RF 0.994, no-phase 0.990,
  phase-only 0.864, RF/LR 0.994/0.988; H1 0.61. **TR (38): RF 0.890, no-phase 0.835, phase-only 0.877, RF/LR 0.882/0.840;
  H1 0.59, H2 0.59** (the earlier 0.977 was the leak). Read-out: on synthetic labels the phase block is neutral for
  SNV/indel (caller features saturate), additive for SV, and the main signal for TR; the classifier-vs-heuristic gap is
  large in every class. The guard dropped 17/27/23 columns, all with 0 % presence on both sides (the never-produced
  annotation/context features) — none of them a real leak this time. τ (0.1 % real pass rate): 0.974 / 0.998 / 0.809;
  τ_rescue (1 %): 0.432 / 0.983 / 0.420. Seeds 1–4 synthetic trios: 140/140 COMPLETED; annotation 135/140.
- **P18 re-issued in rf+phase mode with the corrected seed-0 classifier (2026-09-13, integrate array 54284892, concordance
  54284893; τ at a 0.1 % pass rate on held-out real candidates):** SNV/indel concordant 621, original-only 714 (614 below τ,
  96 demoted, 4 mosaic), module-only 1,722; per proband 38 → **52**; paternal fraction of module YES 0.696; 64 % of YES come
  through the rescue branch (rf ≥ τ_rescue 0.43 with `germline_DNM_phased` and six haplotypes observed), 36 % through
  rf ≥ τ 0.974. SV: 9 / 8,120 / 14, per proband 247 → 0 (23 calls cohort-wide: τ 0.998 lets ~0.6 of ~600 candidates per
  child through — the 0.1 % rule is class-blind, and for SVs the expected true count is that order of magnitude anyway; to be
  set on external truth, P15). TR: 28 / 190 / 2,134, per proband 6 → 57. The SNV/indel per-child YES range is 10–480 — one
  outlier child under investigation. **Bug found in the same pass:** the "5-seed" training jobs were five copies of seed 0
  — `--export=ALL,SEEDS=0,1,2,3,4` is split by sbatch at the commas; `m4_train.sb` now takes the seed list as a positional
  argument and the run is resubmitted into `harness_5seed/`.
- **Per-child YES outlier → per-fold scale (2026-09-13):** the two children with 480 / 220 YES (the SPARK quad) were 449 / 163
  rf ≥ τ calls with rule class `inconclusive`/`unphased`, low child GQ (median 12), not sib-shared (2/480), not recurrent
  (62/480 in any other child). Per fold: fold 0's held-out model put 989 real rows ≥ τ (0.46 %) across all seven fold-0
  families vs 4–81 (≤ 0.04 %) in the other folds, with identical hyperparameters and per-fold AUC 0.9965 — the calibrated
  probability scale, not the ranking, differs between fold models, and one pooled τ then passes 10× more rows in one fold.
  Fix: `rf_q` fold-quantile score (`train/rescore.py`, `phase-dnm rescore`, `integrate --score-column rf_q`), τ as a pass
  rate (0.999 / 0.99) by construction; fold models are now saved per seed/fold (`fold_models/`), which the P27 external arm
  also needs. Rescored integration + concordance queued behind the fold-model training (harness_seed0fm).
- **P27 external-truth arm, spike-ins, seed-0 fold models (2026-09-13, job 54285056, 105 children×class tables):** planted
  germline DNMs (positives) vs the children's raw candidates + planted inherited-missed (negatives), scored by the model that
  held the family out. ROC-AUC RF / RF+P(demote) / **RF+P(demote+rescue)**: SNV/indel 0.986 / 0.989 / **0.998**; SV 0.988 /
  0.996 / **0.999**; TR 0.992 / 0.996 / **0.996** — on real reads with known truth the rescue branch is measurably worth
  having, which synthetic labels could not show. Recall at the probability τ was 0 for SNV/indel and SV (the inflated
  probability scale that `rf_q` replaces; rerun with `rf_q` pending) and 0.83 for TR; mosaic sensitivity at τ: TR
  parental mosaics 0.31, child mosaics 0.02, SNV/SV 0. Heuristic arms are **not evaluable on spike-ins** for SNV/indel and
  SV (planted candidates carry no caller genotype/GQ — the "caller" is the spike); the TR family rule works from allele
  lengths: 0.975 vs RF 0.992. Recorded as such; the WES-confirmed exonic set is the external truth for the slivar arms.
- **External arm rerun on rf_q (2026-09-13, job 54285805):** AUCs unchanged (ranking is scale-free); recall of planted
  germline DNMs at τ_q = 0.999: SNV/indel **0.0**, SV 0.22, TR **0.89**; mosaic sensitivity at τ_q: TR parental 0.38 / child
  0.01, SV child 0.50 / parental 0.04, SNV/indel 0. Why SNV/indel recall is 0 despite AUC 0.986: **planted candidates carry
  no caller features** (GQ, PL, AR, AB are empty — the "caller" is the spike), so the full small-variant model, which leans
  on the caller block, scores them wherever NaN routing sends them; TR is phase-driven and unaffected. Spike-ins are a
  fair external truth for the read/phase arms and a handicap for the full SNV/indel model — the WES-confirmed exonic set is
  that model's external truth (P27). The ablation fold models (phase-only, no-phase) are now saved so the fair spike-in
  comparison can be run; the heuristic arms are reported as not evaluable on spike-ins (no caller genotype).
- **τ_q sweep on the real cohort (2026-09-13, seed-0 rf_q, demoted rows excluded; YES/child median, paternal fraction of
  resolved YES):** SNV/indel 0.999 → 17 / 0.666; 0.998 → 38 / 0.674; 0.997 → 61 / 0.676; 0.995 → 103 / 0.663; 0.99 → 179 /
  0.629. SV: 0.999 → 0 (0–3; 5 resolved); 0.99 → 1 (41 resolved, 0.51). TR: 0.999 → 31 / 0.682; 0.998 → 60 / 0.634; 0.99 →
  277 / 0.590. Reading: the classifier's top of the list is ~55–60 % pure by the parent-of-origin read-out over the whole
  SNV/indel range (the phased-germline six-haplotype subset is purer: 0.71–0.76), so the P15 combination — rescue of phased
  germline rows below τ, demotion above it — is doing real work; the yield-matched point (~60 YES/child ≈ the expected DNM
  count) is τ_q ≈ 0.997. **Provisional τ_q until the WES external truth: SNV/indel 0.997, SV 0.999, TR 0.999; τ_rescue,q 0.99.**
- **P14 attribution, seed-0 fold models (2026-09-13, job 54285775, TreeSHAP |SHAP| share by family on held-out rows):**
  SNV/indel real rows caller 0.74 / context 0.09 / **phase 0.17** (synthetic 0.20; classifier-called real rows 0.21); SV
  0.68 / 0.03 / **0.29** (synthetic 0.30); TR 0.18 / 0.14 / **0.68** (synthetic 0.61; called rows 0.57). Top features:
  SNV/indel `site_qual`, `max_PL0`, `child_PL0`, `cpg_context`, `p_max_alt_hap_frac`; SV `site_qual`, `child_sv_support`,
  `c_alt_hapA`, `max_parent_sv_support`, `p_max_alt_any_hap`; TR `p_max_alt_hap_frac`, `delta_motif_units_nearest_parent`,
  `c_tr_al_hapA_sd`, `child_SD_expanded`, `c_alt_hap_frac`. The read-quality family (C) is absent because those columns are
  empty on real rows until the cohort re-review. Caveats as recorded in P14: the classifier share undercounts phase (the
  transmitted-haplotype test is rf_safe:false, rule layer), and synthetic-label shares can carry leakage — the WES arm checks.
- **P18 rf+phase v2 — rf_q with the provisional per-class τ_q (2026-09-13, integrate 54285849, concordance 54285850,
  seed-0 fold models):** per-fold pass rates are now uniform (SNV/indel 0.26–0.30 % at τ_q 0.997; TR 0.08–0.10 % and SV
  0.02–0.09 % at 0.999) — the fold-0 inflation is gone. SNV/indel: concordant 749, original-only 586 (486 below τ, 96
  demoted, 4 mosaic), module-only 2,658; per proband 38 → **81**; paternal fraction of module YES 0.697; YES = 2,585 via
  rf_q ≥ τ_q + 822 via rescue. SV: 4 / 8,125 / 10, per proband 0 (14 calls cohort-wide). TR: 26 / 192 / 2,044, 6 → 58.
  Per-child SNV/indel YES still has a high tail (19–353) to profile once the 5-seed rf_q exists.
- **GIAB Ashkenazi trio on the filer (2026-09-12, array 54274150, 28–38 min per sample, md5 OK):** unaligned PacBio
  HiFi Revio reads (2023-10-31 release; HG002 48×, HG003 46×, HG004 36×; 78 + 76 + 57 GB) under
  `/expanse/projects/sebat1/jsebat/giab/AshkenazimTrio_PacBio_HiFi-Revio_20231031/`. Coriell LCL DNA — state the
  caveat wherever used. Validation only (P16). Truths to pair with it: HG002 **Q100** assembly-based benchmark (the
  only truth that reaches inside the lab mask; haplotype-resolved), v4.2.1 mapping-based genotypes for all three
  (trio-level inherited/false-positive labels), GIAB stratifications for reporting. Running the lab WDL on the
  three genomes is a separate, larger job to be shown first.
- **M2 cohort review (2026-09-12, array 54274166, 33/33 COMPLETED, 22–79 min per family, max MaxRSS 765 MB): first
  six-haplotype evidence for every raw candidate of all 35 children (numbers below unchanged by the last quad).** Per child (median): SNV `germline_DNM_phased` 192,
  `germline_DNM_unphased` 1,568, `inherited_missed_in_parent` 398, `phase_conflict_artifact` 1,853, `inconclusive`
  8,729; INDEL phased 783, inherited-missed 2,316, conflict 5,706; SV phased 7, conflict 273 of ~627; TR phased 938,
  conflict 17,333 (both-child-haplotype fraction 0.13, was 0.56). **Sanity check (R4 — validates the class, not
  calls): paternal fraction among `germline_DNM_phased` SNVs = 0.571 (95 % CI 0.559–0.583) over all quality,
  0.668 (0.650–0.685) at child GQ ≥ 20;** INDEL 0.511, TR 0.511, SV 0.598 (n = 251). Interpretation (moderate–high
  confidence): the phased-germline SNV class is a mixture of real DNMs (~0.77 paternal) and artefacts (0.50) whose
  real share rises with quality; a two-component estimate gives ~26 % real overall and ~62 % at GQ ≥ 20 → ~50 real
  phased SNV DNMs per child by either route, against ~50 expected for the phased half of ~60–70. INDEL and TR
  phased-germline sets are artefact-dominated (homopolymer errors, stutter), as predicted for the raw set. The
  rule layer therefore enriches but does not purify — which is the P11/P15 design: the classifier carries the
  quality signal, the phase layer the transmission test. Table: `evidence/cohort_phase_class_counts.tsv`.
- **Week 1 (2026-09-12): scaffold, VCF-level minitrio, `orient.py` — done, 12 tests passing.** Package `src/phase_dnm`
  (pure-Python trio VCF reader, `phase-dnm orient` CLI), `tests/make_minitrio_vcf.py` (simulated trio with known
  haplotypes, one crossover per parent per chromosome, HiPhase-shaped per-sample phased VCFs, switch errors,
  planted DNMs, truth tables), `config/thresholds.yaml`, `containers/phase_dnm.def`. Orientation splits blocks at
  located phase-switch errors (P2 addendum). Reader throughput ~50 k sites/s → ~3 min per real family, single core.
  Not yet done from the week-1 list: the read-level minitrio (needs the container; moves to week 2).
- **Two-quad real run (2026-09-12, jobs 54269013/4, `ind-shared` 1 core 4 GB): COMPLETED.** Measured per child:
  **146–152 s**, MaxRSS 240–300 MB (the 4 GB request is 13× too generous → 1 GB from now on); 9.9–16.4 k HiPhase
  blocks, 278–324 split at located switches (313–367 switches per genome ≈ 1 per 8 k phased hets, same order as
  HiPhase's published 1 per 3.3 k variants), 1.9–2.2 M informative sites, Mendelian-inconsistent per informative
  site 0.0006–0.0011 (gate 0.01: pass), ambiguous block-bp fraction 0.094–0.108 (gate 0.10: two of four marginally
  over — see the reason breakdown below before reading this as a problem). Female chrX is 0.20–0.21 ambiguous.
  The first submission (54268976/7) failed on the sbatch spool-path bug noted under Risks.
- **Breakdown of the ambiguous bp (all four children alike):** `LOW_SITES` 8.4–9.9 % (segments with 0–4 votes
  dominate: 2.8–5.6 k of them; 10–19 votes: 0.7–1.0 k), `NO_INFORMATIVE_SITES` 0.6–0.9 %, **`MIXED_VOTES` 0.15–0.19 %**
  (15–19 segments per genome, vote fraction 0.52–0.95, largest 0.45–1.5 Mb) — the reason that would signal a
  pedigree or phasing problem is essentially absent. Female chrX: ~1,000–1,200 segments, 20 % ambiguous, all
  `LOW_SITES`/`NO_INFORMATIVE` with a median of 4–5 informative sites per segment (autosomes 9–13): X blocks are
  short and het-poor, not mis-phased. What-if: orienting 10–19-vote *unanimous* segments recovers 2.6–2.7 % of bp
  and leaves 45–53 non-unanimous small segments per genome ambiguous → adopted as a second tier
  (`small_min_sites: 10`, `small_min_frac: 1.0`; P2). The gate is now stated per reason
  (`max_frac_bp_mixed_votes: 0.01`).
- **Cohort run (2026-09-12, array 54269186, 33 families / 35 complete-trio children, 1 core 1 GB, 20 concurrent):
  33/33 tasks COMPLETED**, 150–312 s per family (quads longest), max MaxRSS 354 MB. Under the two-tier rule:
  block bp oriented median 92.0 % (90.0–97.0), ambiguous 8.0 % (3.0–10.0) = `LOW_SITES` 7.0 % + `NO_INFORMATIVE`
  0.8 % + `MIXED_VOTES` median 0.07 % (0–0.92 %); Mendelian-inconsistent per informative site median 0.0009
  (0.0006–0.0020); 300–470 switches located per genome; 19 male / 16 female children. **All 35 children pass every
  gate** (`cohort_orientation_qc.tsv` on the filer). One male child sits at `MIXED_VOTES` 0.92 % — under the
  1 % gate but 5× the next child. **Separating check:** 0.8 of the 0.92 % is a single 18.2 Mb segment,
  chr1:125.09–143.31 Mb = the 1q12 pericentromeric heterochromatin, holding 384 phased hets and 54 informative
  votes (50:4); the next two are chr9:40.7–41.6 Mb and chr21:10.3–10.8 Mb, also pericentromeric segdup
  regions. The child's switches (405 vs median 379), Mendelian ratio (0.00135 vs cohort max 0.0020) and low-GQ
  share (0.25 vs 0.18) are unremarkable → mismapping-prone regions, not a sample or pedigree problem.
  Consequence: bp-weighted ambiguity is dominated by giant sparse blocks (10 segments ≥ 5 Mb with < 5 informative
  sites/Mb across the cohort, 181.6 Mb, 6 ambiguous), so the gated quantity is now **het-weighted**
  (`frac_het_ambiguous`, `frac_het_mixed_votes`), with bp-weighted values reported alongside.
  **Het-weighted cohort result: 96.9 % of phased hets carry a parent-of-origin label (ambiguous median 3.1 %,
  range 1.9–4.4 %; `MIXED_VOTES` median 0.15 %, max 0.39 %); all 35 children pass.**
- **Week 2 started (2026-09-12): `transmission.py` (M1b) done at the VCF level** — per-parent transmitted/untransmitted
  segments and change-point candidates, `phase-dnm transmission`, 4 tests; 12-seed sweep: 0/858 resolved segment ends
  wrong, 49/49 change points explained by a planted crossover or parental switch, 32/48 crossovers recovered
  (rest between blocks or < 10 votes on a side), 95.9 % parental hets resolved. Crossover vs. parental switch is
  **not** decidable from the VCF (P3) → read-level step M1b2 next, together with `hapdepth`.
- **Cohort run of orient + transmission (2026-09-12, array 54269951, 33/33 COMPLETED, 310–648 s per family, max
  MaxRSS 663 MB).** Transmission per meiosis (35 paternal, 35 maternal): parental hets resolved to
  transmitted/untransmitted **97.6 % (F) / 97.5 % (M)** het-weighted (92.6 / 92.8 % bp-weighted; range
  96.2–98.8 %); Mendelian-inconsistent per informative site 0.0002 (0.0001–0.0003); **change-point candidates
  median 250 per paternal meiosis (193–279) and 273 per maternal (231–315), 18,315 in total** — crossovers plus
  parental switch errors, to be separated by `xo-reads` (M1b2). Resolution intervals: median 17 kb, p90 38 kb
  (1,663 < 1 kb; 214 ≥ 100 kb). Measurement only: the maternal excess of ~23 candidates per meiosis has the sign
  and order of the known maternal excess of crossovers (~1.6×; unverified figure) but could equally be more
  maternal switch errors — the read-level step decides. Table: `cohort_transmission_qc.tsv` on the filer.
- **M1b2 `xo-reads` on the first quad (2026-09-12, job 54270761, 542–557 s per child, 115 MB).** First version
  (spanning-read *count*) passed 137–179 "crossovers" per meiosis — with 15 kb reads ≥ 3 tagged reads span almost
  any gap — and was replaced by **allele concordance** across consecutive phased SNV hets. Result per meiosis:
  **CROSSOVER 32 / 29 (paternal), 59 / 41 (maternal)** against the known ~26 / ~43 (deCODE; ratio ~1.6 reproduced);
  SWITCH_ERROR 5–13; AMBIGUOUS 223–234 (`no_snv_gap` 281 → widened to the nearest flanking SNV hets; `low_reads`
  376 = gaps longer than a read, het deserts; `mixed_disc` 254 = reads disagree with each other, left ambiguous by
  design). Flag: chr9 holds 19 of 161 crossovers across four meioses (~3× its genetic-map share) — pericentromeric
  mismapping; per-chromosome crossover counts vs. genetic-map length become a cohort QC table. Child-side switch
  errors account for 5.7 % of candidates (27× chance) and are flagged, not resolved, by any parent-read test.
- **M1b2 cohort run #1 (2026-09-12, array 54271727, 33/33 COMPLETED, 62–410 s per family, 182 MB).** Per meiosis:
  CROSSOVER median **41 paternal (26–68) / 54 maternal (40–85)**, ratio 1.32; SWITCH_ERROR 9 / 7; ~78 % of candidates
  AMBIGUOUS (`low_reads` 9,385 — a gap longer than a read; `mixed_disc` 4,794; `no_snv_gap` 72). Against ~26 / ~43
  (deCODE; unverified here) and a ratio ~1.6, **paternal is over-called by ~15 and maternal by ~11 — and detection is a
  lower bound, so the excess is false positives.** Separating test: CROSSOVER calls contain a located child switch
  10.9 % of the time vs 4.4 % for AMBIGUOUS (2.5×) — child-side switches explain ~5 per meiosis (41 → 36, 54 → 50),
  not all. Remaining suspect: untested sub-intervals — the concordance span began at the first SNV het *inside*
  the interval, leaving the stretch from an indel-het endpoint untested; the span now always runs from the nearest
  phased SNV het at/before `left` to the nearest at/after `right`, and every row records
  `child_switch_in_interval`. chr19 is 2–3× its physical-length share (its genetic length per Mb is ~1.8× the
  genome mean, so likely genuine); chr9 is no longer flagged. Table: `cohort_crossover_qc.tsv`. Rerun follows.
- **hapdepth cohort run (2026-09-12, array 54271367, 105/105 COMPLETED, 139–895 s per sample, median 7 min,
  max MaxRSS 117 MB).** Autosome-averaged depth per genome: total median **23.1×** (8.3–54.0; offspring 22.35,
  parents 23.65 — matches the report's 22.4/24.4), **hap1 10.0× / hap2 10.0×** (balance 0.99–1.00), untagged 2.7×,
  MAPQ < 20 0.64×; **88 % of depth is haplotagged** (82–96 %). Built-in sex check: chrX/autosome depth ratio median
  0.54 in males (n=54), 1.08 in females (n=51), no exceptions. **Three genomes at 8.3–10.0× (one parent, two
  offspring, ~4 reads per haplotype): the six-haplotype test at k=5 is mostly unreachable there (P7, R11)** — per-sample
  depth becomes a column in every M2 table. The first submission failed 105/105 on the manifest's stale BAM paths
  (see Risks). Table: `hapdepth/cohort_hapdepth_qc.tsv`.
- **M1b2 cohort run #2 (2026-09-12, array 54272721, 33/33, median 84 s per family, 180 MB) — M1b closed at cohort
  level.** With the tested span covering the whole change interval: CROSSOVER median **34 paternal (19–62) / 47
  maternal (36–73)**; excluding intervals that contain a located child switch: **30 / 43**, ratio **1.43**;
  SWITCH_ERROR 9 / 7; ~80 % of candidates AMBIGUOUS (gaps longer than a read 9,798; reads disagreeing 4,696; < 2 SNV
  hets 217). Against ~26 / ~43 (deCODE, unverified here; ratio ~1.6): maternal on target, paternal ~15 % high for a
  method that also *misses* crossovers between parent blocks — the residual is unlocated child switches and
  mismapping. chr17 (2.1×) and chr19 (3.3×) exceed their physical-length share; both have well above-average
  recombination per Mb, so the physical proxy under-predicts them and this is not read as artefact. Reported in the
  paper as: crossover candidates per meiosis with the child-switch flag, the ambiguous fraction, and the
  per-chromosome table (`cohort_crossover_qc.tsv`), never as a recombination map.
- **Read-level steps written (2026-09-12):** `hapdepth` (per-haplotype depth in 1 kb bins over each primary read's
  reference span; pysam, lazily imported) and `xo-reads` (M1b2: weakest-link count of the parent's haplotagged
  reads spanning consecutive phased hets across each change interval → CROSSOVER / SWITCH_ERROR / AMBIGUOUS).
  `hapdepth` smoke on a 5 Mb region of one BAM: 4,632 primary reads in 1.1 s, 56 MB RSS, per-haplotype means
  7.9 / 7.5 + 0.65 untagged + 0.24 low-MAPQ → ~17 min and < 1 GB per ~22× genome; run as a 105-task array
  (`workflow/hapdepth_array.sb`, 1 core, 2 GB, 1 h). pysam 0.24.1 is in the lab's `dnmt3a-py` env, and the WDL's
  hiphase / pbmm2 / sawfish / trgt / pb_wdl_base images are in the miniwdl singularity cache, so no container
  build is needed for week 2.
- **ICR spot check (independent truth from `09_methylation`, nearest-informative-SNV parent of origin at 18
  imprinting control regions, 4 children, 41 ICR rows): 30 concordant, 0 discordant**; 3 ICRs had no local phase
  in the methylation table, 6 fell in gaps between blocks, 2 in `LOW_SITES` segments. Lustre was reachable again
  (backgrounded probe, 20 s kill-timeout).

Three weeks, then a gate. Each week ends with something that runs on `tests/data/minitrio`
**and** on at least two real families (one of each quad, so sib structure is exercised from day
one; the two duos are excluded, P19). Deliverable: `phase/<FAMILY>/` for every complete-trio
family, `phase_qc.json` per trio, a cohort QC table, and the M2 evidence extractor able to
consume it. P# refer to `DESIGN.md`.

What changed after inspecting the filer (2026-09-12): HiPhase already emits the read→haplotag
table, block table and per-chromosome stats for every sample, sawfish emits supporting read
names per SV, and TRGT emits per-read allele lengths. M1 therefore builds **orientation and
transmission** on top of those and does not re-tag BAMs; the only BAM pass is `hapdepth`.

Operating rules that bind every step (from `docs/EXPANSE.md`): nothing heavier than a header
read on the login node; every job shown with partition/cores/memory/walltime **and an ETA**
before submission (`lab_submit`), inspected after (`lab_status`), with a scheduled check;
first run of anything new is small, on `ind-shared` (the `debug` partition is not available);
scripts live in `$HOME` or the filer, never `/tmp`; walltime from measured task times.

Prerequisites: JS's GitHub clean-up (another session) finished — **until then no commit, push or
`git pull` anywhere**; then `git pull` of the on-cluster checkout (at the 09/10 commit; GitHub has
`11_`) and `13_phase_dnm/` copied in from the local staging directory; one small SU slice.

GIAB (Q6: not on Expanse) is pulled during week 3 as a download job so it is ready for week 4:
HG002/HG003/HG004 HiFi aligned BAMs at a coverage matching the cohort (~22–24×; downsample if
the release is deeper), plus the GIAB v4.2.1 benchmark VCFs and beds. Resources, size, source
URLs and ETA are shown before submission (EXPANSE.md rule 3); the data land on the filer, not
Lustre; nothing GIAB touches any fold, calibration or threshold (P16).

---

## Week 1 — scaffold, minitrio, block orientation (M1a)

**Build**
- `13_phase_dnm/` skeleton, `pyproject.toml`, `phase-dnm` CLI stub, Apptainer definition
  (pysam, cyvcf2, polars, numpy, pytest, samtools, bcftools). `config/phase_dnm.env` from the
  example with the real prefixes.
- `tests/make_minitrio.py` → `tests/data/minitrio/`: ~300 kb of GRCh38 with planted het sites
  at realistic density; two parental diploids; a child with **one planted crossover per parent**;
  HiFi-like reads (~15 kb, 1 % indel-dominated error); pbmm2 + HiPhase inside the container so
  the fixture has the *same* output shapes as the WDL (haplotags TSV, blocks TSV, phased VCFs).
  Planted per class (SNV/indel, SV DEL+INS, TR expansion): one germline DNM, one child
  postzygotic mosaic (f_A 0.5), one parental mosaic transmitted (alt_T 0.15), one
  inherited-missed-in-parent, one phase-conflict artefact. Truth in `truth.tsv`. A TRGT run on
  the minitrio gives the spanning-read BAM with `AL:i`; a sawfish run gives
  `supporting_reads.json` — both are what the M2 adapters consume, so they are part of the
  fixture.
- `phasing/orient.py`: family joint phased VCF + manifest → `orientation.tsv` (P2). Parameters
  in `thresholds.yaml`. Sex chromosomes: chrX in a male child is maternal by construction, chrY
  paternal; PAR handled with the PAR BED synthdnm already ships.

**Test / acceptance**
- minitrio: every child block oriented correctly; planted crossovers do not flip the child label.
- two real families: `AMBIGUOUS` fraction (blocks and bp); `n_mendel_err` distribution; spot check
  against the 18 ICR parent-of-origin calls from `09_methylation` (`icr_parent_of_origin_local`).
  Runtime measured — orientation is VCF-only and should be minutes per family on 1–2 cores.

## Week 2 — transmission map, crossovers, haplotype depth (M1b)

**Build**
- `phasing/transmission.py` → `transmission.bed`, `crossovers.tsv` (P3).
- `phasing/hapdepth.py`: one pass per haplotagged BAM → `<SAMPLE>.hapdepth.tsv.gz` (1 kb bins:
  `dp_hap1 dp_hap2 dp_untagged`). This is the only BAM-wide pass in M1; ~60 GB per BAM, so it
  is a job array (rule 4), one task per sample, resources set from a two-sample timing run.
- `phasing/labels.py`: per-sample block-level `XP`/`XT` tables (P4).
- `workflow/m1_family.sb <FAMILY>` (orient → transmission → labels) and `workflow/hapdepth_array.sb`
  (array over samples), both via `lab_submit`; a collector modelled on `01_qc/wait_job.sh`.

**Test / acceptance**
- minitrio: both planted crossovers recovered within planted resolution; no false crossovers;
  segments cover the whole fixture.
- two real families: crossovers per parent within the QC gate; fraction of autosome with the
  transmitted haplotype resolved in both parents; crossover resolution distribution;
  within-block switches failing `xo_min_sites` counted as phase-switch candidates (R6).
  Timing and SU per family logged → ETA for all 33 families is a measurement.

## Week 3 — QC, cohort run, WhatsHap comparison, hand-off (M1c/M1d)

**Build**
- `phasing/qc.py` → `phase_qc.json` (README §3.1) from `phase_stats`, `orientation.tsv`,
  `transmission.bed`, `hapdepth`; cohort table `phase/cohort_phase_qc.tsv` with pass/fail
  against `thresholds.yaml` gates.
- `phase-dnm haplotag --export-bam` (optional IGV export with `XP`/`XT`).
- WhatsHap `--ped` on chr20–22 for 2–3 trios; agreement with oriented HiPhase blocks, listed by
  block → DESIGN.md evidence for P2.
- Cohort run over all 33 complete-trio families (array), with stated ETA and scheduled check.

**Test / acceptance — module gate**
- Every complete-trio family has `phase_qc.json`; no family outside the crossover gate without a
  documented reason; median fraction of autosome with all six haplotypes at k=5 reported (this
  bounds what M2 can ever claim and goes in the paper).
- `hapmatrix.py` builds the six-haplotype matrix for the planted minitrio SNV, SV **and** TR from
  (haplotags TSV, block labels, hapdepth, supporting-reads JSON, spanning-read BAM) — the
  hand-off test. M2 proper starts week 4, on the minitrio and on GIAB (Q6) as the first external
  fixture, before any SPARK candidate is scored.

---

## Deliberately not in Module 1
- No re-phasing (P2) — revisited only if the week-3 WhatsHap comparison says orientation is
  materially worse on the same blocks, with the numbers.
- No SV/TR phasing beyond HiPhase's (they are phased; `PS` present in both VCFs).
- No duos (P19, Q14). No GIAB until week 4 (P16).

## Risks specific to M1
- **sbatch executes a spool copy of the script**, so `BASH_SOURCE`-relative paths resolve to
  `/cm/local/apps/slurm/var/spool/…`. The first two orientation jobs (2026-09-12) died in 1 s on
  `source config/phase_dnm.env` for exactly this reason. Every job script in this module resolves the
  module directory from `PHASE_DNM_HOME`, then `SLURM_SUBMIT_DIR`, and only then `BASH_SOURCE`.
- The Lustre analysis directory was unreachable on 2026-09-11; everything M1 needs is on the
  filer, but the baseline de novo tables for P18 are on Lustre — copy them to the filer first.
- `phase_haplotags` covers phased blocks only; reads in unphased regions have no row and are
  `untagged` by definition (P7). Expect that fraction to be highest exactly where candidates are
  hardest (mask, low mappability) — it is a reported QC number, not a surprise.
- Login-node `bcftools`/`samtools` from the micromamba env need `OPENBLAS_NUM_THREADS=1`
  (they die on allocation otherwise); baked into `phase_dnm.env`.
