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
