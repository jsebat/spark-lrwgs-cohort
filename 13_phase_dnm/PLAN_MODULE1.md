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
