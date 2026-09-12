# Module 1 build plan — phasing orientation, transmission map, QC

## Status
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
