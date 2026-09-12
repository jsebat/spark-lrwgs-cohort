# Module 1 build plan — phasing orientation, transmission map, QC

## Status
- **Week 1 (2026-09-12): scaffold, VCF-level minitrio, `orient.py` — done, 12 tests passing.** Package `src/phase_dnm`
  (pure-Python trio VCF reader, `phase-dnm orient` CLI), `tests/make_minitrio_vcf.py` (simulated trio with known
  haplotypes, one crossover per parent per chromosome, HiPhase-shaped per-sample phased VCFs, switch errors,
  planted DNMs, truth tables), `config/thresholds.yaml`, `containers/phase_dnm.def`. Orientation splits blocks at
  located phase-switch errors (P2 addendum). Reader throughput ~50 k sites/s → ~3 min per real family, single core.
  Not yet done from the week-1 list: the read-level minitrio (needs the container; moves to week 2) and the
  two-real-family run (needs an sbatch line shown first).

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
- The Lustre analysis directory was unreachable on 2026-09-11; everything M1 needs is on the
  filer, but the baseline de novo tables for P18 are on Lustre — copy them to the filer first.
- `phase_haplotags` covers phased blocks only; reads in unphased regions have no row and are
  `untagged` by definition (P7). Expect that fraction to be highest exactly where candidates are
  hardest (mask, low mappability) — it is a reported QC number, not a surprise.
- Login-node `bcftools`/`samtools` from the micromamba env need `OPENBLAS_NUM_THREADS=1`
  (they die on allocation otherwise); baked into `phase_dnm.env`.
