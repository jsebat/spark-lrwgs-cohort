# 05_denovo — phase-aware de novo mutation calling for long-read trios

> **Parent of origin is validated per call (2026-09-16).** An independent read-backed assignment agrees with the
> pipeline on 1,796 of 1,802 tier-1 de novo SNVs (99.7%); the orientation table is right on 1,680 of 1,682 blocks.
> The apparent excess maternal-age effect is leverage from one quad with a 42-year-old mother, not misassignment.
> Validation step: `workflow/m3_poo_readcheck.sb` + `tools/poo_eval.py`. Details: [DESIGN.md P35](DESIGN.md).



Module of `spark-lrwgs-cohort`. Turns trio phase from a by-product into primary evidence for
de novo mutation (DNM) calling, uniformly for SNV/indel, SV and TR. Sub-modules, built in
this order: **M2 six-haplotype review (all classes) → feature registry / likelihood / spike-in →
M3 integration → M4 SynthDNM retraining (nested CV).**

**M1 — phasing, orientation and the transmission map — is no longer part of this module.** It is
`../02_phasing` (package `trio_phase`, CLI `trio-phase`), which has no dependency on de novo calling
and is useful on its own. This module CONSUMES its tables **by path** under `$PHASE_DIR` (§2.1) and
imports nothing from it. Run `02_phasing` first; `workflow/Snakefile` has
`$PHASE_DIR/<family>/m1_xoreads.status` as an input of `m2_review`, so a family whose phasing has not
been run fails the DAG at build time.

Conventions inherited from the repo: scripts take family / sample identifiers as **arguments**
and never embed them; configuration in `config/cohort.env` (this module adds
`config/phase_dnm.env`); no identifiers, names or data in the tree; `scripts/phi_scan.sh` before
every push. Design decisions with rationale are in `DESIGN.md` (numbered P#), so the Methods
section can cite them the way `docs/METHODS.md` cites `D#`.

Status: **proposal — nothing here runs yet.** Open questions are at the top of `DESIGN.md`;
assumptions made in their absence are marked `ASSUMED`.

---

## 1. Layout

```
05_denovo/
├── README.md                  this file: layout, interfaces, file formats
├── DESIGN.md                  open questions, decisions P1…, SynthDNM one-pager, CV design,
│                              overclaim register
├── PLAN_MODULE1.md            running build/results log (M1 origins, then M2-M4)
├── config/
│   ├── phase_dnm.env.example  module config (sourced after config/cohort.env)
│   ├── features.yaml          FEATURE REGISTRY — every feature: classes, source, rf_safe
│   └── thresholds.yaml        k (min reads/haplotype), margins, class thresholds; versioned
├── src/phase_dnm/             installable package, one CLI: `phase-dnm <cmd>`
│   ├── records.py             CandidateRecord — the shared normalised candidate
│   ├── hapmatrix.py           six-haplotype evidence matrix — the shared extractor
│   ├── adapters/              ONE class-specific hook each: which reads support ALT?
│   │   ├── snv_indel.py       base / indel at position from aligned pairs
│   │   ├── sv.py              junction (SA / large CIGAR D,I) + per-hap depth in interval
│   │   └── tr.py              per-read allele length in the TRGT locus (spanning reads)
│   ├── classify/
│   │   ├── rules.py           phase_class + transparent rule score
│   │   └── likelihood.py      per-haplotype read-count likelihood, posterior over hypotheses
│   ├── features/
│   │   ├── registry.py        loads features.yaml; REFUSES rf_safe:false columns in RF matrices
│   │   └── extract.py         evidence matrix + caller/context/read features → one vector
│   ├── sim/spike.py           spike-in harness: edit haplotagged reads, all three classes
│   ├── eval/                  ONE harness for every arm (DESIGN P22): heuristics.py reproduces the original
│   │                          pipeline's per-class filters verbatim (pass + sweep_score); harness.py scores
│   │                          H, H+P, RF, RF+P on identical outer folds and external truth
│   ├── train/                 M4: folds.py (swap-closed), swap.py, nested_cv.py, calibrate.py
│   └── io/                    tables.py (core columns), vcfinfo.py (INFO equivalents)
├── workflow/                  Snakefile + profiles/expanse (ASSUMED, Q9) + *.sb wrappers
├── tests/
│   ├── data/minitrio/         ~300 kb synthetic trio: 3 BAMs, phased VCF, sawfish- and
│   │                          TRGT-style records; one planted DNM per class + one of each
│   │                          failure mode per class
│   └── test_*.py              the three adapters are exercised by the SAME test suite
└── containers/                phase_dnm.def (Apptainer): pysam, cyvcf2, polars, xgboost, sklearn
```

**Python ≥ 3.10 is required** (`pyproject.toml`). On Expanse the login-node `python3` is 3.6 and is never
used for this module; run inside `containers/phase_dnm.def` or with the lab micromamba environment named in
`config/phase_dnm.env` (`PHASE_DNM_PYTHON`). Module 1 has no compiled dependencies (pure-Python VCF reader),
so the same code runs unchanged on a laptop, the login node's env and inside a job.

CLI surface — one command per sub-module. Every command takes `--family --proband --father
--mother` (or `--ped`) and resolves paths from the env, never from embedded identifiers. The phasing
commands (`orient`, `transmission`, `xo-reads`, `hapdepth`, `phase-qc`) are now `trio-phase <cmd>` in
`../02_phasing`:

```
phase-dnm candidates    phase-dnm review         phase-dnm reclassify  phase-dnm likelihood   phase-dnm features
phase-dnm spike         phase-dnm integrate      phase-dnm concordance   phase-dnm train       phase-dnm classify
```

## 2. Interfaces

### 2.1 Phasing inputs — produced by `../02_phasing`

Not this module. `02_phasing/README.md` documents the inputs (per-sample HiPhase-phased VCFs,
haplotagged BAMs, the cohort manifest) and the full schema of every table it writes to
`$PHASE_DIR/<FAMILY>/`. What this module reads from there, by path:

| file | read by |
|---|---|
| `<CHILD>.orientation.tsv` | `review`, `reclassify`, `features`, `spike`, `annotate` — `LabelTables` (child HP → parent of origin) and `Geometry` (block length, distance to a block edge) |
| `<CHILD>.transmission.tsv` | `review`, `spike` — `LabelTables` (parent HP → transmitted / untransmitted) |
| `<CHILD>.changepoints.resolved.tsv` | `review`, `reclassify`, `features`, `spike`, `annotate` — distance to a crossover, `NEAR_CHILD_SWITCH` |
| `$PHASE_DIR/hapdepth/<SAMPLE>.hapdepth.tsv.gz` | the SV adapter's per-haplotype interval depth |

`config/phase_dnm.env` must name the same `PHASE_DIR` as `02_phasing/config/phasing.env`; on the
cluster the two env files may simply be the same file. There is no Python import in either direction.

### 2.2 Module 2 — six-haplotype review

`review` (the only BAM pass) writes the immutable `<child>.<class>.evidence.review.tsv`; `reclassify` re-runs the rule layer from its count columns under the current `thresholds.yaml` version and writes the working `<child>.<class>.evidence.tsv` that `likelihood` (adds `lik_post_*`, `phase_score`) and `features` consume (`workflow/m2_reclassify_family.sb` chains the three). Rule changes therefore never re-touch the BAMs.

**Input: `candidates.tsv`** — the shared candidate record. Produced by
`phase-dnm candidates --class {snv_indel,sv,tr}` from the joint callsets (GLnexus BCF, sawfish
VCF, TRGT VCF) **and** from every existing candidate list (`*.denovo.hiconf.tsv`,
`denovo_sv_priority` output, `tr_outliers` output), so that high- and low-quality candidates
enter unfiltered (P5). Columns:

```
family_id sample_id variant_id chrom start end ref alt variant_class{SNV,INDEL,SV,TR}
caller caller_gt caller_gq caller_dp caller_qual caller_filter
source_tier{HIGH,LOW,UNFILTERED} source_list mask_overlap{0,1} class_payload
```

`class_payload` is JSON and is the only class-specific field at this stage:
- SNV/INDEL: `{}`
- SV: `{"svtype","svlen","cipos","ciend","end2","mateid","caller_id"}`
- TR: `{"trid","motifs","struc","child_AL","child_SD","child_ALLR","father_AL","mother_AL","expanded_allele_idx"}`

**Output: `evidence/<FAMILY>.<CLASS>.evidence.parquet`** — one row per candidate. The
**six-haplotype matrix** is flattened with row prefixes `M1_ M2_ F1_ F2_ CM_ CP_`; the
transmission labelling is carried in *separate* columns so the matrix itself stays symmetric in
the parents (required for `rf_safe`, P11):

| per-haplotype field | meaning |
|---|---|
| `dp` | spanning reads assigned to this haplotype (SV/TR: reads spanning the event with anchors both sides) |
| `alt`, `ref`, `amb` | reads supporting ALT / REF / neither — the adapter decides |
| `mapq_mean`, `mapq0_frac`, `nm_alt_mean`, `nm_ref_mean`, `clip_alt_frac` | read-level quality, alt vs ref |
| `al_mean`, `al_sd`, `al_n` (TR only) | per-read allele length on this haplotype |
| sample-level `untagged_dp`, `untagged_alt` | reads with no HP tag — still counted, as the existing DNM logic does |

Plus: `F_transmitted_hap{1,2,NA}`, `M_transmitted_hap{1,2,NA}`, `dist_crossover_bp`,
`dist_block_edge_bp`, `n_informative_1kb`, `local_rephase_agreement` (the `pofo_local.py`
logic — do ALT reads agree with each other on flanking hets, independent of HP tags), the
derived phase features of `features.yaml` section D, `phase_class`, `rule_score`,
`lik_post_*`, `parent_of_origin`, `poo_confidence`, `poo_reason`, `flags`.

**Output: `evidence/<FAMILY>.<CLASS>.reads.jsonl.gz`** — one JSON object per candidate, per
read `{rid, sample_role, HP, PS, XP|XT, support{ALT,REF,AMB}, mapq, nm, clip, (TR) al}`.
`rid` is a salted blake2b hash of the read name so nothing sample-derived leaves the cluster.

### 2.3 Module 3 — integration / final tables

One **unfiltered** table per class group and child, `final/<FAMILY>.<child>.<snv_indel|sv|tr>.dnm.tsv` (+ `.vcf` sites file; Parquet when pyarrow is present): common core
in this order, then class-specific columns, then the full registered feature vector:

```
family_id sample_id chrom start end ref alt variant_class caller caller_gt caller_qual
rf_prob dnm_call{YES,NO} call_mode{rf+phase,phase_only} decision_reason mosaic_flag
parent_of_origin{paternal,maternal,undetermined} poo_reason poo_confidence
phase_class hap_obs_k3 hap_obs_k5 child_alt_hap_frac child_alt_other_hap
transmitted_parent_alt_reads untransmitted_parent_alt_reads phase_score flags
[SV: svtype svlen bp_precision]
[TR: trid motif child_AL_pat child_AL_mat father_AL_T father_AL_U mother_AL_T mother_AL_U]
<features.yaml columns…>
```

`rf_prob` comes from M4 (`--rf-probs`); until then `dnm_call` is the documented provisional phase-only decision (`call_mode`, DESIGN P15). Everything else is from M2 and is **filled for every row,
including NO calls.** `phase-dnm integrate` builds the table (`workflow/m3_integrate_family.sb`); `phase-dnm concordance` produces the P18 tables against `baselines/` (`workflow/m3_concordance.sb`). VCF INFO equivalents (`io/vcfinfo.py`), one header block per class VCF:

```
PDNM_PROB  PDNM_CALL  PDNM_POO  PDNM_POOR(reason)  PDNM_POOC  PDNM_CLASS  PDNM_HAPOBS
PDNM_CHF(child alt-hap frac)  PDNM_CAO(child alt other hap)  PDNM_TALT  PDNM_UALT
PDNM_SCORE  PDNM_FLAGS
```

### 2.4 Module 4 — training artefacts

`train/<RUN>/folds.json` (swap-closed family groups + seed), `pairings.json` (whose offspring
were paired with whose parents, per fold), `training_manifest.json` (families, classes,
feature set, rf_safe check hash, software versions), `cv_report.json` (outer-fold PR-AUC,
ROC-AUC, Brier, reliability bins; per class, per coverage stratum, per ablation),
`model.<class>.json`, `calibration.<class>.json`, `thresholds.<class>.json`.

## 3. QC and acceptance

### 3.1 `phase_qc.json`
Written by `../02_phasing` (`trio-phase phase-qc`); see `02_phasing/README.md`. Its gates are the
entry condition for everything below: a child that fails them should not be reviewed.

### 3.2 Definition of done (per sub-module)
Nothing is done for one class until the other two have a working path **and** a passing test on
`tests/data/minitrio`. The minitrio carries, **per class**, one planted germline DNM, one
postzygotic mosaic, one parental mosaic transmitted, one inherited-missed-in-parent and one
phase-conflict artefact.
