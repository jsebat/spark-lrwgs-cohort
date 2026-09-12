# 13_phase_dnm — phase-aware de novo mutation calling for long-read trios

Module of `spark-lrwgs-cohort`. Turns trio phase from a by-product into primary evidence for
de novo mutation (DNM) calling, uniformly for SNV/indel, SV and TR. Four sub-modules, built in
this order: **M1 phasing + transmission map → M2 six-haplotype review (all classes) → feature
registry / likelihood / spike-in → M3 integration → M4 SynthDNM retraining (nested CV).**

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
13_phase_dnm/
├── README.md                  this file: layout, interfaces, file formats
├── DESIGN.md                  open questions, decisions P1…, SynthDNM one-pager, CV design,
│                              overclaim register
├── PLAN_MODULE1.md            week-by-week build plan for M1
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
│   ├── phasing/
│   │   ├── orient.py          M1a  orient child HiPhase blocks pat/mat by Mendelian vote
│   │   ├── transmission.py    M1b  transmitted/untransmitted map per parent + crossovers
│   │   ├── haplotag.py        M1c  write XP (child) / XT (parent) tags into BAMs
│   │   └── qc.py              M1d  phase_qc.json
│   ├── classify/
│   │   ├── rules.py           phase_class + transparent rule score
│   │   └── likelihood.py      per-haplotype read-count likelihood, posterior over hypotheses
│   ├── features/
│   │   ├── registry.py        loads features.yaml; REFUSES rf_safe:false columns in RF matrices
│   │   └── extract.py         evidence matrix + caller/context/read features → one vector
│   ├── sim/spike.py           spike-in harness: edit haplotagged reads, all three classes
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

CLI surface — one command per sub-module. Every command takes `--family --proband --father
--mother` (or `--ped`) and resolves paths from the env, never from embedded identifiers:

```
phase-dnm orient        phase-dnm transmission   phase-dnm haplotag    phase-dnm phase-qc
phase-dnm candidates    phase-dnm review         phase-dnm features
phase-dnm spike         phase-dnm integrate      phase-dnm train       phase-dnm classify
```

## 2. Interfaces

### 2.1 Module 1 — phasing, orientation, transmission map

**Inputs** (all already produced by HiFi-human-WGS-WDL v3.3.1 for every family; paths and
globs in `config/phase_dnm.env.example`, confirmed on the filer 2026-09-12):
- per-sample haplotagged BAMs (`HP:i`, `PS:i`, `MM/ML`, `rq` verified on reads)
- HiPhase 1.6.0 outputs per sample: `phase_haplotags` (read_name → haplotag per block, ~2.8 M
  rows), `phase_blocks`, `phase_stats` (per-chromosome NG50 etc.)
- family joint small-variant VCF (all members) and HiPhase-phased small-variant, **SV (sawfish,
  phased — `PS` present)** and TRGT VCFs
- for the M2 adapters: `sv_supporting_reads/*.json.gz` (sawfish id → sample → read names) and
  `trgt_spanning_reads/*.bam` (per-read `AL:i`, `HP`, `PS`)
- the cohort manifest (`sample_id family_id father_id mother_id sex affected role …`) — the
  pedigree source; the two mother–child duos are excluded (DESIGN P19)
- reference FASTA (GRCh38 no-alt)

**Outputs** → `$PHASE_DIR/<FAMILY>/`

| file | format | content |
|---|---|---|
| `<SAMPLE>.phased.blocks.bed` | BED | one row per HiPhase block: `chrom start end PS n_het_phased block_len` |
| `<CHILD>.orientation.tsv` | TSV | one row per **child block segment** (implemented, `phase-dnm orient`): `chrom phase_block_id segment n_segments start end switch_pos n_het_phased n_inf_pat n_inf_mat n_informative vote_frac orientation{HAP1_PAT,HAP1_MAT,AMBIGUOUS} reason{OK,SPLIT_AT_SWITCH,LOW_SITES,MIXED_VOTES,NO_INFORMATIVE_SITES} n_dissent`. A HiPhase block whose votes run one sign then the other is a phase-switch error; it is split at the located change point into segments oriented separately (P2), so a read's label is looked up by `(phase_block_id, pos)`, not by `phase_block_id` alone. |
| `<CHILD>.orientation.dissent.tsv` | TSV | positions voting against their block's orientation (switch-error / genotype-error candidates): `chrom pos phase_block_id block_orientation` |
| `<CHILD>.orientation.summary.json` | JSON | per-chromosome and total counters (`n_het`, `n_het_phased`, `n_informative`, `n_uninformative`, `n_mendel_inconsistent`, `n_low_gq`, `n_parent_missing`, `n_skipped_sex_chrom`), block counts, `frac_bp_ambiguous`, `mendel_inconsistent_per_informative`, the parameters and `thresholds_version` used |
| `<FAMILY>.transmission.bed` | BED | one row per contiguous segment per parent: `chrom start end parent{F,M} transmitted_hap{1,2} parent_PS n_informative confidence boundary_type{CROSSOVER,BLOCK_EDGE,CHROM_END}` |
| `<FAMILY>.crossovers.tsv` | TSV | `chrom left_bound right_bound parent resolution_bp n_inf_left n_inf_right in_parent_block{Y,N}` — a within-block switch (`Y`) is a crossover **or** a phase-switch error; see P3 |
| `<SAMPLE>.hapdepth.tsv.gz` | TSV | per-haplotype depth in fixed bins (default 1 kb): `chrom start end dp_hap1 dp_hap2 dp_untagged` — feeds `hap_obs` and the SV adapter without re-reading BAMs |
| `<SAMPLE>.tagged.bam` (+.bai) | BAM | HiPhase tags preserved. **Child** reads gain `XP:A:{P,M,U}` (parent of origin; U = block unresolved). **Parent** reads gain `XT:A:{T,U,N}` (transmitted / untransmitted / not resolved). Both tags are in the SAM range reserved for local use. |
| `<FAMILY>.phase_qc.json` | JSON | §3.1 |

### 2.2 Module 2 — six-haplotype review

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

One **unfiltered** table per class, `final/<FAMILY>.<CLASS>.dnm.tsv` (+ Parquet): common core
in this order, then class-specific columns, then the full registered feature vector:

```
family_id sample_id chrom start end ref alt variant_class caller caller_gt caller_qual
rf_prob dnm_call{YES,NO} parent_of_origin{paternal,maternal,undetermined} poo_reason poo_confidence
phase_class hap_obs_k3 hap_obs_k5 child_alt_hap_frac child_alt_other_hap
transmitted_parent_alt_reads untransmitted_parent_alt_reads phase_score flags
[SV: svtype svlen bp_precision]
[TR: trid motif child_AL_pat child_AL_mat father_AL_T father_AL_U mother_AL_T mother_AL_U]
<features.yaml columns…>
```

`rf_prob`, `dnm_call` come from M4; everything else from M2 and is **filled for every row,
including NO calls.** VCF INFO equivalents (`io/vcfinfo.py`), one header block per class VCF:

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

### 3.1 `phase_qc.json` (M1)
Per sample: block NG50, fraction of het SNVs phased, fraction of reads haplotagged,
per-haplotype mean depth. Per trio: Mendelian-error rate at informative sites, block
orientation ambiguity rate (fraction of child blocks `AMBIGUOUS`), switch/flip rate estimated
from trio consistency, crossovers per parent (expect roughly 25–45; far outside is a QC
failure, not a discovery), fraction of the autosomal genome with the transmitted haplotype
resolved in both parents, fraction of the genome with all six haplotypes at ≥k reads (k=3,5).

### 3.2 Definition of done (per sub-module)
Nothing is done for one class until the other two have a working path **and** a passing test on
`tests/data/minitrio`. The minitrio carries, **per class**, one planted germline DNM, one
postzygotic mosaic, one parental mosaic transmitted, one inherited-missed-in-parent and one
phase-conflict artefact.
