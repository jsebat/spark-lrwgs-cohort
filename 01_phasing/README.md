# 01_phasing — trio-resolved and physical phasing for long-read families

Module of `spark-lrwgs-cohort`. Turns a family's read-based (HiPhase) phase blocks into **biological**
phase: which of the child's two haplotypes came from which parent, which of each parent's two
haplotypes was transmitted, and where along a chromosome the transmitted haplotype changes — with the
crossovers separated from the parental switch errors using the parents' own reads.

It answers questions that stand on their own — parent of origin for any inherited variant, imprinting,
recombination maps, phase-aware QC of a cohort — and it is also the first stage of `13_phase_dnm`,
which was its first consumer and where this code lived until it was extracted.

> **Numbering.** `01_qc` already exists, so `01_phasing` shares the 01 slot. The intent (JS) is that
> phasing runs early, straight after `00_upstream`; if the duplicate prefix is unwelcome, renaming this
> directory is safe — nothing outside it refers to it by path except the two README pointers and
> `13_phase_dnm/workflow/Snakefile`'s comment.

**Provenance.** This was the `M1` sub-module of `13_phase_dnm`. The extraction was a move with **no
behaviour change**: the six source files are byte-identical to the versions that produced the cohort
tables now on the filer, and the thresholds are copied verbatim, version string included. What changed
is the package name (`phase_dnm` → `trio_phase`), the CLI name (`phase-dnm` → `trio-phase`), the
environment file (`config/phase_dnm.env` → `config/phasing.env`) and the addition of `audit`.
Design decisions P2 (block orientation), P3 (change points and their read-level resolution), P7
(per-haplotype depth) and P19 (duos excluded) are recorded in `13_phase_dnm/DESIGN.md`; decisions taken
during and since the extraction are in §7 below.

---

## 1. Layout

```
01_phasing/
├── README.md                  this file
├── config/
│   ├── phasing.env.example    module config (copy to config/phasing.env; gitignored)
│   └── thresholds.yaml        every threshold, versioned; written into each summary JSON
├── src/trio_phase/            installable package, one CLI: `trio-phase <cmd>`
│   ├── cli.py                 orient | transmission | xo-reads | hapdepth | phase-qc | audit
│   ├── audit.py               declaration-vs-artefact invariants; exits non-zero on any FAIL
│   ├── io/vcf.py              pure-Python VCF reader, trio site merge, manifest/pedigree
│   └── phasing/
│       ├── orient.py          orient the child's phase blocks pat/mat by Mendelian vote
│       ├── transmission.py    transmitted/untransmitted map per parent + change points
│       ├── xo_reads.py        change point -> CROSSOVER / SWITCH_ERROR, from parent reads (pysam)
│       ├── hapdepth.py        per-haplotype depth in fixed bins (pysam)
│       └── qc.py              phase_qc.json + the cohort QC table with gates
├── workflow/                  Snakefile + the five sbatch scripts + profiles/expanse
└── tests/                     synthetic mini-trio fixture; the same suite that covered M1, plus the audit
```

**Python ≥ 3.10.** `orient`, `transmission`, `phase-qc` and `audit` are pure Python (stdlib + PyYAML) and
run anywhere. `xo-reads` and `hapdepth` need `pysam` (`pip install -e '.[reads]'`). On Expanse the
login-node `python3` is 3.6 and is never used; `PHASING_PYTHON` (or the older `PHASE_DNM_PYTHON`) names
a Python ≥ 3.10.

`io/vcf.py` is a **copy** of `13_phase_dnm/src/phase_dnm/io/vcf.py`, not a move: `candidates.py`,
`annotate.py` and `features/extract.py` over there still use it. The two copies are byte-identical
today; if one changes, the other should be looked at.

## 2. Inputs

All produced per family by HiFi-human-WGS-WDL v3.3.1; paths and globs live in `config/phasing.env`.

| input | used by | note |
|---|---|---|
| per-sample **HiPhase-phased small-variant VCF** (`*.small_variants.phased.vcf.gz`, `PS` present) | `orient`, `transmission`, `xo-reads` | the per-sample split of the family joint VCF |
| per-sample **haplotagged BAM** + index (`HP:i`, `PS:i`) | `xo-reads` (parents), `hapdepth` (every sample) | indexes live in the parallel `*_index/` directory, **not** beside the BAM |
| the **cohort manifest** TSV (`sample_id family_id father_id mother_id sex affected role …`) | everything | the pedigree source; identifiers are never embedded in code |
| reference FASTA (GRCh38 no-alt) | — | for context; nothing here reads it directly |

Only **complete trios** are processed: an offspring whose father and mother are both samples of the
manifest. Duos are skipped with a note (DESIGN P19) and have no outputs, which the audit expects.

## 3. Outputs → `$PHASE_DIR/<FAMILY>/`

Column lists below are the writers' own (`orient.ORIENTATION_COLUMNS`,
`transmission.SEGMENT_COLUMNS` / `CHANGE_COLUMNS`, `xo_reads.COLUMNS`); the audit compares them against
the headers on disk rather than restating them, so this table cannot drift silently.

### `<CHILD>.orientation.tsv`
One row per **segment of a child phase block**. A block whose votes run one sign and then the other
carries a phase-switch error; it is split at the located change point and each segment is oriented
separately (P2), so a position's label is looked up by `(phase_block_id, pos)`, never by block alone.

```
chrom  phase_block_id  segment  n_segments  start  end  switch_pos  n_het_phased
n_inf_pat  n_inf_mat  n_informative  vote_frac  orientation  reason  n_dissent
```
* `orientation` ∈ `HAP1_PAT` | `HAP1_MAT` | `AMBIGUOUS`
* `reason` ∈ `OK` | `OK_SMALL_UNANIMOUS` | `SPLIT_AT_SWITCH` | `LOW_SITES` | `MIXED_VOTES` | `NO_INFORMATIVE_SITES`

### `<CHILD>.orientation.dissent.tsv`
Positions that voted against their block's orientation — switch-error and genotype-error candidates.
`chrom pos phase_block_id block_orientation`. May legitimately be empty.

### `<CHILD>.transmission.tsv`
One row per segment of a **parent** phase block.

```
chrom  start  end  parent  transmitted  parent_phase_block_id  segment  n_segments
n_het_phased  n_hap1  n_hap2  vote_frac  reason  left_boundary  right_boundary  change_pos
```
* `parent` ∈ `F` | `M`; `transmitted` ∈ `HAP1` | `HAP2` | `UNRESOLVED`
* `reason` ∈ `OK` | `LOW_SITES` | `MIXED_VOTES` | `NO_INFORMATIVE_SITES`
* `left_boundary`, `right_boundary` ∈ `BLOCK_EDGE` | `CHANGE_POINT`

### `<CHILD>.changepoints.tsv`
One row per within-block change of the transmitted haplotype. At the VCF level a crossover and a
parental switch error are indistinguishable (P3), so every row is emitted as a candidate.

```
chrom  parent  parent_phase_block_id  left_pos  right_pos  resolution_bp
n_left  n_right  left_hap  right_hap  status
```
* `status` is always `CANDIDATE`.

### `<CHILD>.changepoints.resolved.tsv`
The same rows after `xo-reads` has read the parent's haplotagged reads across every gap between
consecutive phased hets in the interval and measured allele concordance.

```
… the eleven columns above, then …
n_parent_hets_in_interval  n_gaps  weakest_gap_start  weakest_gap_end
weakest_gap_informative_reads  weakest_gap_discordant  max_gap_disc_frac
min_spanning_reads_required  child_switch_in_interval
```
* `status` ∈ `CROSSOVER` | `SWITCH_ERROR` | `AMBIGUOUS` | `UNTESTABLE`
  — every gap informative and concordant → `CROSSOVER` (the parent's phasing is intact across the
  interval, so the haplotype really changed); any gap with discordance ≥ `switch_min_disc` →
  `SWITCH_ERROR` located at that gap; otherwise `AMBIGUOUS`. `UNTESTABLE` is declared by the
  classifier and has not occurred on this cohort.
* `child_switch_in_interval` ∈ `Y` | `N` — a located switch in the *child's* phasing inside the same
  interval, which inflates crossover calls; the QC table reports crossovers with and without these.

### Summaries and per-sample files
| file | content |
|---|---|
| `<CHILD>.orientation.summary.json` | per-chromosome and total counters, block counts, `frac_bp_ambiguous`, `mendel_inconsistent_per_informative`, the parameters used, `thresholds_version` |
| `<CHILD>.transmission.summary.json` | per parent: blocks, resolved/unresolved segments and hets, change points, `frac_het_resolved`, `thresholds_version` |
| `<CHILD>.phase_qc.json` | written by `phase-qc`: block counts, crossovers per parent (raw and excluding child-switch overlaps), depth and haplotype balance per member, gate verdicts |
| `$PHASE_DIR/hapdepth/<SAMPLE>.hapdepth.tsv.gz` | `chrom start end dp_hap1 dp_hap2 dp_untagged dp_lowmapq` in 1 kb bins (+ `.summary.json`) |
| `$PHASE_DIR/cohort_phase_qc.tsv` | one row per child with the gate verdict and flags (+ `.summary.json`) |

## 4. Running it

```bash
cp config/phasing.env.example config/phasing.env      # then edit; it is gitignored
# on a machine that already has 13_phase_dnm configured, this reproduces the setup exactly:
#   ln -s ../../13_phase_dnm/config/phase_dnm.env config/phasing.env
```

One family, in order:

```bash
sbatch workflow/m1_orient_family.sb   <FAMILY_ID>    # orient + transmission (VCF level, ~150 s/child, 300 MB)
sbatch workflow/m1_xoreads_family.sb  <FAMILY_ID>    # resolve the change points from parent reads (~9 min/child)
```

The whole cohort, as job arrays (EXPANSE.md rule 4) — the list files are data and live under
`$PHASE_DIR`, never in the repo:

```bash
sbatch --array=1-$(wc -l < "$LIST")%20 --export=ALL,LIST="$LIST" workflow/m1_orient_array.sb
sbatch --array=1-$(wc -l < "$LIST")%20 --export=ALL,LIST="$LIST" workflow/m1_xoreads_array.sb
sbatch --array=1-$(wc -l < "$S")%25   --export=ALL,LIST="$S"    workflow/hapdepth_array.sb   # one task per SAMPLE
trio-phase phase-qc --manifest "$MANIFEST" --phase-dir "$PHASE_DIR" --hapdepth-dir "$PHASE_DIR/hapdepth"
trio-phase audit    --manifest "$MANIFEST" --phase-dir "$PHASE_DIR" --hapdepth-dir "$PHASE_DIR/hapdepth"
```

or the whole DAG: `snakemake all --profile workflow/profiles/expanse -j 40`.

The scripts read `PHASING_HOME` / `PHASING_PYTHON` and fall back to the older `PHASE_DNM_HOME` /
`PHASE_DNM_PYTHON`, so an environment set up for `13_phase_dnm` works unchanged.

Tests: `python -m pytest` (32 tests, a ~300 kb synthetic trio, no cluster needed).

## 5. The audit — what it asserts

`trio-phase audit` exists for the reason `13_phase_dnm/src/phase_dnm/audit.py` exists: every expensive
defect in this pipeline has been something **declared and never produced**, with nothing comparing the
declaration against the artefact. It prints its evidence either way and **exits non-zero on any FAIL**,
so it belongs in the DAG (it is the last rule of the Snakefile) and before any number leaves the cluster.

| check | assertion |
|---|---|
| `trios` | every complete trio of the manifest has all seven artefacts; tables for a child the manifest does not make a complete trio are a WARN |
| `nonempty` | every table has rows (except `orientation.dissent.tsv`, where empty is the ideal result and only a WARN); `changepoints.resolved` has exactly as many rows as `changepoints`; every child has at least one `CROSSOVER` |
| `schema` | the header on disk **is** the writing module's declared column list — imported from the writer, not restated |
| `vocabulary` | every categorical column's values lie inside the declared vocabulary; declared values that never occur are reported as INFO, values outside it are a FAIL |
| `overlap` | orientation segments do not overlap **within one phase block**, and transmission segments do not overlap within one `(parent, phase block)` |
| `summaries` | both summaries parse and carry a `thresholds_version`, and the cohort was produced by exactly one |
| `hapdepth` | with `--hapdepth-dir`: every manifest sample has a readable depth table and a summary |

**Overlap: within a block, not across blocks.** The invariant that holds is *within* a phase block —
a block is split only at a located switch, so its segments are disjoint by construction and an overlap
would give a position two orientations. Across blocks it does **not** hold, and asserting it would be a
false alarm: HiPhase emits small blocks nested inside larger ones. Measured on the cohort of
2026-09-16 (35 complete trios, 491,982 orientation segments): **0** within-block overlaps, **193**
across-block overlaps, only **2** of them between two oriented (non-`AMBIGUOUS`) segments. The audit
therefore FAILs on the first and WARNs on the second; `--strict-across-block` promotes the WARN to a
FAIL for a caller that needs an exact cover.

Cohort baseline the audit reproduced on 2026-09-16 (0 FAIL, 1 WARN, 40 checks):

```
35 complete trios, all 7 artefacts each
orientation          7,102- 22,080 rows/child      AMBIGUOUS 263,469  HAP1_PAT 114,765  HAP1_MAT 113,748
orientation.dissent    214-    693 rows/child
transmission        11,765- 36,436 rows/child      UNRESOLVED 462,019  HAP1 225,397  HAP2 207,304
changepoints           452-    570 rows/child      status CANDIDATE 18,315 (all)
changepoints.resolved  452-    570 rows/child      AMBIGUOUS 14,711  CROSSOVER 2,987  SWITCH_ERROR 617
hapdepth               105/105 samples
```

## 6. What consumes this module

`13_phase_dnm` (phase-aware de novo calling) reads the tables **by path** under `$PHASE_DIR` — its
`config/phase_dnm.env` names the same directory, and `workflow/m2_*.sb` and `m4_*.sb` pass
`--orientation`, `--transmission` and `--changepoints` to the review, feature, spike and annotate
steps. There is deliberately **no Python import** in either direction: the interface is the file format
documented in §3, and `13_phase_dnm`'s own `LabelTables` (in `evidence/hapmatrix.py`) is the reader.
`13_phase_dnm/workflow/Snakefile` lists `$PHASE_DIR/<family>/m1_xoreads.status` as an input of
`m2_review`, so a family whose phasing has not been run fails that DAG at build time rather than
reviewing against absent labels.

Anything else that wants parent of origin, transmission or recombination positions should read the
same tables the same way.

## 7. Decisions taken in the extraction

1. **Nothing in `evidence/hapmatrix.py` moved.** The only phasing-related thing there is `LabelTables`,
   which *reads* the orientation / transmission / changepoint tables for the six-haplotype matrix. It is
   consumer-side code, and moving it would have created exactly the import dependency from
   `13_phase_dnm` into this module that the split is meant to avoid. It stays where it is.
2. **`PLAN_MODULE1.md` stayed in `13_phase_dnm`.** Despite the name it is a running log covering M2–M4
   as well, and `DESIGN.md` — which must not be edited — cites it by bare filename.
3. **The sbatch scripts kept their `m1_` filenames.** Renaming them on top of moving them would have
   turned a reviewable rename into an add/delete pair for no functional gain. The prefix is historical.
4. **`config/thresholds.yaml` is a verbatim copy, version string and all.** The `version` is written
   into every summary JSON, so bumping it would change output bytes. From now on the two files version
   independently, which is the point: a change to the de novo decision layer no longer re-versions the
   phasing tables. *(Observation from the cluster, 2026-09-16: the tables on the filer carry
   `thresholds_version` `0.1.0` while `13_phase_dnm` has since moved to `0.3.0` for M2/M3 reasons — a
   pre-existing drift that this split prevents recurring.)*
5. **`PAR_GRCH38` was restated, not imported.** `13_phase_dnm/src/phase_dnm/features/extract.py` used to
   import the GRCh38 pseudoautosomal coordinates from `phasing/orient.py`; it now defines the same
   two-interval constant locally, with a comment pointing here.
6. **`io/vcf.py` was copied, not moved** (see §1).
