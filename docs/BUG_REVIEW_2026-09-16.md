# Pipeline-wide bug review — 2026-09-16

A read-only review of every module for concrete defects, requested by JS. Four reviewers each covered a disjoint
module set; every HIGH item and the most consequential MEDIUM items below were then **verified at source by a second
reader** before being listed, and several were confirmed by reproducing them (marked *repro*). Items the second
reader could not confirm from the code are marked *reported*. Nothing here is a style comment.

Severity: **high** = wrong results or data loss; **medium** = crash or missed data on realistic input; **low** =
robustness / documentation. Status: **FIXED** (commit named), **OPEN**, or **DECISION** (fixing it changes a result
and is JS's call).

---

## 05_denovo — phase-aware de novo calling

| # | sev | where | defect | status |
|---|---|---|---|---|
| D1 | high | `workflow/m2_review_family.sb:51` | A literal `\n` (backslash, n) sat inside the `review` command line, so bash passed the word `n` as a stray positional argument → argparse exit 2 → `set -e` aborts for every family. Introduced in `6110a3a`; the cohort review predates it. *verified by bytes* | **FIXED** (this commit) |
| D2 | high | `workflow/Snakefile` rule `m2_reclassify` / `m3_integrate` | Nothing in the `all → m3_integrate → m2_reclassify → m2_review` chain depended on `annotate_real`. `m2_reclassify_family.sb` passes `--annot` only if the per-child table already exists, and `integrate` applies the population rules (gnomAD, leave-one-family-out cohort count, founder panel, sib-shared) only when those columns are filled — so a fresh `snakemake all` produced tier-1 calls with the mask rule alone and a normal summary. `m4_annotate.sb` and reclassify also both wrote `features.tsv` with no ordering. *verified* | **FIXED**: `m2_reclassify` now requires `annotate_real` for every class group |
| D3 | high | `evidence/hapmatrix.py::reclassify_row` | (a) `sv_row` values came straight from `csv.DictReader` as strings; `classify_sv_interval`/`_rule_score_sv` compare them with floats → `TypeError`. (b) The filter kept the raw `C_sv_*` columns but dropped the derived block (`c_sv_junc_both_ends`, …) that `classify_sv_interval` reads, so every reclassified deletion fell to `inconclusive / NO_JUNCTION_BOTH_ENDS`. Neither fired in production because the one cohort reclassify of SV tables ran while the `sv_` columns were blank — i.e. the "carry SV evidence through reclassify" intent had never executed. *repro (both)* | **FIXED** + regression test `tests/test_reclassify_sv_strings.py` |
| D4 | high | `evidence/readers.py::_sv_signature_matches` (SA branch) | Only `read.reference_start` is compared to the breakpoints, never `reference_end`, so a junction read whose primary segment lies to the LEFT of a breakpoint (clipped at its 3′ end) can never match; and the inner `any(... for bp in bps)` re-binds `bp`, making the "other breakpoint" constraint a no-op. Affects every junction read not named in sawfish's supporting-reads list — by design all of them in the spike arm and in synthetic trios. The planted-deletion result "36 of 47 land as inconclusive at 20/50 kb" that motivated routing large SVs to the clinical arm (P31/P32) is the expected signature of this bug. *verified by reading* | **OPEN** — fix, then re-measure planted recall ≥ 20 kb before P32's conclusion is repeated |
| D5 | high | `features/extract.py` (`multiallelic`, `rf_safe: true`, in the frozen SNV/indel model) | Positives (synthetic trios) read `n_alts` from the 105-sample cohort BCF; negatives from the 3–4-sample family VCF. Measured: synthetic 0.34–0.36 vs real 0.38–0.55 (family-dependent). A systematic, non-variant difference the classifier can learn — the mechanism that retired `site_qual`. *measured on the cohort* | **DECISION** — mark `rf_safe: false` or compute from one source; either changes the frozen model |
| D6 | high | `workflow/Snakefile` (`SCORE_MODE=frozen` default) | The frozen models were refit on all 33 families with every real candidate (true DNMs included) as label 0, and `score.py` builds `rf_q` from the scored rows' own ranks; scoring the training cohort with them is in-sample. Correct for a new cohort (P21), wrong for this one — the reported call set used out-of-fold rescoring (`SCORE_MODE=train`). *verified* | **FIXED (documentation)**: stated in the Snakefile; consider auto-detecting manifest overlap with the training manifest |
| D7 | medium | `workflow/Snakefile` rule `m4` synthetic loop | `bash a && bash b` inside a `for` under `set -e`: a failing first command does not abort, `.synthetic.done` is touched, `m4_train` runs with missing positives. *reported, consistent with bash semantics* | OPEN |
| D8 | medium | `workflow/m4_annotate_synthetic.sb` (TR) | `CV=""` for `tr` → synthetic TR annotation has empty `cohort_AC_loo` → `fillna(0)` reads missing as private → the P24 rarity gate keeps every TR positive; missing annotation is indistinguishable from "private". *reported* | OPEN |
| D9 | medium | `annotate.py::sib_shared_sites` | Prefilter compares raw `rec.pos` with NORMALISED positions (anchor-base indels shift +1), so `sib_shared` is never 1 for indels. *reported* | OPEN |
| D10 | medium | `cli.py::cmd_external` | Reads `tau_q` from the harness `tau.<cls>.json` (hard-coded 0.999 in `rescore.py`/`score.py`) not from `thresholds.yaml final.tau_q`; external-truth recall is reported at an operating point the pipeline does not run at (`tools/sv_size_recall.py` was patched for this; `external.py` was not). *reported* | OPEN |
| D11 | medium | `tests/test_hapmatrix.py::test_reclassify_from_evidence_row_reproduces_the_rule_layer` | Failed at HEAD: it appended `NEAR_CHILD_SWITCH` and still expected the original parent of origin, but today's P34 change deliberately clears PoO inside a switch window. Test never re-run after the change. *repro* | **FIXED** (test now expects `undetermined` + `PHASE_SWITCH_RISK`, and checks the flag-free round trip) |
| D12 | low | `sim/edit.py::apply_breakpoint` | Leading soft clip takes `continue` before `q += qry_len`; planted junction-read bases are read from the wrong query offset (lengths stay consistent, so the all-"C" test cannot see it). | OPEN |
| D13 | low | `classify/likelihood.py:139` | UNT row is `(alt, dp − alt)`, counting untagged AMB reads as REF, unlike tagged rows. | OPEN |
| D14 | low | docs | `models/FROZEN.md`/P21 record a different sha, feature count and row count from the shipped `snv_indel.training_manifest.json`; `integrate.py` docstring describes the pre-0.3.0 rescue rule; README says `dnm_call ∈ {YES,NO}` but `CANDIDATE` is emitted; `rf_prob` column carries `rf_q`. | OPEN |
| D15 | low | `tests/test_annotate_tr.py` | Asserts on `inspect.getsource` substrings; cannot fail on behaviour. | OPEN |

## 02_phasing / 01_qc / 00_upstream / scripts

| # | sev | where | defect | status |
|---|---|---|---|---|
| P1 | high | `01_qc/family_downstream.sb:316-320` | TRGT dump `bcftools query … > "$DUMP" 2>/dev/null`: a failed query yields an empty dump, "loci dumped: 0", every child "0 candidate expansions", exit 0. The SNV section of the same script guards this; the TR section does not. *verified* | OPEN |
| P2 | medium | `02_phasing/src/trio_phase/phasing/xo_reads.py:161` | `hets = hets[::step] + [hets[-1]]` duplicates the last het whenever `(len−1) % step == 0`; `pos_index` then maps that position to the last index only, het `n−2` never receives support, the two gaps around it read 0 informative reads and the change point is forced to `AMBIGUOUS`. Any interval with > 400 parental hets and the wrong parity cannot be called `CROSSOVER`. *verified by arithmetic* | OPEN |
| P3 | medium | `02_phasing/…/cli.py:66` vs `qc.py:127` | `cmd_orient` writes the raw manifest `sex`; the depth gate fires only for `"M"/"F"`. A PED-coded manifest (`1/2`) silently disables `SEX_DEPTH_MISMATCH`. *reported* | OPEN |
| P4 | medium | `02_phasing/workflow/Snakefile` vs `audit.py:421` | hapdepth computed only for complete-trio samples but the audit demands every manifest sample → `rule audit` fails on any duo. *reported* | OPEN |
| P5 | medium | `02_phasing/workflow/hapdepth_array.sb:33`, `m1_orient_family.sb:42` | `X=$(ls … | head -1)` under `set -euo pipefail` exits silently (status 2) when the glob matches nothing; the documented fallbacks and error messages are unreachable. *reported, consistent with bash semantics* | OPEN |
| P6 | medium | `01_qc/fam_downstream_array.sb:4-9` | Resolves the worker under `$DATA_ROOT` (the output root), not the repo. *reported* | OPEN |
| P7 | medium | `00_upstream/cohort/glnexus_full.sb:99-107` | Publishing the cohort BCF sits inside `if [ -n "$BASE" ]`; with no container found the callset is never published and the job prints "done", exit 0. *reported* | OPEN |
| P8 | medium | `.gitignore` `*.tsv` vs `00_upstream/README.md` | The two `*.template.tsv` map files the README tells users to fill in are git-ignored and absent from a clone. *verified with git check-ignore by the reviewer* | OPEN |
| P9 | low | multiple | Hard-coded account/partition/personal paths in `01_qc/*.sb` despite `config/cohort.env.example`; README/`pyproject`/error strings still say `01_phasing`/`13_phase_dnm`/`03_panel`; `thresholds.yaml hapdepth:` section never read; `test_qc.py:66` right-hand side always true. | OPEN |

## 03_tiering / 04_panel / 07_inheritance / 08_burden

| # | sev | where | defect | status |
|---|---|---|---|---|
| T1 | high | `03_tiering/lof_table.sb:23` | `sort -u -k1,1 -k2,2n` de-duplicates on chrom+pos only: distinct HC-LoF alleles at one position (multiallelic, or SNV + indel) collapse to one row in `lof_sites.tsv` and everything downstream. *verified* | OPEN |
| T2 | high | `03_tiering/lr_vep.sb:62` + `08_burden/lof_persample.sb:31`, `07_inheritance/inherited_lof_table.py:93` | gnomAD AF comes from the dbNSFP plugin, an SNV-only database; every indel LoF (frameshifts — the bulk of HC LoF) has an empty AF, which the rare filter maps to −1 = "absent = rare". Common frameshifts enter `lof_tiered.rare.tsv`, the burden counts, the TDT and the inherited tables. Contradicts METHODS §4. *verified* | OPEN — needs a gnomAD genomes VCF annotation for indels |
| T3 | high | `03_tiering/tr_cohort.sb:13` | `I="$WORKFLOW_ROOT/02_tiering"` — no such directory; `known_repeats.py` and `tr_outliers.py` cannot be found, rc echoed only. *verified* | OPEN (one-line fix) |
| T4 | high | `07_inheritance/inherited_report.py:145` | `bcftools query -i 'INFO/gnomad_af<0.005'` is false for a missing tag: variants absent from gnomAD — the rarest class — never reach any model; the `af=-1 → absent` branch is unreachable. *verified (bcftools missing-value semantics)* | OPEN |
| T5 | medium | `04_panel/apply_pon_loo.sb:32-36` | Panel-building pipeline inside `bash -c` has no pipefail; `&& index` binds only to the last `view -G`. Whether a header-only panel results depends on bcftools' handling of a closed upstream; the "sites=0" line is not checked, unlike `build_pon.sb:33`. *reported; plausible* | OPEN (add the zero-sites abort) |
| T6 | medium | `08_burden/lof_persample.sb`, `tdt.py`, `07_inheritance/inherited_lof_table.py` | Sites keyed by (chrom,pos), not allele; at co-located records one allele's tier is applied to carriers of another. *reported* | OPEN |
| T7 | medium | `08_burden/tdt.py:108-129`, `03_tiering/tr_outliers.py:154` | chrX: a father's hemizygous X allele (`count("1")==1`) is treated as a heterozygous informative parent; PED sex unused. *verified* | OPEN |
| T8 | medium | `07_inheritance/prioritised_consequence.py:63,75` | Insertion frame uses `end − start` of the POS–END interval (0 or 1 bp), not the insertion length carried in `detail`; every coding insertion's frameshift/in-frame call is from a 0 or 1. *verified* | OPEN |
| T9 | medium | `07_inheritance/inherited_tr_table.py:103`, `03_tiering/known_repeats.py:200` | Spanning-read QC uses `max(SD)` over alleles, so the expanded allele's own support is never tested. *reported* | OPEN |
| T10 | medium | `04_panel/apply_pon_loo.sb:49,57` | SVs annotated by exact CHROM/POS/REF/ALT while METHODS §3 says reciprocal-overlap matching; breakpoint drift → "absent from panel". *reported* | OPEN |
| T11 | medium | `03_tiering/tr_cohort.sb` | Required env (`SAMPLE_QC`, `SIF_CACHE`, `PON_INVARIANCE`) not passed/validated; the D55 gate silently OFF with a WARNING while METHODS says it is applied. *reported* | OPEN |
| T12 | medium | several `.sb` | `2>/dev/null` + `wc -l` swallow split-vep failures; mask build without pipefail can write a valid empty `.gz`; resume guards accept truncated files; `--force-samples` silently drops absent panel members. *reported* | OPEN |
| T13 | low | docs/config | COMPOUND_HET documented, not implemented; DDG2P claimed in `inherited_lof_table.py`, code uses SFARI ∪ LOEUF; "exact binomial" is a normal approximation with p > 1 at T = NT; `cohort.env.example` variables defined but paths hard-coded module-wide; `gencode_cds.bed` col 4 is `GENEA,GENEB` at overlaps and is compared as one name. | OPEN |

## 09_transmission / 10_ascertainment / 11_phenotype / 12_methylation / 13_episignatures / 14_x_inactivation

| # | sev | where | defect | status |
|---|---|---|---|---|
| X1 | high | `14_x_inactivation/05_escape_status.py::load_islands` (l.62-86) | Never reads the chromosome column: pools `(pos, beta)` from every row of the combined pileup and bins chrX intervals against it. The pileup regions are chrX islands **plus autosomal islands** (README l.63), so every chrX island mean is contaminated by autosomal CpGs at the same numeric coordinate. `04_xci_skew.island_means` keys by chromosome; 05 does not. *verified at source; consequence depends on the pileups actually spanning autosomes, which the README specifies* | OPEN — every escape call from 05 must be recomputed |
| X2 | medium | `14_x_inactivation/02_hap_pileups.sb:60` | Passes `--regions` to `aligned_bam_to_cpg_scores`; the installed pb-CpG-tools (`--help` checked) has no such option, so every array task dies on the argument parser under `set -e`. *verified against the binary's help* | OPEN — the module's documented pileup step cannot have run as written |
| X3 | medium | `14_x_inactivation/02_hap_pileups.sb:49` | `samtools view | head -2000 | grep -q` under `pipefail`: `head` exits early, samtools gets SIGPIPE, the `if !` fires and a valid tagged BAM is rejected with "no MM/ML tags". *verified* | OPEN |
| X4 | medium | `14_x_inactivation/04_xci_skew.py:143` vs `03_trio_phase_x.py:279` | Trio mode requires `block_start`, which 03 emits only with an optional `--blocks-dir`; run as documented every female falls silently to the folded estimator. *reported* | OPEN |
| X5 | medium | `14_x_inactivation/07_denovo_x_screen.py:274` | Father GT must equal the string `"0/0"`; haploid male X genotypes (`"0"`) — the module's own premise elsewhere — can never pass, so `denovo_candidates.tsv` is empty by construction. *reported* | OPEN |
| X6 | medium | `12_methylation/meth_stage2b.py:306-313` | Per island the minimum p over all SNVs within ±10 kb is taken (multiplicity uncounted) and then labelled `fdr_bh` with m = islands, without step-up monotonicity; the printed FDR < 0.05 count is anti-conservative. *reported* | OPEN |
| X7 | medium | `12_methylation/meth_stage1.py:77-87` | A profile with a different row order/count is written by positional index against the first sample's islands after only a WARNING. *reported* | OPEN |
| X8 | medium | `12_methylation/*`, `deconv.R:23` | DNA source inferred from a sample-ID prefix rather than the manifest `cohort` column; on another cohort the epithelial-fraction zero point becomes NaN and every founder covariate is imputed. *reported* | OPEN |
| X9 | medium | `11_phenotype/severity.py:145-150` | `frac()` counts SPARK's "not asked" code as a negative response in every categorical percentage; `valid()` in the same script excludes it for the Fisher tables. *reported* | OPEN |
| X10 | medium | `09_transmission/imprinted_tdt_all.py:141` | Output directory never created; on a fresh `DATA_ROOT` the script does all its tabix work then dies with `FileNotFoundError`. *reported* | OPEN |
| X11 | low | several | Duplicate-family resolution differs between `gene_tdt.py` and `imprinted_tdt_all.py`; hard-coded personal paths in `10_ascertainment/ascertain.py`; `probands.txt` read without skipping `#`; Mann-Whitney without tie correction on a 4-level ordinal; `.sb` files without `set -e`; `13_episignatures` report text says ≥15× while config sets 10; a typed-in annotation string; substring version detection; `rssh() … \|\| true` hides remote failures; README numbering drift (`12_x_inactivation`, "methylation (09)"). | OPEN |

## 06_clinical (built today)

| # | sev | where | defect | status |
|---|---|---|---|---|
| C1 | medium | `src/clinical/smallvar_annot.py` | 43,053 of 93,432 panel-span SNV/indel candidates (46%) are absent from the freeze-1 VEP VCFs and therefore have no consequence in arm B. Most are low-quality raw candidates the freeze excluded, but the gap is not characterised. The recorded VEP/LOFTEE resources of `03_tiering` are absent at their paths (2026-09-16), so re-annotation is blocked upstream. | OPEN — list written to `smallvar_annotation.tsv.missing` |
| C2 | low | `src/clinical/quality.py::tr_quality` | REVIEW gate is permissive (CANDIDATE, or rule ≥ 3, or rf ≥ 0.3): 141 TR REVIEW rows across 34 children in the first run. Tunable; left as is because arm B's job is to surface. | OPEN (tunable) |

---

## What was NOT reviewed

Runtime data was not inspected except where a finding says *measured*. External-tool semantics were verified only where stated (bcftools missing-value comparison, GNU sort `-u` with keys, pb-CpG-tools `--help`, bash `set -e` with `&&` lists). The parent-of-origin orientation logic was excluded because it was validated per call the same day (05_denovo DESIGN P35). Design choices (parents as clogit controls, cohort-wide references, exact-match SNV panel) were treated as intended. Nothing was executed except read-only checks.
