# DESIGN — 13_phase_dnm

Decisions are numbered **P#** and each carries its rationale, so `docs/METHODS.md` can cite
them. `ASSUMED` marks a decision made without an answer to an open question; it is revisited
when the question is answered. Started 2026-09-12 from the public repo
(`jsebat/spark-lrwgs-cohort`), `docs/METHODS.md`, the current `james-guevara/synthdnm` code,
`docs/EXPANSE.md`, the cohort report draft (2026-09-12), and a direct inspection of the cohort
directory on the lab filer (via the WSL multiplexed ssh; all identifiers redacted before they
left the cluster).

---

## 0. What the repo, code and data actually say (corrections to the brief)

| brief assumed | what is actually there | consequence |
|---|---|---|
| SVs from pbsv | **sawfish 2.2.1** joint calling (347,631 records) | HiPhase's documented inputs are DeepVariant/pbsv/TRGT, but **the WDL did run HiPhase 1.6.0 on the sawfish joint VCF** — `phased_sv_vcf/` exists with `PS` per sample (Q3 closed). SV allele haplotypes are therefore available from the VCF *and* from reads. |
| small variants from DeepTrio | **DeepVariant + GLnexus v1.4.3** (no trio-aware caller) | the "DeepTrio-only" comparison requires running DeepTrio, at least on GIAB (P16). |
| SynthDNM = random forest | paper (verified) is RF; the **current repo trains XGBoost** (`train.py`; feature sets universal 21 / gatk 28 / ssc 29) | M4 targets the XGBoost code path; RF is a baseline on the same folds. |
| swap = child paired with random surrogate hom-ref parents | `swap_pedigree.py` **pairs families A↔B and exchanges offspring** | swap-closure unit is the *pair* (P12). |
| SynthDNM CV is family-grouped | `train.py`: `StratifiedKFold` on rows + `train_test_split`; no grouping, nesting or calibration | M4 is new training code. |
| 33 trios | **31 SPARK trios + 2 SPARK mother–child duos (father not sequenced) + 1 SPARK quad (2 affected) + 1 REACH quad (affected + unaffected)** = 37 offspring, 35 probands (report §2) | duos have no paternal haplotypes: transmission on the paternal side is undefined (P19). Complete-trio children for M1/M4: 33. |
| multi-sib families as a mosaic truth set | **2 quads** | illustration only (R2). |
| candidates include sib-shared variants | report §3.5 / METHODS §7: slivar de novo model **excludes calls shared by both siblings** | M2 ingests from the family joint VCF, before that exclusion (P5, P6). |
| coverage [X] | cases 22.4×, controls 24.4× (METHODS §8); ~11× per haplotype | sets k (P7) and the mosaic floor (P10). |

Citation check: SynthDNM = Lian A, Guevara J, Xia K, Sebat J. *Bioinformatics* 2021;37(20):3640–3641, doi:10.1093/bioinformatics/btab225, PMID 33821956 — **verified** (Europe PMC). HiPhase = Holt et al., *Bioinformatics* 2024;40(2):btae042 — **verified** (title/journal). Anything else below marked "unverified" I could not check.

### 0.1 Data inventory (confirmed on the filer, 2026-09-12; no identifiers)

`$LRC` = the long-read cohort directory on the persistent lab filer (`/expanse/projects/sebat1/...`; exact prefix in `config/phase_dnm.env`).

| what | where | notes |
|---|---|---|
| manifest, 105 rows | `$LRC/MANIFEST_lrWGS_105.tsv` | `sample_id family_id father_id mother_id sex affected role lr_haplotagged_bam lr_family_sv_vcf lr_family_smallvar_vcf lr_family_trgt_vcf`. 36 affected + 1 unaffected offspring; 65 unaffected + 3 affected parents. This is the `manifest.tsv` the module takes. |
| joint small variants | `$LRC/freeze1/freeze1.cohort.bcf` (11 GB, 105 samples) | FORMAT `GT RNC DP AD GQ PL`; INFO `AF AQ AC AN`. SynthDNM's *universal* set is fully computable (Q8 closed). |
| joint SVs | `$LRC/freeze1/freeze1.cohort.sv.vcf.gz` | FORMAT `GT GQ PL AD CN CNQ PS`; INFO `SVTYPE SVLEN END HOMLEN HOMSEQ INSLEN INSSEQ SVCLAIM EVENT MATEID IMPRECISE` |
| joint TRs | `$LRC/freeze1/freeze1.cohort.trgt.vcf.gz` (TRGT 5.0.0) | FORMAT `GT AL ALLR SD MC MS AP AM PS`; INFO `TRID END MOTIFS STRUC` |
| panel of normals | `$LRC/freeze1/pon/freeze1.pon65.{snv,sv,trgt}.sites.vcf.gz` + `trgt.invariance.tsv` | 65 unaffected founders |
| per-family WDL outputs | `$LRC/callsets/<FAMILY>/out/<item>/<0..n>/` — 35 families | indexes in parallel `<item>_index/` dirs |
| — haplotagged BAMs | `merged_haplotagged_bam/` | verified: `HP`, `PS`, `MM/ML`, `rq` on every read in a probed window |
| — HiPhase blocks | `phase_blocks/*.hiphase.blocks.tsv` | `phase_block_id chrom start end num_variants` |
| — **HiPhase haplotags** | `phase_haplotags/*.hiphase.haplotags.tsv.gz` (~2.8 M rows/sample) | `chrom phase_block_id read_name haplotag` — **this is the read→HP sidecar; M1 does not have to build it** |
| — HiPhase stats | `phase_stats/*.hiphase.stats.tsv` | per chromosome: `num_phased`, `num_blocks`, `block_ng50`, … → `phase_qc.json` input |
| — phased VCFs | `phased_{small_variant,sv,trgt}_vcf/` | per-sample split of the family joint VCF; `PS`, `PF` |
| — family joint small-variant VCF | `joint_small_variants_vcf/` | all members; **the UNFILTERED candidate source** (P5) |
| — **TRGT spanning reads** | `trgt_spanning_reads/*.trgt.spanning.sorted.bam` (~1.8 GB/sample) | per read: `AL:i` (allele length), `MC`, `HP`, `PS` — **the TR adapter's alt-support hook is already computed per read** |
| — **sawfish supporting reads** | `sv_supporting_reads/*.supporting_reads.json.gz` | `{sawfish_id: {sample: [read_name, …]}}` — **the SV adapter's alt-support hook, joinable to haplotags by read name** |
| — gVCFs, `sv_discover_tars`, `sv_depth_bw`, `sv_copynum_bedgraph`, `mosdepth_*` | same tree | per-haplotype depth still has to come from the BAM (`hapdepth`) |
| short-read WES | SPARK iWES v3 DeepVariant pVCF (cohort-wide); `$LRC/shortread/wes_crams/` (98 samples), `xhmm/` (101) | 101/105 lrWGS samples are in the WES ped (REACH family absent). **Exonic SNV/indel concordance only.** GraphTyper on WES CRAMs → exonic SV candidates only. |
| short-read WGS | 3 CRAMs locally; SPARK iWGS v1.1 release exists cohort-wide | not usable as a per-trio comparator |
| short-read TR | none | TR has **no** orthogonal comparator; validation = spike-ins + family segregation + methylation-on-the-same-molecule where applicable |
| existing de novo set (baseline for the concordance table, P18) | `L/tiering/denovo_tiered.tsv` (1,332 SNV/indel across 35 probands, median 38), `L/tiering/denovo_sv/denovo_sv_PRIORITIZED_v2.tsv` (coding, disease-gene SVs), `L/tr/tr_denovo_expansions.tsv` (191 loci) | `L` = the Lustre analysis directory; **was unreachable 2026-09-11** — copy to the filer before use |
| lab operating rules | `docs/EXPANSE.md`; wrappers `longread-autism/cohort/lab_commands.sh` (`lab_submit`, `lab_status`, `lab_usage`, …) | every job: show resources first, inspect after submitting, state an ETA |

### 0.2 The P13 number, measured (chr22, first callset family, two children)

Child het + both parents hom-ref in the family joint VCF, no other filter:

| filter | chr22, child 1 | chr22, child 2 | ×75 ≈ genome |
|---|---|---|---|
| none (raw) | 378 | 324 | **~24–28 k** |
| min trio GQ ≥ 20 | 37 | 23 | ~1.7–2.8 k |
| min trio GQ ≥ 30 | 4 | 6 | ~300–450 |
| GQ ≥ 20 and DP ≥ 10 in all three | 30 | 22 | ~1.6–2.3 k |
| indel fraction of raw | 0.18 | 0.18 | |

Against ~70 true germline DNMs per genome (**unverified** exact figure; order of magnitude): raw negatives carry ~0.3 % true DNMs — SynthDNM's negative-set assumption **holds** on HiFi if negatives are the raw set, as Jon confirmed is SynthDNM's practice. At GQ ≥ 30 the contamination would be 15–25 % and the scheme would fail. Decision P13 follows. The current pipeline's filtered set (median 38/proband) is ~half the expected rate, so rescue is plausible in principle; how much of the gap is the 11 % mask vs. caller filtering is one of the first things M2 measures.

---

## 1. Open questions

**Closed 2026-09-12** (answers recorded above): Q1 negatives = the full unfiltered putative set, no quality stratification up front; low-quality calls stay in the all-rows RF output table and drop out only from the final list. Q2 measured (0.2). Q3 HiPhase phased sawfish SVs. Q4 BAMs tagged, on the filer. Q5 WES only (exonic SNV/indel; WES CRAMs for exonic SV re-genotyping; nothing for TR). Q8 all FORMAT fields present.

**Closed 2026-09-12 (second round, JS)**
- **Q6** GIAB HG002/HG003/HG004 HiFi is **not** on Expanse. It has to be pulled (a download job on `ind-shared`, resources and ETA shown first per EXPANSE.md rule 3; release/coverage to be chosen to match ~22–24× — proposal in `PLAN_MODULE1.md`). Until then GIAB is a week-4 item, and P16 holds.
- **Q7** SynthDNM has **never** been run on the long-read Freeze-1 callset (only short-read SPARK WES / SSC). We start from `main`; the "SynthDNM as published" baseline (P14) is therefore our own run of the unmodified *universal* model on our folds — there is no prior lrWGS run to reproduce or to worry about leakage from.
- **Q9** Snakemake with an Expanse profile — agreed.
- **Q10** one Apptainer image — `ASSUMED` yes (`singularitypro/4.1.2`); the WDL's own `pb_wdl_base` image (in the miniwdl singularity cache) already carries samtools/bcftools and is what `05_denovo/sv_table.py` uses — reuse it for the read-level steps until the module image exists.
- **Q11** closed by inspection: HiPhase already wrote the read→haplotag table (P4).
- **Q12** **One paper, the pipeline paper**, with the SPARK lrWGS cohort as its proof of concept. There is no standalone DNM-method paper. Consequences: this module is a *section* of that paper — the "phase as evidence" framing is kept, but validation is scoped to what the pipeline paper needs (P20); the GIAB comparison becomes the reproducibility anchor for the pipeline section rather than a method-paper benchmark; nothing in R1–R10 relaxes.
- **Q13** classic swap first; haplotype swap as ablation with the chimera test — `ASSUMED`, unchallenged.
- **Q14** **Exclude the duos** (P19 default). Optional half-trio M2 rows only if useful for a specific candidate.
- **Q15** the fuller de novo SV list exists as an intermediate, never printed: `denovo_sv_priority.sb` step 3 aggregates every family's `run_<FAM>/analysis/denovo_sv/<FAM>.<PROBAND>.denovo_sv.tsv` into `$T/denovo_sv/denovo_sv_all.bed` **before** the CDS intersect. P18 uses that (and the per-family TSVs) as the SV baseline, not the coding-prioritised table. Both are on Lustre.

**Process constraint (JS, 2026-09-12):** a GitHub clean-up is in progress in another session. **No commits, pushes, or `git pull` on the on-cluster checkout until that is finished.** The module stays staged locally until then.

---

## 2. Decisions

### P1 — Module numbering, scope, and reuse
New directory `13_phase_dnm/`. Reuses, not duplicates: `09_methylation/pofo_local.py` (nearest-informative-site rephasing with phase-switch flag → `local_rephase_agreement`; fallback parent-of-origin for untagged reads); `05_denovo/sv_read_review.sh` (junction detection by `SA` tag / CIGAR `D ≥ MINBIG`, per-haplotype depth inside vs flanks, het-SNV persistence → SV adapter and two SV features); `02_tiering/tr_outliers.py` (TRGT `AL`/`SD`, spanning-read QC, margin rule → TR adapter and TR candidate generator). GraphTyper's three-way outcome (METHODS §5) is an orthogonal *label* for exonic SVs, never a feature.

### P2 — Phasing stack: keep HiPhase 1.6.0, add orientation, do not re-phase
Read-based phasing stays HiPhase (already run per family; jointly phased small variants, sawfish SVs and TRGT). Chromosome-scale, parent-of-origin-labelled child haplotypes come from **orienting each existing child block** with a Mendelian vote at informative sites (one parent het, other hom): the block is `HAP1_PAT`/`HAP1_MAT` when the majority is ≥ `orient_min_frac` (0.95) with ≥ `orient_min_sites` (20), else `AMBIGUOUS`; dissenting sites are `n_mendel_err` with positions kept.
*Phase-switch errors inside a block* (added after the week-1 fixture showed them): a block whose informative votes run one sign and then the other is not "mixed", it is two correctly phased halves joined by a switch error. `orient.py` segments the position-ordered votes by an exact dynamic programme over run boundaries (`segment_votes`): at most `max_splits` (3) cuts, each of which must raise the total majority count by ≥ `split_min_sites` (2) — so a lone dissenting vote is a genotype error and never splits a block, a run of two at a block edge does. Each segment is then oriented on its own (`min_sites`, `min_frac`), so a switch near a block end gives one oriented segment plus a short `AMBIGUOUS: LOW_SITES` tail, and a switch-and-back gives three segments. A block whose segmentation orients nothing is reported once as `AMBIGUOUS: MIXED_VOTES` — the pedigree-error / sample-swap signature — and counts no switches. The fixture that forced this: a 60-vote block with two paternal votes then 58 maternal (vote fraction 0.967, so it had passed as "clean") whose first 13 kb were mislabelled. Consequence for every downstream lookup: a read's parent-of-origin label is keyed by `(phase_block_id, position)`, not by block id alone. The number of located switches per sample is an M1 QC number to compare with HiPhase's published rate (R6).
*Measured on the two quads (2026-09-12, four children):* 9.9–16.4 k HiPhase blocks each; 1.9–2.2 M informative sites; Mendelian-inconsistent per informative site 0.0006–0.0011; 313–367 switches located per genome (≈ 1 per 8 k phased hets — same order as HiPhase's published 1 per 3.3 k variants, R6); ambiguous block-bp 9.4–10.8 %, of which `MIXED_VOTES` (the pedigree/phasing-problem signature) is only 0.15–0.19 %, the rest short blocks. **Independent check:** at the 18 imprinting control regions, `09_methylation`'s nearest-informative-SNV parent-of-origin calls agree with the oriented segment in 30 of 30 comparable cases (0 discordant). *Second tier:* a segment with 10–19 informative votes is oriented only if unanimous (`small_min_sites`, `small_min_frac`); with a per-site inconsistency rate of ~0.1 %, ten independent errors agreeing is not a realistic failure mode, whereas the realistic one — a mis-phased block — produces mixed votes and stays `AMBIGUOUS`. Recovers ~2.6 % of bp per genome.
*Tradeoffs.* WhatsHap `--ped` (PedMEC) would re-phase all three samples jointly and give longer blocks, but is slow at 22× genome-wide, SNV-only, and would create a second phasing disagreeing with the `HP` tags already on 2.8 M reads/sample. Orientation is one pass over the VCF, keeps one phasing, and its failure mode is measurable (`AMBIGUOUS` rate, `n_mendel_err`). The 18 imprinting-control regions already phased to parent of origin by `09_methylation` are free ground truth for orientation. WhatsHap `--ped` is an M1 QC comparison on 2–3 trios, chr20–22 only.

### P3 — Transmission map and crossovers
For each parent, inside each *parent* phase block, compare the parent's two haplotypes with the child's oriented haplotype from that parent at informative sites; the matching haplotype is `transmitted`. Segments are emitted per contiguous run with typed boundaries: `BLOCK_EDGE` (parent block ends; no crossover information) or `CROSSOVER` (switch inside a parent block) — called only with ≥ `xo_min_sites` (10) informative sites on each side, else flagged as a possible phase-switch error with lowered `confidence`. Crossovers per parent per meiosis is a QC gate (expect roughly 25–45 autosomal; **unverified** exact figure). Crossover resolution is recorded; a candidate inside it gets `NEAR_CROSSOVER`. Note the child's `XP` label does *not* change at a crossover — only the parent's transmitted haplotype does.

### P4 — Tags, and why M1 does not rewrite BAMs
HiPhase already emits `read_name → haplotag` per block for every read (`phase_haplotags/`). M1 therefore emits, per sample, a **block-level** table — `phase_block_id → {XP: P|M|U}` for the child, `phase_block_id → {XT: T|U|N, segment}` for parents — and M2 resolves any read's parent-of-origin / transmission label by `read_name → phase_block_id → label`. No BAM copy, no 105-BAM pass. `phase-dnm haplotag --export-bam` writes `XP:A`/`XT:A` (SAM local-use range) into a new BAM for IGV only. `HP`/`PS` are never overwritten.

### P5 — Candidates enter unfiltered; the mask is a flag, not a filter
Small variants: every site in the **family joint VCF** where a child is het and both parents hom-ref, any GQ (~25 k per child, 0.2). SVs: every sawfish record where the child carries an allele absent from both parents' GTs. TRs: every TRGT locus where a child allele exceeds both parental alleles by ≥ 1 motif unit (no p99+3 margin; that is tiering, not calling). Every existing candidate list is also ingested and cross-referenced (`source_list`), so rescue/demotion can be reported against the current pipeline (P18). `source_tier` and `mask_overlap` are columns. Rationale: the point of the module is to rescue and reject, and the 338 Mb mask (METHODS §2) is where phase evidence is both most valuable and least reliable — it must be *measured* there. The mask remains a filter for burden analyses.

### P6 — Sib-shared variants are kept
Report §3.5 excludes calls shared by both siblings; here they are kept and flagged `SIB_SHARED` — a variant in both sibs and absent from both parents' genotypes is the parental-mosaic signature. Two quads: illustration, not a truth set (R2).

### P7 — The six-haplotype matrix
Rows `M1 M2 F1 F2 CM CP`. A haplotype is *observed* when `dp ≥ k`; `hap_obs_k` counts over six; k=3 and k=5 are always both computed. At ~11× per haplotype, k=5 is the working threshold in unique sequence; k=3 is a fallback that carries `LOW_HAP_DEPTH`. Untagged reads count in `untagged_dp/alt` and in child-level allele support (as the existing DNM logic does) but never in a haplotype row. A read whose `HP` disagrees with `local_rephase_agreement` moves to `amb` and the candidate gets `HP_REPHASE_CONFLICT`. TR rows carry per-read allele-length distributions (from `AL:i` in the spanning-read BAM); "alt" for TR is the expanded child allele (`class_payload.expanded_allele_idx`).

### P8 — `phase_class`: rule layer (transparent)
Evaluated in order; first match wins. `T` transmitted parental haplotype, `U` untransmitted, `A` the child's alt-carrying haplotype, `O` the other; `f_A = alt/(alt+ref)` on A; `e` = per-read error allowance (`thresholds.yaml`: max(1 read, 5 %)).

| class | rule |
|---|---|
| `phase_conflict_artifact` | alt > e on **both** child haplotypes, **or** alt > e on ≥ 2 parental haplotypes with no parental het GT |
| `inherited_missed_in_parent` | T observed and alt_T ≥ max(3, 0.3·dp_T) |
| `parental_mosaic_transmitted` | T observed, e < alt_T < 0.3·dp_T, `dp_T ≥ mosaic_min_dp` (P10); else `inconclusive` + `PARENTAL_ALT_LOW_DEPTH` |
| `child_postzygotic_mosaic` | alt confined to A, `f_A ≤ 0.7`, `dp_A ≥ mosaic_min_dp`, alt_O ≤ e, all four parental haplotypes observed with alt ≤ e |
| `germline_DNM_phased` | alt confined to A, `f_A ≥ 0.8`, alt_O ≤ e, T observed with alt_T ≤ e, U observed |
| `germline_DNM_unphased` | alt confined to A (or child unphased with AB 0.25–0.75), parental alt ≤ e on every *observed* haplotype, but `hap_obs_k5 < 6` or T unresolved |
| `inconclusive` | everything else, with a reason flag |

`rule_score` = number of satisfied germline criteria (0–6). All thresholds live in `thresholds.yaml` with a version string written into every output row.

### P9 — Likelihood layer
Per-haplotype alt/ref counts are binomial with hypothesis-specific expected alt fractions, H ∈ {germline DNM on A; inherited from T; parental mosaic on T at fraction m; child postzygotic on A at fraction c; artefact at rate ε on all six}. ε per class from cohort hom-ref sites in the same mask/mappability stratum; m, c integrated over Uniform(0.02, 0.5). Priors from expected rates vs. the measured candidate counts (0.2). Output: posterior per hypothesis; `lik_post_germline` is `phase_score`. Both layers always emitted.

### P10 — Mosaic evidence floor (overclaim guard)
A mosaic label requires ≥ `mosaic_min_dp` (15) reads on the relevant haplotype **and** ≥ 3 minority-allele reads **and** those reads pass MAPQ ≥ 20, not soft-clipped at the site, NM within the haplotype's ref-read distribution. Otherwise `inconclusive` + `MOSAIC_UNDERPOWERED`. At ~11× per haplotype, f_A 1.0 vs 0.8 is two reads. Mosaics from this cohort are reported as *candidates with read counts*, never as a rate, unless confirmed by targeted deep sequencing (R1).

### P11 — `rf_safe` policy and positive-generation schemes
**P11a (classic swap, primary).** A synthetic trio has real surrogate parents with real haplotagged reads; what is missing is only the transmitted/untransmitted label. A feature is `rf_safe: true` iff computable from (child evidence) ∪ (parent evidence symmetric under swapping the parent's two haplotypes and under swapping the two parents). `features/registry.py` refuses to build an RF matrix with `rf_safe: false` columns; the check hash goes into the manifest. Transmission-dependent features live in P8/P9 only, applied after the RF.
A second asymmetry for Methods: in positives the child variant is truly present (inherited); in negatives most child calls are artefacts. The RF learns largely "is the child genotype real"; only P8/P9 test "is it absent from the transmitted haplotype". Hence the combined final rule (P15).
**P11b (haplotype swap, ablation; Q13).** For a variant inherited on T, delete the alt-bearing reads on T in the parent (downsampling U to keep the T:U ratio); the transmission map survives. Guard: a *chimera test* — a classifier on read-level features alone must not distinguish edited from unedited parental windows; if it does, P11b is reported as leaky.

### P12 — Folds are closed under the swap; pairings are drawn inside folds
Because `swap_pedigree.py` exchanges offspring between paired families, the unit that must not cross a fold boundary is the **pair**. Adopted: fold over families first (all members of a quad together), then draw pairings *within* each fold. With 33 complete-trio children and 5 outer folds: 6–7 families → 3 pairs per fold (+ one unpaired family contributing negatives only, or paired in a second within-fold seed). Extra seeds only within-fold. Inner folds built identically from outer-training families. Stratify by variant class (separate models per class by default; pooled model with class feature as ablation) and by DNA source (saliva vs blood; the REACH quad is the only blood family, so it must not be alone in a fold).

### P13 — The negative set is the raw putative set (measured)
From 0.2: raw negatives ≈ 25 k/child, true DNMs ≈ 70 → ~0.3 % contamination; the SynthDNM assumption holds on HiFi **only** for the raw set. Negatives are therefore the P5 unfiltered candidates. `cv_report.json` reports `negative_contamination_estimate` per class. Iterative relabelling (dropping negatives that a fold-external model calls DNM) runs as an ablation and is labelled partially circular wherever it appears.

### P14 — Metrics, calibration, baselines
Outer folds: PR-AUC, ROC-AUC, Brier, reliability curves; isotonic calibration fitted on inner out-of-sample predictions. Repeated nested CV (≥ 5 outer-fold seeds) — 33 trios give ~6 families per fold; report mean ± range. Baselines on identical folds: XGBoost without phase features; phase features only; RF; logistic regression; SynthDNM *universal* 21-feature model as published. Permutation importance on outer held-out folds, grouped by feature family. Headline performance on external sets (WES-confirmed exonic DNMs, GraphTyper-confirmed exonic SVs, GIAB), never on synthetic held-out labels alone.

### P15 — Final call: RF + phase layer, thresholds from outer folds
`dnm_call = YES` iff `rf_prob ≥ τ_class`, **or** (`rf_prob ≥ τ_rescue,class` and `phase_class == germline_DNM_phased` and `hap_obs_k5 == 6`) — rescue. `dnm_call = NO` regardless of `rf_prob` when `phase_class ∈ {phase_conflict_artifact, inherited_missed_in_parent}` — demotion. Mosaic classes are never `YES` and are reported separately. τ values chosen on outer-fold predictions at a target FDR from the external sets, recorded with the fold seed. Rescued/demoted counts per class go in the cohort summary.

### P16 — GIAB is held out of everything
HG002/HG003/HG004 HiFi (Q6) never enter folds, pairing, ε estimation, calibration or thresholds. It is the only fully external test; DeepTrio-only runs there.

### P17 — Adapter inputs come from the WDL outputs that already exist
- **SNV/indel**: pysam pileup at the site in the three haplotagged BAMs; `HP` from the read, cross-checked against `phase_haplotags`.
- **SV**: alt-supporting read names from `sv_supporting_reads.json.gz` (sawfish `--report-supporting-reads`), joined to `phase_haplotags` by read name for haplotype; *plus* the `sv_read_review.sh` junction scan around both breakpoints in **all three** BAMs, because sawfish's per-sample list is at sawfish's own sensitivity and a hom-ref parent's 1–2 junction reads (the parental-mosaic evidence) must be counted independently. Per-haplotype depth inside vs flanks from `hapdepth`.
- **TR**: per-read `AL:i` from `trgt_spanning_reads/*.bam` (already `HP`-tagged), so the per-haplotype allele-length distribution is a direct read of tags; no re-genotyping.
Rationale: three adapters, one join key (`read_name`), no new callers.

### P18 — Concordance with the original pipeline is an M3 deliverable
For every class, a table of the original pipeline's de novo set (the filtered `denovo_tiered.tsv` — 1,332 calls, median 38/proband — plus the per-family hiconf lists, the prioritised SV table, `tr_denovo_expansions.tsv`; Q15) against the module's final `dnm_call`, with: concordant YES, original-only (demoted: by `phase_class`), module-only (rescued: by `source_tier`/`mask_overlap`), and per-proband rate before/after. This is both the requested comparison and the R9 guard: a rescue that moves the per-proband rate from ~38 towards ~70 with the same parent-of-origin ratio is credible; one that overshoots is not.

### P19 — Duos
The two mother–child duos have no paternal reads; `F1/F2` rows are unobservable, paternal transmission is undefined. Default (Q14): excluded from M1 transmission, from M4 folds and from cohort rates; optionally run in M2 as half-trios with `poo = undetermined:NO_FATHER` and `hap_obs ≤ 4`, clearly separated in every table.

---

## 3. SynthDNM one-pager: what the code does, and what it implies for phase features

**Construction.** `swap_pedigree.py` shuffles complete families, pairs them, and writes a PED in which family A's parents get family B's offspring and vice versa. `extract_dnm_features.py` extracts FORMAT-level trio features at sites that look de novo *under the swapped pedigree* (child het, surrogate parents hom-ref; an `AC=2` condition appears in the extractor) → `truth=1`; and at putative DNM sites under the *real* pedigree → `truth=0`. `preprocess_features.py` derives `child_AR/min_AR/max_AR` (ref/(alt+1)), min/max over parents of `GQ, DP, PL0–2`, `child_AB`, `indel_flag`, sex-aware `haploid_flag` (PSAM + PAR BED), and samples `--sample_n` rows. `train.py` fits XGBoost (1000 trees, depth 6, lr 0.1, early stopping) on a stratified row-level split; `classify.py` emits `synthdnm_prob`. Published: ~96 % recall of denovo-db SSC DNMs.

**Implications.**
1. Positives use a real child and real surrogate parents: any feature that is a function of (child) or (a parent on its own, haplotypes unlabelled) exists for positives. Anything needing the parent–child relationship does not → `rf_safe: false`.
2. Pairing is an exchange → folds must contain whole pairs (P12). The current code has no grouping.
3. Negatives inherit the caller's label noise; on HiFi this is tolerable only for the raw set (0.2, P13).
4. Positives are inherited het variants at ~50 % VAF, clonal on one haplotype; mosaics are out of distribution (P10, P15).
5. The `universal` set is FORMAT-only and fully available in the GLnexus BCF — the natural baseline and the first 21 rows of `features.yaml`. `gatk`/`ssc` INFO features do not exist here.

---

## 4. Overclaim register

| R# | risk | what settles it |
|---|---|---|
| R1 | child or parental mosaics at ~11× per haplotype | P10 floors; amplicon ≥ 1000× or ddPCR on a random subset of mosaic calls **and** of `germline_DNM_phased` controls; report confirmation rate, not a mosaic rate |
| R2 | "multi-sib shared DNMs validate the parental-mosaic class" from 2 quads | "illustrated in two families"; SPARK WES multi-sib families give a population shared-DNM rate but on another platform |
| R3 | per-class thresholds for SV and TR from ~33 trios | true de novo SVs cohort-wide are single digits; SV/TR thresholds provisional; SV evaluated by GraphTyper-confirmed exonic events and spike-ins; TR has no orthogonal comparator at all |
| R4 | parent-of-origin ratio "~75–80 % paternal" as validation | sanity check on the phased subset only; binomial CI; validates the cohort, not individual calls |
| R5 | paternal-age slope from 33 trios | expected ~1–2 DNMs/yr (**unverified**); wide CI; report, do not lean on |
| R6 | switch/flip error from trio consistency | lower bound (informative sites only); compare with HiPhase's published rate |
| R7 | RF performance on synthetic held-out folds | ceiling effect; headline numbers from external sets (P14) |
| R8 | read-level haplotype swap as leak-free | chimera test (P11b) |
| R9 | rescue inflating DNM counts | P18 table; rate before/after vs expected ~60–80 (**unverified**); rescued calls must show the same parent-of-origin ratio |
| R10 | WES concordance as "validation" of the small-variant set | exonic only (~1–2 % of DNMs); state the denominator; non-detection in WES is never refutation (same rule as METHODS §5) |

---

## 5. Paper framing (P20) — one pipeline paper, cohort as proof of concept (Q12)

The module is reported as the de novo calling section of the SPARK lrWGS pipeline paper, not as
a method paper. What that fixes:
- **Claims** are pipeline-level: "phase-aware review changes the de novo call set in these
  measurable ways on this cohort" (P18 table: rescued, demoted, mosaic-flagged, per class, with
  parent-of-origin ratio and per-proband rate before/after). No claim of general superiority over
  other long-read DNM callers is made or needed.
- **Validation scope** is what supports those claims: WES concordance for exonic SNV/indel (R10),
  GraphTyper on WES CRAMs for exonic SVs, spike-ins for all three classes (the only TR
  validation), family segregation, parent-of-origin and rate sanity checks (R4, R5, R9), and GIAB
  as the reproducibility anchor once pulled (Q6). The DeepTrio comparison is kept only on GIAB.
- **Methods text** is assembled from P1–P19 exactly as `docs/METHODS.md` is assembled from D#;
  each P# is one paragraph.
- **What is deliberately left for later**: mosaic *rates*, TR de novo *spectrum*, paternal-age
  *slope* as findings. They appear as descriptive tables with the overclaim guards, not as results
  the paper rests on.
