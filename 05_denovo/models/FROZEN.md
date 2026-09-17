# Frozen classifiers (P21) — release 2026-09-17

Three pre-trained XGBoost classifiers, one per variant class group, refit on the whole SPARK long-read cohort (33 complete-trio
families, 35 children, 105 HiFi genomes at ~22–24×) with the hyperparameters chosen by swap-closed, family-grouped nested
cross-validation over five outer-fold seeds (DESIGN P12–P14, P23–P24). Training run `harness_v7`, 2026-09-17, on the
re-reviewed evidence tables (thresholds 0.2.0 observability; read-quality block present on both sides).

| class group | file | sha256 (first 16) | features | training rows | positives (synthetic) |
|---|---|---|---|---|---|
| SNV/indel | `snv_indel.xgb.json` | `efc3d3c0d51dc227` | 52 | 749,095 | 50,210 |
| SV | `sv.xgb.json` | `8e5635e70f4a51c9` | 40 | 52,084 | 29,650 |
| TR | `tr.xgb.json` | `c9bfeedcab78b3d6` | 38 | 696,654 | 9,498 |

Full checksums, feature lists, hyperparameters (`max_depth 6, learning_rate 0.05, n_estimators 600, subsample 0.8,
colsample_bytree 0.8, min_child_weight 10`), XGBoost 3.2.0 and the feature-registry checksum are in each
`<class>.training_manifest.json`. Positives are synthetic de novo variants from within-fold pedigree swaps; negatives are the
raw, unfiltered candidate sets of the real trios (SynthDNM construction). No population-frequency, cohort-recurrence or
parent-of-origin quantity is a feature (`rf_safe: false` in `config/features.yaml`, P24): those enter as rules after the score.

## Operating points (thresholds.yaml 0.3.0, DESIGN P15 amendment 2026-09-14)
Decision on `rf_q` = 1 − pass rate among the cohort's candidates of the same variant class. Tier 1 (`dnm_call YES`) and tier 2
(`CANDIDATE`, for validation), both behind the six-haplotype germline review gate and the rule layer (gnomAD AF < 0.001 or
long-read SV catalogue AF < 0.001; leave-one-family-out founder-panel recurrence 0; leave-one-family-out cohort allele count 0;
segdup/simple-repeat mask — a filter for SNV/indel and TR, a flag for SV):

| class | tier 1 τ_q | tier 2 τ_q | cohort read-out (per proband, tier 1) |
|---|---|---|---|
| SNV/indel | 0.99 | 0.95 | 85 (73 SNV + 12 indel; 2,964 cohort-wide, 35 children), paternal fraction 0.77 (n 2,478 determined) |
| SV | 0.99 | 0.97 | 0.9 (30 cohort-wide), paternal fraction 0.52 (n 25) |
| TR | 0.999 | 0.997 | 1.1 (37 cohort-wide), paternal fraction 0.93 (n 29) |

## Using the frozen models on another cohort (transfer path)
```
phase-dnm score --class-group snv_indel --model models/snv_indel.xgb.json --manifest MANIFEST.tsv \
                --evidence-dir $EVIDENCE_DIR --out-dir $TRAIN_DIR/frozen          # rf_probs/ + tau.snv_indel.json
RF_PROBS_DIR=$TRAIN_DIR/frozen/rf_probs TAU_DIR=$TRAIN_DIR/frozen bash workflow/m3_integrate_family.sb <FAMILY>
```
`score` verifies the model checksum against its manifest, scores every child's `features.rf.tsv` (Module 2 output), and
computes `rf_q` against the scored cohort's own candidates per variant class, so the operating points above keep their
pass-rate meaning. **Transfer caveat (P21):** valid for the same callers and versions (DeepVariant + GLnexus, sawfish 2.2.1,
TRGT 5.0.0, HiPhase 1.6.0; HiFi-human-WGS-WDL v3.3.1) at similar coverage. With different callers or depth, retrain:

## Retraining recipe (the workflow is the deliverable, DESIGN P21)
```
phase-dnm swap       --manifest MANIFEST.tsv --out-dir $TRAIN_DIR/folds --n-folds 5 --seeds 0,1,2,3,4      # swap-closed family folds
workflow/m4_synthetic_trio.sb + m4_annotate_synthetic.sb   (every synthetic trio through the same M2 chain)
workflow/m4_train.sb <class> 0,1,2,3,4   (OUT_DIR=…)        # nested CV, arms, ablations, frozen refit -> models/
phase-dnm rescore    (per-variant-class rf_q from the fold models)  ->  m3_integrate_family.sb  ->  phase-dnm concordance
```
or `snakemake m4 --profile workflow/profiles/expanse`.
