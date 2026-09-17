# 06_clinical — the targeted clinical arm

The genome-wide de novo caller (`05_denovo`) answers one question with one set of rules for every variant: what is
the de novo mutation rate, and which calls survive a uniform threshold. This module answers a different question —
**is there a reportable variant in this child** — and is allowed to be biased towards clinical yield, because nothing
it produces ever enters a recall, FDR or rate estimate (05_denovo DESIGN **P31**). Every row it writes carries
`discovery_mode = targeted_clinical`.

It is a stage of the pipeline, documented as one, because this is how the cohort's two known pathogenic events were
actually found: a 35 kb *MECP2* deletion the genome-wide caller does call (tier 1, depth rule 6/6), and a 302 bp
*DNMT3A* coding deletion it leaves **below threshold** (rule 5/6, rf 0.72, `BELOW_TAU`). A workflow restricted to
called variants would not see the second one. The scripts that made the original discovery are preserved unchanged in
[`legacy/`](legacy/) as provenance; this module re-implements the same logic over `05_denovo`'s evidence tables, which
already contain every measurement the original read-level review made by hand.

## Two arms, opposite order (P32, P33)

| | arm A — rare + large | arm B — clinically led |
|---|---|---|
| what selects a candidate | **size and rarity**: SV ≥ 5 kb, cohort AC 0 (leave-one-family-out), founder-panel recurrence 0 | **gene panel and impact**: exon/UTR of a SFARI ∪ DDG2P gene; LoF / missense / splice for SNV-indel; STRchive or panel exon for TR |
| which candidates | every SV `05_denovo` scored, any `dnm_call` | every candidate, any `dnm_call` — including below threshold |
| order | **quality first**, then clinical annotation | **clinical relevance first**, then quality |
| what the quality verdicts may be used for | validation truth for the evidence (an unbiased sample) | reporting that variant only — never as a truth set (conditioned on being interesting) |

Both arms record **two verdicts independently**: `quality_verdict` (PASS / REVIEW / FAIL, with reasons and a
criteria count) and `clinical_relevance` (HIGH / MODERATE / LOW / NONE with the reason). A quality-PASS,
clinically-NONE row is a complete result, not a null one.

## Gene criterion: panel membership, not constraint

The gene criterion is SFARI ∪ DDG2P membership. `s_het` is reported on every row and **never gates**. Constraint
metrics fail systematically for clonal-haematopoiesis driver genes: *DNMT3A* scores s_het 0.0063 despite being a
haploinsufficiency gene (germline LoF → Tatton-Brown-Rahman syndrome), because somatic CHIP variants in blood-derived
population exomes inflate its apparent LoF tolerance — gnomAD flags its pLI unreliable for exactly this reason, as for
*TET2*, *ASXL1*, *PPM1D*. A constraint filter would discard one of the two pathogenic events this cohort contains,
and would do so for a whole class of genes.

UTRs count. A deletion of a 5′ UTR or of a last exon's 3′ UTR can abolish expression without touching a coding base;
a CDS-only rule cannot see it. The SV impact rule is the original one (DEL over CDS → coding LoF; INS breakpoint in
CDS → coding LoF; DUP over CDS → coding, flagged as gain; INV breakpoint in gene → coding; BND excluded), extended
with a labelled UTR class.

## Quality rubric (`src/clinical/quality.py`)

Decided entirely from what `05_denovo` measured; no BAM is opened here.

* **Any class → FAIL** when the phase class says not de novo (`phase_conflict_artifact`, `inherited_missed_in_parent`,
  `parental_mosaic_transmitted`), when the variant is seen in another family or in the founder panel, or when it
  fails the site filter. Mask overlap ≥ 0.5 of an SV interval fails; a touch is noted.
* **DEL / DUP with an interval** (the P8 deletion path): child depth halved (≤ 0.7) or gained (≥ 1.3); neither parent
  depleted (≥ 0.85) and both measurable; junction reads at both breakpoints; ≥ 2 child junction reads; loss of
  heterozygosity inside (≥ 5 kb only). PASS = depth + parents + (junctions or LOH). REVIEW = depth without
  junctions, or junctions without depth. FAIL = no depth change and heterozygous sites persist. A parent with junction
  reads or lost depth is FAIL (inherited).
* **INS / INV / short SVs**: ≥ 3 child supporting reads, no parental support, alt confined to one child haplotype,
  rule score ≥ 5, germline phase class → PASS; rule ≥ 3 with clean parents, or rf ≥ 0.3 → REVIEW.
* **SNV / indel**: ≥ 3 alt reads, confined to one haplotype, parents alt-free (max 1 read on any parental haplotype),
  rule ≥ 5, germline phase class, child depth ≥ 8 → PASS; gnomAD AF ≥ 0.001 or parental alt reads → FAIL;
  postzygotic mosaic → REVIEW (clinically relevant if real).
* **TR**: called tier 1 with a germline phase class → PASS; candidate or rule ≥ 3 → REVIEW.

The deletion-path constants are the ones P31 records as having been chosen with *MECP2* and *DNMT3A* in view. That is
why those two events are this module's **acceptance test** and are not evidence that the rubric works.

## Running it

```
sbatch workflow/run_clinical.sb          # annotate-smallvar → run → acceptance, on the cohort's 05_denovo outputs
```

The job sources `05_denovo/config/phase_dnm.env` for `FINAL_DIR`, `EVIDENCE_DIR`, `PHASE_DNM_PYTHON`, `LRC`, and takes
the clinical resources from the environment (`GFF3`, `GENE_SETS`, `SHET`, `VEP_GLOB`, `DBNSFP_PARQUET`,
`CLINICAL_OUT`) with the cohort's paths as defaults. Consequences for SNV/indel candidates come from the per-chromosome
VEP + LOFTEE VCFs written by `03_tiering`; only candidates inside a panel-gene span are looked up. Missense tiers
(ClinPred / AlphaMissense / popEVE / MPC rankscore thresholds, `legacy/miss_tier.py`) attach when dbNSFP parquet can be
read and are otherwise left blank with a note.

Outputs in `$CLINICAL_OUT`:

| file | content |
|---|---|
| `arm_A.rare_large_sv.tsv` | one row per rare ≥ 5 kb SV, quality verdict, then gene / impact / panel / s_het |
| `arm_B.clinically_led.tsv` | one row per (candidate, panel gene): SV, SNV/indel, TR; clinical relevance and quality |
| `clinical_report.md` | per child: PASS/REVIEW rows with a clinical relevance above NONE, clinically led rows first |
| `smallvar_annotation.tsv` (+ `.missing`) | the consequences looked up, and the panel-span candidates the VEP VCFs did not contain |
| `acceptance_expect.tsv` | generated by the job from coordinates — sample identifiers never enter this repository |
| `summary.json` | counts per arm and verdict |

### Second VEP pass (C1)
`run_clinical.sb` pass 1 writes `smallvar_annotation.tsv.missing`: panel-span SNV/indel candidates the freeze-1 VEP VCFs do not
carry (raw candidates the freeze excluded; 46 % in the first run). `sbatch workflow/vep_missing.sb` annotates them with the
identical VEP 115 + LOFTEE + dbNSFP invocation of `03_tiering/lr_vep.sb` into `$CLINICAL_OUT/vep_missing/<chrom>.vep.vcf.gz`, and
a second `run_clinical.sb` pass consults that source for whatever the first lacked (`--vep-glob` may be repeated).

## Acceptance test

`clinical acceptance` requires arm B to contain the *MECP2* deletion (chrX:154062405–154097859, 35,454 bp) at quality
**PASS** and the *DNMT3A* deletion (chr2:25244326–25244628, 302 bp) at **REVIEW or better**, each in the child that
carries it. The job derives the sample identifiers from `05_denovo`'s tables at run time. Expected reasons: *MECP2*
passes on depth 0.59, LOH 0.0, both parents intact; *DNMT3A* is depth 0.38 and phased but has no junction read at both
ends of a 302 bp event and no heterozygosity measurement below 5 kb, which is exactly the REVIEW case.

## Runs on the cohort (2026-09-16 blind run; 2026-09-17 after the pipeline-wide fixes and the second VEP pass)

Blind — no gene named in advance — the module's own shortlist (quality PASS or REVIEW, relevance HIGH or MODERATE) is
**108 rows** over 33 children on the retrained call set (harness_v7). Its only HIGH + PASS rows are the 35 kb *MECP2*
deletion (PASS 6/6, tier 1), the 302 bp *DNMT3A* deletion (PASS 5/5, **tier 1** — below the genome-wide threshold in the
first run, tier 1 after the multiallelic feature was dropped and the SV model refit), a frameshift `LoF_HC` indel in
*DLX6* (below threshold) and a start-lost SNV in *ETFA* (tier 1). Thirteen further `LoF_HC` variants in panel genes sit at
REVIEW (12 indels, 1 SNV; 9 of them below the genome-wide threshold) — the list this arm exists to hand an analyst. The
rest of the shortlist is TR REVIEW rows. Acceptance: PASS for both known events in both runs.

| | rows | PASS | REVIEW | FAIL |
|---|---|---|---|---|
| arm A (rare + large SV) | 68 | 1 | 1 | 66 |
| arm B (clinically led): SV 121, SNV 509, indel 112, TR 5,807 | 6,549 | 7 | 147 | 6,395 |

The first run had a known gap: 43,053 of the 93,432 SNV/indel candidates inside a panel-gene span were absent from the
freeze-1 VEP VCFs (raw candidates the freeze excluded) and carried no consequence. `workflow/vep_missing.sb` (see *Second
VEP pass* above) re-annotated all of them with the 03_tiering invocation — 65 `LoF_HC` among them — and the second pass
annotates **93,432 of 93,432** (0 missing; missense tiers attached for 428 via dbNSFP parquet). No HIGH + PASS row came
from the re-annotated set: the gap held low-quality candidates, as suspected, but that is now measured rather than assumed.

## What this module must never be used for

Any performance number, any threshold fit, any claim about the genome-wide classifier. Arm A's quality verdicts may
validate the *evidence* (do the depth and junction criteria agree with a careful reader?) but may not re-fit the
selector that chose them. Arm B's verdicts describe the panel, not the caller. GIAB stays held out of both.
