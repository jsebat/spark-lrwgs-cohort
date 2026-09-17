# 06_clinical — design decisions

The governing decisions are recorded in `05_denovo/DESIGN.md` P31 (the targeted arm is a pipeline stage; MECP2 is an
acceptance case, not validation), P32 (one classifier; large SVs go to the targeted arm; two verdicts recorded
independently) and P33 (two arms in opposite order). This file records what is specific to the implementation.

### C1 — Inputs are 05_denovo's tables, all candidates, not a re-review of BAMs (2026-09-16)

The original discovery (`legacy/`) collected de novo SV calls from the WDL's per-family tables, intersected CDS, then
ran a read-level review script by hand for the survivors. The de novo module now measures for EVERY candidate what
that script measured for a handful — per-haplotype depth inside versus flanks for all six haplotypes, junction reads
by haplotype at both breakpoints, heterozygous-SNV persistence, parental support — and writes it to the evidence
tables. The clinical arm therefore opens no BAM. The raw six-haplotype SV columns live in
`EVIDENCE_DIR/<family>/evidence/<child>.sv.evidence.tsv` (the final table keeps only derived features), joined on
`variant_id`.

### C2 — The gene panel is data; identifiers are data (2026-09-16)

SFARI and DDG2P gene lists are read from files; nothing is hard-coded. The acceptance expectations name genes and
coordinates only; the job derives the carrying sample from the tables at run time, so the public repository never
contains a sample identifier.

### C3 — Consequence annotation reuses 03_tiering's VEP output and is restricted to panel-gene spans (2026-09-16)

The de novo module scores ~34,000 SNV/indel candidates per child. Arm B needs a consequence only for candidates inside
a panel gene, so candidates are first intersected with panel-gene spans (first to last exon) and only those are
fetched from the per-chromosome VEP+LOFTEE VCFs. Candidates absent from those VCFs are written to
`smallvar_annotation.tsv.missing` rather than silently dropped. The VEP resources used by 03_tiering
(`g2mh/scripts_for_rare_pipeline/...`) were found to be absent at their recorded paths on 2026-09-16, so the module
consumes the annotated VCFs rather than re-running VEP; re-annotation is a 03_tiering concern.

### C4 — The quality rubric is deterministic and explains itself (2026-09-16)

Every verdict carries `+`/`-` criteria and named rejections, so a reviewer disagrees with a rule, not a label. The
verdict is a function of the row alone, which makes it reproducible without the cluster and testable
(`tests/test_quality.py`). Its constants are the P8 deletion-path constants; changing them is a 05_denovo decision.

### C5 — REVIEW is a first-class outcome (2026-09-16)

The DNMT3A case is the reason: a 302 bp coding deletion with depth halved and the alt phased to one haplotype, but a
read cannot present a junction at both ends of an event shorter than itself in the way the rule expects, and there is
no heterozygosity measurement below 5 kb. That is not a FAIL and not a PASS; it is the case where an analyst looks.
The module says so rather than forcing a binary.

### C6 — Arm A keeps its FAILs (2026-09-16)

Arm A's value is that its candidates were chosen by size and rarity alone, so its quality verdicts are an unbiased
sample of what the caller saw at ≥ 5 kb. Inherited and artefact phase classes are therefore kept in the table and
FAIL there, rather than being filtered out upstream, so the FAIL fraction is itself a measurement.

### C7 — Missense tiers are optional and say when they are absent (2026-09-16)

James's four-rankscore tiers need dbNSFP parquet; the module's Python environment lacks `duckdb`, and `pyarrow` is used
when present. Without it the impact is still `missense_variant` and the tier column is blank, which the README states.
