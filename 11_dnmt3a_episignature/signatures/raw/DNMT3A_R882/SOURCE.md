# DNMT3A_R882 — raw sources (DNMT3A R882 dominant-negative, DNMT3A Overgrowth Syndrome)

Download date: 2026-09-09. Files are unmodified copies of the publisher ESM files.

## Smith AM et al. 2021, Nature Communications 12:4549
"Functional and epigenetic phenotypes of humans and mice with DNMT3A Overgrowth Syndrome"
DOI 10.1038/s41467-021-24800-7 — PMID 34315901 — PMC8316576

Design (human WGBS, peripheral blood): 3 R882 patients (UPN 624400 R882H, 154605 R882H, 894912 R882C in AML remission)
vs 15 healthy donors (13 unrelated normal donors, columns `TWGB.NM*` / `TWGB.ND*`, plus 2 unaffected siblings
`H_KA.867535...sibling`, `H_KA.978897...sib`). The 8 non-R882 patients are carried passively in the same table
(columns with a UPN in {411168, 511909, 228211, 295041, 930075, 723972, 786396, 518693} or "TBRS" in the name).
Mapping: biscuit (MGI analysis-workflows bisulfite CWL v1.5.0). DMRs: metilene; >10 CpGs, mean group difference >0.2, FDR-filtered.
Raw human sequence data: controlled access, dbGaP phs000159.

### Smith2021_NatCommun_41467_2021_24800_MOESM4_ESM.xlsx — Supplementary Data 1  **PRIMARY**
- Download URL: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-021-24800-7/MediaObjects/41467_2021_24800_MOESM4_ESM.xlsx
  (847,113 bytes; identical file also in the Europe PMC bundle https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8316576/supplementaryFiles)
- Publisher description (MOESM3): "Differentially methylated regions (DMRs) identified in peripheral blood using WGBS data
  from DNMT3A Overgrowth Syndrome (DOS) patients with R882. Genomic coordinates in columns A-D, functional regions F-L,
  mean methylation values per sample N-AM. Size of DMR in bp."
- Sheets: `Figures` (Fig 1c, 1d, 1e, 1f, 1g, Supplemental 2h), `dmrs.annot.summarized`
- Rows: 2,209 DMRs; chromosomes chr1-22 and chrX
- Columns: dmr.id (chr_start_end), chromosome, start, end, closest_gene, islands, shores, shelves, genes, enhancers,
  promoters, tss (0/1 flags), size (bp; median 546, range 42-4,683), then 26 per-sample mean-methylation columns
  (15 control, 3 R882, 8 non-R882).
- Genome build: hg38 (not written in the sheet; inferred from the paper quoting the UPN 518693 deletion in hg38, the
  GRCh38-based biscuit pipeline, and the first DMR chr1:1,034,563-1,034,983 being annotated to AGRN, which only holds in hg38).
  0- vs 1-based convention not stated; check against CpG positions before converting to BED.
- Direction / effect size: no explicit column. Derived here as mean(R882) minus mean(controls): negative for all 2,209 DMRs
  (median -0.378, range -0.896..-0.191), i.e. uniformly hypomethylated.
- Usable for a DMR list: YES (hg38 regions + derived delta beta; no lift-over needed).

### Smith2021_NatCommun_41467_2021_24800_MOESM3_ESM.pdf — Description of Additional Supplementary Files
- https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-021-24800-7/MediaObjects/41467_2021_24800_MOESM3_ESM.pdf
- Text source for the column descriptions above.

Related: Supplementary Data 2 (non-R882 DMRs, 332 regions) is the TBRS_LOF companion table and lives in ../TBRS_LOF/.
No separate R882 table exists in Jeffries 2019, Aref-Eshghi 2020, Levy 2022 or Heyn 2019.
