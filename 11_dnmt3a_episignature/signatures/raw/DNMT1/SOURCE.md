# DNMT1 — ADCA-DN / HSAN1E episignature sources

Download date: 2026-09-09.

## Source A (primary probe list): EpiSign ADCADN — Aref-Eshghi et al. 2020 AJHG, Table S2/S3

File lives in `../EPISIGN_SHARED/ArefEshghi2020_AJHG_mmc2.xlsx` (see `../EPISIGN_SHARED/SOURCE.md`).
DOI 10.1016/j.ajhg.2020.01.019. Sheet `S2- Probes`, column `ADCADN == True`: **104 probes**,
with `Chr` and `Position (hg19)`. Sheet `S3-Methylation levels` gives mean beta for `ADCADN` and
`Control` on the same 3643 probes, so delta-beta = `ADCADN - Control` and direction are derivable.
Genome build hg19. Usable: yes.

## Source B (region/DMR list): Kernohan KD et al. 2016, Clin Epigenetics 8:91

"Identification of a methylation profile for DNMT1-associated autosomal dominant cerebellar ataxia,
deafness, and narcolepsy". DOI 10.1186/s13148-016-0254-x. PMID 27602171. PMCID PMC5011850.

| Item | Value |
|---|---|
| Table | Additional file 1 = "Supplementary Table 1: Significant regions detected by methylation array in ADCA-DN versus control individuals" |
| Download URL | https://static-content.springer.com/esm/art%3A10.1186%2Fs13148-016-0254-x/MediaObjects/13148_2016_254_MOESM1_ESM.pdf |
| Local file | `Kernohan2016_ClinEpigenetics_AdditionalFile1_MOESM1.pdf` (178 KB, 4 pages) |
| Format | PDF table only (no xlsx/csv offered; this is the only additional file of the paper) |
| Rows | 82 regions (regex `^chr\w+\s+\d+\s+\d+\s+[\d.]+\s+\d+\s+[\d.-]+` over pypdf text recovers all 82) |
| Columns | Chr, Region start, Region end, Meth. average (mean beta in ADCA-DN), # probes (5–18 consecutive 450K probes), Estimate (effect size, ADCA-DN minus control mean beta; all positive, 0.20–0.34), Nearest gene (with strand), Distance to nearest gene, Distance to nearest CpG island |
| Probe IDs | **no** (region-level only; individual cg IDs not listed) |
| Genome build | hg19 (stated in main-text Methods; array = Illumina Infinium HumanMethylation450). Not restated inside the PDF. |
| Direction | all 82 regions hypermethylated in ADCA-DN (Estimate > 0) |
| Array | Illumina 450K |

Usable: **yes** as a DMR (region) list after hg19 -> hg38 liftOver; must be parsed from PDF text
(pypdf extraction is clean and column-aligned). No probe IDs, so probe-level harmonisation must go
through Source A.

## Not applicable

Levy et al. 2022 HGG Advances Table S3 (`../EPISIGN_SHARED/Levy2022_HGGAdv_mmc2.xlsx`) does not
contain an ADCADN probe list (only the paper's 19 new episignatures).
