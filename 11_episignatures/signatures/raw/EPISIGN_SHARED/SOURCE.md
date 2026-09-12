# EPISIGN_SHARED — multi-disorder EpiSign supplementary tables

Shared source files referenced by DNMT1, DNMT3B, NSD1_SOTOS, NEG_KMT2D and NEG_CHD7.
Download date: 2026-09-09. Downloaded with curl from Elsevier's CDN (`ars.els-cdn.com`);
the cell.com article pages themselves return HTTP 403 to non-browser clients.

---

## 1. Aref-Eshghi E et al. 2020, Am J Hum Genet 106(3):356-370

"Evaluation of DNA Methylation Episignatures for Diagnosis and Phenotype Correlations in
42 Mendelian Neurodevelopmental Disorders". DOI 10.1016/j.ajhg.2020.01.019. PMID 32109418.
PMCID PMC7058829. Elsevier PII S0002-9297(20)30019-7.

| File | Label | URL |
|---|---|---|
| `ArefEshghi2020_AJHG_mmc1.pdf` (544 KB, 2 pp) | Document S1 (supplemental figures/notes) | https://ars.els-cdn.com/content/image/1-s2.0-S0002929720300197-mmc1.pdf |
| `ArefEshghi2020_AJHG_mmc2.xlsx` (2.0 MB) | Tables S1–S5 | https://ars.els-cdn.com/content/image/1-s2.0-S0002929720300197-mmc2.xlsx |
| `ArefEshghi2020_AJHG_mmc3.pdf` (3.4 MB) | Document S2, article plus supplemental data | https://ars.els-cdn.com/content/image/1-s2.0-S0002929720300197-mmc3.pdf |

### `ArefEshghi2020_AJHG_mmc2.xlsx` — sheets

Every sheet has a title in spreadsheet row 1 and the column header in row 2
(`pandas.read_excel(..., header=1)`).

| Sheet | Data rows | Columns |
|---|---|---|
| `S1- Samples` | 686 | Syndrome, Subset (Training/Testing), id, Genetic change |
| `S2- Probes` | 3643 (3643 unique cg IDs) | `Probes` (Illumina cg ID), `Chr` (chr1–chr22, autosomes only), `Position (hg19)`, then 34 boolean columns, one per episignature: ADCADN, ADNP_C, ADNP_T, ATRX, AUTS18, BAFopathy2, BFLS, CdLS, CHARGE2, CJS, Down, Dup7, EEOC, FHS, GTPTS, HMA, ICF1, ICF2_3_4, Kabuki, KDVS, Kleefstra1, MRD51, MRX93, MRX97, MRXSN, MRXSSR, RMNS, RSTS, SBBYSS, SETD1B, Sotos, TBRS, WDSTS, Williams |
| `S3-Methylation levels` | 3643 (same probe set as S2) | `Probe`, then mean beta (0–1) per condition: ADCADN, ATRX, AUTS18, BAFopathy, BFLS, CdLS, CHARGE, **Control**, Down, Dup7, EEOC, FLHS, GTPTS, HMA, HVDAS_C, HVDAS_T, ICF1, ICF2_3_4, Kabuki, KDVS, Kleefstra, MRD51, MRX93, MRX97, MRXSCJ, MRXSN, MRXSSR, RMNS, RSTS, SBBYSS, SETD1B, Sotos, TBRS, WDSTS, Williams |
| `S4- Unresolved` | 9 | idat, Predicted phenotype, id, sex, age, Phenotype |
| `S5- Uncertain` | 12 | id, Queried phenotype, Sequence variant, Additional information, DNA methylation classification |

Note the column-name mismatch between S2 and S3 for the same disorder: S2 `CHARGE2` = S3 `CHARGE`;
S2 `BAFopathy2` = S3 `BAFopathy`; S2 `FHS` = S3 `FLHS`; S2 `Kleefstra1` = S3 `Kleefstra`;
S2 `ADNP_C/ADNP_T` = S3 `HVDAS_C/HVDAS_T`.

Probe counts (S2 boolean == True) for the signatures used in this project:

| Project ID | S2 column | n probes |
|---|---|---|
| DNMT1 | ADCADN | 104 |
| DNMT3B | ICF1 | 113 |
| NSD1_SOTOS | Sotos | 112 |
| NEG_KMT2D | Kabuki | 153 |
| NEG_CHD7 | CHARGE2 | 148 |
| (TBRS_LOF) | TBRS | 139 |

Content per probe: probe ID yes; chr/pos yes (hg19, explicit in header); direction and effect size
not given directly but derivable as `S3[<disorder>] - S3[Control]` (mean-beta difference, case minus
control) for every probe in the table. Genome build: hg19 (GRCh37). Array: Illumina 450K / EPIC
(probes common to both).

Usable for building probe lists: **yes** — this is the primary EpiSign probe source for all five IDs.
Coordinates must be re-derived from the Illumina manifest / Zhou-lab hg38 annotation, not taken
from this table, per project rules.

---

## 2. Levy MA et al. 2022, HGG Advances 3(1):100075

"Novel diagnostic DNA methylation episignatures expand and refine the epigenetic landscapes of
Mendelian disorders". DOI 10.1016/j.xhgg.2021.100075. PMID 35047860. PMCID PMC8756545.
Elsevier PII S2666-2477(21)00056-7.

| File | Label | URL |
|---|---|---|
| `Levy2022_HGGAdv_mmc1.pdf` (4.4 MB, 15 pp) | Document S1 (Figures S1–S2 and supplemental notes) | https://ars.els-cdn.com/content/image/1-s2.0-S2666247721000567-mmc1.pdf |
| `Levy2022_HGGAdv_mmc2.xlsx` (283 KB) | Tables S1–S3 | https://ars.els-cdn.com/content/image/1-s2.0-S2666247721000567-mmc2.xlsx |
| `Levy2022_HGGAdv_mmc3.pdf` (8.9 MB) | Document S2, article plus supplemental data | https://ars.els-cdn.com/content/image/1-s2.0-S2666247721000567-mmc3.pdf |

### `Levy2022_HGGAdv_mmc2.xlsx` — sheets

Title in row 1, blank row 2, header in row 3 (`header=2`).

| Sheet | Data rows | Columns |
|---|---|---|
| `Table S1 samples` | 235 | Episignature, Sex, Age, Variant |
| `Table S2 controls` | 19 | Episignature, Case samples, Control samples, Ratio, Notes |
| `Table S3 probes` | 7076 | Episignature, Probe (cg ID), Mean beta value difference (percentage-point scale, range −57.9 to +61.3; case minus control), p value, Adjusted p value |

Table S3 covers **only the 19 new episignatures** of that paper (ARTHS, BEFAHRS, CSS4_c.2650,
CSS9, CSS_c.6200, Chr16p11.2del, DYT28, GADEVS, KDM2B, KDM4B, LLS, MKHK_IDR4, MLASA2, MRXSA, PHMDS,
RENS1, RSTS1, RSTS2, VCFS). It does **not** contain probe lists for ADCADN, ICF, Sotos, Kabuki, CHARGE
or TBRS, even though those disorders are in the paper's 57-episignature classifier. No chr/pos columns
and no genome build stated in the table; probe IDs only.

Usable for the five IDs in this task: **no** (disorders absent). Kept as a shared reference file.
