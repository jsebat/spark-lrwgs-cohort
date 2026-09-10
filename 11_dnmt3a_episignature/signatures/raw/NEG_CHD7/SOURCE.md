# NEG_CHD7 — CHARGE syndrome (CHD7 LOF) episignature sources (negative control)

Download date: 2026-09-09.

## Source A (disorder-specific paper): Butcher DT et al. 2017, Am J Hum Genet 100(5):773-788

"CHARGE and Kabuki Syndromes: Gene-Specific DNA Methylation Signatures Identify Epigenetic Mechanisms
Linking These Clinically Overlapping Conditions". DOI 10.1016/j.ajhg.2017.04.004. PMID 28475860.
PMCID PMC5420353. Elsevier PII S0002-9297(17)30148-9.

Files (identical copies also in `../NEG_KMT2D/`), fetched from Elsevier's CDN because cell.com returns 403:

| File | Label | URL |
|---|---|---|
| `Butcher2017_AJHG_mmc1.pdf` (602 KB, 10 pp) | Document S1: Figures S1–S4 and Tables S1, S2, S5, S14 | https://ars.els-cdn.com/content/image/1-s2.0-S0002929717301489-mmc1.pdf |
| `Butcher2017_AJHG_mmc2.xlsx` (153 KB) | Document S2: Tables S3, S4, S6–S13 | https://ars.els-cdn.com/content/image/1-s2.0-S0002929717301489-mmc2.xlsx |
| `Butcher2017_AJHG_mmc3.pdf` (2.5 MB) | Document S3: article plus supplemental data | https://ars.els-cdn.com/content/image/1-s2.0-S0002929717301489-mmc3.pdf |

### `Butcher2017_AJHG_mmc2.xlsx` sheets relevant to CHD7

All sheets: title in row 1, blank row 2, header in row 3 (`header=2`). Sheet names carry a "TS" prefix
but the in-sheet titles say "Table S..." (sheet `TS3.CHD7phenotype` is titled Table S2, etc.).

| Sheet | In-sheet title | Data rows | Columns |
|---|---|---|---|
| `TS8.CHD7DNAmsig` | Table S8. DNAm data for CpG sites comprising the CHD7LOF DNAm signature | **163 cg probes** (166 rows; last 3 are blank/"Abbreviations" footnotes to drop) | TargetID (cg ID), p_limma_raw, p_limma_fdr, p_MannWhitneyU_raw, p_MannWhitneyU_fdr, deltaBeta (CHD7LOF minus control, beta scale), absDeltaBeta, DNAm_Effect ("loss" 117 / "GAIN" 46), Mean_Disease, Mean_Ctrl, UCSC_RefGene_Name, UCSC_RefGene_Accession, UCSC_RefGene_Group, UCSC_CpG_Islands_Name, Relation_to_UCSC_CpG_Island, CHR (1–21 autosomes; no chrX), Position, Classifcation [sic] ("Yes" = 75 probes retained in the reduced classification model, "-" otherwise) |
| `TS10.CHD7DMR` | Table S10. DMRs for the CHD7LOF DNAm signature | 13 | chr, start, end, avg_deltaBeta, p.value, length, genes, num_CpGs, cpgs (semicolon list), num_signature_cpgs, signature_cpgs |
| `TS12.CHD7GREATbioterms` | Table S12 GO terms | 134 | not needed |
| `TS3.CHD7phenotype` | Table S2 phenotype data | — | not needed |

Other sheets (`TS4.KMT2Dphenotype`, `TS6.Controls`, `TS7.GEODNAmcontrols`, `TS9.KMT2DDNAmsig`,
`TS11.KMT2DDMR`, `TS13. KMT2DGREATbioterms`) are documented in `../NEG_KMT2D/SOURCE.md`.

Probe IDs: yes. chr/pos: yes (CHR + Position, plus UCSC CpG-island strings such as
`chr1:3511709-3512118`). Genome build: **not stated in the xlsx or Document S1**; coordinates and
annotations are the Illumina HumanMethylation450 manifest columns, which are hg19 (GRCh37) — treat
as hg19 and verify against the manifest rather than trusting the table. Direction and effect size:
yes (signed deltaBeta, DNAm_Effect). Array: Illumina 450K (363,979 autosomal probes analysed).

Usable: **yes** — 163-probe list (or 75-probe classifier subset) with signed delta-beta, plus 13 DMRs.

## Source B (EpiSign probe list): Aref-Eshghi et al. 2020 AJHG, Table S2/S3

File lives in `../EPISIGN_SHARED/ArefEshghi2020_AJHG_mmc2.xlsx` (see `../EPISIGN_SHARED/SOURCE.md`).
DOI 10.1016/j.ajhg.2020.01.019. Sheet `S2- Probes`, column **`CHARGE2` == True: 148 probes**, with
`Chr` and `Position (hg19)`. The matching mean-beta column in sheet `S3-Methylation levels` is named
`CHARGE` (not CHARGE2); delta-beta = `CHARGE - Control`. Usable: yes.

## Not applicable

Levy et al. 2022 HGG Advances Table S3 (`../EPISIGN_SHARED/Levy2022_HGGAdv_mmc2.xlsx`) has no CHARGE
probe list.
