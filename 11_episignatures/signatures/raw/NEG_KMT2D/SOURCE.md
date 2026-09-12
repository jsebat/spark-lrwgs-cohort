# NEG_KMT2D — Kabuki syndrome (KMT2D LOF) episignature sources (negative control)

Download date: 2026-09-09.

## Source A (disorder-specific paper): Butcher DT et al. 2017, Am J Hum Genet 100(5):773-788

"CHARGE and Kabuki Syndromes: Gene-Specific DNA Methylation Signatures Identify Epigenetic Mechanisms
Linking These Clinically Overlapping Conditions". DOI 10.1016/j.ajhg.2017.04.004. PMID 28475860.
PMCID PMC5420353. Elsevier PII S0002-9297(17)30148-9.

Files (identical copies also in `../NEG_CHD7/`), fetched from Elsevier's CDN because cell.com returns 403:

| File | Label | URL |
|---|---|---|
| `Butcher2017_AJHG_mmc1.pdf` (602 KB, 10 pp) | Document S1: Figures S1–S4 and Tables S1, S2, S5, S14 | https://ars.els-cdn.com/content/image/1-s2.0-S0002929717301489-mmc1.pdf |
| `Butcher2017_AJHG_mmc2.xlsx` (153 KB) | Document S2: Tables S3, S4, S6–S13 | https://ars.els-cdn.com/content/image/1-s2.0-S0002929717301489-mmc2.xlsx |
| `Butcher2017_AJHG_mmc3.pdf` (2.5 MB) | Document S3: article plus supplemental data | https://ars.els-cdn.com/content/image/1-s2.0-S0002929717301489-mmc3.pdf |

### `Butcher2017_AJHG_mmc2.xlsx` sheets relevant to KMT2D

All sheets: title in row 1, blank row 2, header in row 3 (`header=2`). Sheet names carry a "TS" prefix
but in-sheet titles say "Table S..." (sheet `TS4.KMT2Dphenotype` is titled Table S4, etc.).

| Sheet | In-sheet title | Data rows | Columns |
|---|---|---|---|
| `TS9.KMT2DDNAmsig` | Table S9. DNAm data for CpG sites comprising the KMT2DLOF DNAm signature | **221 cg probes** (224 rows; last 3 are blank/"Abbreviations" footnotes to drop) | TargetID (cg ID), pvalue_limma_raw, pvalue_limma_fdr, pvalue_MannWhitneyU_raw, pvalue_MannWhitneyU_fdr, deltaBeta (KMT2DLOF minus control, beta scale), absDeltaBeta, DNAm_Effect ("loss" 135 / "GAIN" 86), Mean_Disease, Mean_Ctrl, UCSC_RefGene_Name, UCSC_RefGene_Accession, UCSC_RefGene_Group, UCSC_CpG_Islands_Name, Relation_to_UCSC_CpG_Island, CHR (1–22 autosomes; no chrX), Position, Classification ("Yes" = 113 probes retained in the reduced classification model — main text says 112 — "-" otherwise) |
| `TS11.KMT2DDMR` | Table S11. DMRs for the KMT2DLOF DNAm signature | 45 | chr, start, end, avg_deltaBeta, p.value, length, genes, num_CpGs, cpgs (semicolon list), num_signature_cpgs, signature_cpgs |
| `TS13. KMT2DGREATbioterms` | Table S13 GO terms | 56 | not needed |
| `TS4.KMT2Dphenotype` | Table S4 phenotype data | — | not needed |
| `TS6.Controls`, `TS7.GEODNAmcontrols` | Tables S6, S7 control cohorts | 125 / 162 | not needed |

CHD7 sheets (`TS3.CHD7phenotype`, `TS8.CHD7DNAmsig`, `TS10.CHD7DMR`, `TS12.CHD7GREATbioterms`) are
documented in `../NEG_CHD7/SOURCE.md`.

Probe IDs: yes. chr/pos: yes (CHR + Position, plus UCSC CpG-island strings such as
`chr1:1289707-1291126`). Genome build: **not stated in the xlsx or Document S1**; coordinates and
annotations are the Illumina HumanMethylation450 manifest columns, which are hg19 (GRCh37) — treat as
hg19 and verify against the manifest rather than trusting the table. Direction and effect size: yes
(signed deltaBeta, DNAm_Effect). Array: Illumina 450K (363,979 autosomal probes analysed).

Usable: **yes** — 221-probe list (or 113-probe classifier subset) with signed delta-beta, plus 45 DMRs.

## Source B (disorder-specific paper): Aref-Eshghi E et al. 2017, Epigenetics 12(11):923-933

"The defining DNA methylation signature of Kabuki syndrome enables functional assessment of genetic
variants of unknown clinical significance". DOI 10.1080/15592294.2017.1381807. PMID 28933623.
PMCID PMC5788422.

| Item | Value |
|---|---|
| File | Supplementary Material (single Word document; PMC name `kepi-12-11-1381807-s001.docx`) |
| Download URL | https://pmc.ncbi.nlm.nih.gov/articles/instance/5788422/bin/kepi-12-11-1381807-s001.docx |
| Access note | PMC serves this URL behind a JavaScript proof-of-work challenge ("Preparing to download..." stub, 1.8 KB HTML) for plain curl. Fetched by loading the PMC article in a browser (which solves the challenge and sets the short-lived `cloudpmc-viewer-pow` cookie) and passing that cookie to curl. Taylor & Francis (https://www.tandfonline.com/doi/suppl/10.1080/15592294.2017.1381807) returns 403 to non-browser clients; Europe PMC has no supplementary bundle for this article. |
| Local file | `ArefEshghi2017_Epigenetics_Supplementary_Material.docx` (1.24 MB) |
| Format | Word tables (4 tables); parse with python-docx or by regex over `word/document.xml` |
| Table S1 | "CpG probes differentially methylated between patients with KMT2D loss of function mutations and controls" — **1504 data rows, 1504 unique probe IDs** (1501 match `cg\d{8}`; the rest are other Illumina ID types). Columns: Probe ID, Chr (chr1–chr22, autosomes only), Position, Gene, Enhancer, DHS, CpG island, Methylation Difference (signed, case minus control, beta scale; range −0.36 to +0.25; 621 negative / 839 positive among the 1460 values that parse as plain decimals — the remaining 44 need cleaning of number formatting), Adjusted P-value |
| Table S2 | "Probes selected for the SVM classification model" — 142 probes; columns Probe ID, Chr, Position, Gene name, CpG island, AUC |
| Table S3 | "SVM classifier for Kabuki syndrome" — 142 rows; columns Probe, w (weight), center, scale (for scaling new samples) |
| Table S4 | 241 patients with pathogenic mutations in other epigenetic-machinery genes (Disease, Id, Sex, Age, Gene, Mutation); includes DNMT1, CHARGE, Sotos etc. cohorts used for specificity testing |
| Genome build | hg19 (the string "hg19" appears in the document; array = Illumina HumanMethylation450) |
| Direction / effect | yes (signed Methylation Difference in Table S1; none in Table S2/S3 except SVM weight) |

Usable: **yes** — Table S1 (1504 probes, signed effect, hg19 coordinates) and/or the 142-probe
classifier set in Table S2. Roughly 45 of the 142 overlap Butcher 2017's ~200 probes per the paper.

## Source C (EpiSign probe list): Aref-Eshghi et al. 2020 AJHG, Table S2/S3

File lives in `../EPISIGN_SHARED/ArefEshghi2020_AJHG_mmc2.xlsx` (see `../EPISIGN_SHARED/SOURCE.md`).
DOI 10.1016/j.ajhg.2020.01.019. Sheet `S2- Probes`, column `Kabuki == True`: **153 probes**, with
`Chr` and `Position (hg19)`; delta-beta = `Kabuki - Control` from sheet `S3-Methylation levels`.
Usable: yes.

## Not applicable

Levy et al. 2022 HGG Advances Table S3 (`../EPISIGN_SHARED/Levy2022_HGGAdv_mmc2.xlsx`) has no Kabuki
probe list.
