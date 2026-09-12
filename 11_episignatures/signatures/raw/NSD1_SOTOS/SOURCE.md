# NSD1_SOTOS — Sotos syndrome (NSD1) episignature sources

Download date: 2026-09-09.

## Source A (large probe list with effect sizes): Choufani S et al. 2015, Nat Commun 6:10207

"NSD1 mutations generate a genome-wide DNA methylation signature". DOI 10.1038/ncomms10207.
PMID 26690673. PMCID PMC4703864.

| Item | Value |
|---|---|
| Table | Supplementary Data 3 — "List of NSD1+/- significant CpG sites" |
| Publisher file name | `ncomms10207-s4.xlsx` (Nature numbers the SI PDF as s1, so Supplementary Data 3 = s4) |
| Download URL used | https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4703864/supplementaryFiles (Europe PMC zip of all PMC4703864 supplementary files, 10.5 MB; `ncomms10207-s4.xlsx` and `ncomms10207-s1.pdf` extracted from it) |
| URLs that did not work | https://pmc.ncbi.nlm.nih.gov/articles/instance/4703864/bin/ncomms10207-s4.xlsx (PMC returns a JavaScript proof-of-work "Preparing to download" stub to curl); https://www.nature.com/articles/ncomms10207 (redirects to idp.nature.com for non-browser clients) |
| Local files | `Choufani2015_NatCommun_SupplementaryData3_ncomms10207-s4.xlsx` (1.5 MB); `Choufani2015_NatCommun_SupplementaryInformation_ncomms10207-s1.pdf` (8.7 MB; Supplementary Figures 1-2, Methods, References) |
| Sheet | `Suppl. Data 3` (single sheet); title in row 1, blank row 2, header in row 3 (`header=2`) |
| Rows | 7087 data rows, 7087 unique cg IDs (main text quotes 7,085) |
| Columns | Illumina ID; p-value; Bonferroni corrected p-value; deltaBeta (Sotos minus control, beta scale); Absolute deltaBeta; DNA methylation effect ("Loss of methylation" 7038 / "Gain of methylation" 47); Mean not-SS; Mean SS; Regression absolute deltaBeta; Regression p-value; Bonferroni regression p-value; Genome_Build (37); Chromosome (1-22, X); Genomic Location (NCBI, hg19); Strand; Relation_to_UCSC_CpG_Island; Gene Symbol (Genomic features) |
| Probe IDs | yes |
| chr/pos | yes, hg19 / GRCh37 (explicit `Genome_Build = 37` column and header "NCBI, hg19"). Includes chrX probes. |
| Direction / effect size | yes: signed `deltaBeta` and categorical effect; overwhelmingly hypomethylation in Sotos |
| Array | Illumina HumanMethylation450 |

Usable: **yes** (probe-level list with signed delta-beta; selection criterion was >20% absolute
mean difference plus Bonferroni significance). chrX probes should be handled per project sex policy.

## Source B (EpiSign probe list): Aref-Eshghi et al. 2020 AJHG, Table S2/S3

File lives in `../EPISIGN_SHARED/ArefEshghi2020_AJHG_mmc2.xlsx` (see `../EPISIGN_SHARED/SOURCE.md`).
DOI 10.1016/j.ajhg.2020.01.019. Sheet `S2- Probes`, column `Sotos == True`: **112 probes**, with
`Chr` and `Position (hg19)`; delta-beta derivable from sheet `S3-Methylation levels` as
`Sotos - Control`. Usable: yes.

## Not applicable

Levy et al. 2022 HGG Advances Table S3 (`../EPISIGN_SHARED/Levy2022_HGGAdv_mmc2.xlsx`) has no Sotos
probe list.
