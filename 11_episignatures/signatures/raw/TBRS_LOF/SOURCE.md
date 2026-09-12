# TBRS_LOF — raw sources (DNMT3A loss-of-function / Tatton-Brown-Rahman syndrome)

Download date for every file in this folder: 2026-09-09.
Nothing here has been edited; files are byte-for-byte as served. Row counts below exclude header/legend rows.
"Direction" was derived by me from the tables themselves (case mean minus control mean) where the table has no
explicit sign column; those derivations are labelled as such.

Access notes (why some files come from mirrors rather than the publisher page):
- genome.cshlp.org supplement page (https://genome.cshlp.org/content/29/7/1057/suppl/DC1) returns 404 / SSO login
  redirect for anonymous curl; files were taken from the Europe PMC supplementary-file bundle for PMC6633263 (open access).
- cell.com and sciencedirect.com return HTTP 403 to non-browser clients. The NCBI PMC "bin" download for PMC7058829
  sits behind a JavaScript proof-of-work bot challenge, which I did not attempt to bypass. The Elsevier CDN
  (ars.els-cdn.com) served the files directly under the article PII S0002929720300197.
- Levy 2022 (open access) came from the Europe PMC supplementary-file bundle for PMC8756545.

---------------------------------------------------------------------------------------------------

## 1. Jeffries AR et al. 2019, Genome Research 29:1057-1066
"Growth disrupting mutations in epigenetic regulatory molecules are associated with abnormalities of epigenetic aging"
DOI 10.1101/gr.243584.118 — PMID 31160375 — PMC6633263
Cohort for the DNMT3A tables: Amish family carriers of DNMT3A c.2312G>A p.(Arg771Gln) vs wild-type relatives,
peripheral blood, Illumina HumanMethylation450 (450K) array, limma (DMPs) and DMRcate (DMRs).
Download URL (all Jeffries files, one zip):
https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6633263/supplementaryFiles?includeInlineImage=false
(zip members were prefixed `supp_`; renamed here with prefix `Jeffries2019_GenomeRes_`.)

### Jeffries2019_GenomeRes_gr.243584.118_Supplemental_Table_S1.xlsx  — Supplemental Table S1 (DMPs)  **PRIMARY**
- Sheets: `Legend`, `Significant probes`
- Rows: 2,606 probes (all unique cg IDs)
- Columns: TargetID (Illumina probe ID), FC (= delta beta carriers minus WT, per legend), t, P.Value, adj.P.Val, B,
  CHR_37, Position, Strand, UCSC_RefGene_Name, UCSC_RefGene_Accession, UCSC_RefGene_Group, UCSC_CpG_Islands_Name,
  Relation_to_UCSC_CpG_Island, Phantom, DMR, Enhancer, Regulatory_Feature_Name, Regulatory_Feature_Group,
  Closest_TSS_gene_name, Closest_TSS_Transcript
- Genome build: GRCh37/hg19 (column literally named CHR_37; Illumina 450K annotation)
- Direction/effect size: yes — FC is delta beta. 2,576 probes FC<0 (hypomethylated in carriers), 30 FC>0; range -0.478..+0.454.
- Usable for a probe list: YES (probe ID + hg19 position + signed delta beta). Map to hg38 via Illumina/Zhou manifest, not by lift-over of the hg19 position.

### Jeffries2019_GenomeRes_gr.243584.118_Supplemental_Table_S3.xlsx  — Supplemental Table S3 (DMRs, DMRcate)  **PRIMARY (region-level)**
- Sheets: `Legend`, `DMRs`
- Rows: 388 DMRs
- Columns: (unnamed index), coord (chrN:start-end), no.cpgs, minfdr, Stouffer, maxbetafc, meanbetafc
- Genome build: hg19 (coordinates derived from the 450K build-37 annotation; not restated in the sheet)
- Direction/effect size: yes — meanbetafc / maxbetafc are signed delta beta. 387 of 388 DMRs meanbetafc<0.
- Usable: YES for a DMR list (needs hg19 -> hg38 liftOver; coordinates are 1-based probe-to-probe spans).

### Other Jeffries files (retained, not signature tables)
- ..._Supplemental_Table_S4.xlsx — GO (21,681 rows) and KEGG (321 rows) enrichment; not usable as a probe list.
- ..._Supplemental_Material.pdf — figure legends and Tables S2, S5, S6, S7 (cell proportions, NSD1/KMT2D variants, primers).
- ..._Supplemental_Data_S1.pdf, ..._Supplemental_Data_S2.pdf, ..._Supplemental_Code_S1.txt, ..._29_7_1057__index.html.
Caveat: the cohort carries a single missense allele (R771Q) shown by the authors to behave as loss-of-function; it is
not a truncating/deletion cohort. Record this in INVENTORY.md when merging with other LOF sources.

---------------------------------------------------------------------------------------------------

## 2. Aref-Eshghi E et al. 2020, Am J Hum Genet 106:356-370
"Evaluation of DNA Methylation Episignatures for Diagnosis and Phenotype Correlations in 42 Mendelian Neurodevelopmental Disorders"
DOI 10.1016/j.ajhg.2020.01.019 — PMID 32109418 — PMC7058829 — PII S0002-9297(20)30019-7
Platform: Illumina 450K / EPIC, peripheral blood (EpiSign v2 probe sets).
Download URLs (Elsevier CDN):
- https://ars.els-cdn.com/content/image/1-s2.0-S0002929720300197-mmc1.pdf   -> ArefEshghi2020_AJHG_mmc1.pdf  (544,262 B; 2-page cover/TOC of Supplemental Data)
- https://ars.els-cdn.com/content/image/1-s2.0-S0002929720300197-mmc2.xlsx  -> ArefEshghi2020_AJHG_mmc2.xlsx (1,995,980 B)  **PRIMARY**
- https://ars.els-cdn.com/content/image/1-s2.0-S0002929720300197-mmc3.pdf   -> ArefEshghi2020_AJHG_mmc3.pdf  (3,426,625 B; 24 pages of supplemental figures)
Tried and blocked: https://www.cell.com/ajhg/fulltext/S0002-9297(20)30019-7 (403),
https://www.sciencedirect.com/science/article/pii/S0002929720300197 (403),
https://pmc.ncbi.nlm.nih.gov/articles/instance/7058829/bin/mmc2.xlsx (returns a proof-of-work challenge page, not bypassed),
Europe PMC supplementaryFiles for PMC7058829 ("not open access").

### ArefEshghi2020_AJHG_mmc2.xlsx
- Sheets: `S1- Samples`, `S2- Probes`, `S3-Methylation levels`, `S4- Unresolved`, `S5- Uncertain`  (each has a title row, header on row 2)
- `S2- Probes` ("Probes selected for all episignatures"): 3,643 probes x 37 columns:
  Probes (cg ID), Chr, Position (hg19), then one boolean column per episignature:
  ADCADN, ADNP_C, ADNP_T, ATRX, AUTS18, BAFopathy2, BFLS, CdLS, CHARGE2, CJS, Down, Dup7, EEOC, FHS, GTPTS, HMA, ICF1,
  ICF2_3_4, Kabuki, KDVS, Kleefstra1, MRD51, MRX93, MRX97, MRXSN, MRXSSR, RMNS, RSTS, SBBYSS, SETD1B, Sotos, TBRS, WDSTS, Williams.
  **TBRS = 139 probes** (autosomal, chr1-20 represented). Also relevant for other panel members: ADCADN (DNMT1) 104,
  ICF1 (DNMT3B) 113, Sotos 112, Kabuki 153, CHARGE2 148.
- `S3-Methylation levels` ("Mean methylation levels for probes in various conditions"): 3,643 probes x 36 columns:
  Probe, then mean beta per group incl. `Control` and `TBRS` (group names differ slightly from S2: BAFopathy, CHARGE, FLHS, HVDAS_C/T, Kleefstra, MRXSCJ).
- `S1- Samples`: 686 samples (Syndrome, Subset train/test, id, Genetic change). TBRS: 10 training + 4 testing.
- Genome build: hg19 (explicit in the column header). No hg38 column.
- Direction/effect size: no signed column in S2. Derived here from S3 as TBRS mean minus Control mean for the 139 TBRS
  probes: all 139 negative (hypomethylated), median -0.179, range -0.384..-0.109.
- Usable: YES — probe IDs + derived delta beta. Overlap with Jeffries Table S1 DMPs: 45 probes.

---------------------------------------------------------------------------------------------------

## 3. Levy MA et al. 2022, HGG Advances 3:100075
"Novel diagnostic DNA methylation episignatures expand and refine the epigenetic landscapes of Mendelian disorders"
DOI 10.1016/j.xhgg.2021.100075 — PMID 35047860 — PMC8756545
Download URL: https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8756545/supplementaryFiles?includeInlineImage=false
(members mmc1.pdf, mmc2.xlsx, mmc3.pdf; renamed with prefix `Levy2022_HGGAdv_`). cell.com page returned 403.

### Levy2022_HGGAdv_mmc2.xlsx
- Sheets: `Table S1 samples` (235 rows), `Table S2 controls` (19 rows), `Table S3 probes` (7,076 rows)
- `Table S3 probes` = "List of probes used for the 19 new episignatures": columns Episignature, Probe, Mean beta value difference
  (values on a percent-like scale, e.g. -7.24), p value, Adjusted p value. No chromosome/position, no genome build stated.
- Episignatures present: ARTHS, BEFAHRS, Chr16p11.2del, CSS4_c.2650, CSS9, CSS_c.6200, DYT28, GADEVS, KDM2B, KDM4B, LLS,
  MKHK_IDR4, MLASA2, MRXSA, PHMDS, RENS1, RSTS1, RSTS2, VCFS.
- **TBRS: not present (0 rows).** This paper only publishes probes for its 19 new signatures.
- Usable for TBRS_LOF: NO. Retained for reference only (no DNMT3A/TBRS content). mmc1.pdf/mmc3.pdf = supplemental methods/figures.

---------------------------------------------------------------------------------------------------

## 4. Smith AM et al. 2021, Nature Communications 12:4549
"Functional and epigenetic phenotypes of humans and mice with DNMT3A Overgrowth Syndrome"
DOI 10.1038/s41467-021-24800-7 — PMID 34315901 — PMC8316576
Human data: WGBS of peripheral blood, 11 DOS patients (3 R882: UPN 624400 R882H, 154605 R882H, 894912 R882C;
8 non-R882: 411168 C583Y, 511909 R301W, 228211 F414fsTer7, 295041 R736H, 930075 R688H, 723972 Y660H, 518693 135-kb
whole-gene deletion, 786396 I310N) vs 15 healthy donors (13 unrelated normal donors NM*/ND* + 2 unaffected siblings
867535, 978897). Reads mapped with biscuit (MGI analysis-workflows bisulfite CWL v1.5.0); DMRs called with metilene,
required >10 CpGs, mean group difference >0.2, FDR-filtered. Raw human data controlled-access (dbGaP phs000159).
Download URL base: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-021-24800-7/MediaObjects/
(each file = base + original file name; the same files are also in the Europe PMC bundle for PMC8316576).

### Smith2021_NatCommun_41467_2021_24800_MOESM5_ESM.xlsx — Supplementary Data 2 (non-R882 DMRs)  **PRIMARY for TBRS_LOF**
- Sheets: `Figures` (lists Fig 1c-g, Supp 2h), `dmrs.annot.summarized`
- Rows: 332 DMRs, chr1-22 (no X/Y)
- Columns: dmr.id (chr_start_end), chromosome, start, end, closest_gene, islands, shores, shelves, genes, enhancers,
  promoters, tss (0/1 annotation flags), size (bp; median 542, range 133-3,309), then 26 per-sample mean-methylation
  columns (15 controls, 8 non-R882, 3 R882). Sample group is only recoverable from the column name (UPN / "TBRS" / "NM"/"ND"/"sib").
- Genome build: hg38. Not stated in the sheet; inferred from (i) the paper reports the UPN 518693 deletion in hg38,
  (ii) the biscuit CWL pipeline is GRCh38-based, (iii) the first DMR chr1:1,034,534-1,034,983 is annotated to AGRN,
  which is true in hg38 (AGRN chr1:1,020,120-1,056,118) but not hg19. Coordinate convention (0- vs 1-based) not stated;
  verify against CpG positions before use.
- Direction/effect size: no explicit column. Derived here: non-R882 mean minus control mean is negative for all 332 DMRs
  (median -0.246, range -0.446..-0.180).
- Usable: YES as an hg38 DMR list with derived delta beta. Caveat: the "non-R882" group is mostly missense; only
  UPN 228211 (frameshift) and 518693 (deletion) are unambiguous haploinsufficiency alleles.

### Smith2021_NatCommun_41467_2021_24800_MOESM4_ESM.xlsx — Supplementary Data 1 (R882 DMRs)
- Same layout, 2,209 DMRs, chr1-22 + X. All hypomethylated in R882 vs controls (median delta -0.378). This is the
  DNMT3A_R882 signature; the copy in ../DNMT3A_R882/ is the one to build from (documented there).

### Other Smith files (retained, not human signature tables)
- ..._MOESM1_ESM.pdf Supplementary Information (figures/legends); ..._MOESM3_ESM.pdf Description of Additional Supplementary Files.
- ..._MOESM6/7_ESM.xlsx Supplementary Data 3/4: human scRNA-seq / bulk RNA-seq DEGs.
- ..._MOESM8/9/10_ESM.xlsx Supplementary Data 5/6/7: mouse bone-marrow WGBS DMRs (Dnmt3a R878H/+, -/-, +/-; mm10).
- ..._MOESM11_ESM.xlsx Supplementary Data 8: mouse scRNA-seq DEGs. ..._MOESM13_ESM.zip: Source Data for figures.
