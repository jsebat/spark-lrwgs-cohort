# DNMT3A_GOF — raw sources (DNMT3A PWWP gain-of-function, Heyn-Sproul-Jackson microcephalic dwarfism)

Download date: 2026-09-09. Files are unmodified publisher ESM files.

## Heyn P et al. 2019, Nature Genetics 51:96-105
"Gain-of-function DNMT3A mutations cause microcephalic dwarfism and hypermethylation of Polycomb-regulated regions"
DOI 10.1038/s41588-018-0274-x — PMID 30478443 — PMC6520989
Patients: P1, P2 DNMT3A c.988T>C p.W330R (de novo), P3 c.997G>A p.D333N. Human DMR tables come from reduced-representation
bisulfite sequencing (RRBS) of patient vs control primary fibroblasts (SI Fig. legends: fibroblasts n = 3 controls,
n = 2 DNMT3A patients P1/P2). Blood was profiled by Infinium MethylationEPIC (incl. 2 TBRS overgrowth patients as
comparators) but NO probe-level supplementary table was published; processed data are in GEO GSE120558. Raw human
sequencing: EGA EGAS00001003231 (exome), EGAS00001003232 (RNA-seq/RRBS/ChIP-seq), controlled access.
Download URL base: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41588-018-0274-x/MediaObjects/
(each file = base + original name, e.g. .../41588_2018_274_MOESM3_ESM.txt). Europe PMC has no supplementary bundle for PMC6520989 (404).

### Heyn2019_NatGenet_41588_2018_274_MOESM3_ESM.txt — Supplementary Table 2  **PRIMARY (hypermethylated)**
- Publisher legend: "Human patient fibroblasts hypermethylated DMRs. Supplied as BED formatted tab-delimited text file. All co-ordinates in hg19."
- Format: 4-column BED, no header, tab-delimited: chrom, start, end, name (DMR_up_1 .. DMR_up_1140)
- Rows: 1,140 DMRs; chromosomes chr1-22 only (no X/Y); median width 1,288 bp; DMR coordinate range 112,617-244,007,281
- Genome build: hg19 (explicit). BED => 0-based half-open starts.
- Direction: hypermethylated in patients (encoded in table identity / "DMR_up" names). No effect size, no p-value, no CpG count.
- Usable: YES as a region list after hg19 -> hg38 liftOver; direction = +1 for every region; no per-region weight available.

### Heyn2019_NatGenet_41588_2018_274_MOESM4_ESM.txt — Supplementary Table 3 (hypomethylated)
- Legend: "Human patient fibroblasts hypomethylated DMRs. Supplied as BED formatted tab-delimited text file. All co-ordinates in hg19."
- 4-column BED, 738 DMRs (DMR_down_1 .. DMR_down_738), chr1-22, median width 1,130 bp, hg19.
- Usable: YES (same caveats). Include as the negative-direction half of the GOF signature if a two-sided shape vector is wanted.

### Heyn2019_NatGenet_41588_2018_274_MOESM5_ESM.txt — Supplementary Table 5 (mouse)
- "Mouse neural differentiation hypermethylated DMRs ... All co-ordinates in mm10." 342 rows, chr1-19/X/Y, median width 76 bp.
- NOT usable for the human panel (mouse genome).

### Heyn2019_NatGenet_41588_2018_274_MOESM6_ESM.txt — Supplementary Table 10 (mouse)
- "Mouse neural differentiation hypomethylated DMRs ... mm10." 222 rows. NOT usable (mouse).

### Heyn2019_NatGenet_41588_2018_274_MOESM1_ESM.pdf — Supplementary Text and Figures
- Supplementary Figs 1-11, Supplementary Tables 1, 4, 6-9 (patient summary, GO terms, oligos, sequencing stats) and
  Supplementary Note (methods: RRBS aligned to hg19; EPIC arrays annotated with IlluminaHumanMethylationEPICanno.ilm10b2.hg19).

### Heyn2019_NatGenet_41588_2018_274_MOESM2_ESM.pdf — Reporting Summary (not data).

Caveats for INVENTORY.md: tissue is fibroblast, not blood/saliva; regions are RRBS-covered CpG-dense loci (Polycomb
targets), so the CpG-density-matched permutation null matters here; no effect sizes means shape_r can only use +/-1 signs.
