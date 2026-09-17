# DNMT3B — ICF1 syndrome episignature sources

Download date: 2026-09-09.

## Source A (primary probe list): EpiSign ICF1 — Aref-Eshghi et al. 2020 AJHG, Table S2/S3

File lives in `../EPISIGN_SHARED/ArefEshghi2020_AJHG_mmc2.xlsx` (see `../EPISIGN_SHARED/SOURCE.md`).
DOI 10.1016/j.ajhg.2020.01.019. Sheet `S2- Probes`, column `ICF1 == True`: **113 probes**, with
`Chr` and `Position (hg19)`. Sheet `S3-Methylation levels` gives mean beta for `ICF1` and `Control`
on the same probes, so delta-beta = `ICF1 - Control` and direction are derivable. (A separate
`ICF2_3_4` signature, 95 probes, is also present but is not DNMT3B and is not used.)
Genome build hg19. Usable: yes.

## Source B (checked, NOT usable as a probe/DMR list): Velasco G et al. 2018, Hum Mol Genet 27(14):2409-2424

"Comparative methylome analysis of ICF patients identifies heterochromatin loci that require ZBTB24,
CDCA7 and HELLS for their methylated state". DOI 10.1093/hmg/ddy130. PMID 29659838. Not in PMC
(subscription article; Europe PMC reports isOpenAccess=N, hasSuppl=N).

| Item | Value |
|---|---|
| Supplementary file offered by OUP | a single PDF, "Supplementary Information" (Tables S1–S8, Figures S1–S18, Supplementary Methods) |
| Download URL (base) | https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/hmg/27/14/10.1093_hmg_ddy130/1/ddy130_supplementary_information.pdf |
| Access note | the base URL returns HTTP 403; OUP serves it only through a time-limited signed URL (`?Expires=...&Signature=...&Key-Pair-Id=...`) embedded in the article page https://academic.oup.com/hmg/article/27/14/2409/4975537. The signed link was taken from that page in a browser session and fetched with curl. |
| Local file | `Velasco2018_HMG_ddy130_supplementary_information.pdf` (4.8 MB, 44 pages) |
| Tables inside | S1 genetic status of ICF patients; S2 GO analysis of HypoMPs common to all ICF; **S3 list of DMPs experimentally validated in the study (35 distinct cg IDs)**; S4 GO of ICF1 HypoMPs; S5 germline-enriched genes hypomethylated in promoter in ICF; S6 GO of ICF2/3/4 HypoMPs; S7 hypomethylated gene clusters common to ICF2/3/4; S8 primers |
| Probe IDs | only the 35 validated DMPs in Table S3 (mixed ICF1 and ICF2/3/4 loci); the genome-wide ICF1 HypoMP/DMP lists are **not** provided in the supplement |
| chr/pos | no (one `chr:pos`-style string in the whole PDF) |
| Genome build | not stated in the supplement (array = Illumina HM450K) |
| Direction / effect | Table S3 reports validation beta values, not a genome-wide effect column |

Usable: **no** — no downloadable genome-wide ICF1 DMP/DMR table; only 35 hand-picked validation loci.
Do not approximate an ICF1 list from this file. The DNMT3B signature therefore relies on Source A only.
Any genome-wide list would have to come from the raw data deposit referenced in the main text (not
part of this supplement) and would need a separate analysis step.

## Not applicable

Levy et al. 2022 HGG Advances Table S3 (`../EPISIGN_SHARED/Levy2022_HGGAdv_mmc2.xlsx`) has no ICF
probe list.
