# Signature inventory (PLAN step 02)

Generated 2026-09-09T21:58:36 by workflow/scripts/02_harmonize_signatures.py. region_pad = 250 bp; regions merged within signature; hg38; BED 0-based half-open. Source provenance: signatures/raw/<ID>/SOURCE.md and signatures/sources.yaml.

| ID | Sources | n array probes (mapped) | n regions (tier1 / tier2) | Direction | Status | Notes |
|---|---|---|---|---|---|---|
| TBRS_LOF | jeffries2019; arefeshghi2020; smith2021_nonR882 | 2745 | 165 / 2255 | 28 hyper / 2226 hypo | OK | tier1 165 regions (165 hypo); 1 merged regions with mixed direction (direction '.'); median region 502 bp |
| DNMT3A_R882 | smith2021_R882 | 0 | 2173 | 0 hyper / 2173 hypo | OK | median region 1055 bp |
| DNMT3A_GOF | heyn2019 | 0 | 1810 | 1078 hyper / 725 hypo | OK | 7 merged regions with mixed direction (direction '.'); median region 1717 bp |
| DNMT1 | arefeshghi2020; kernohan2016 | 104 | 169 | 168 hyper / 1 hypo | OK | median region 502 bp |
| DNMT3B | arefeshghi2020 | 113 | 102 | 0 hyper / 102 hypo | OK | median region 502 bp |
| NSD1_SOTOS | choufani2015; arefeshghi2020 | 7197 | 4254 | 31 hyper / 4223 hypo | OK | median region 502 bp |
| NEG_KMT2D | butcher2017; arefeshghi2020 | 374 | 265 | 116 hyper / 149 hypo | OK | median region 502 bp |
| NEG_CHD7 | butcher2017; arefeshghi2020 | 311 | 234 | 101 hyper / 133 hypo | OK | median region 502 bp |

Direction = sign of case-minus-control delta beta; 'hyper' = more methylated in cases. TBRS_LOF tier1 = regions supported by >=2 independent sources with concordant direction and no within-region conflict; tier2 = union of all sources. Probe coordinates come from the Zhou-lab hg38 manifests (masked probes excluded); hg19 DMRs were lifted with hg19ToHg38.over.chain (see resources/*/SOURCE.md).

## Source notes and caveats

**TBRS_LOF**
- jeffries2019: 450K, Amish DNMT3A R771Q (LOF-behaving missense) carriers vs WT relatives; limma DMPs (Table S1)
- jeffries2019: same cohort, DMRcate DMRs (Table S3), hg19 probe-to-probe spans
- arefeshghi2020: EpiSign v2 TBRS probe set (139 probes); delta derived as S3 mean(TBRS) - mean(Control)
- smith2021_nonR882: WGBS blood, 8 non-R882 DOS patients vs 15 controls (metilene DMRs, Supp Data 2); hg38 1-based (verified: start = C of first CpG, end = G of last CpG)

**DNMT3A_R882**
- smith2021_R882: WGBS blood, 3 R882H/C DOS patients vs 15 controls (Supp Data 1)

**DNMT3A_GOF**
- heyn2019: HESJAS patient fibroblast RRBS, hypermethylated DMRs (Supp Table 2), hg19 BED; no effect size published
- heyn2019: same study, hypomethylated DMRs (Supp Table 3)

**DNMT1**
- arefeshghi2020: EpiSign v2 ADCADN probe set (104 probes)
- kernohan2016: 450K blood, ADCA-DN vs controls; 82 regions parsed from Additional file 1 PDF; Estimate = case-control difference

**DNMT3B**
- arefeshghi2020: EpiSign v2 ICF1 probe set (113 probes). Velasco 2018 has no genome-wide list (35 validation DMPs only) - not used

**NSD1_SOTOS**
- choufani2015: 450K blood, NSD1+/- Sotos vs controls, 7087 significant CpGs with signed deltaBeta (Supp Data 3)
- arefeshghi2020: EpiSign v2 Sotos probe set (112 probes)

**NEG_KMT2D**
- butcher2017: 450K blood, KMT2D LOF Kabuki signature, 221 probes with signed deltaBeta (Table S9)
- arefeshghi2020: EpiSign v2 Kabuki probe set (153 probes)

**NEG_CHD7**
- butcher2017: 450K blood, CHD7 LOF CHARGE signature, 163 probes with signed deltaBeta (Table S8)
- arefeshghi2020: EpiSign v2 CHARGE probe set (148 probes; S2 column CHARGE2, S3 column CHARGE)

