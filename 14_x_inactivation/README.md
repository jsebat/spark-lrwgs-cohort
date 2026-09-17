# X-chromosome inactivation (14_x_inactivation)

Measures X-inactivation skew in every female in a long-read cohort, from the same HiFi reads used for variant
calling. No array, no expression assay, no extra library.

On the inactive X the promoter CpG islands are methylated and on the active X they are not. Each read comes
from one cell and one haplotype, so the methylation difference between a female's two haplotypes at those
islands estimates the fraction of cells that silenced each X directly.

## The problem this workflow exists to solve

Read-based phasers (HiPhase, WhatsHap without a pedigree) phase from reads alone, so a single sample's X comes
back in on the order of a thousand independent blocks. Haplotype 1 in one block has no relationship to
haplotype 1 in the next. Skew is a property of the whole chromosome, so:

- the **signed** hap1-minus-hap2 difference cancels across blocks and collapses toward zero;
- the **absolute** difference survives but is a folded statistic, biased upward by noise, and carries no
  direction — it cannot say *which* X is silenced.

The pedigree fixes this and the standard pipelines do not use it. The father is hemizygous across the
non-pseudoautosomal X, so at every heterozygous site in a daughter the paternal allele is known outright: it is
whatever single allele the father carries. Comparing that with the child's phased genotype orients each block,
and orienting every block stitches them onto one chromosome-wide pair of labels, paternal and maternal.

`03_trio_phase_x.py` does this. Within-block vote consistency is the built-in check: an incorrect rule produces
votes near 50:50, so consistency below about 0.9 means something is wrong with the pedigree or the call set.

## Two estimators, and when to use which

| Estimator | Needs | Gives | Script |
|---|---|---|---|
| Transmission-phased | Both parents sequenced | Signed skew — names the silenced X | `03` + `04 --mode trio` |
| Folded, noise-corrected | Nothing beyond the sample | Magnitude only | `04 --mode folded` |

The folded estimator deconvolves the observed median \|hap1 − hap2\| in quadrature against a noise floor built
from autosomal CpG islands matched on methylation level. Autosomes are diploid and not subject to inactivation,
so their haplotype asymmetry in the same reads measures assay noise with no true signal. Matching on
methylation level matters because sampling noise peaks near 0.5 and an unmethylated control would understate
the floor.

Run both where you can and check they agree before trusting the folded one on samples that cannot be
trio-phased (in practice, the parents).

## Assigning an unphased structural variant to a parental haplotype

A hemizygous deletion retains no heterozygous site, so phasers leave the SV itself unphased and the deleted
haplotype is not directly observable. `06_sv_parental_haplotype.py` recovers it from reads that cross a
breakpoint: such a read carries the junction *and* extends into flanking phased sequence, where it picks up a
haplotype tag. Each junction read is oriented in its own phase block, and the script requires the reads to
agree before it will make a call. This is what connects a de novo SV to the direction of skew.

## Escape status from the cohort itself

Whether skew direction matters for a given gene depends on that gene being subject to inactivation. Rather than
take an escape list from the literature, `05_escape_status.py` derives it from the same data: a gene subject to
inactivation has a promoter island unmethylated in males (one X, active) and near half-methylated in females
(one active, one inactive), while an escaping gene is unmethylated in both. Recovering known escapees such as
*KDM6A* and *DDX3X* is a useful positive control on any new cohort.

## Pipeline

```bash
# 1. island sets: chrX islands plus autosomal islands for the noise floor
bash 01_make_island_sets.sh --cgi <cpg_islands.bed> --out islands/

# 2. per-haplotype CpG pileups (SLURM array over samples)
sbatch 02_hap_pileups.sb            # edit the header for your account/partition

# 3. orient every X phase block by transmission (trio children only)
python 03_trio_phase_x.py --ped cohort.ped --vcf-dir phased_vcf/ \
       --blocks-dir phase_blocks/ --out orientation.tsv

# 4. skew for every female
python 04_xci_skew.py --ped cohort.ped --pileup-dir pileups/ --islands islands/ \
       --orientation orientation.tsv --out xci_skew.tsv

# 5. which X-linked genes are subject to inactivation here
python 05_escape_status.py --ped cohort.ped --pileup-dir pileups/ \
       --islands islands/chrX_cgi.bed --genes genes.bed --out escape.tsv

# 6. parental haplotype of an unphased SV
python 06_sv_parental_haplotype.py --ped cohort.ped --sample <ID> --bam <bam> \
       --vcf-dir phased_vcf/ --region chrX:START-END --out sv_haplotype.tsv

# 7. X-linked de novo and rare-variant screen against skew
python 07_denovo_x_screen.py --ped cohort.ped --vep chrX.vep.vcf.gz --calls cohort.bcf \
       --skew xci_skew.tsv --gene-sets sets/ --out screen/
```

Every script takes `--help`. No sample identifiers, cohort paths or data files are committed here; supply them
through the arguments above.

## Interpreting the output

`04_xci_skew.py` reports skew as the fraction of cells silencing the more often inactivated X, so 0.5 is
balanced and 1.0 is complete. The conventional clinical bands are 70:30 (skewed), 80:20 (marked) and 90:10
(extreme).

Two cautions that apply to any cohort. Skewing is acquired with age, so mothers and daughters are not one
population and should not be pooled when defining outliers — in the pilot cohort 37 per cent of mothers versus
12 per cent of daughters exceeded 70:30. And the measurement is whatever tissue was sequenced; saliva or blood
inactivation need not represent brain.
