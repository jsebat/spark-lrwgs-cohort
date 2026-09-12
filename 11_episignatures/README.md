# Methylation episignatures - SPARK lrWGS pilot

Scores every genome in the cohort against a panel of published DNA-methylation episignatures and reports
which of them can carry interpretation. The panel is not specific to any one gene: it harmonises 55
published signatures to a common reference build, derives concordance-filtered cores wherever two or more
independent sources describe the same syndrome, and scores each sample by shape correlation against a
leave-one-out null. A gene-specific run (for example the Tatton-Brown-Rahman signature used to confirm a
DNMT3A lesion in this cohort) is one configuration of that panel, not a separate workflow.

See `CLAUDE.md` (study design, hard rules) and `PLAN.md` (implementation steps).
Steps 03-10 run on Expanse (SDSC) via Snakemake + SLURM next to the SPARK lrWGS data; nothing sample-level is committed.

```bash
conda env create -f envs/py.yaml && conda env create -f envs/r.yaml
snakemake -n                                                      # dry run
snakemake --use-conda --cores 2 results/01_variant/annotation.md  # step 01 (public data only)
```
