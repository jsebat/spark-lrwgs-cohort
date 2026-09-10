# DNMT3A episignature panel - SPARK lrWGS pilot

See `CLAUDE.md` (study design, hard rules) and `PLAN.md` (implementation steps).
Steps 03-10 run on Expanse (SDSC) via Snakemake + SLURM next to the SPARK lrWGS data; nothing sample-level is committed.

```bash
conda env create -f envs/py.yaml && conda env create -f envs/r.yaml
snakemake -n                                                      # dry run
snakemake --use-conda --cores 2 results/01_variant/annotation.md  # step 01 (public data only)
```
