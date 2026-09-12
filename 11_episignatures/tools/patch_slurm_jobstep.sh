#!/usr/bin/env bash
# snakemake-executor-plugin-slurm-jobstep (0.6.x) strips all SLURM_* env vars except a keep-list before its nested
# `srun`. Expanse keeps slurm.conf at a non-default path advertised only via SLURM_CONF, so the nested srun fails with
# "Could not establish a configuration source". This adds SLURM_CONF to the keep-list. Idempotent. Run after env creation:
#   tools/patch_slurm_jobstep.sh
set -euo pipefail
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
F=$(ls "$MAMBA_ROOT_PREFIX"/envs/dnmt3a-py/lib/python3*/site-packages/snakemake_executor_plugin_slurm_jobstep/__init__.py)
if grep -q '"SLURM_CONF",' "$F"; then echo "already patched: $F"; exit 0; fi
cp "$F" "$F.orig"
sed -i 's/^\(\s*\)"SLURM_JOB_ID",$/\1"SLURM_CONF",\n\1"SLURM_JOB_ID",/' "$F"
grep -q '"SLURM_CONF",' "$F" && echo "patched: $F (backup $F.orig)" || { echo "patch failed"; exit 1; }
