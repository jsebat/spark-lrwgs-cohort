#!/bin/bash -l
# Launch a detached Snakemake run on an Expanse login node with the SLURM executor profile.
#   tools/launch_snakemake.sh <target> [extra snakemake args...]
# Login shell (-l) so SLURM_CONF and the module environment are set; jobs inherit them via sbatch --export.
# Logs to logs/snakemake_<target>.log. Refuses to start if another run is already active.
set -euo pipefail
cd "$(dirname "$0")/.."
TARGET="${1:?target rule/file}"; shift || true
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
eval "$(micromamba shell hook -s bash)"; micromamba activate dnmt3a-py
if pgrep -u "$USER" -f "[s]nakemake --profile" >/dev/null; then
  echo "a snakemake run is already active:"; pgrep -u "$USER" -af "[s]nakemake --profile" | cut -c1-120; exit 1
fi
snakemake --unlock --snakefile workflow/Snakefile >/dev/null 2>&1 || true
mkdir -p logs
LOG="logs/snakemake_${TARGET//\//_}.log"
setsid nohup snakemake --profile workflow/profiles/expanse "$TARGET" "$@" > "$LOG" 2>&1 < /dev/null &
echo "launched snakemake -> $TARGET (pid $!), log $LOG; SLURM_CONF=${SLURM_CONF:-unset}"
