#!/bin/bash
#SBATCH -J dnmt3a_smoke04
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH -N 1 -n 1 -c 8 --mem=32G -t 02:00:00
#SBATCH -o logs/smoke04.%j.out
# Step 04 smoke test: chr22 slice of one non-proband child's merged BAM -> 04_extract_methylation.py -> parquet + QC.
# Submit from the repo root on Expanse:  sbatch tools/smoke04.sbatch
set -euo pipefail
export MAMBA_ROOT_PREFIX=$HOME/micromamba OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
eval "$(micromamba shell hook -s bash)"; micromamba activate dnmt3a-py
S=/expanse/lustre/scratch/jsebat/temp_project/dnmt3a-episignature/smoke; mkdir -p "$S"
# columns: 1 sample_id, 3 role, 13 bam_source, 14 bam_paths, 15 is_proband   (sample ID is never printed)
BAM=$(awk -F'\t' 'NR>1 && $3=="child" && $13=="merged" && $15=="False" {print $14; exit}' results/03_manifest/samples.frozen.tsv)
echo "picked one child sample; bam size $(stat -L -c %s "$BAM") bytes"
samtools view -@8 -b -o "$S/slice.chr22.bam" "$BAM" chr22 && samtools index "$S/slice.chr22.bam"
echo "slice: $(stat -c %s "$S/slice.chr22.bam") bytes, reads $(samtools view -c "$S/slice.chr22.bam")"
/usr/bin/time -v python workflow/scripts/04_extract_methylation.py --config config/config.yaml --sample SMOKE \
  --bams "$S/slice.chr22.bam" --out "$S/SMOKE.cpg.parquet" --threads 8 --workdir "$S/work" 2>&1 | grep -E '^\{|Elapsed|Maximum resident'
python - <<'PY'
import pandas as pd
d = pd.read_parquet("/expanse/lustre/scratch/jsebat/temp_project/dnmt3a-episignature/smoke/SMOKE.cpg.parquet")
print(d.head(3).to_string()); print("rows", len(d), "chroms", d.chrom.unique().tolist()); print(d[["beta", "depth"]].describe().round(3).to_string())
PY
echo SMOKE_DONE
