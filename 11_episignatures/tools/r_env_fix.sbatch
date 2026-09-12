#!/bin/bash
#SBATCH -J dnmt3a_renv
#SBATCH -A ddp195
#SBATCH -p ind-shared
#SBATCH -N 1 -n 1 -c 4 --mem=16G -t 03:00:00
#SBATCH -o logs/r_env_fix.%j.out
# methylclock is not on bioconda and its BiocManager install failed on compiled deps (png, XML, xml2, RCurl, Rhtslib,
# Rhdf5lib, nloptr, systemfonts ...). Install those and the heavy Bioconductor deps from conda, then let BiocManager
# fetch the pure-R remainder (methylclockData, methylclock). Submit from the repo root: sbatch tools/r_env_fix.sbatch
set -uo pipefail
export MAMBA_ROOT_PREFIX=$HOME/micromamba OPENBLAS_NUM_THREADS=4
micromamba install -y -n dnmt3a-r -c conda-forge -c bioconda \
  r-xml r-xml2 r-rcurl r-png r-nloptr r-systemfonts r-fs r-tidyverse r-ggpubr r-quadprog r-rpmm r-performanceanalytics \
  r-gridextra r-devtools r-tibble r-dplyr bioconductor-minfi bioconductor-impute bioconductor-experimenthub \
  bioconductor-rhtslib bioconductor-rhdf5lib bioconductor-rhdf5filters bioconductor-preprocesscore bioconductor-planet \
  bioconductor-biobase bioconductor-genomicranges || echo "CONDA_INSTALL_HAD_ERRORS"
micromamba run -n dnmt3a-r Rscript -e 'options(Ncpus=4); BiocManager::install(c("methylclockData","methylclock"), ask=FALSE, update=FALSE); cat("methylclock", requireNamespace("methylclock", quietly=TRUE), "\n")'
micromamba run -n dnmt3a-r Rscript -e 'suppressMessages(library(methylclock)); cat("methylclock loads\n")' || echo "METHYLCLOCK_LOAD_FAILED"
echo R_ENV_FIX_DONE
