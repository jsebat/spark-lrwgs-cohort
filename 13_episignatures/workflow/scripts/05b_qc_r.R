#!/usr/bin/env Rscript
# Step 05b - HEpiDISH cell composition and methylclock DNAm age on the imputed probe matrix.
# Usage: Rscript 05b_qc_r.R <results/05_qc dir>
# Inputs : probe_betas.imputed.csv (ProbeID + one column per sample), probe_missing.csv (0/1 mask), pheno.csv
# Outputs: r_composition.csv (sample, epi_frac, fib_frac, immune_frac, + blood subtypes),
#          r_clocks.csv (sample, dnam_age_<clock>..., pct_clock_cpgs_imputed, n_clock_cpgs_missing_<clock>)
suppressPackageStartupMessages({ library(EpiDISH); library(methylclock) })
args <- commandArgs(trailingOnly = TRUE); dir <- args[1]
B <- read.csv(file.path(dir, "probe_betas.imputed.csv"), check.names = FALSE)
M <- read.csv(file.path(dir, "probe_missing.csv"), check.names = FALSE)
rownames(B) <- B$ProbeID; B$ProbeID <- NULL; beta <- as.matrix(B)
rownames(M) <- M$ProbeID; M$ProbeID <- NULL; miss <- as.matrix(M) == 1
pheno <- read.csv(file.path(dir, "pheno.csv"))
cat("matrix:", nrow(beta), "probes x", ncol(beta), "samples\n")

## ---- HEpiDISH: epithelial / fibroblast / immune, then immune subtypes (RPC)
data(centEpiFibIC.m); data(centBloodSub.m)
ref_hit <- intersect(rownames(centEpiFibIC.m), rownames(beta))
cat("HEpiDISH reference CpGs present:", length(ref_hit), "/", nrow(centEpiFibIC.m), "\n")
hep <- hepidish(beta.m = beta, ref1.m = centEpiFibIC.m, ref2.m = centBloodSub.m, h.CT.idx = 3, method = "RPC")
comp <- data.frame(sample = rownames(hep), epi_frac = hep[, "Epi"], fib_frac = hep[, "Fib"],
                   immune_frac = rowSums(hep[, setdiff(colnames(hep), c("Epi", "Fib")), drop = FALSE]),
                   hep[, setdiff(colnames(hep), c("Epi", "Fib")), drop = FALSE], check.names = FALSE)
comp$n_hepidish_ref_cpgs_present <- length(ref_hit)
write.csv(comp, file.path(dir, "r_composition.csv"), row.names = FALSE)

## ---- methylclock: Horvath (pan-tissue), Hannum, skinHorvath, PedBE (buccal/saliva paediatric), Levine
clocks <- c("Horvath", "Hannum", "Levine", "skinHorvath", "PedBE")
df <- data.frame(ProbeID = rownames(beta), beta, check.names = FALSE)
chk <- tryCatch(checkClocks(df), error = function(e) NULL)
age <- DNAmAge(df, clocks = clocks, age = pheno$age[match(colnames(beta), pheno$sample_id)], fastImputation = FALSE, normalize = FALSE)
age <- as.data.frame(age)
out <- data.frame(sample = age$id)
for (cl in clocks) if (cl %in% colnames(age)) out[[paste0("dnam_age_", cl)]] <- age[[cl]]
## % of clock CpGs that were cohort-mean imputed, per sample, using the union of clock CpGs found in the matrix
clock_cpgs <- unique(unlist(lapply(clocks, function(cl) {
  cf <- tryCatch(get(paste0("coef", cl), envir = asNamespace("methylclock")), error = function(e) NULL)
  if (is.null(cf)) cf <- tryCatch(methylclockData::get_coefs(cl), error = function(e) NULL)
  if (is.null(cf)) character(0) else as.character(cf$CpGmarker)
})))
clock_cpgs <- intersect(clock_cpgs, rownames(miss))
cat("clock CpGs in matrix:", length(clock_cpgs), "\n")
if (length(clock_cpgs) > 0) {
  pct <- colMeans(miss[clock_cpgs, , drop = FALSE]) * 100
  out$pct_clock_cpgs_imputed <- round(pct[out$sample], 2)
} else out$pct_clock_cpgs_imputed <- NA
if (!is.null(chk)) writeLines(capture.output(print(chk)), file.path(dir, "r_checkClocks.txt"))
write.csv(out, file.path(dir, "r_clocks.csv"), row.names = FALSE)
cat("done\n")
