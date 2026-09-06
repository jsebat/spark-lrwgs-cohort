# Rare tandem-repeat expansion burden: affected offspring vs unaffected founders, family-stratified, depth-adjusted.
suppressMessages(library(survival))
O <- Sys.getenv("OUT", Sys.getenv("OUT"))
ps <- read.delim(file.path(O, "tr_persample.tsv"), stringsAsFactors = FALSE)
qc <- read.delim(Sys.getenv("SAMPLE_QC"), stringsAsFactors = FALSE)
ps$depth <- qc$depth_mean[match(ps$sample, qc$SAMPLE)]
ps <- ps[ps$loci_called > 0 & !is.na(ps$depth), ]
ps$case <- ps$affected
cat("samples", nrow(ps), " cases", sum(ps$case), " families", length(unique(ps$family)), "\n")
cat("outlier loci per sample: cases median", median(ps$outlier_loci[ps$case == 1]), " controls median", median(ps$outlier_loci[ps$case == 0]), "\n")
cat("cor(depth, outlier_loci) =", round(cor(ps$depth, ps$outlier_loci), 3), "  cor(depth, loci_called) =", round(cor(ps$depth, ps$loci_called), 3), "\n")
ps$out_std <- ps$outlier_loci / sd(ps$outlier_loci)
ps$outcds_std <- ps$outlier_cds_loci / max(1e-9, sd(ps$outlier_cds_loci))
show <- function(m, label) { co <- summary(m)$coefficients; cat("\n==", label, "==\n"); print(round(cbind(OR = exp(co[, 1]), lo = exp(co[, 1] - 1.96 * co[, 3]), hi = exp(co[, 1] + 1.96 * co[, 3]), p = co[, 5]), 4)); cat("n =", m$n, " events =", m$nevent, "\n") }
show(clogit(case ~ out_std + strata(family), data = ps), "M1 outlier loci (per SD) only")
show(clogit(case ~ out_std + depth + strata(family), data = ps), "M2 outlier loci (per SD) + depth")
show(clogit(case ~ out_std + depth + loci_called + strata(family), data = ps), "M3 + loci called")
show(clogit(case ~ outcds_std + depth + strata(family), data = ps), "M4 CDS outlier loci (per SD) + depth")
w <- wilcox.test(outlier_per_10k ~ case, data = ps); cat("\nWilcoxon outlier rate per 10k loci, cases vs controls: p =", signif(w$p.value, 3), "\n")
print(aggregate(cbind(outlier_loci, outlier_cds_loci, outlier_per_10k, loci_called, depth) ~ case, data = ps, FUN = function(x) round(median(x), 2)))
