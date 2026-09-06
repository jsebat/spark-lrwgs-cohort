# Conditional logistic regression, inherited LoF tier burden -- with GENOME-WIDE
# COVERAGE as the primary nuisance covariate (JS 2026-09-06).
#
# Rationale: cases average 22.37x and controls 24.41x (8.3% deficit), and the cohort
# spans 8.1x-52.5x. Lower depth means fewer detected variants, which deflates case
# counts technically rather than biologically. Coverage is therefore the most
# important nuisance covariate, ahead of the site-level lowqual_skipped I used first
# (which is partly circular: it is measured at the same sites as the outcome).
suppressPackageStartupMessages({library(survival); library(data.table)})

base <- "/expanse/lustre/projects/ddp195/jsebat/longread-autism/tiering"
d  <- fread(file.path(base, "lof_counts_persample.tsv"))
qc <- fread(file.path(base, "sample_qc.tsv"))
d  <- merge(d, qc, by = "SAMPLE", all.x = TRUE)
cat("merged samples:", nrow(d), " with depth:", sum(!is.na(d$depth_mean)), "\n")

inf <- d[, .(nca = sum(case_binary == 1), nco = sum(case_binary == 0)), by = FID][nca > 0 & nco > 0]
d <- d[FID %in% inf$FID & !is.na(depth_mean)]
cat("informative families:", nrow(inf), " samples:", nrow(d), "\n\n")

cat("=== covariate means by status ===\n")
print(d[, .(n = .N, t1 = mean(tier1), allLoF = mean(all_lof),
            depth = mean(depth_mean), lowq = mean(lowqual_skipped),
            readq = mean(read_q_mean), mapped = mean(mapped_pct)),
        by = .(case = ifelse(case_binary == 1, "case", "control"))])

cat("\n=== ARE THE COVARIATES ORTHOGONAL? (Pearson r) ===\n")
cv <- d[, .(depth_mean, lowqual_skipped, read_q_mean, mapped_pct, all_lof, tier1)]
cm <- cor(cv, use = "pairwise.complete.obs")
print(round(cm, 3))
r_dl <- cm["depth_mean", "lowqual_skipped"]
cat(sprintf("\n  depth vs lowqual_skipped r = %.3f -> %s\n", r_dl,
            ifelse(abs(r_dl) < 0.5, "sufficiently orthogonal to include both",
                   "COLLINEAR: do not include both")))
cat(sprintf("  depth vs all_lof        r = %.3f (coverage drives detected burden)\n",
            cm["depth_mean", "all_lof"]))

show <- function(label, fit) {
  cat("=====", label, "=====\n")
  if (is.null(fit)) { cat("  failed to converge\n\n"); return(invisible()) }
  s <- summary(fit); co <- s$coefficients
  cat(sprintf("  %-18s %9s %8s %8s %8s %9s\n", "term", "beta", "OR", "SE", "z", "p"))
  for (i in seq_len(nrow(co)))
    cat(sprintf("  %-18s %9.4f %8.4f %8.4f %8.3f %9.4f\n", rownames(co)[i],
                co[i, "coef"], co[i, "exp(coef)"], co[i, "se(coef)"],
                co[i, "z"], co[i, "Pr(>|z|)"]))
  cat(sprintf("  n=%d events=%d loglik=%.3f\n\n", s$n, s$nevent, fit$loglik[2]))
}

fit <- function(f) tryCatch(clogit(f, data = d), error = function(e) NULL)

show("MODEL A: tiers + genome-wide depth (primary)",
     fit(case_binary ~ tier1 + tier2 + tier3 + depth_mean + strata(FID)))
if (abs(r_dl) < 0.5) {
} else {
  cat("===== M3 skipped: depth and lowqual are collinear =====\n\n")
}
show("MODEL B: tier 1 alone + depth",
     fit(case_binary ~ tier1 + depth_mean + strata(FID)))

cat("OR > 1 means the tier is enriched in cases. Read the depth coefficient as a check:\n")
cat("if depth carries a large positive effect, detected burden is coverage-driven and\n")
cat("the unadjusted tier estimates were biased downward for the cases.\n")
