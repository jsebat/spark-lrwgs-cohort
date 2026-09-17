# Deconvolution v2: vectorized projected-gradient NNLS on Loyfer U25 markers + two-compartment estimate.
A <- read.delim(file.path(Sys.getenv("METH_REFS"), "Atlas.U25.l4.hg38.tsv"), check.names = FALSE)
cts <- setdiff(colnames(A), c("chr","start","end","startCpG","endCpG","target","name","direction"))
X <- as.matrix(A[, cts]); X[is.na(X)] <- NA
M <- read.delim(file.path(Sys.getenv("DATA_ROOT"), "meth", "marker_means_U25.tsv"), check.names = FALSE)
mcols <- grep("^m[0-9]+$", colnames(M))
# atlas values are the UNMETHYLATED fraction (U markers: target ~0.9, others ~0.02); convert sample methylation to 1 - m
for (k in mcols) M[[k]] <- 1 - as.numeric(M[[k]])
nnls_pg <- function(X, y, iters = 20000) {
  ok <- !is.na(y) & rowSums(is.na(X)) == 0; Xo <- X[ok, , drop = FALSE]; yo <- y[ok]
  w <- rep(1 / ncol(Xo), ncol(Xo)); L <- 2 * max(eigen(crossprod(Xo), symmetric = TRUE, only.values = TRUE)$values); lr <- 1 / L
  for (i in 1:iters) { g <- 2 * crossprod(Xo, Xo %*% w - yo); w <- pmax(0, w - lr * g) }
  s <- sum(w); list(w = if (s > 0) w / s else w, rss = sum((Xo %*% (w / max(s, 1e-9)) - yo)^2), n = sum(ok), raw_sum = s)
}
out <- data.frame(sample = M$sample, family = M$family, role = M$role, dna_blood = M$dna_blood, n_markers = NA, rss = NA, raw_sum = NA, stringsAsFactors = FALSE)
for (ct in cts) out[[ct]] <- NA
for (i in seq_len(nrow(M))) { y <- as.numeric(M[i, mcols]); f <- nnls_pg(X, y); out$n_markers[i] <- f$n; out$rss[i] <- round(f$rss, 3); out$raw_sum[i] <- round(f$raw_sum, 3); out[i, cts] <- round(as.numeric(f$w), 4) }
out$epithelial <- rowSums(out[, intersect(c("Head-Neck-Ep", "Epid-Kerat"), cts), drop = FALSE]); out$blood <- rowSums(out[, grep("^Blood", cts, value = TRUE), drop = FALSE])
# two-compartment check from marker means at Head-Neck-Ep-target and blood-target markers
tg <- A$target; hn <- which(tg == "Head-Neck-Ep"); bl <- which(tg %in% grep("^Blood", cts, value = TRUE))
epi_ref <- mean(X[hn, "Head-Neck-Ep"]); blood_ref_at_hn <- mean(rowMeans(X[hn, grep("^Blood", cts, value = TRUE)]))
u_hn <- apply(M[, mcols], 1, function(v) mean(as.numeric(v[hn]), na.rm = TRUE))
zero <- mean(u_hn[M$dna_blood == 1])   # empirical zero: blood genomes measured on this platform
out$epi_2comp <- pmax(0, pmin(1, (u_hn - zero) / (epi_ref - zero)))
out$u_HeadNeckEp <- round(u_hn, 4)
gr <- which(tg == "Blood-Granul"); u_gr <- apply(M[, mcols], 1, function(v) mean(as.numeric(v[gr]), na.rm = TRUE)); out$granul_2comp <- round(pmax(0, pmin(1, u_gr / mean(X[gr, "Blood-Granul"]))), 3)
write.table(out, file.path(Sys.getenv("DATA_ROOT"), "meth", "deconvolution_U25.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("atlas refs at Head-Neck-Ep markers: epithelium", round(epi_ref, 3), " blood mean", round(blood_ref_at_hn, 3), "\n")
cat("by DNA source: median NNLS epithelial, NNLS blood, 2-compartment epithelial, raw_sum, rss\n")
print(aggregate(cbind(epithelial, blood, epi_2comp, granul_2comp, raw_sum, rss) ~ dna_blood, data = out, FUN = function(x) round(median(x), 3)))
cat("blood genomes:\n"); print(out[out$dna_blood == 1, c("sample", "epithelial", "blood", "epi_2comp", "Blood-Granul", "Blood-T", "Blood-Mono+Macro", "raw_sum")])
cat("saliva epithelial quantiles (NNLS / 2comp):\n"); print(round(quantile(out$epithelial[out$dna_blood == 0], c(.05,.25,.5,.75,.95)), 3)); print(round(quantile(out$epi_2comp[out$dna_blood == 0], c(.05,.25,.5,.75,.95)), 3))
cat("cor(NNLS epithelial, 2comp) =", round(cor(out$epithelial, out$epi_2comp), 3), "\n")
