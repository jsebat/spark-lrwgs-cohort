#!/usr/bin/env python3
"""Stage 2b (D74 addendum 6). Stdlib only; bcftools via container.
2b-1 LOI v2: cohort-relative screen at robust imprinted-like islands.
2b-2 Gene annotation (GENCODE v44, +/- 50 kb) of imprinted-like islands not near catalogued imprinted genes.
2b-3 cis genotype vs methylation in founders (n = 68): SNVs +/- 10 kb (MAF >= 0.05), SVs overlapping +/- 10 kb,
     TR summed allele length for TRGT loci overlapping the island +/- 1 kb; OLS with epithelial fraction and DNA source.
2b-4 composition vs phenotype.
"""
import os, csv, gzip, re, math, subprocess, collections, bisect
ME = os.environ["DATA_ROOT"] + "/meth/"; MAN = os.environ["MANIFEST"]
BCF = os.environ["COHORT_BCF"]; SVV = os.environ["COHORT_SV_VCF"]; TRV = os.environ["COHORT_TRGT_VCF"]
GTF = os.environ["GENCODE_GTF"]
SIF = os.environ["SIF"]; H = os.path.expanduser("~"); sx = ["singularity", "exec", "-B", "/expanse:/expanse", "-B", H + ":" + H, SIF]
TAB = "\t"; NL = "\n"
man = {r["sample_id"]: r for r in csv.DictReader(open(MAN), delimiter=TAB)}
def run(cmd): return subprocess.run(cmd, stdout=subprocess.PIPE, universal_newlines=True).stdout
def num(v):
    try: return float(v)
    except (TypeError, ValueError): return None
def load(name):
    with open(ME + "island_%s.tsv" % name) as fh:
        rd = csv.reader(fh, delimiter=TAB); hdr = next(rd); return hdr[4:], [r for r in rd]
samples, CB = load("combined"); _, H1 = load("hap1"); _, H2 = load("hap2")
islands = [(r[0], int(r[1]), int(r[2]), r[3]) for r in CB]
si = {s: i for i, s in enumerate(samples)}
cat = {(r["chrom"], int(r["start"]), int(r["end"])): r for r in csv.DictReader(open(ME + "asm_island_catalog.tsv"), delimiter=TAB)}
dec = {r["sample"]: r for r in csv.DictReader(open(ME + "deconvolution_U25.tsv"), delimiter=TAB)}
epi = {s: num(dec[s]["epi_2comp"]) if s in dec else None for s in samples}
BLOOD_PREFIXES = tuple(x for x in os.environ.get("BLOOD_SAMPLE_PREFIX", "REACH").split(",") if x)   # see meth_stage1.py (X8)
blood = {s: 1.0 if s.startswith(BLOOD_PREFIXES) else 0.0 for s in samples}
founders = [s for s in samples if man[s]["father_id"] == "0" and man[s]["mother_id"] == "0"]
print("samples", len(samples), "founders", len(founders), flush=True)

def median(xs): xs = sorted(xs); n = len(xs); return xs[n // 2] if n % 2 else 0.5 * (xs[n // 2 - 1] + xs[n // 2])
def ols(y, X):
    """OLS with intercept; X = list of covariate lists (per sample). returns (beta1, se1, p1, n) for the first covariate; two-sided t via normal approx if n>30 else t (df) approx"""
    n = len(y); k = len(X[0]) + 1
    A = [[1.0] + list(x) for x in X]
    keep = [0] + [a for a in range(1, len(A[0])) if len(set(round(r[a], 9) for r in A)) > 1]
    A = [[r[a] for a in keep] for r in A]; k = len(keep)
    # normal equations
    XtX = [[sum(A[i][a] * A[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
    Xty = [sum(A[i][a] * y[i] for i in range(n)) for a in range(k)]
    # solve via Gaussian elimination
    M = [row[:] + [Xty[r]] for r, row in enumerate(XtX)]
    for c in range(k):
        piv = max(range(c, k), key=lambda r: abs(M[r][c]))
        if abs(M[piv][c]) < 1e-12: return None
        M[c], M[piv] = M[piv], M[c]
        for r in range(k):
            if r != c:
                f = M[r][c] / M[c][c]
                for cc in range(c, k + 1): M[r][cc] -= f * M[c][cc]
    beta = [M[r][k] / M[r][r] for r in range(k)]
    resid = [y[i] - sum(A[i][a] * beta[a] for a in range(k)) for i in range(n)]
    df = n - k
    if df < 5: return None
    s2 = sum(r * r for r in resid) / df
    # inverse of XtX for var(beta1): re-solve with identity columns (cheap for small k)
    def inv_col(j):
        Mi = [row[:] + [1.0 if r == j else 0.0] for r, row in enumerate(XtX)]
        for c in range(k):
            piv = max(range(c, k), key=lambda r: abs(Mi[r][c])); Mi[c], Mi[piv] = Mi[piv], Mi[c]
            for r in range(k):
                if r != c:
                    f = Mi[r][c] / Mi[c][c]
                    for cc in range(c, k + 1): Mi[r][cc] -= f * Mi[c][cc]
        return [Mi[r][k] / Mi[r][r] for r in range(k)]
    v11 = inv_col(1)[1] * s2
    if v11 <= 0: return None
    se = math.sqrt(v11); t = beta[1] / se
    # two-sided p from t with df (Student t via regularized incomplete beta, simple numeric)
    def t_p(t, df):
        x = df / (df + t * t)
        # regularized incomplete beta I_x(df/2, 1/2) via continued fraction (Lentz)
        a, b = df / 2.0, 0.5
        def betacf(a, b, x):
            MAXIT, EPS = 200, 3e-14; qab = a + b; qap = a + 1; qam = a - 1; c = 1.0; d = 1 - qab * x / qap
            d = 1.0 / (d if abs(d) > 1e-300 else 1e-300); h = d
            for m in range(1, MAXIT + 1):
                m2 = 2 * m; aa = m * (b - m) * x / ((qam + m2) * (a + m2)); d = 1 + aa * d; d = 1.0 / (d if abs(d) > 1e-300 else 1e-300); c = 1 + aa / (c if abs(c) > 1e-300 else 1e-300); h *= d * c
                aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2)); d = 1 + aa * d; d = 1.0 / (d if abs(d) > 1e-300 else 1e-300); c = 1 + aa / (c if abs(c) > 1e-300 else 1e-300); de = d * c; h *= de
                if abs(de - 1) < EPS: break
            return h
        lbeta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
        if x <= 0: return 0.0
        if x >= 1: return 1.0
        bt = math.exp(math.log(x) * a + math.log(1 - x) * b - lbeta)
        if x < (a + 1) / (a + b + 2): return bt * betacf(a, b, x) / a
        return 1 - bt * betacf(b, a, 1 - x) / b
    p = t_p(t, df)
    return (beta[1], se, p, n)

# ---------------- 2b-1 LOI v2
robust = [k for k, r in cat.items() if r["class"].startswith("imprinted_like") and int(r["n_children_ASM"]) >= 15 and int(r["n_parents_with_data"]) > 0 and int(r["n_parents_ASM"]) / int(r["n_parents_with_data"]) >= 0.6]
idx_of = {(c, a, b): i for i, (c, a, b, _) in enumerate(islands)}
with open(ME + "loi_v2.tsv", "w") as fh:
    fh.write(TAB.join(["chrom", "start", "end", "class", "genes", "sample", "role", "combined", "cohort_median", "robust_z", "hap1", "hap2", "pattern"]) + NL); n = 0
    for k in robust:
        i = idx_of[k]; vals = [num(CB[i][4 + si[s]]) for s in samples]; v = [x for x in vals if x is not None]
        if len(v) < 30: continue
        med = median(v); mad = median([abs(x - med) for x in v]) * 1.4826
        if mad < 0.02: mad = 0.02
        for s in samples:
            cb = num(CB[i][4 + si[s]]); h1 = num(H1[i][4 + si[s]]); h2 = num(H2[i][4 + si[s]])
            if cb is None or h1 is None or h2 is None: continue
            if h1 == 0 and h2 == 0 and cb > 0.3: continue           # no phased CpGs
            z = (cb - med) / mad
            if abs(z) > 3 and abs(h1 - h2) < 0.15:
                fh.write(TAB.join(map(str, [k[0], k[1], k[2], cat[k]["class"], cat[k]["imprinted_genes_within_50kb"], s, man[s]["role"], round(cb, 3), round(med, 3), round(z, 2), round(h1, 3), round(h2, 3), "gain" if z > 0 else "loss"])) + NL); n += 1
print("LOI v2: robust islands", len(robust), " candidates", n, flush=True)

# ---------------- 2b-2 gene annotation of imprinted-like islands
genes_by_chr = collections.defaultdict(list)
with gzip.open(GTF, "rt") as fh:
    for ln in fh:
        if ln.startswith("#"): continue
        f = ln.split(TAB)
        if f[2] != "gene" or 'gene_type "protein_coding"' not in f[8] and 'gene_type "lncRNA"' not in f[8]: continue
        m = re.search(r'gene_name "([^"]+)"', f[8]); genes_by_chr[f[0]].append((int(f[3]), int(f[4]), m.group(1) if m else ""))
with open(ME + "imprinted_like_annotation.tsv", "w") as fh:
    fh.write(TAB.join(["chrom", "start", "end", "class", "n_children_ASM", "n_mat", "n_pat", "catalogued_imprinted_genes_50kb", "genes_within_50kb", "nearest_gene"]) + NL)
    for k, r in cat.items():
        if not r["class"].startswith("imprinted_like"): continue
        c, a, b = k; near = [(max(0, max(gs - b, a - ge)), gn) for gs, ge, gn in genes_by_chr.get(c, []) if gs - 50000 <= b and ge + 50000 >= a]
        near.sort(); fh.write(TAB.join(map(str, [c, a, b, r["class"], r["n_children_ASM"], r["n_maternal_methylated"], r["n_paternal_methylated"], r["imprinted_genes_within_50kb"], ";".join(gn for _, gn in near[:8]), near[0][1] if near else ""])) + NL)
print("annotation written", flush=True)

# ---------------- 2b-3 cis tests in founders
fidx = [si[s] for s in founders]
cov = [[epi[s] if epi[s] is not None else 0.1, blood[s]] for s in founders]
# islands to test: cohort SD > 0.05 among founders OR in ASM catalogue with class imprinted/frequent/sequence_dependent
test_idx = []
for i, k4 in enumerate(islands):
    c, a, b, _ = k4
    if c in ("chrX", "chrY", "chrM"): continue
    y = [num(CB[i][4 + j]) for j in fidx]
    yy = [v for v in y if v is not None]
    if len(yy) < 50: continue
    m = sum(yy) / len(yy); sd = math.sqrt(sum((v - m) ** 2 for v in yy) / (len(yy) - 1))
    klass = cat.get((c, a, b), {}).get("class", "")
    if sd > 0.05 or klass.startswith(("imprinted_like", "sequence_dependent", "frequent")): test_idx.append(i)
print("islands to test:", len(test_idx), flush=True)
def write_bed(idxs, pad, path):
    with open(path, "w") as fh:
        for i in idxs: c, a, b, _ = islands[i]; fh.write(TAB.join([c, str(max(0, a - pad)), str(b + pad)]) + NL)
by_chr = collections.defaultdict(list)
for i in test_idx: c, a, b, _ = islands[i]; by_chr[c].append((a, b, i))
for c in by_chr: by_chr[c].sort()
starts = {c: [a for a, b, i in v] for c, v in by_chr.items()}
def islands_near(c, pos, pad):
    lst = by_chr.get(c, []); out = []
    if not lst: return out
    k = bisect.bisect_right(starts[c], pos + pad)
    j = k - 1
    while j >= 0 and lst[j][0] >= pos - pad - 60000:
        a, b, i = lst[j]
        if a - pad <= pos <= b + pad: out.append(i)
        j -= 1
    return out
def gt_add(g):
    g = g.replace("|", "/")
    if "." in g: return None
    return sum(1 for x in g.split("/") if x != "0")
results = {"snv": [], "sv": [], "tr": []}
# SNVs
write_bed(test_idx, 10000, ME + "_cis_snv.bed")
out = subprocess.run(sx + ["bcftools", "query", "-R", ME + "_cis_snv.bed", "-s", ",".join(founders), "-i", 'TYPE="snp" && QUAL>30', "-f", "%CHROM" + TAB + "%POS[" + TAB + "%GT]" + NL, BCF], stdout=subprocess.PIPE).stdout
best = {}
ntest = {}   # per-island count of SNVs tested, for the within-island Bonferroni (X6)
n_tests = 0
for ln in out.decode("latin-1").splitlines():
    f = ln.split(TAB)
    if len(f) < 2 + len(founders): continue
    c, pos = f[0], int(f[1]); g = [gt_add(x) for x in f[2:2 + len(founders)]]
    ok = [j for j, v in enumerate(g) if v is not None]
    if len(ok) < 50: continue
    af = sum(g[j] for j in ok) / (2 * len(ok))
    if af < 0.05 or af > 0.95: continue
    for i in islands_near(c, pos, 10000):
        y = []; X = []
        for j in ok:
            v = num(CB[i][4 + fidx[j]])
            if v is None: continue
            y.append(v); X.append([g[j]] + cov[j])
        if len(y) < 50: continue
        r = ols(y, X); n_tests += 1
        if r: ntest[i] = ntest.get(i, 0) + 1
        if r and (i not in best or r[2] < best[i][2]): best[i] = (pos, r[2], r[0], af, len(y))
with open(ME + "cis_snv_best.tsv", "w") as fh:
    fh.write(TAB.join(["chrom", "start", "end", "class", "best_snv_pos", "beta_per_allele", "p", "maf", "n", "fdr_bh"]) + NL)
    # the minimum p over the SNVs tested for an island is not a p-value: Bonferroni within the island first, then
    # Benjamini-Hochberg across islands with the step-up monotonicity (the TR block already had it) (X6)
    padj = {i: min(1.0, v[1] * ntest.get(i, 1)) for i, v in best.items()}
    items = sorted(best.items(), key=lambda kv: padj[kv[0]]); m = len(items)
    q = [0.0] * m; running = 1.0
    for rank in range(m, 0, -1):
        running = min(running, padj[items[rank - 1][0]] * m / rank); q[rank - 1] = running
    for rank, (i, (pos, p, beta, af, n)) in enumerate(items, 1):
        c, a, b, _ = islands[i]; fh.write(TAB.join(map(str, [c, a, b, cat.get((c, a, b), {}).get("class", ""), pos, round(beta, 4), "%.3g" % padj[i], round(af, 3), n, "%.3g" % q[rank - 1]])) + NL)
print("cis SNV: islands with a tested SNV", len(best), " tests", n_tests, " FDR<0.05 (island-Bonferroni then BH):", sum(1 for x in q if x < 0.05), flush=True)
# SVs
write_bed(test_idx, 10000, ME + "_cis_sv.bed")
out = run(sx + ["bcftools", "query", "-R", ME + "_cis_sv.bed", "-s", ",".join(founders), "-i", 'INFO/SVTYPE!="BND"', "-f", "%CHROM" + TAB + "%POS" + TAB + "%INFO/END" + TAB + "%INFO/SVTYPE" + TAB + "%INFO/SVLEN[" + TAB + "%GT]" + NL, SVV])
with open(ME + "cis_sv.tsv", "w") as fh:
    fh.write(TAB.join(["chrom", "start", "end", "class", "sv_pos", "sv_end", "svtype", "svlen", "n_carriers", "beta_carrier", "p", "n"]) + NL); nsv = 0
    for ln in out.splitlines():
        f = ln.split(TAB)
        if len(f) < 5 + len(founders): continue
        c, pos = f[0], int(f[1]); end = int(f[2]) if f[2] not in (".", "") else pos; g = [gt_add(x) for x in f[5:5 + len(founders)]]
        ok = [j for j, v in enumerate(g) if v is not None]; carriers = sum(1 for j in ok if g[j] > 0)
        if carriers < 3 or carriers > len(ok) - 3: continue
        for i in sorted(set(islands_near(c, pos, 10000) + islands_near(c, end, 10000))):
            y = []; X = []
            for j in ok:
                v = num(CB[i][4 + fidx[j]])
                if v is None: continue
                y.append(v); X.append([1.0 if g[j] > 0 else 0.0] + cov[j])
            if len(y) < 50: continue
            r = ols(y, X)
            if r and r[2] < 0.01:
                a, b = islands[i][1], islands[i][2]; fh.write(TAB.join(map(str, [c, a, b, cat.get((c, a, b), {}).get("class", ""), pos, end, f[3], f[4], carriers, round(r[0], 4), "%.3g" % r[2], r[3]])) + NL); nsv += 1
print("cis SV associations p<0.01:", nsv, flush=True)
# TRs: loci overlapping island +/- 1 kb
write_bed(test_idx, 1000, ME + "_cis_tr.bed")
out = run(sx + ["bcftools", "query", "-R", ME + "_cis_tr.bed", "-s", ",".join(founders), "-f", "%CHROM" + TAB + "%POS" + TAB + "%INFO/END" + TAB + "%INFO/TRID" + TAB + "%INFO/MOTIFS[" + TAB + "%AL:%SD]" + NL, TRV])
def rank(v):
    s = sorted(range(len(v)), key=lambda i: v[i]); r = [0] * len(v); i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and v[s[j + 1]] == v[s[i]]: j += 1
        for k in range(i, j + 1): r[s[k]] = (i + j) / 2 + 1
        i = j + 1
    return r
def spearman(x, y):
    rx, ry = rank(x), rank(y); n = len(x); mx = sum(rx) / n; my = sum(ry) / n
    sxx = sum((a - mx) ** 2 for a in rx); syy = sum((b - my) ** 2 for b in ry)
    return sum((a - mx) * (b - my) for a, b in zip(rx, ry)) / math.sqrt(sxx * syy) if sxx and syy else 0.0
with open(ME + "cis_tr_all.tsv", "w") as fh:
    fh.write(TAB.join(["chrom", "start", "end", "class", "trid", "motifs", "n", "spearman_rho_sumlen", "beta_per_bp_sumlen", "p_sumlen", "beta_per_bp_maxlen", "p_maxlen", "sumlen_range"]) + NL); ntr = 0; nall = 0
    for ln in out.splitlines():
        f = ln.split(TAB)
        if len(f) < 5 + len(founders): continue
        c, pos, end, trid, motifs = f[0], int(f[1]), int(f[2]) if f[2] not in (".", "") else int(f[1]), f[3], f[4]
        cells = f[5:5 + len(founders)]
        for i in sorted(set(islands_near(c, pos, 1000) + islands_near(c, end, 1000))):
            y = []; Xs = []; Xm = []; sl = []
            for j, cell in enumerate(cells):
                al, sd = (cell.split(":") + ["", ""])[:2]
                if al in (".", ""): continue
                try:
                    L = [int(float(x)) for x in al.split(",")]; S = [int(float(x)) for x in sd.split(",")]
                except ValueError: continue
                if min(S) < 5: continue
                v = num(CB[i][4 + fidx[j]])
                if v is None: continue
                y.append(v); Xs.append([float(sum(L))] + cov[j]); Xm.append([float(max(L))] + cov[j]); sl.append(sum(L))
            if len(y) < 40 or len(set(sl)) < 3: continue
            nall += 1
            rs = ols(y, Xs); rm = ols(y, Xm); rho = spearman(sl, y)
            if rs and (rs[2] < 0.01 or abs(rho) > 0.4):     # the reported set is what the message says it is (X19)
                a, b = islands[i][1], islands[i][2]; fh.write(TAB.join(map(str, [c, a, b, cat.get((c, a, b), {}).get("class", ""), trid, motifs, len(y), round(rho, 3), round(rs[0], 5) if rs else "", "%.3g" % rs[2] if rs else "", round(rm[0], 5) if rm else "", "%.3g" % rm[2] if rm else "", "%d-%d" % (min(sl), max(sl))])) + NL); ntr += 1
print("TR x island pairs tested:", nall, " reported (p<0.01 or |rho|>0.4):", ntr, flush=True)

_rows = [l.rstrip(NL).split(TAB) for l in open(ME + "cis_tr_all.tsv")][1:]
_seen = set(); _u = []
for r in _rows:
    key = (r[0], r[1], r[2], r[4])
    if key in _seen: continue
    _seen.add(key); _u.append(r)
_pv = [(float(r[9]), i) for i, r in enumerate(_u) if r[9] not in ("", "None")]; _pv.sort(); m = len(_pv); q = [""] * len(_u); prev = 1.0
for rank_, (p, i) in reversed(list(enumerate(_pv, 1))):
    prev = min(prev, p * m / rank_); q[i] = "%.3g" % prev
with open(ME + "cis_tr.tsv", "w") as fh:
    fh.write(TAB.join(["chrom", "start", "end", "class", "trid", "motifs", "n", "spearman_rho_sumlen", "beta_per_bp_sumlen", "p_sumlen", "q_sumlen_BH", "beta_per_bp_maxlen", "p_maxlen", "sumlen_range"]) + NL)
    sel = [(float(r[9]) if r[9] not in ("", "None") else 1.0, i) for i, r in enumerate(_u) if (r[9] not in ("", "None") and float(r[9]) < 0.01) or abs(float(r[7])) > 0.4]
    for p, i in sorted(sel): r = _u[i]; fh.write(TAB.join(r[:10] + [q[i]] + r[10:]) + NL)
print("TR pairs unique:", len(_u), " with p:", m, " BH q<0.05:", sum(1 for v in q if v not in ("",) and float(v) < 0.05), " q<0.10:", sum(1 for v in q if v not in ("",) and float(v) < 0.10), " reported:", len(sel), flush=True)
# ---------------- 2b-4 composition vs phenotype
with open(ME + "composition_vs_phenotype.tsv", "w") as fh:
    fh.write("group\tn\tmedian_epithelial\n")
    for lab, sel in (("affected_offspring", [s for s in samples if man[s]["role"] == "offspring" and man[s]["affected"] == "affected"]), ("unaffected_offspring", [s for s in samples if man[s]["role"] == "offspring" and man[s]["affected"] != "affected"]), ("parents", [s for s in samples if man[s]["role"] == "parent"]), ("blood", [s for s in samples if s.startswith(BLOOD_PREFIXES)]), ("saliva", [s for s in samples if not s.startswith(BLOOD_PREFIXES)])):
        v = [epi[s] for s in sel if epi[s] is not None]; fh.write("%s\t%d\t%s\n" % (lab, len(v), round(median(v), 3) if v else ""))
print(open(ME + "composition_vs_phenotype.tsv").read()); print("stage 2b done")
