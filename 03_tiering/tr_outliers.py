#!/usr/bin/env python3
"""Cohort-level tandem-repeat analysis on the joint TRGT callset (streamed; stdlib + numpy).

Rules are fixed before running (D70):
  allele length      = total motif count of the allele (sum of MC parts), in motif units
  QC per genotype    = every allele has spanning reads SD >= MIN_SD and purity AP >= MIN_AP
  founder reference  = the 65 unaffected founders (130 alleles) at each locus, QC-passed only;
                       a locus needs >= MIN_FOUNDER_ALLELES founder alleles to be scored
  OUTLIER allele     = length > founder 99th percentile AND > founder p99 + MIN_EXCESS units
  DE NOVO expansion  = affected offspring allele > max(parental alleles) + MIN_EXCESS,
                       both parents QC-passed, and > founder maximum (so it is also an outlier)
  TDT                = informative meiosis: exactly one parent carries an outlier allele
                       (heterozygous: the other allele is not an outlier), other parent none;
                       transmitted if the affected child carries an outlier allele.
Outputs (OUT dir): tr_locus_stats.tsv, tr_outlier_calls.tsv, tr_denovo_expansions.tsv,
tr_persample.tsv, tr_tdt.tsv, tr_summary.txt.
"""
import os, sys, subprocess, collections, math

def pct(sorted_vals, q):
    """nearest-rank percentile on a sorted list"""
    k = max(0, min(len(sorted_vals) - 1, int(math.ceil(q / 100.0 * len(sorted_vals))) - 1))
    return sorted_vals[k]

VCF = os.environ["TRGT_VCF"]
PED = os.environ["PED"]
OUT = os.environ["OUT"]
CDS = os.environ.get("CDS_BED", "")            # optional: chrom start end gene  (0-based BED)
SIF = os.environ["SIF"]
MIN_SD, MIN_AP = int(os.environ.get("TR_MIN_SPANNING", "5")), float(os.environ.get("TR_MIN_PURITY", "0.80"))
MIN_EXCESS = int(os.environ.get("TR_MIN_EXCESS", "3"))
MIN_FOUNDER_ALLELES = int(os.environ.get("TR_MIN_FOUNDER_ALLELES", "100"))
os.makedirs(OUT, exist_ok=True)

# ---------------- pedigree
aff, pat, mat, fam = {}, {}, {}, {}
for ln in open(PED):
    f = ln.rstrip("\n").split("\t")
    if len(f) < 6 or f[0].startswith("#"): continue
    fam[f[1]], pat[f[1]], mat[f[1]], aff[f[1]] = f[0], f[2], f[3], (f[5] == "2")
sx = ["singularity", "exec", "-B", "/expanse:/expanse", SIF]
samples = subprocess.run(sx + ["bcftools", "query", "-l", VCF], stdout=subprocess.PIPE, universal_newlines=True, check=True).stdout.split()
si = {s: i for i, s in enumerate(samples)}
founders_unaff = [s for s in samples if pat.get(s, "0") == "0" and mat.get(s, "0") == "0" and s in aff and not aff[s]]
offspring = [s for s in samples if pat.get(s, "0") != "0" and mat.get(s, "0") != "0" and pat[s] in si and mat[s] in si]
aff_off = [s for s in offspring if aff.get(s)]
print("samples %d; unaffected founders %d; offspring with both parents %d (affected %d)" % (len(samples), len(founders_unaff), len(offspring), len(aff_off)), flush=True)
fidx = [si[s] for s in founders_unaff]

# ---------------- optional CDS annotation (interval index per chrom)
cds = collections.defaultdict(list)
if CDS and os.path.exists(CDS):
    for ln in open(CDS):
        f = ln.split("\t")
        if len(f) >= 3: cds[f[0]].append((int(f[1]), int(f[2]), f[3].strip() if len(f) > 3 else ""))
    for c in cds: cds[c].sort()
    print("CDS intervals loaded:", sum(len(v) for v in cds.values()), flush=True)
import bisect
def cds_hit(chrom, s, e):
    iv = cds.get(chrom)
    if not iv: return ""
    i = bisect.bisect_left(iv, (e, 0, "")); genes = set()
    j = i
    while j >= 0 and j < len(iv):
        if iv[j][0] < e and iv[j][1] > s: genes.add(iv[j][2])
        j -= 1
        if j >= 0 and iv[j][1] < s - 200000: break
    return ",".join(sorted(g for g in genes if g))

def parse_alleles(cell):
    """cell = 'AL|SD' -> list of (length_bp, ok). Length = AL (whole-region allele length in bp).
    QC = spanning reads only: in this catalog a locus is a broad multi-motif region, so allele purity (AP)
    is ~0.4-0.5 by construction and cannot be used as a filter (D70 addendum). MC undercounts region length."""
    al, sd = (cell.split("|") + ["", ""])[:2]
    if al in ("", "."): return []
    out = []
    sds = sd.split(",")
    for k, a in enumerate(al.split(",")):
        try: L = int(float(a))
        except ValueError: continue
        try: s_ok = int(float(sds[k])) >= MIN_SD
        except (ValueError, IndexError): s_ok = False
        out.append((L, s_ok))
    return out

def excess_bp(motifs):
    """minimum excess in bp = 3 motif units of the shortest motif, at least 6 bp"""
    ml = [len(m) for m in motifs.split(",") if m and m != "."]
    return max(6, MIN_EXCESS * (min(ml) if ml else 2))

fmt = "%INFO/TRID\\t%CHROM\\t%POS\\t%INFO/END\\t%INFO/MOTIFS[\\t%AL|%SD]\\n"
proc = subprocess.Popen(sx + ["bcftools", "query", "-f", fmt, VCF], stdout=subprocess.PIPE, universal_newlines=True, bufsize=1 << 20)

per_called = collections.Counter(); per_out = collections.Counter(); per_out_cds = collections.Counter()
locus_f = open(os.path.join(OUT, "tr_locus_stats.tsv"), "w"); locus_f.write("trid\tchrom\tpos\tend\tmotifs\tcds_genes\tn_founder_alleles\tfounder_p50\tfounder_p99\tfounder_max\tthreshold\tn_outlier_alleles_cohort\tn_outlier_carriers\n")
out_f = open(os.path.join(OUT, "tr_outlier_calls.tsv"), "w"); out_f.write("trid\tchrom\tpos\tcds_genes\tsample\trole\taffected\tallele_len\tthreshold\tfounder_max\tin_founders\n")
dn_f = open(os.path.join(OUT, "tr_denovo_expansions.tsv"), "w"); dn_f.write("trid\tchrom\tpos\tend\tmotifs\tcds_genes\tchild\tfamily\taffected\tchild_alleles\tfather_alleles\tmother_alleles\tfounder_max\texcess_over_parents\n")
tdt = collections.Counter()   # (cds_flag) -> T / NT
n_loci = n_scored = 0
for line in proc.stdout:
    f = line.rstrip("\n").split("\t")
    if len(f) < 5 + len(samples): continue
    n_loci += 1
    trid, chrom, pos, end, motifs = f[0], f[1], int(f[2]), f[3], f[4]
    cells = f[5:5 + len(samples)]
    G = [parse_alleles(c) for c in cells]
    # founder distribution (QC-passed alleles only)
    fa = [L for i in fidx for (L, ok) in G[i] if ok]
    if len(fa) < MIN_FOUNDER_ALLELES: continue
    n_scored += 1
    fa.sort(); p50 = float(pct(fa, 50)); p99 = float(pct(fa, 99)); fmax = int(fa[-1])
    exc = excess_bp(motifs)
    thr = p99 + exc
    # leave-one-family-out thresholds: a sample is never compared with a reference containing its own family's
    # founders (otherwise founders sit inside their own reference and offspring do not -> biased offspring counts)
    fam_alleles = collections.defaultdict(list)
    for i in fidx:
        for (L, ok) in G[i]:
            if ok: fam_alleles[fam.get(samples[i], "")].append(L)
    loo_thr = {}; loo_max = {}
    for fm, vals in fam_alleles.items():
        rest = list(fa)
        for v in vals: rest.remove(v)
        if len(rest) >= MIN_FOUNDER_ALLELES - 10:
            loo_thr[fm] = pct(rest, 99) + exc; loo_max[fm] = rest[-1]
    def thr_for(s): return loo_thr.get(fam.get(s, ""), thr)
    def fmax_for(s): return loo_max.get(fam.get(s, ""), fmax)
    genes = cds_hit(chrom, pos - 1, int(end) if end not in ("", ".") else pos) if cds else ""
    # per-sample outlier status
    carriers = []; n_out_alleles = 0
    is_out = {}
    for s in samples:
        al = G[si[s]]
        if not al or not all(ok for _, ok in al):
            is_out[s] = None; continue
        per_called[s] += 1
        t_s = thr_for(s)
        o = [L for L, _ in al if L > t_s]
        is_out[s] = bool(o)
        if o:
            n_out_alleles += len(o); carriers.append(s); per_out[s] += 1
            if genes: per_out_cds[s] += 1
            role = "founder" if pat.get(s, "0") == "0" else "offspring"
            out_f.write("\t".join(map(str, [trid, chrom, pos, genes, s, role, "yes" if aff.get(s) else "no", max(o), round(t_s, 1), fmax_for(s), "yes" if s in founders_unaff else "no"])) + "\n")
    locus_f.write("\t".join(map(str, [trid, chrom, pos, end, motifs, genes, len(fa), round(p50, 1), round(p99, 1), fmax, round(thr, 1), n_out_alleles, len(carriers)])) + "\n")
    # de novo expansions and TDT in complete trios with an affected child
    for c in aff_off:
        ca, pa, ma = G[si[c]], G[si[pat[c]]], G[si[mat[c]]]
        if not (ca and pa and ma) or not all(ok for _, ok in ca + pa + ma): continue
        cl = [L for L, _ in ca]; pl = [L for L, _ in pa]; ml = [L for L, _ in ma]
        pmax = max(pl + ml); t_c = thr_for(c); fm_c = fmax_for(c)
        if max(cl) > pmax + exc and max(cl) > fm_c:
            dn_f.write("\t".join(map(str, [trid, chrom, pos, end, motifs, genes, c, fam.get(c, ""), "yes", ",".join(map(str, cl)), ",".join(map(str, pl)), ",".join(map(str, ml)), fm_c, max(cl) - pmax])) + "\n")
        # TDT on outlier alleles (family-LOO threshold)
        if chrom in ("chrX", "chrY"):
            continue                          # a father is hemizygous there: his single allele is not a het transmission test (T7)
        po = [L > t_c for L in pl]; mo = [L > t_c for L in ml]
        if sum(po) + sum(mo) == 1:           # exactly one outlier allele among the four parental alleles
            transmitted = any(L > t_c for L in cl)
            tdt[("cds" if genes else "noncds", "T" if transmitted else "NT")] += 1
            tdt[("all", "T" if transmitted else "NT")] += 1
proc.wait()
if proc.returncode != 0:
    sys.exit("FATAL: bcftools query exited %d -- the outlier tables above are truncated, do not read them as results" % proc.returncode)
for fh in (locus_f, out_f, dn_f): fh.close()

# ---------------- per-sample table and TDT
with open(os.path.join(OUT, "tr_persample.tsv"), "w") as fh:
    fh.write("sample\tfamily\trole\taffected\tloci_called\toutlier_loci\toutlier_cds_loci\toutlier_per_10k\n")
    for s in samples:
        role = "founder" if pat.get(s, "0") == "0" else "offspring"
        n = per_called[s]
        fh.write("\t".join(map(str, [s, fam.get(s, ""), role, 1 if aff.get(s) else 0, n, per_out[s], per_out_cds[s], round(1e4 * per_out[s] / n, 3) if n else ""])) + "\n")
with open(os.path.join(OUT, "tr_tdt.tsv"), "w") as fh:
    fh.write("set\tT\tNT\tT_frac\tp_binom\n")
    for k in ("all", "cds", "noncds"):
        T, NT = tdt[(k, "T")], tdt[(k, "NT")]; n = T + NT
        if n:
            # two-sided exact binomial
            from math import comb
            lo = sum(comb(n, i) for i in range(0, min(T, NT) + 1)) / 2 ** n
            p = min(1.0, 2 * lo)
            fh.write("%s\t%d\t%d\t%.4f\t%.3g\n" % (k, T, NT, T / n, p))
with open(os.path.join(OUT, "tr_summary.txt"), "w") as fh:
    fh.write("loci in callset %d; loci scored (>= %d founder alleles) %d\n" % (n_loci, MIN_FOUNDER_ALLELES, n_scored))
    fh.write("rules: length = AL (bp); QC = SD>=%d per allele (AP not used: catalog regions are multi-motif, AP ~0.5 by construction); outlier > founder p99 + max(6 bp, %d x shortest motif); de novo > max(parents) + same excess and > founder max\n" % (MIN_SD, MIN_EXCESS))
    fh.write("TDT: %s\n" % dict(tdt))
print(open(os.path.join(OUT, "tr_summary.txt")).read(), flush=True)
