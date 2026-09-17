"""Transmission test on inherited LoF tier burden, long-read cohort.

Companion to the clogit. For each rare tiered LoF site, in each complete trio with
an AFFECTED offspring, count transmissions from HETEROZYGOUS parents:
    T  = het parent transmitted the alt allele to the affected child
    NT = het parent did not transmit
Under the null T/(T+NT) = 0.5. Over-transmission of a tier means burden enrichment.

Why this and not only clogit: the clogit used the parents as controls, so the control
values are not independent of the case (a child's genome is half of each parent's) and
the contrast is confounded with generation. The TDT compares transmitted against
non-transmitted alleles WITHIN the same parent, which removes both problems.

Only informative meioses count: parent het, genotypes called and QC-passing in
parent and child. Homozygous-reference parents contribute nothing.
"""
import os
import subprocess
import sys
import collections
import math

T = "/expanse/lustre/projects/ddp195/jsebat/longread-autism/tiering"
R = "/expanse/lustre/projects/ddp195/jsebat/longread-autism"
BCF = "/expanse/projects/sebat1/longread_cohort_2026/freeze1/freeze1.cohort.bcf"
TAB, NL = chr(9), chr(10)
GQ_MIN, DP_MIN = 20, 10
SIF = subprocess.Popen("ls /expanse/projects/sebat1/jsebat/singularity_cache/miniwdl/*pb_wdl_base*.sif",
                       shell=True, stdout=subprocess.PIPE).communicate()[0].decode().split()[0]

# ---- tiered rare sites
tier = {}
with open(os.path.join(T, "lof_tiered.rare.tsv")) as fh:
    hdr = fh.readline().rstrip(NL).split(TAB)
    ix = dict((h, i) for i, h in enumerate(hdr))
    for ln in fh:
        f = ln.rstrip(NL).split(TAB)
        tier[(f[ix["chrom"]], int(f[ix["pos"]]), f[ix["ref"]], f[ix["alt"]])] = f[ix["lof_tier"]]   # keyed by allele (T6)
print("rare tiered LoF sites: %d" % len(tier))

# ---- trios with an affected offspring
trios = []
import glob
for ped in sorted(glob.glob(os.path.join(R, "run_*", "analysis", "*.ped"))):
    if os.environ.get("EXCLUDE_FAMILY", "__none__") in ped:
        continue
    fam = os.path.basename(ped)[:-4]
    rows = []
    for ln in open(ped):
        f = ln.split()
        if len(f) >= 6:
            rows.append(f)
    for f in rows:
        kid, fa, mo, aff = f[1], f[2], f[3], f[5]
        if fa == "0" or mo == "0" or aff != "2":
            continue
        trios.append((fam, kid, fa, mo))
print("complete trios with an affected offspring: %d" % len(trios))

samples_needed = set()
for _, k, fa, mo in trios:
    samples_needed.update([k, fa, mo])

# ---- genotypes at the tiered sites
tgt = os.path.join(T, "rare_targets.txt")
cmd = ("singularity exec -B /expanse:/expanse " + SIF +
       " bcftools query -T " + tgt +
       " -f '%CHROM" + TAB + "%POS" + TAB + "%REF" + TAB + "%ALT[" + TAB + "%SAMPLE|%GT|%GQ|%DP]" + NL + "' " + BCF)
p = subprocess.Popen(["bash", "-c", cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = p.communicate()
if p.returncode != 0:
    sys.exit("bcftools failed: " + err.decode()[:400])
lines = out.decode().splitlines()
print("genotype rows: %d" % len(lines))


def parse(cell):
    q = cell.split("|")
    if len(q) < 4:
        return None
    try:
        gq, dp = int(q[2]), int(q[3])
    except ValueError:
        return None
    g = q[1]
    if g in ("./.", ".|.", "."):
        return None
    if gq < GQ_MIN or dp < DP_MIN:
        return None
    a = g.replace("|", "/").split("/")
    return a


cnt = collections.defaultdict(lambda: [0, 0])   # tier -> [T, NT]
for ln in lines:
    f = ln.split(TAB)
    if len(f) < 3:
        continue
    key = (f[0], int(f[1]), f[2], f[3])
    t = tier.get(key)
    if t is None:
        continue
    on_x = f[0] in ("chrX", "chrY")
    gt = {}
    for cell in f[4:]:
        q = cell.split("|")
        if q:
            gt[q[0]] = parse(cell)
    for fam, kid, fa, mo in trios:
        gk = gt.get(kid)
        if gk is None:
            continue
        nk = gk.count("1")
        for par in (fa, mo):
            if on_x and par == fa:
                continue                  # a father is hemizygous on X: his single allele is not a heterozygous transmission test (T7)
            gp = gt.get(par)
            if gp is None:
                continue
            if gp.count("1") != 1:      # only heterozygous parents are informative
                continue
            # with one het parent, whether THIS parent transmitted is inferable only
            # when the other parent is hom-ref; otherwise the child's alt could come
            # from either. Restrict to that unambiguous case.
            other = mo if par == fa else fa
            go = gt.get(other)
            if go is None or go.count("1") != 0:
                continue
            if nk >= 1:
                cnt[t][0] += 1
            else:
                cnt[t][1] += 1

print("")
print("=" * 66)
print("  TRANSMISSION TEST -- inherited LoF, affected offspring")
print("=" * 66)
print("  %-9s %8s %8s %8s %9s %10s" % ("tier", "T", "NT", "total", "T_frac", "p(normal)"))   # normal approximation with continuity correction, not an exact binomial
tot_t = tot_nt = 0
for t in ("lof_t1", "lof_t2", "lof_t3"):
    a, b = cnt.get(t, [0, 0])
    n = a + b
    tot_t += a
    tot_nt += b
    if n == 0:
        print("  %-9s %8d %8d %8d %9s %10s" % (t, a, b, n, "-", "-"))
        continue
    frac = 1.0 * a / n
    # two-sided exact binomial via normal approx with continuity correction
    se = math.sqrt(0.25 / n)
    z = (abs(frac - 0.5) - 0.5 / n) / se if se > 0 else 0
    pv = min(1.0, math.erfc(max(0.0, z) / math.sqrt(2)))   # z < 0 when |frac-0.5| < the correction gave p > 1 (T21)
    print("  %-9s %8d %8d %8d %9.4f %10.4f" % (t, a, b, n, frac, pv))
n = tot_t + tot_nt
if n:
    frac = 1.0 * tot_t / n
    se = math.sqrt(0.25 / n)
    z = (abs(frac - 0.5) - 0.5 / n) / se
    print("  %-9s %8d %8d %8d %9.4f %10.4f" % ("ALL", tot_t, tot_nt, n, frac, min(1.0, math.erfc(max(0.0, z) / math.sqrt(2)))))
print("=" * 66)
print("  T_frac > 0.5 = over-transmission to affected offspring (burden enrichment).")
print("  Only unambiguous meioses counted: one het parent, other parent hom-ref.")
