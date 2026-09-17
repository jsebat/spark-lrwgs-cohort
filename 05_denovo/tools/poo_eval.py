"""Score the pipeline's parent-of-origin against the independent read-backed assignment, then apply JS's test:
re-split each child's de novo SNVs by the read-backed parent and regress paternal-origin ~ paternal age and
maternal-origin ~ maternal age. Also simulate PERFECT assignment for these trios' ages, so the size of difference
that is achievable at n = 33 is known before the pipeline is judged against it.

usage: poo_eval.py <poo_readcheck_dir> <milad_ages.tsv> <manifest.tsv>
"""
import csv
import glob
import math
import os
import random
import sys
from collections import Counter, defaultdict

TAB = chr(9)
D, AGES, MAN = sys.argv[1:4]

rows = []
for f in sorted(glob.glob(os.path.join(D, "*.poo_readcheck.tsv"))):
    rows.extend(csv.DictReader(open(f, newline=""), delimiter=TAB))
sex = []
for f in sorted(glob.glob(os.path.join(D, "*.sex.tsv"))):
    sex.extend(csv.DictReader(open(f, newline=""), delimiter=TAB))
print("children with output: %d   SNV rows: %d" % (len(set(r["child"] for r in rows)), len(rows)))

# ------------------------------------------------------------------ genetic sex vs manifest
print()
print("=== genetic sex (chrX 20-60 Mb het fraction): males should be ~0, females ~0.2-0.3")
bad = 0
for r in sex:
    hf = float(r["chrX_het_frac"]) if r["chrX_het_frac"] not in ("", "None") else None
    gsex = None if hf is None else ("F" if hf > 0.1 else "M")
    ok = gsex == r["manifest_sex"]
    bad += not ok
    if not ok or r["role"] != "child":
        pass
    if not ok:
        print("  *** %s %s role=%s manifest=%s chrX_het_frac=%s -> genetic %s" % (r["family"], r["sample"], r["role"], r["manifest_sex"], r["chrX_het_frac"], gsex))
print("  samples checked: %d   mismatches: %d" % (len(sex), bad))

# ------------------------------------------------------------------ per-call agreement
t1 = [r for r in rows if r["dnm_call"] == "YES"]
print()
print("=== tier-1 SNVs: %d" % len(t1))
rb_ok = [r for r in t1 if r["readback_poo"] in ("paternal", "maternal")]
print("read-backed assignment available (>=2 informative votes, >=80%% one parent): %d (%.1f%%)" % (len(rb_ok), 100.0 * len(rb_ok) / max(1, len(t1))))
print("read-backed verdicts:", Counter(r["readback_poo"] for r in t1))
ctrl = [r for r in rb_ok if r["ref_control_ok"] != ""]
print("REF-read control (other haplotype votes the OTHER parent): %d of %d pass (%.1f%%)" % (
    sum(int(r["ref_control_ok"]) for r in ctrl), len(ctrl), 100.0 * sum(int(r["ref_control_ok"]) for r in ctrl) / max(1, len(ctrl))))

both = [r for r in rb_ok if r["pipeline_poo"] in ("paternal", "maternal")]
cm = Counter((r["pipeline_poo"], r["readback_poo"]) for r in both)
print()
print("confusion, PIPELINE (rows) vs READ-BACKED (cols), n=%d:" % len(both))
print("%-20s %10s %10s" % ("", "rb:paternal", "rb:maternal"))
for p in ("paternal", "maternal"):
    print("%-20s %10d %10d" % ("pipeline:" + p, cm[(p, "paternal")], cm[(p, "maternal")]))
agree = cm[("paternal", "paternal")] + cm[("maternal", "maternal")]
print("agreement: %d / %d = %.1f%%" % (agree, len(both), 100.0 * agree / max(1, len(both))))
if cm[("paternal", "paternal")] + cm[("paternal", "maternal")]:
    print("of pipeline-PATERNAL calls, read-back says maternal: %.1f%%" % (100.0 * cm[("paternal", "maternal")] / (cm[("paternal", "paternal")] + cm[("paternal", "maternal")])))
if cm[("maternal", "paternal")] + cm[("maternal", "maternal")]:
    print("of pipeline-MATERNAL calls, read-back says paternal: %.1f%%" % (100.0 * cm[("maternal", "paternal")] / (cm[("maternal", "paternal")] + cm[("maternal", "maternal")])))

# where does the disagreement come from?  orientation table vs reads, and HP tags vs reads
ot = [r for r in rb_ok if r["agree_orient_table_readback"] != ""]
print()
print("orientation TABLE (HAP1_PAT/MAT) vs what the reads say hap1 is: agree %d / %d = %.1f%%" % (
    sum(int(r["agree_orient_table_readback"]) for r in ot), len(ot), 100.0 * sum(int(r["agree_orient_table_readback"]) for r in ot) / max(1, len(ot))))
ho = [r for r in rb_ok if r["agree_hporient_readback"] != ""]
print("HP-tag + table recomputation vs read-back: agree %d / %d = %.1f%%" % (
    sum(int(r["agree_hporient_readback"]) for r in ho), len(ho), 100.0 * sum(int(r["agree_hporient_readback"]) for r in ho) / max(1, len(ho))))
# undetermined-by-pipeline but resolvable by reads
und = [r for r in t1 if r["pipeline_poo"] == "undetermined" and r["readback_poo"] in ("paternal", "maternal")]
print("pipeline undetermined but read-back resolves: %d  (%s)" % (len(und), Counter(r["poo_reason"].split(":")[0] for r in und).most_common(5)))
# disagreements by flag / reason
dis = [r for r in both if r["pipeline_poo"] != r["readback_poo"]]
print("disagreements by pipeline poo_reason:", Counter(r["poo_reason"] for r in dis).most_common(6))
print("disagreements: alt reads HP split (hp1,hp2,untagged) examples:", [(r["alt_hp1"], r["alt_hp2"], r["alt_hp0"]) for r in dis[:12]])
print("disagreements: per child:", Counter(r["child"] for r in dis).most_common(8))

# ------------------------------------------------------------------ JS's test: age regressions after re-splitting
ages = {}
for r in csv.DictReader(open(AGES, newline=""), delimiter=TAB):
    try:
        ages[r["child_sp_id"]] = (float(r["father_age_at_child_birth"]), float(r["mother_age_at_child_birth"]))
    except ValueError:
        pass


def pear(x, y):
    n = len(x); mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x); syy = sum((b - my) ** 2 for b in y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    r = sxy / math.sqrt(sxx * syy) if sxx and syy else 0.0
    return r, sxy / sxx if sxx else 0.0


per = defaultdict(Counter)
for r in t1:
    per[r["child"]]["total"] += 1
    per[r["child"]]["pipe_" + r["pipeline_poo"]] += 1
    per[r["child"]]["rb_" + r["readback_poo"]] += 1
kids = [k for k in per if k in ages]
pa = [ages[k][0] for k in kids]; ma = [ages[k][1] for k in kids]
print()
print("=== JS's test, n=%d trios with ages" % len(kids))
print("%-42s %7s %9s" % ("", "r", "slope/yr"))
for lab, key, x in (("total SNV ~ paternal age", "total", pa),
                    ("PIPELINE paternal-origin ~ paternal age", "pipe_paternal", pa),
                    ("PIPELINE maternal-origin ~ maternal age", "pipe_maternal", ma),
                    ("READ-BACKED paternal-origin ~ paternal age", "rb_paternal", pa),
                    ("READ-BACKED maternal-origin ~ maternal age", "rb_maternal", ma),
                    ("READ-BACKED maternal-origin ~ PATERNAL age", "rb_maternal", pa)):
    r, s = pear(x, [per[k][key] for k in kids])
    print("%-42s %7.3f %9.3f" % (lab, r, s))
print("read-backed fraction maternal among assigned: %.3f" % (
    sum(per[k]["rb_maternal"] for k in kids) / max(1, sum(per[k]["rb_maternal"] + per[k]["rb_paternal"] for k in kids))))

# ------------------------------------------------------------------ what PERFECT assignment would look like at this n
# model: paternal ~ Poisson(a_p + 1.5*pat_age adjusted so the mean matches the cohort), maternal ~ Poisson(a_m + 0.37*mat_age),
# assigned fraction as observed; 2000 replicates; report the distribution of (r, slope) for both regressions.
random.seed(1)
tot_mean = sum(per[k]["total"] for k in kids) / len(kids)
mp, mm = sum(pa) / len(pa), sum(ma) / len(ma)
frac_assigned = sum(per[k]["rb_maternal"] + per[k]["rb_paternal"] for k in kids) / max(1, sum(per[k]["total"] for k in kids))
b_p, b_m = 1.5, 0.37
# intercepts so that E[pat + mat] at the mean ages equals the observed total, split ~80/20
a_p = 0.80 * tot_mean - b_p * mp
a_m = 0.20 * tot_mean - b_m * mm
rp, sp, rm, sm = [], [], [], []
for _ in range(2000):
    P = [max(0, int(random.gauss(a_p + b_p * x, math.sqrt(max(1, a_p + b_p * x)))) ) for x in pa]
    M = [max(0, int(random.gauss(a_m + b_m * x, math.sqrt(max(1, a_m + b_m * x)))) ) for x in ma]
    # thin to the assigned fraction
    P = [sum(1 for _ in range(v) if random.random() < frac_assigned) for v in P]
    M = [sum(1 for _ in range(v) if random.random() < frac_assigned) for v in M]
    r1, s1 = pear(pa, P); r2, s2 = pear(ma, M)
    rp.append(r1); sp.append(s1); rm.append(r2); sm.append(s2)


def q(v):
    v = sorted(v); n = len(v)
    return v[n // 2], v[int(0.05 * n)], v[int(0.95 * n)]


print()
print("=== simulated PERFECT assignment for these %d trios (paternal 1.5/yr, maternal 0.37/yr, assigned frac %.2f), 2000 reps" % (len(kids), frac_assigned))
print("paternal-origin ~ paternal age:  r median %.2f (90%% %.2f-%.2f)   slope median %.2f (%.2f-%.2f)" % (q(rp) + q(sp)))
print("maternal-origin ~ maternal age:  r median %.2f (90%% %.2f-%.2f)   slope median %.2f (%.2f-%.2f)" % (q(rm) + q(sm)))
print("-> this is the size of difference a perfect assignment can show at n=%d; judge the observed rows above against it." % len(kids))
