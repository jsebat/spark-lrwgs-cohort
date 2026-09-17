# -*- coding: utf-8 -*-
"""SPARK phenotype: severity / cognition predictors for the lrWGS cohort.

Q1  which predictors are available for most lrWGS probands
Q2  are lrWGS probands more severe than SPARK ASD probands overall
Q3  are the 4 de novo SV carriers more severe than the other lrWGS probands

Uses core_descriptive_variables (consolidated, one row per person) plus the
approximated cognitive impairment table. Aggregate statistics only, except the
four SV carriers, which JS asked to see individually. Standard library only.
"""
import csv
import sys
import os
import math
import statistics as st
from collections import Counter, defaultdict

P = os.environ.get("PHENO_DIR") or sys.exit("set PHENO_DIR to the SPARK phenotype release directory")
L = os.path.dirname(os.path.abspath(__file__))
CORE = os.path.join(P, "core_descriptive_variables-2026-06-25.csv")
ACI = os.path.join(P, "approximated_cognitive_impairment-2026-06-25.csv")

ours = set(x.strip() for x in open(os.environ.get("PROBANDS", "probands.txt")) if x.strip() and not x.startswith("#"))
sv = dict(l.split()[:2] for l in open(os.environ.get("CARRIERS", "carriers.tsv")) if l.strip() and not l.startswith("#"))

NUMERIC = ["scq_total_final_score", "rbsr_total_final_score", "fsiq", "viq", "nviq",
           "vineland_abc_ss_latest", "reported_cog_test_score_latest",
           "used_words_age_mos", "walked_age_mos", "age_onset_mos", "diagnosis_age",
           "age_at_registration_years"]
CATEG = ["cognitive_impairment_latest", "language_level_latest", "id_confirmed",
         "regress_lang_y_n", "regress_other_y_n", "dcdq_dcd", "neuro_sz", "sex"]


def num(x):
    try:
        v = float(x)
        return v if not math.isnan(v) else None
    except (TypeError, ValueError):
        return None


core = {}
with open(CORE, newline="", encoding="utf-8", errors="replace") as fh:
    for r in csv.DictReader(fh):
        core[r["subject_sp_id"]] = r
print("core rows: %d" % len(core))
aci = {}
with open(ACI, newline="", encoding="utf-8", errors="replace") as fh:
    for r in csv.DictReader(fh):
        aci[r["subject_sp_id"]] = r

found = ours & set(core)
print("lrWGS probands in core table: %d of %d" % (len(found), len(ours)))
missing = ours - set(core)
if missing:
    print("  NOT in phenotype release: %s" % ", ".join(sorted(missing)))

# SPARK comparison pool: all ASD-positive individuals in core (children + adults)
asd_vals = Counter(r.get("asd", "").strip() for r in core.values())
print("asd column values in core: %s" % dict(asd_vals.most_common(5)))
spark_asd = [s for s, r in core.items() if r.get("asd", "").strip().lower() in ("true", "1", "yes")]
if not spark_asd:
    spark_asd = list(core)
    print("  asd filter matched nothing -> using ALL core rows; check whether core is ASD-only")
print("SPARK ASD individuals (comparison pool): %d" % len(spark_asd))

# ---------------------------------------------------------------- Q1 coverage
print("")
print("=" * 78)
print("Q1  COVERAGE in lrWGS probands (n=%d)" % len(found))
print("=" * 78)
cov = []
for v in NUMERIC:
    n = sum(1 for s in found if num(core[s].get(v)) is not None)
    cov.append((v, n))
for v in CATEG:
    n = sum(1 for s in found if core[s].get(v, "").strip() not in ("", "NA", "NaN"))
    cov.append((v, n))
n_aci = sum(1 for s in found if s in aci and aci[s].get("derived_cog_impair", "").strip() not in ("", "NA"))
cov.append(("derived_cog_impair (approx table)", n_aci))
n_ml = sum(1 for s in found if s in aci and aci[s].get("ml_predicted_cog_impair", "").strip() not in ("", "NA"))
cov.append(("ml_predicted_cog_impair (approx table)", n_ml))
for v, n in sorted(cov, key=lambda x: -x[1]):
    print("  %-40s %3d / %d  (%3.0f%%)" % (v, n, len(found), 100.0 * n / len(found)))

# ---------------------------------------------------------------- Q2 ours vs SPARK
print("")
print("=" * 78)
print("Q2  lrWGS probands vs ALL SPARK ASD  (numeric: median [IQR]; higher SCQ/RBSR = more severe; lower IQ/Vineland = more impaired)")
print("=" * 78)


def summ(vals):
    vals = sorted(v for v in vals if v is not None)
    if not vals:
        return None
    q = lambda p: vals[min(len(vals) - 1, int(p * len(vals)))]
    return (len(vals), st.median(vals), q(0.25), q(0.75), st.mean(vals))


def mw_u_p(a, b):
    """Mann-Whitney U, normal approximation, two-sided."""
    a = [v for v in a if v is not None]
    b = [v for v in b if v is not None]
    if len(a) < 3 or len(b) < 3:
        return None
    allv = sorted([(v, 0) for v in a] + [(v, 1) for v in b])
    ranks = {}
    i = 0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1][0] == allv[i][0]:
            j += 1
        r = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            ranks[k] = r
        i = j + 1
    ra = sum(ranks[k] for k, (v, g) in enumerate(allv) if g == 0)
    na, nb = len(a), len(b)
    u = ra - na * (na + 1) / 2.0
    mu = na * nb / 2.0
    # tie-corrected variance: on a 4-level ordinal nearly every observation is tied and the uncorrected sd is too large,
    # so p was systematically too big (X15)
    ntot = na + nb
    tie_groups = {}
    for v, _ in allv:
        tie_groups[v] = tie_groups.get(v, 0) + 1
    tie_term = sum(t ** 3 - t for t in tie_groups.values())
    sd = math.sqrt(na * nb / 12.0 * ((ntot + 1) - tie_term / (ntot * (ntot - 1)))) if ntot > 1 else 0.0
    if sd == 0:
        return None
    z = (u - mu) / sd
    return 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))


for v in NUMERIC:
    a = [num(core[s].get(v)) for s in found]
    b = [num(core[s].get(v)) for s in spark_asd]
    sa, sb = summ(a), summ(b)
    if not sa or not sb:
        continue
    p = mw_u_p(a, b)
    print("  %-32s lrWGS n=%2d med=%6.1f [%5.1f-%5.1f]   SPARK n=%6d med=%6.1f [%5.1f-%5.1f]   MW p=%s"
          % (v, sa[0], sa[1], sa[2], sa[3], sb[0], sb[1], sb[2], sb[3],
             ("%.3g" % p) if p is not None else "-"))

print("")
print("  categorical (fraction with the impairment / flag):")


def frac(vals, pos):
    vals = [x.strip() for x in vals if x.strip() not in ("", "NA", "na_survey_logic")]   # not asked is unknown, not "no" (X9); matches valid() below
    if not vals:
        return (0, 0.0)
    k = sum(1 for x in vals if x in pos)
    return (len(vals), 100.0 * k / len(vals))


CAT_POS = {"cognitive_impairment_latest": ("TRUE", "True", "true", "1", "Yes", "yes"),
           "id_confirmed": ("TRUE", "True", "true", "1", "Yes", "yes"),
           "regress_lang_y_n": ("TRUE", "True", "true", "1", "Yes", "yes"),
           "regress_other_y_n": ("TRUE", "True", "true", "1", "Yes", "yes"),
           "neuro_sz": ("TRUE", "True", "true", "1", "Yes", "yes"),
           "dcdq_dcd": ("TRUE", "True", "true", "1", "Yes", "yes")}
for v, pos in CAT_POS.items():
    na, fa = frac([core[s].get(v, "") for s in found], pos)
    nb, fb = frac([core[s].get(v, "") for s in spark_asd], pos)
    print("  %-32s lrWGS %5.1f%% (n=%2d)    SPARK %5.1f%% (n=%6d)" % (v, fa, na, fb, nb))
# approx cognitive impairment
for v in ("derived_cog_impair", "ml_predicted_cog_impair"):
    a = [aci[s].get(v, "") for s in found if s in aci]
    b = [aci[s].get(v, "") for s in spark_asd if s in aci]
    na, fa = frac(a, ("TRUE", "True", "true", "1", "Yes", "yes"))
    nb, fb = frac(b, ("TRUE", "True", "true", "1", "Yes", "yes"))
    print("  %-32s lrWGS %5.1f%% (n=%2d)    SPARK %5.1f%% (n=%6d)" % (v, fa, na, fb, nb))
print("")
print("  raw value distributions (to check coding):")
for v in ("cognitive_impairment_latest", "language_level_latest", "id_confirmed"):
    c = Counter(core[s].get(v, "").strip() for s in found)
    print("    %-30s lrWGS: %s" % (v, dict(c)))
    c2 = Counter(core[s].get(v, "").strip() for s in spark_asd)
    print("    %-30s SPARK: %s" % ("", dict(c2.most_common(6))))

# ---------------------------------------------------------------- Q3 SV carriers
print("")
print("=" * 78)
print("Q3  the 4 de novo SV carriers vs the other lrWGS probands")
print("=" * 78)
others = [s for s in found if s not in sv]
keyv = ["scq_total_final_score", "rbsr_total_final_score", "fsiq", "vineland_abc_ss_latest",
        "used_words_age_mos", "walked_age_mos", "age_onset_mos"]
print("  %-10s %-8s %-4s %-6s %-6s %-6s %-6s %-6s %-6s %-6s %-9s %-12s %s"
      % ("proband", "gene", "sex", "SCQ", "RBSR", "FSIQ", "VABS", "words", "walk", "onset",
         "cogimp", "lang_level", "id_conf"))
for s, g in sv.items():
    r = core.get(s)
    if not r:
        print("  %-10s %-8s  NOT IN PHENOTYPE RELEASE" % (s, g))
        continue
    f = lambda k: (r.get(k, "") or "").strip() or "."
    print("  %-10s %-8s %-4s %-6s %-6s %-6s %-6s %-6s %-6s %-6s %-9s %-12s %s"
          % (s, g, f("sex"), f("scq_total_final_score"), f("rbsr_total_final_score"), f("fsiq"),
             f("vineland_abc_ss_latest"), f("used_words_age_mos"), f("walked_age_mos"),
             f("age_onset_mos"), f("cognitive_impairment_latest")[:9], f("language_level_latest")[:12],
             f("id_confirmed")))
print("")
print("  other lrWGS probands (n=%d), median [IQR]:" % len(others))
for v in keyv:
    so = summ([num(core[s].get(v)) for s in others])
    sc = [num(core[s].get(v)) for s in sv if s in core]
    sc = [x for x in sc if x is not None]
    if so:
        print("  %-28s others n=%2d med=%6.1f [%5.1f-%5.1f]   SV carriers: %s"
              % (v, so[0], so[1], so[2], so[3], ", ".join("%.0f" % x for x in sc) or "-"))
for v, pos in (("cognitive_impairment_latest", CAT_POS["cognitive_impairment_latest"]),):
    no, fo = frac([core[s].get(v, "") for s in others], pos)
    nc, fc = frac([core[s].get(v, "") for s in sv if s in core], pos)
    print("  %-28s others %5.1f%% (n=%d)    SV carriers %5.1f%% (n=%d)" % (v, fo, no, fc, nc))


# ---------------------------------------------------------------- Fisher / ordinal
def fisher_2x2(a, b, c, d):
    """two-sided Fisher exact p for [[a,b],[c,d]] via hypergeometric enumeration."""
    from math import lgamma, exp
    def lchoose(n, k):
        return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
    n = a + b + c + d
    r1, c1 = a + b, a + c
    def pmf(x):
        return exp(lchoose(r1, x) + lchoose(n - r1, c1 - x) - lchoose(n, c1))
    p_obs = pmf(a)
    lo, hi = max(0, c1 - (n - r1)), min(r1, c1)
    return sum(pmf(x) for x in range(lo, hi + 1) if pmf(x) <= p_obs + 1e-12)

def yes(x):
    return x.strip().lower() in ("true", "1", "yes")
def valid(x):
    return x.strip() not in ("", "NA", "na_survey_logic")

print("")
print("=" * 78)
print("FISHER: cognitive_impairment_latest (parent-reported; 100% coverage)")
print("=" * 78)
def ci_counts(ids):
    v = [core[s].get("cognitive_impairment_latest", "") for s in ids if s in core]
    v = [x for x in v if valid(x)]
    return sum(1 for x in v if yes(x)), sum(1 for x in v if not yes(x))
a, b = ci_counts(found); c, d = ci_counts(spark_asd)
print("  lrWGS  impaired %d / %d  (%.1f%%)" % (a, a + b, 100.0 * a / (a + b)))
print("  SPARK  impaired %d / %d  (%.1f%%)" % (c, c + d, 100.0 * c / (c + d)))
print("  Fisher two-sided p = %.3g" % fisher_2x2(a, b, c, d))
a2, b2 = ci_counts(list(sv)); c2, d2 = ci_counts(others)
print("  SV carriers impaired %d / %d ;  other lrWGS %d / %d  (%.1f%%)   Fisher p = %.3g"
      % (a2, a2 + b2, c2, c2 + d2, 100.0 * c2 / max(1, c2 + d2), fisher_2x2(a2, b2, c2, d2)))

print("")
print("=" * 78)
print("ORDINAL: language_level_latest (0_no_words worst ... 3_long_sentences best)")
print("=" * 78)
LANG = {"0_no_words": 0, "1_single_words": 1, "2_short_sentences": 2, "3_long_sentences": 3}
def lang(ids):
    return [LANG[core[s]["language_level_latest"].strip()] for s in ids
            if s in core and core[s].get("language_level_latest", "").strip() in LANG]
la, lb, lo = lang(found), lang(spark_asd), lang(others)
def dist(v):
    c = Counter(v); n = len(v)
    return "  ".join("%s:%5.1f%%" % (k, 100.0 * c[k] / n) for k in (0, 1, 2, 3)) if n else "-"
print("  lrWGS  n=%3d  %s   mean=%.2f" % (len(la), dist(la), st.mean(la)))
print("  SPARK  n=%3d  %s   mean=%.2f" % (len(lb), dist(lb), st.mean(lb)))
print("  MW p (lrWGS vs SPARK) = %.3g" % mw_u_p(la, lb))
print("  minimally verbal (0-1) lrWGS %.1f%%  vs SPARK %.1f%%"
      % (100.0 * sum(1 for x in la if x <= 1) / len(la), 100.0 * sum(1 for x in lb if x <= 1) / len(lb)))
lc = lang(list(sv))
print("  SV carriers: %s   others n=%d mean=%.2f   MW p = %s"
      % (lc, len(lo), st.mean(lo), ("%.3g" % mw_u_p(lc, lo)) if len(lc) >= 3 else "n<3"))

print("")
print("=" * 78)
print("derived_cog_impair (approximated table, 88% coverage)")
print("=" * 78)
def dci(ids):
    v = [aci[s].get("derived_cog_impair", "") for s in ids if s in aci]
    v = [x for x in v if valid(x)]
    return sum(1 for x in v if yes(x)), sum(1 for x in v if not yes(x))
a, b = dci(found); c, d = dci(spark_asd)
print("  lrWGS %d/%d (%.1f%%)   SPARK %d/%d (%.1f%%)   Fisher p = %.3g"
      % (a, a + b, 100.0 * a / (a + b), c, c + d, 100.0 * c / (c + d), fisher_2x2(a, b, c, d)))
a2, b2 = dci(list(sv)); c2, d2 = dci(others)
print("  SV carriers %d/%d   others %d/%d   Fisher p = %.3g" % (a2, a2 + b2, c2, c2 + d2, fisher_2x2(a2, b2, c2, d2)))
