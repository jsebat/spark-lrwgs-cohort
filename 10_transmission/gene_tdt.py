#!/usr/bin/env python3
"""Stage 3 (D75): parent-of-origin TDT of SVOPL loss-of-function alleles in SPARK WES, SPARK WGS and SSC.
Source: James Guevara's rare-variant pipeline annotated family genotype tables (output/indexed/chr7.merged.tsv.gz).
Strata: mother->proband, father->proband, mother->unaffected sibling, father->unaffected sibling.
Informative meiosis: exactly one parent heterozygous (other homozygous reference), child called; all three pass
GQ >= 20, DP >= 10; heterozygous calls need allele balance 0.25-0.75. WES/WGS duplicate families resolved to WGS.
Variant sets: (A) recurrent stop-gain chr7:138656474 G>A; (B) all other LOFTEE-HC LoF with gnomAD/cohort AF < 0.001;
(C) all HC LoF. Per-cohort counts, then pooled counts with a two-sided exact binomial P (H0: T = NT).
"""
import os, subprocess, csv, collections, math, sys
D = os.environ["SHORTREAD_ROOT"]; OUT = os.environ["DATA_ROOT"] + "/transmission"; os.makedirs(OUT, exist_ok=True)
SIF = os.environ["SIF"]; sx = ["singularity", "exec", "-B", "/expanse:/expanse", SIF]
GENE = os.environ.get("GENE", "SVOPL"); REGION = os.environ.get("REGION", "chr7:138594284-138701362"); CHR = REGION.split(":")[0]
REC = os.environ.get("RECURRENT", "chr7:138656474:G:A")
COHORTS = [("SPARK_WGS", D + "/nf_rare_spark_wgs/output/indexed/%s.merged.tsv.gz" % CHR, os.environ["PED_SPARK_WGS"]),
           ("SPARK_WES", D + "/nf_rare_spark_wes/output/indexed/%s.merged.tsv.gz" % CHR, os.environ["PED_SPARK_WES"]),
           ("SSC", D + "/nf_rare_ssc/output/indexed/%s.merged.tsv.gz" % CHR, os.environ["PED_SSC"])]
TAB = "\t"
def binom_two_sided(k, n):
    if n == 0: return float("nan")
    lp = [math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) - n * math.log(2.0) for i in range(n + 1)]
    pk = lp[k]
    return min(1.0, sum(math.exp(v) for v in lp if v <= pk + 1e-12))
def read_ped(p):
    ped = {}
    for ln in open(p):
        f = ln.split()
        if len(f) >= 6: ped[f[1]] = dict(fam=f[0], pat=f[2], mat=f[3], sex=f[4], aff=f[5])
    return ped
def gt_alt_count(gt, allele_num):
    g = gt.replace("|", "/")
    if "." in g: return None
    return sum(1 for a in g.split("/") if a == str(allele_num))
results = collections.defaultdict(lambda: collections.Counter())   # (set, cohort, stratum) -> Counter(T/NT)
variant_info = {}; carriers_by_cohort = collections.defaultdict(set); fam_seen = collections.defaultdict(set); events = []
wgs_fams = set()
for name, table, pedf in COHORTS:
    ped = read_ped(pedf)
    hdr = subprocess.run(["bash", "-c", "zcat %s | head -1" % table], stdout=subprocess.PIPE, universal_newlines=True).stdout.rstrip("\n").split(TAB)
    ci = {c: i for i, c in enumerate(hdr)}
    txt = subprocess.run(sx + ["tabix", table, REGION], stdout=subprocess.PIPE, universal_newlines=True).stdout
    # collect per (variant, allele) -> sample -> (GT alt count, GQ, DP, AB)
    per_var = collections.defaultdict(dict); info = {}
    for ln in txt.splitlines():
        f = ln.split(TAB)
        if len(f) < len(hdr): continue
        if f[ci["SYMBOL"]] != GENE or f[ci["LoF"]] != "HC": continue
        cons = f[ci["Consequence"]]
        if not any(k in cons for k in ("stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant")): continue
        try: an = int(float(f[ci["allele_num"]]))
        except ValueError: an = 1
        alt = f[ci["ALT_var"]] if f[ci["ALT_var"]] not in ("", ".") else f[ci["ALT"]].split(",")[an - 1] if an - 1 < len(f[ci["ALT"]].split(",")) else f[ci["ALT"]]
        key = "%s:%s:%s:%s" % (f[ci["#CHROM"]], f[ci["POS"]], f[ci["REF"]], alt)
        s = f[ci["SAMPLE"]]; gt = f[ci["GT"]]
        try: gq = float(f[ci["GQ"]]); dp = float(f[ci["DP"]])
        except ValueError: gq = dp = 0.0
        try: ar, aa = float(f[ci["AD_ref"]]), float(f[ci["AD_alt"]]); ab = aa / (ar + aa) if ar + aa > 0 else None
        except ValueError: ab = None
        n_alt = gt_alt_count(gt, an)
        per_var[key][s] = (n_alt, gq, dp, ab)
        gaf = f[ci["gnomAD4.1_joint_AF"]]; caf = f[ci["AF"]]
        info[key] = (cons, gaf, caf, f[ci["CANONICAL"]], f[ci["HGVSp"]])
    print("%s: %s HC LoF variants %d; sample rows %d" % (name, GENE, len(per_var), sum(len(v) for v in per_var.values())), flush=True)
    for key, samp in per_var.items():
        cons, gaf, caf, canon, hgvsp = info[key]; variant_info[key] = info[key]
        try: g = float(gaf) if gaf not in ("", ".") else 0.0
        except ValueError: g = 0.0
        try: c = float(caf) if caf not in ("", ".") else 0.0
        except ValueError: c = 0.0
        sets = ["C_all_HC_LoF"]
        if key == REC: sets.append("A_recurrent_stopgain")
        elif g < 0.001 and c < 0.001: sets.append("B_other_rare_HC_LoF")
        # families with a carrier parent
        fams = collections.defaultdict(list)
        for s, v in samp.items():
            if s in ped: fams[ped[s]["fam"]].append(s)
        for fam, members in fams.items():
            if name == "SPARK_WES" and fam in wgs_fams: continue           # WGS takes precedence for duplicated families
            children = [s for s in ped if ped[s]["fam"] == fam and ped[s]["pat"] != "0" and ped[s]["mat"] != "0"]
            for ch in children:
                pa, ma = ped[ch]["pat"], ped[ch]["mat"]
                if pa not in samp or ma not in samp or ch not in samp: continue     # need all three called in this variant's family rows
                gp, gm, gc = samp[pa], samp[ma], samp[ch]
                def ok(v, need_het=False):
                    n, gq, dp, ab = v
                    if n is None or gq < 20 or dp < 10: return False
                    if n == 1 and (ab is None or ab < 0.25 or ab > 0.75): return False
                    return True
                if not (ok(gp) and ok(gm) and ok(gc)): continue
                if gp[0] == 1 and gm[0] == 0: parent = "father"
                elif gm[0] == 1 and gp[0] == 0: parent = "mother"
                else: continue
                role = "proband" if ped[ch]["aff"] == "2" else ("sibling" if ped[ch]["aff"] == "1" else None)
                if role is None or gc[0] is None or gc[0] > 1: continue
                tr = "T" if gc[0] == 1 else "NT"
                for st in sets:
                    results[(st, name, parent + "->" + role)][tr] += 1
                events.append((name, key, fam, ch, role, parent, tr))
                if name == "SPARK_WGS": wgs_fams.add(fam)
                carriers_by_cohort[name].add(pa if parent == "father" else ma)
# ---------------- output
strata = ["mother->proband", "father->proband", "mother->sibling", "father->sibling"]
sets = ["A_recurrent_stopgain", "B_other_rare_HC_LoF", "C_all_HC_LoF"]
with open(OUT + "/%s_tdt.tsv" % GENE, "w") as fh:
    fh.write(TAB.join(["variant_set", "cohort", "stratum", "T", "NT", "T_frac", "p_binomial"]) + "\n")
    for st in sets:
        for name, _, _ in COHORTS + [("POOLED", None, None)]:
            for s in strata:
                if name == "POOLED": T = sum(results[(st, c, s)]["T"] for c, _, _ in COHORTS); NT = sum(results[(st, c, s)]["NT"] for c, _, _ in COHORTS)
                else: T, NT = results[(st, name, s)]["T"], results[(st, name, s)]["NT"]
                n = T + NT; fh.write(TAB.join(map(str, [st, name, s, T, NT, "%.3f" % (T / n) if n else "", "%.3g" % binom_two_sided(T, n) if n else ""])) + "\n")
with open(OUT + "/%s_lof_variants.tsv" % GENE, "w") as fh:
    fh.write(TAB.join(["variant", "consequence", "gnomAD_AF", "cohort_AF_first_seen", "canonical", "HGVSp"]) + "\n")
    for k, v in sorted(variant_info.items()): fh.write(TAB.join([k] + list(v)) + "\n")
with open(OUT + "/%s_transmission_events.tsv" % GENE, "w") as fh:
    fh.write(TAB.join(["cohort", "variant", "family", "child", "child_role", "transmitting_parent", "transmitted"]) + "\n")
    for e in events: fh.write(TAB.join(e) + "\n")
print(open(OUT + "/%s_tdt.tsv" % GENE).read())
# maternal vs paternal transmission contrast in probands (pooled), set A and C: 2x2 Fisher-like via chi-square approx
for st in sets:
    Tm = sum(results[(st, c, "mother->proband")]["T"] for c, _, _ in COHORTS); Nm = sum(results[(st, c, "mother->proband")]["NT"] for c, _, _ in COHORTS)
    Tp = sum(results[(st, c, "father->proband")]["T"] for c, _, _ in COHORTS); Np = sum(results[(st, c, "father->proband")]["NT"] for c, _, _ in COHORTS)
    print("%s pooled probands: maternal T/NT %d/%d (%.3f) vs paternal %d/%d (%.3f)" % (st, Tm, Nm, Tm / max(1, Tm + Nm), Tp, Np, Tp / max(1, Tp + Np)))
