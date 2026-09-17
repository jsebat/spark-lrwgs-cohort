#!/usr/bin/env python3
"""Stage 3 extension (D76): expressed-allele TDT of LoF alleles in ALL catalogued imprinted genes, SPARK WES + WGS + SSC,
plus collapsed burden by s_het tier with clonal-hematopoiesis (CHIP) genes flagged and excluded from the tier tests.
Per gene and cohort: transmissions of LOFTEE-HC LoF alleles from mother / father to affected probands / unaffected
siblings; test allele = expressed-allele parent (maternal for maternally expressed genes, paternal for paternally
expressed); isoform-dependent / unknown / random genes tested both ways and flagged. Variant sets: rare (gnomAD and
cohort AF < 0.001), common (>= 0.001), all HC LoF. Filters as in svopl_tdt.py (D75)."""
import os, subprocess, collections, math, csv
D = os.environ["SHORTREAD_ROOT"]; ME = os.environ["DATA_ROOT"] + "/meth/"; OUT = os.environ["DATA_ROOT"] + "/transmission/"
SIF = os.environ["SIF"]; sx = ["singularity", "exec", "-B", "/expanse:/expanse", SIF]; TAB = "\t"
COHORTS = [("SPARK_WGS", D + "/nf_rare_spark_wgs/output/indexed/%s.merged.tsv.gz", os.environ["PED_SPARK_WGS"]),
           ("SPARK_WES", D + "/nf_rare_spark_wes/output/indexed/%s.merged.tsv.gz", os.environ["PED_SPARK_WES"]),
           ("SSC", D + "/nf_rare_ssc/output/indexed/%s.merged.tsv.gz", os.environ["PED_SSC"])]
# CHIP driver genes (Jaiswal 2014, Genovese 2014, Bick 2020, Kessler 2022, Beauchamp 2021 consensus set)
CHIP = set("DNMT3A TET2 ASXL1 JAK2 TP53 PPM1D SF3B1 SRSF2 U2AF1 CBL GNB1 GNAS IDH1 IDH2 KRAS NRAS CALR MPL BCOR BCORL1 STAG2 ZRSR2 CUX1 RUNX1 ETV6 PHF6 EZH2 KMT2D CREBBP NOTCH1 MYD88 CHEK2 ATM ZNF318 YLPM1 SRCAP ZBTB33 MTA2 SPRED2 BRCC3 KDM6A PRPF8 SETDB1 SETD2 SUZ12 CTCF WT1 NF1 RAD21 SMC1A SMC3 PTPN11 FLT3 NPM1 CEBPA GATA2 LUC7L2 CSF3R".split())
def binom_two_sided(k, n):
    if n == 0: return float("nan")
    lp = [math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) - n * math.log(2.0) for i in range(n + 1)]
    return min(1.0, sum(math.exp(v) for v in lp if v <= lp[k] + 1e-12))
def read_ped(p):
    ped = {}
    for ln in open(p):
        f = ln.split()
        if len(f) >= 6: ped[f[1]] = dict(fam=f[0], pat=f[2], mat=f[3], aff=f[5])
    return ped
peds = {name: read_ped(pf) for name, _, pf in COHORTS}
children_by_fam = {}
for name, ped in peds.items():
    d = collections.defaultdict(list)
    for s, r in ped.items():
        if r["pat"] != "0" and r["mat"] != "0": d[r["fam"]].append(s)
    children_by_fam[name] = d
headers = {}
def header(table):
    if table not in headers:
        headers[table] = subprocess.run(["bash", "-c", "zcat %s | head -1" % table], stdout=subprocess.PIPE, universal_newlines=True).stdout.rstrip("\n").split(TAB)
    return headers[table]
# genes
genes = []
for ln in open(ME + "imprinted_genes.bed"):
    c, s, e, g, ea, strand = ln.rstrip("\n").split(TAB); genes.append((g, c, int(s) + 1, int(e), ea))
shet = {}
for ln in open(os.environ["SHET_TABLE"]):
    f = ln.rstrip("\n").split(TAB)
    if len(f) >= 2:
        try: shet[f[0]] = float(f[1])
        except ValueError: pass
gb = {}
for r in csv.DictReader(open(os.environ["GENEBAYES_TABLE"]), delimiter=TAB):
    gb[r["ensg"].split(".")[0]] = r
ensg_of = {}
import gzip, re
with gzip.open(os.environ["GENCODE_GTF"], "rt") as fh:
    want = set(g for g, *_ in genes)
    for ln in fh:
        if ln.startswith("#"): continue
        f = ln.split(TAB)
        if f[2] != "gene": continue
        m = re.search(r'gene_name "([^"]+)"', f[8]); gi = re.search(r'gene_id "([^".]+)', f[8])
        if m and gi and m.group(1) in want: ensg_of[m.group(1)] = gi.group(1)
def tier(g):
    s = shet.get(g)
    if s is None: return "lof_t3"
    return "lof_t1" if s >= 0.18 else ("lof_t2" if s >= 0.03 else "lof_t3")
# gnomAD v4.1 constraint flags: outlier_lof = observed LoF exceeds expectation (the browser's clonal-hematopoiesis warning); JS: use gnomAD AND the curated driver list
gflags = collections.defaultdict(set)
with open(os.environ["GNOMAD_CONSTRAINT"]) as fh:
    rd = csv.DictReader(fh, delimiter=TAB)
    for r in rd:
        if "outlier_lof" in (r.get("constraint_flags") or ""): gflags[r["gene"]].add("gnomAD_outlier_lof")
def chip_flag(g):
    fl = set(gflags.get(g, []))
    if g in CHIP: fl.add("CHIP_driver_curated")
    return ";".join(sorted(fl))
strata = ["mother->proband", "father->proband", "mother->sibling", "father->sibling"]
sets = ["rare", "common", "all"]
res = collections.defaultdict(collections.Counter)   # (gene, set, cohort, stratum) -> T/NT
nvar = collections.defaultdict(set); events = []
wgs_fams = set()
for g, c, s0, e0, ea in genes:
    region = "%s:%d-%d" % (c, s0, e0)
    for name, tpl, _ in COHORTS:
        table = tpl % c
        if not os.path.exists(table): continue
        hdr = header(table); ci = {k: i for i, k in enumerate(hdr)}
        ped = peds[name]
        txt = subprocess.run(sx + ["tabix", table, region], stdout=subprocess.PIPE, universal_newlines=True).stdout
        per_var = collections.defaultdict(dict); info = {}
        for ln in txt.splitlines():
            f = ln.split(TAB)
            if len(f) < len(hdr) or f[ci["SYMBOL"]] != g or f[ci["LoF"]] != "HC": continue
            cons = f[ci["Consequence"]]
            if not any(k in cons for k in ("stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant")): continue
            try: an = int(float(f[ci["allele_num"]]))
            except ValueError: an = 1
            alts = f[ci["ALT"]].split(","); alt = f[ci["ALT_var"]] if f[ci["ALT_var"]] not in ("", ".") else (alts[an - 1] if an - 1 < len(alts) else alts[0])
            key = "%s:%s:%s:%s" % (c, f[ci["POS"]], f[ci["REF"]], alt)
            gt = f[ci["GT"]].replace("|", "/")
            n_alt = None if "." in gt else sum(1 for a in gt.split("/") if a == str(an))
            try: gq = float(f[ci["GQ"]]); dp = float(f[ci["DP"]])
            except ValueError: gq = dp = 0.0
            try: ar, aa = float(f[ci["AD_ref"]]), float(f[ci["AD_alt"]]); ab = aa / (ar + aa) if ar + aa > 0 else None
            except ValueError: ab = None
            per_var[key][f[ci["SAMPLE"]]] = (n_alt, gq, dp, ab)
            info[key] = (f[ci["gnomAD4.1_joint_AF"]], f[ci["AF"]])
        for key, samp in per_var.items():
            gaf, caf = info[key]
            def fl(v):
                try: return float(v) if v not in ("", ".") else 0.0
                except ValueError: return 0.0
            rare = fl(gaf) < 0.001 and fl(caf) < 0.001
            vsets = ["all", "rare" if rare else "common"]; nvar[(g, name)].add(key)
            fams = collections.defaultdict(list)
            for s in samp:
                if s in ped: fams[ped[s]["fam"]].append(s)
            for fam in fams:
                if name == "SPARK_WES" and (fam, key) in wgs_fams: continue
                for ch in children_by_fam[name].get(fam, []):
                    pa, ma = ped[ch]["pat"], ped[ch]["mat"]
                    if pa not in samp or ma not in samp or ch not in samp: continue
                    gp, gm, gc = samp[pa], samp[ma], samp[ch]
                    def ok(v):
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
                    for st in vsets: res[(g, st, name, parent + "->" + role)][tr] += 1
                    events.append((g, name, key, fam, ch, role, parent, tr, "rare" if rare else "common"))
                    if name == "SPARK_WGS": wgs_fams.add((fam, key))
    print("%s done: variants WGS/WES/SSC = %d/%d/%d" % (g, len(nvar[(g, "SPARK_WGS")]), len(nvar[(g, "SPARK_WES")]), len(nvar[(g, "SSC")])), flush=True)
# ---------------- per-gene table
def pooled(g, st, s): return sum(res[(g, st, n, s)]["T"] for n, _, _ in COHORTS), sum(res[(g, st, n, s)]["NT"] for n, _, _ in COHORTS)
rows = []
with open(OUT + "imprinted_tdt_per_gene.tsv", "w") as fh:
    fh.write(TAB.join(["gene", "expressed_allele", "test_parent", "tier", "s_het", "chip_flag", "variant_set", "n_variants_WGS/WES/SSC", "test_T", "test_NT", "test_frac", "p_test", "silenced_T", "silenced_NT", "sib_test_T", "sib_test_NT", "WGS_T/NT", "WES_T/NT", "SSC_T/NT"]) + "\n")
    for g, c, s0, e0, ea in genes:
        tp = "mother" if ea == "Maternal" else ("father" if ea == "Paternal" else "both")
        for st in sets:
            for parent in (["mother", "father"] if tp == "both" else [tp]):
                other = "father" if parent == "mother" else "mother"
                T, NT = pooled(g, st, parent + "->proband"); sT, sNT = pooled(g, st, other + "->proband"); bT, bNT = pooled(g, st, parent + "->sibling")
                if T + NT == 0: continue
                p = binom_two_sided(T, T + NT)
                row = [g, ea, parent + ("(flagged:both-way)" if tp == "both" else ""), tier(g), "%.4f" % shet[g] if g in shet else "NA", chip_flag(g), st, "%d/%d/%d" % (len(nvar[(g, "SPARK_WGS")]), len(nvar[(g, "SPARK_WES")]), len(nvar[(g, "SSC")])), T, NT, "%.3f" % (T / (T + NT)), "%.3g" % p, sT, sNT, bT, bNT] + ["%d/%d" % (res[(g, st, n, parent + "->proband")]["T"], res[(g, st, n, parent + "->proband")]["NT"]) for n, _, _ in COHORTS]
                fh.write(TAB.join(map(str, row)) + "\n"); rows.append((p, st, g, T, NT, tp))
# ranked list (rare set, single-direction genes first)
ranked = sorted([r for r in rows if r[1] == "rare"], key=lambda r: r[0])
ntests = len(ranked)
with open(OUT + "imprinted_tdt_ranked_rare.tsv", "w") as fh:
    fh.write(TAB.join(["rank", "gene", "test_parent_dir", "T", "NT", "frac", "p", "bonferroni_p", "chip_flag", "tier"]) + "\n")
    for i, (p, st, g, T, NT, tp) in enumerate(ranked, 1): fh.write(TAB.join(map(str, [i, g, tp, T, NT, "%.3f" % (T / (T + NT)), "%.3g" % p, "%.3g" % min(1.0, p * ntests), chip_flag(g), tier(g)])) + "\n")
# ---------------- collapsed burden by tier, expressed-allele transmissions only, single-direction genes, CHIP-flagged excluded
with open(OUT + "imprinted_tdt_tier_burden.tsv", "w") as fh:
    fh.write(TAB.join(["variant_set", "tier", "n_genes", "genes_excluded_CHIP", "test_T", "test_NT", "test_frac", "p", "silenced_T", "silenced_NT", "sib_test_T", "sib_test_NT", "WGS_T/NT", "WES_T/NT", "SSC_T/NT"]) + "\n")
    for st in sets:
        for tr_ in ("lof_t1", "lof_t2", "lof_t3", "all_tiers"):
            gs = [x for x in genes if (tier(x[0]) == tr_ or tr_ == "all_tiers") and x[4] in ("Maternal", "Paternal")]
            excl = [x[0] for x in gs if chip_flag(x[0])]; gs = [x for x in gs if not chip_flag(x[0])]
            T = NT = sT = sNT = bT = bNT = 0; per = {n: [0, 0] for n, _, _ in COHORTS}
            for g, c, s0, e0, ea in gs:
                parent = "mother" if ea == "Maternal" else "father"; other = "father" if parent == "mother" else "mother"
                a, b = pooled(g, st, parent + "->proband"); T += a; NT += b
                a, b = pooled(g, st, other + "->proband"); sT += a; sNT += b
                a, b = pooled(g, st, parent + "->sibling"); bT += a; bNT += b
                for n, _, _ in COHORTS: per[n][0] += res[(g, st, n, parent + "->proband")]["T"]; per[n][1] += res[(g, st, n, parent + "->proband")]["NT"]
            fh.write(TAB.join(map(str, [st, tr_, len(gs), ";".join(excl), T, NT, "%.3f" % (T / (T + NT)) if T + NT else "", "%.3g" % binom_two_sided(T, T + NT) if T + NT else "", sT, sNT, bT, bNT] + ["%d/%d" % tuple(per[n]) for n, _, _ in COHORTS])) + "\n")
with open(OUT + "imprinted_tdt_events.tsv", "w") as fh:
    fh.write(TAB.join(["gene", "cohort", "variant", "family", "child", "role", "parent", "transmitted", "set"]) + "\n")
    for e in events: fh.write(TAB.join(e) + "\n")
with open(OUT + "imprinted_gene_flags.tsv", "w") as fh:
    fh.write("gene\texpressed_allele\ttier\ts_het\tchip_flag\n")
    for g, c, s0, e0, ea in genes: fh.write(TAB.join([g, ea, tier(g), "%.4f" % shet[g] if g in shet else "NA", chip_flag(g)]) + "\n")
print(open(OUT + "imprinted_tdt_tier_burden.tsv").read())
print("top 15 ranked (rare):"); print("".join(open(OUT + "imprinted_tdt_ranked_rare.tsv").readlines()[:16]))
print("CHIP-flagged imprinted genes:", [ (g, chip_flag(g)) for g, *_ in genes if chip_flag(g)])
