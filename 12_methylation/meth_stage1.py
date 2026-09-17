#!/usr/bin/env python3
"""Methylation stage 1 (D74). Stdlib only; tabix/bcftools via container.
A. Imprinted-gene coordinates (GENCODE v44) and promoter CpG islands.
B. Cohort CpG-island methylation matrix from methbat profiles (combined, hap1, hap2, ASM p) for all samples.
C. Haplotype-convention and parent-of-origin test at known imprinted control regions in every trio child.
D. LoF variants in imprinted genes (LOFTEE HC) across the cohort.
E. Marker means on Loyfer U25 atlas regions per sample (NNLS deconvolution done in deconv.R).
"""
import os, sys, glob, gzip, csv, subprocess, collections, math, re, bisect
L = os.environ["DATA_ROOT"]; OUT = L + "/meth"; os.makedirs(OUT, exist_ok=True)
REF = os.environ["METH_REFS"]
GTF = os.environ["GENCODE_GTF"]
MAN = os.environ["MANIFEST"]
SIF = os.environ["SIF"]; sx = ["singularity", "exec", "-B", "/expanse:/expanse", SIF]
TAB = "\t"; NL = "\n"
man = {r["sample_id"]: r for r in csv.DictReader(open(MAN), delimiter=TAB)}
def run(cmd): return subprocess.run(cmd, stdout=subprocess.PIPE, universal_newlines=True).stdout

# ---------------- A. imprinted genes -> coordinates
imp = {r["gene"]: r for r in csv.DictReader(open(REF + "/geneimprint_human_imprinted.tsv"), delimiter=TAB)}
alias = {"AIM1": "CRYBG1", "IGF2AS": "IGF2-AS", "TCEB3C": "ELOA3", "PWCR1": "SNORD116-1", "SNORD115": "SNORD115-1", "ZNF127AS": "MKRN3"}
want = set(imp) | set(alias.values())
genes = {}
with gzip.open(GTF, "rt") as fh:
    for ln in fh:
        if ln.startswith("#"): continue
        f = ln.split(TAB)
        if f[2] != "gene": continue
        m = re.search(r'gene_name "([^"]+)"', f[8]); g = m.group(1) if m else ""
        if g in want: genes[g] = (f[0], int(f[3]), int(f[4]), f[6])
found = {g: genes.get(g) or genes.get(alias.get(g, "")) for g in imp}
with open(OUT + "/imprinted_genes.bed", "w") as fh:
    for g, v in sorted(found.items(), key=lambda x: (x[1] is None, x[0])):
        if v: fh.write(TAB.join([v[0], str(v[1] - 1), str(v[2]), g, imp[g]["expressed_allele"], v[3]]) + NL)
print("imprinted genes with coordinates: %d of %d; missing: %s" % (sum(1 for v in found.values() if v), len(imp), [g for g, v in found.items() if not v]), flush=True)

# ---------------- B. methbat island matrix
prof = {}
for p in glob.glob(L + "/run_*/*/out/methbat_profile/*/*.methbat.profile.tsv"):
    if "_LAST" in p: continue
    sid = os.path.basename(p).split(".")[0]
    if sid in man: prof[sid] = p          # manifest samples only (excludes unrelated runs)
print("methbat profiles:", len(prof), flush=True)
islands = None; M = {}
for s, p in sorted(prof.items()):
    rows = []; keys = []
    with open(p) as fh:
        rd = csv.DictReader((l for l in fh if not l.startswith("#")), delimiter=TAB)
        for r in rd:
            keys.append((r["chrom"], int(r["start"]), int(r["end"]), r["cpg_label"]))
            rows.append((r.get("mean_combined_methyl", ""), r.get("mean_hap1_methyl", ""), r.get("mean_hap2_methyl", ""), r.get("asm_fishers_pvalue", ""), r.get("summary_label", "")))
    if islands is None: islands = keys
    elif keys != islands: print("WARNING island order differs for", s)
    M[s] = rows
samples = sorted(M)
with open(OUT + "/islands.bed", "w") as fh:
    for k in islands: fh.write(TAB.join([k[0], str(k[1]), str(k[2]), k[3]]) + NL)
for col, name in ((0, "combined"), (1, "hap1"), (2, "hap2"), (3, "asm_p"), (4, "label")):
    with open(OUT + "/island_%s.tsv" % name, "w") as fh:
        fh.write(TAB.join(["chrom", "start", "end", "cpg_label"] + samples) + NL)
        for i, k in enumerate(islands):
            fh.write(TAB.join([k[0], str(k[1]), str(k[2]), k[3]] + [M[s][i][col] for s in samples]) + NL)
print("island matrix written: %d islands x %d samples" % (len(islands), len(samples)), flush=True)
isl_by_chr = collections.defaultdict(list)
for i, k in enumerate(islands): isl_by_chr[k[0]].append((k[1], k[2], i))
prom = {}
for g, v in found.items():
    if not v: continue
    c, s0, e0, strand = v; tss = s0 if strand == "+" else e0
    prom[g] = [i for (a, b, i) in isl_by_chr[c] if a <= tss + 2000 and b >= tss - 2000]
with open(OUT + "/imprinted_promoter_islands.tsv", "w") as fh:
    fh.write(TAB.join(["gene", "expressed_allele", "n_islands", "islands"]) + NL)
    for g, h in prom.items(): fh.write(TAB.join([g, imp[g]["expressed_allele"], str(len(h)), ";".join("%s:%d-%d" % islands[i][:3] for i in h)]) + NL)
print("imprinted genes with a promoter island: %d" % sum(1 for h in prom.values() if h), flush=True)

# ---------------- C. haplotype convention + parent of origin at known imprinted control regions
ICR = [("SNRPN/SNURF", "chr15", 24954000, 24956500, "maternal"), ("H19 ICR", "chr11", 1998000, 2003000, "paternal"), ("KCNQ1OT1", "chr11", 2697000, 2701000, "maternal"),
       ("MEG3/DLK1 IG-DMR", "chr14", 100810000, 100830000, "paternal"), ("PEG3", "chr19", 56837000, 56842000, "maternal"), ("MEST", "chr7", 130485000, 130494000, "maternal"),
       ("PLAGL1", "chr6", 144006000, 144010000, "maternal"), ("GRB10", "chr7", 50780000, 50784000, "maternal"), ("PEG10/SGCE", "chr7", 94655000, 94660000, "maternal"),
       ("NNAT", "chr20", 37518000, 37522000, "maternal"), ("GNAS NESP/XL", "chr20", 58838000, 58856000, "maternal"), ("ZDBF2/GPR1", "chr2", 206249000, 206260000, "paternal"),
       ("L3MBTL1", "chr20", 43513000, 43520000, "maternal"), ("NAP1L5", "chr4", 88696000, 88699000, "maternal"), ("INPP5F", "chr10", 119811000, 119815000, "maternal"),
       ("FAM50B", "chr6", 3849000, 3852000, "maternal"), ("ZNF597/NAA60", "chr16", 3442000, 3446000, "paternal"), ("DIRAS3", "chr1", 67787000, 67793000, "maternal")]
children = [s for s in samples if man.get(s, {}).get("father_id", "0") != "0" and man[s]["father_id"] in samples and man[s]["mother_id"] in samples]
vcf_of = {}
for fam in set(man[s]["family_id"] for s in children):
    v = os.path.join(os.environ["PHASED_VCF_DIR"], "%s.phased.vcf.gz" % fam)
    if os.path.exists(v) and os.path.exists(v + ".tbi"): vcf_of[fam] = v
print("phased VCFs with index:", len(vcf_of), "of", len(set(man[s]["family_id"] for s in children)), flush=True)
ICR_BED = OUT + "/icr_regions.bed"
with open(ICR_BED, "w") as fh:
    for name, c, a, b, exp in sorted(ICR, key=lambda x: (x[1], x[2])): fh.write(TAB.join([c, str(a), str(b), name]) + NL)
_hap_cache = {}
def hap_means_all(s):
    if s in _hap_cache: return _hap_cache[s]
    res = collections.defaultdict(dict)
    for h in ("hap1", "hap2"):
        p = [x for x in glob.glob(os.path.join(os.environ["CPG_DIR"], "%s.%s.bed.gz" % (s, h))) if "_LAST" not in x]
        if not p: _hap_cache[s] = None; return None
        txt = run(sx + ["tabix", "-R", ICR_BED, sorted(p)[-1]])
        acc = collections.defaultdict(list)
        for l in txt.splitlines():
            f = l.split(TAB)
            if len(f) < 6 or int(f[5]) < 4: continue
            pos = int(f[1])
            for name, c, a, b, exp in ICR:
                if f[0] == c and a <= pos <= b: acc[name].append(float(f[3]))
        for name, c, a, b, exp in ICR:
            v = acc.get(name, []); res[name][h] = (sum(v) / len(v) / 100.0, len(v)) if v else (None, 0)
    _hap_cache[s] = res; return res
def maternal_hap(child, c, a, b):
    """fraction of informative phased het sites (same phase block as the region) whose hap1 allele is maternal"""
    fam = man[child]["family_id"]; vcf = vcf_of.get(fam)
    if not vcf: return (None, 0)
    txt = run(sx + ["bcftools", "query", "-s", ",".join([child, man[child]["father_id"], man[child]["mother_id"]]), "-r", "%s:%d-%d" % (c, max(1, a - 300000), b + 300000), "-i", 'TYPE="snp"', "-f", "%POS[" + TAB + "%GT:%PS]" + NL, vcf])
    ps_in = set(); rows = []
    for ln in txt.splitlines():
        f = ln.split(TAB)
        if len(f) < 4: continue
        pos = int(f[0]); cg = f[1].split(":")[0]; fg = f[2].split(":")[0]; mg = f[3].split(":")[0]; ps = f[1].split(":")[1] if ":" in f[1] else "."
        rows.append((pos, cg, fg, mg, ps))
        if a <= pos <= b and "|" in cg and cg[0] != cg[2]: ps_in.add(ps)
    if not ps_in: return (None, 0)
    n = 0; k = 0
    for pos, cg, fg, mg, ps in rows:
        if ps not in ps_in or "|" not in cg or cg[0] == cg[2]: continue
        h1, h2 = cg[0], cg[2]
        if "." in fg or "." in mg: continue
        fa = set(fg.replace("|", "/").split("/")); mo = set(mg.replace("|", "/").split("/"))
        if h1 in mo and h1 not in fa and h2 in fa: k += 1; n += 1
        elif h2 in mo and h2 not in fa and h1 in fa: n += 1
        elif h1 in fa and h1 not in mo and h2 in mo: n += 1
        elif h2 in fa and h2 not in mo and h1 in mo: k += 1; n += 1
    return ((k / n) if n else None, n)
tally = collections.Counter()
with open(OUT + "/icr_parent_of_origin_check.tsv", "w") as fh:
    fh.write(TAB.join(["child", "family", "icr", "expected_methylated", "hap1_meth", "hap2_meth", "n_cpg_h1", "n_cpg_h2", "frac_hap1_maternal", "n_informative", "methylated_hap", "methylated_parent", "concordant"]) + NL)
    for ch in children:
        allm = hap_means_all(ch)
        if allm is None: continue
        for name, c, a, b, exp in ICR:
            hm = allm.get(name)
            if not hm or hm["hap1"][0] is None or hm["hap2"][0] is None: tally["no_hap_data"] += 1; continue
            h1, h2 = hm["hap1"][0], hm["hap2"][0]
            methhap = "hap1" if h1 - h2 > 0.4 else ("hap2" if h2 - h1 > 0.4 else "none")
            fm, n = maternal_hap(ch, c, a, b)
            par = ""
            if methhap != "none" and fm is not None and n >= 3 and (fm >= 0.9 or fm <= 0.1):
                hap1_is_mat = fm >= 0.9
                par = ("maternal" if hap1_is_mat else "paternal") if methhap == "hap1" else ("paternal" if hap1_is_mat else "maternal")
            conc = "" if not par else ("yes" if par == exp else "NO")
            tally[conc or ("no_ASM" if methhap == "none" else "no_phase")] += 1
            fh.write(TAB.join(map(str, [ch, man[ch]["family_id"], name, exp, round(h1, 3), round(h2, 3), hm["hap1"][1], hm["hap2"][1], "" if fm is None else round(fm, 2), n, methhap, par, conc])) + NL)
print("ICR parent-of-origin check:", dict(tally), flush=True)

# ---------------- D. LoF variants in imprinted genes
T = L + "/tiering"; gset = set(found) | set(alias.values())
def scan(path, label):
    rows = []
    if not os.path.exists(path): print("missing", path); return rows
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter=TAB):
            sym = r.get("symbol") or r.get("gene") or ""
            if sym in gset: rows.append(r)
    print("%s: %d rows in imprinted genes" % (label, len(rows)), flush=True)
    return rows
lof_rare = scan(T + "/lof_tiered.tsv", "lof_tiered (rare)"); lof_all = scan(T + "/lof_sites.tsv", "lof_sites (all)")
with open(OUT + "/imprinted_lof_sites.tsv", "w") as fh:
    cols = sorted(set(k for r in lof_rare + lof_all for k in r.keys()))
    fh.write(TAB.join(["source"] + cols) + NL)
    for src, rows in (("rare", lof_rare), ("all", lof_all)):
        for r in rows: fh.write(TAB.join([src] + [str(r.get(c, "")) for c in cols]) + NL)

# ---------------- E. marker means on Loyfer U25 regions (one tabix -R per sample)
atlas = list(csv.DictReader(open(REF + "/Atlas.U25.l4.hg38.tsv"), delimiter=TAB))
ATLAS_BED = OUT + "/atlas_U25.bed"
with open(ATLAS_BED, "w") as fh:
    for i, r in enumerate(atlas): fh.write(TAB.join([r["chr"], r["start"], r["end"], str(i)]) + NL)
by_chr = collections.defaultdict(list)
for i, r in enumerate(atlas): by_chr[r["chr"]].append((int(r["start"]), int(r["end"]), i))
for c in by_chr: by_chr[c].sort()
starts = {c: [a for a, b, i in v] for c, v in by_chr.items()}
with open(OUT + "/marker_means_U25.tsv", "w") as fh:
    fh.write(TAB.join(["sample", "family", "role", "dna_blood"] + ["m%d" % i for i in range(len(atlas))]) + NL)
    for s in samples:
        p = [x for x in glob.glob(os.path.join(os.environ["CPG_DIR"], "%s.combined.bed.gz" % s)) if "_LAST" not in x]
        if not p: continue
        txt = run(sx + ["tabix", "-R", ATLAS_BED, sorted(p)[-1]])
        acc = collections.defaultdict(list)
        for l in txt.splitlines():
            f = l.split(TAB)
            if len(f) < 6 or int(f[5]) < 4: continue
            c = f[0]; pos = int(f[1]); lst = by_chr.get(c)
            if not lst: continue
            k = bisect.bisect_right(starts[c], pos) - 1
            if k >= 0 and lst[k][0] <= pos <= lst[k][1]: acc[lst[k][2]].append(float(f[3]) / 100.0)
        vals = []
        for i in range(len(atlas)):
            v = acc.get(i); vals.append(("%.4f" % (sum(v) / len(v))) if v else "NA")
        fh.write(TAB.join([s, man[s]["family_id"], man[s]["role"], "1" if s.startswith("REACH") else "0"] + vals) + NL)
        print("markers", s, "regions with data:", len(acc), flush=True)
print("stage 1 done")
