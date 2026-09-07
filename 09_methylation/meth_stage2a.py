#!/usr/bin/env python3
"""Stage 2a (D74): cohort allele-specific-methylation (ASM) island catalogue with parent of origin.
Per child: ASM islands = methbat asm_fishers_pvalue < 1e-3 and |hap1 - hap2| >= 0.30 (autosomes).
Parent of origin of the methylated haplotype from the NEAREST informative phased SNVs (+/- 20 kb, >= 2 sites,
>= 90% agreement), per D74 addendum 2. Aggregate across children: imprinted-like islands (>= 5 assigned children,
>= 90% same parent), sequence-dependent ASM (frequent ASM, parent random), and a loss-of-imprinting screen at
imprinted-like islands in every genome with data."""
import os, csv, subprocess, collections
ME = os.environ["DATA_ROOT"] + "/meth/"; MAN = os.environ["MANIFEST"]
SIF = os.environ["SIF"]; H = os.path.expanduser("~"); sx = ["singularity", "exec", "-B", "/expanse:/expanse", "-B", H + ":" + H, SIF]
TAB = "\t"; NL = "\n"
man = {r["sample_id"]: r for r in csv.DictReader(open(MAN), delimiter=TAB)}
def load(name):
    with open(ME + "island_%s.tsv" % name) as fh:
        rd = csv.reader(fh, delimiter=TAB); hdr = next(rd); samples = hdr[4:]; rows = [r for r in rd]
    return samples, rows
samples, H1 = load("hap1"); _, H2 = load("hap2"); _, AP = load("asm_p"); _, CB = load("combined")
islands = [(r[0], int(r[1]), int(r[2]), r[3]) for r in H1]
def num(v):
    try: return float(v)
    except ValueError: return None
children = [s for s in samples if man[s]["father_id"] in samples and man[s]["mother_id"] in samples]
print("islands", len(islands), "samples", len(samples), "children", len(children), flush=True)
# imprinted genes for annotation
imp = []
for ln in open(ME + "imprinted_genes.bed"):
    f = ln.rstrip(NL).split(TAB); imp.append((f[0], int(f[1]), int(f[2]), f[3], f[4]))
def near_imprinted(c, a, b, w=50000):
    return ";".join("%s(%s)" % (g, e) for cc, s, e_, g, e in imp if cc == c and s - w <= b and e_ + w >= a)

# ---------------- per-child ASM islands and parent of origin
asm_calls = collections.defaultdict(dict)   # island idx -> child -> (methylated_hap, parent or "", n_local, switch)
per_child = {}
for ch in children:
    ci = samples.index(ch); fam = man[ch]["family_id"]
    idx = []
    for i, (c, a, b, lab) in enumerate(islands):
        if c in ("chrX", "chrY", "chrM"): continue
        h1, h2, p = num(H1[i][4 + ci]), num(H2[i][4 + ci]), num(AP[i][4 + ci])
        if h1 is None or h2 is None or p is None: continue
        if p < 1e-3 and abs(h1 - h2) >= 0.30: idx.append(i)
    per_child[ch] = len(idx)
    if not idx: continue
    bed = ME + "_tmp_%s.bed" % ch
    with open(bed, "w") as fh:
        for i in idx: c, a, b, _ = islands[i]; fh.write(TAB.join([c, str(max(0, a - 20000)), str(b + 20000)]) + NL)
    vcf = os.path.join(os.environ["PHASED_VCF_DIR"], "%s.phased.vcf.gz" % fam)
    txt = subprocess.run(sx + ["bcftools", "query", "-s", ",".join([ch, man[ch]["father_id"], man[ch]["mother_id"]]), "-R", bed, "-i", 'TYPE="snp"', "-f", "%CHROM" + TAB + "%POS[" + TAB + "%GT]" + NL, vcf], stdout=subprocess.PIPE, universal_newlines=True).stdout
    os.remove(bed)
    inf = collections.defaultdict(list)   # chrom -> [(pos, 'mat'|'pat')]
    for ln in txt.splitlines():
        f = ln.split(TAB)
        if len(f) < 5: continue
        c, pos, cg, fg, mg = f[0], int(f[1]), f[2], f[3].replace("|", "/"), f[4].replace("|", "/")
        if "|" not in cg or cg[0] == cg[2] or "." in fg or "." in mg: continue
        h1a, h2a = cg[0], cg[2]; fa = set(fg.split("/")); mo = set(mg.split("/"))
        if h1a in mo and h1a not in fa and h2a in fa: inf[c].append((pos, "mat"))
        elif h2a in mo and h2a not in fa and h1a in fa: inf[c].append((pos, "pat"))
        elif h1a in fa and h1a not in mo and h2a in mo: inf[c].append((pos, "pat"))
        elif h2a in fa and h2a not in mo and h1a in mo: inf[c].append((pos, "mat"))
    for c in inf: inf[c].sort()
    n_assigned = 0
    for i in idx:
        c, a, b, _ = islands[i]; h1, h2 = num(H1[i][4 + ci]), num(H2[i][4 + ci])
        mh = "hap1" if h1 > h2 else "hap2"
        loc = [k for p, k in inf.get(c, []) if a - 20000 <= p <= b + 20000]
        par = ""
        if len(loc) >= 2:
            fm = sum(1 for k in loc if k == "mat") / len(loc)
            if fm >= 0.9 or fm <= 0.1:
                hap1_mat = fm >= 0.9; par = ("maternal" if hap1_mat else "paternal") if mh == "hap1" else ("paternal" if hap1_mat else "maternal"); n_assigned += 1
        asm_calls[i][ch] = (mh, par, len(loc))
    print("child %s: ASM islands %d, parent assigned %d" % (ch, len(idx), n_assigned), flush=True)

# ---------------- aggregate
cat = open(ME + "asm_island_catalog.tsv", "w")
cat.write(TAB.join(["chrom", "start", "end", "cpg_label", "n_children_with_data", "n_children_ASM", "n_maternal_methylated", "n_paternal_methylated", "n_unassigned", "class", "imprinted_genes_within_50kb", "n_parents_ASM", "n_parents_with_data"]) + NL)
imprinted_like = []
for i, (c, a, b, lab) in enumerate(islands):
    if c in ("chrX", "chrY", "chrM"): continue
    calls = asm_calls.get(i, {})
    nd = sum(1 for ch in children if num(H1[i][4 + samples.index(ch)]) is not None and num(H2[i][4 + samples.index(ch)]) is not None)
    nm = sum(1 for v in calls.values() if v[1] == "maternal"); npat = sum(1 for v in calls.values() if v[1] == "paternal"); nu = sum(1 for v in calls.values() if v[1] == "")
    parents_asm = 0; parents_data = 0
    for s in samples:
        if s in children or man[s]["role"] != "parent": continue
        si = samples.index(s); h1, h2, p = num(H1[i][4 + si]), num(H2[i][4 + si]), num(AP[i][4 + si])
        if h1 is None or h2 is None: continue
        parents_data += 1
        if p is not None and p < 1e-3 and abs(h1 - h2) >= 0.30: parents_asm += 1
    nasg = nm + npat; klass = ""
    if nasg >= 5 and max(nm, npat) / nasg >= 0.9: klass = "imprinted_like_" + ("maternal_methylated" if nm > npat else "paternal_methylated")
    elif len(calls) >= 5 and nasg >= 4 and 0.25 <= nm / nasg <= 0.75: klass = "sequence_dependent_ASM"
    elif len(calls) >= 5: klass = "frequent_ASM_unclassified"
    elif len(calls) >= 1: klass = "sporadic_ASM"
    if not calls and parents_asm == 0: continue
    cat.write(TAB.join(map(str, [c, a, b, lab, nd, len(calls), nm, npat, nu, klass, near_imprinted(c, a, b), parents_asm, parents_data])) + NL)
    if klass.startswith("imprinted_like"): imprinted_like.append((i, klass))
cat.close()
# ---------------- loss-of-imprinting screen at imprinted-like islands: genomes with data but no ASM and extreme combined methylation
with open(ME + "loi_candidates.tsv", "w") as fh:
    fh.write(TAB.join(["chrom", "start", "end", "cpg_label", "class", "sample", "role", "hap1", "hap2", "combined", "pattern"]) + NL)
    n = 0
    for i, klass in imprinted_like:
        c, a, b, lab = islands[i]
        for s in samples:
            si = samples.index(s); h1, h2, cb = num(H1[i][4 + si]), num(H2[i][4 + si]), num(CB[i][4 + si])
            if h1 is None or h2 is None or cb is None: continue
            if abs(h1 - h2) < 0.15 and (cb <= 0.2 or cb >= 0.8):
                fh.write(TAB.join(map(str, [c, a, b, lab, klass, s, man[s]["role"], round(h1, 3), round(h2, 3), round(cb, 3), "both_unmethylated" if cb <= 0.2 else "both_methylated"])) + NL); n += 1
print("imprinted-like islands:", len(imprinted_like), " LOI candidate genome x island rows:", n, flush=True)
with open(ME + "asm_per_child.tsv", "w") as fh:
    fh.write("child\tn_ASM_islands\n"); [fh.write("%s\t%d\n" % (ch, n)) for ch, n in per_child.items()]
cls = collections.Counter(l.split(TAB)[9] for l in open(ME + "asm_island_catalog.tsv") if not l.startswith("chrom"))
print("catalog classes:", dict(cls)); print("stage 2a done")
