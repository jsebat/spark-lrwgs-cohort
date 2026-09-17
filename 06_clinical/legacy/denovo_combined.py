"""Combined de novo priority list (JS): de novo in NDD/SFARI genes UNION
de novo in Tier 1 constrained genes (s_het >= 0.18), across SNV, indel, SV, TR.

"Tier 1 gene" here means the GENE is highly constrained (s_het >= 0.18), so a
de novo missense in such a gene qualifies even though it is not itself a LoF
tier 1 call. The two criteria are a union, not an intersection.
"""
import os
import csv
import glob
import collections

T = "/expanse/lustre/projects/ddp195/jsebat/longread-autism/tiering"
R = "/expanse/lustre/projects/ddp195/jsebat/longread-autism"
GS = "/expanse/projects/sebat1/s3/data/sebat/nf_rare_spark_wes/resources/gene_sets"
TAB, NL = chr(9), chr(10)
SHET_T1 = 0.18


def load_set(fn):
    p = os.path.join(GS, fn)
    if not os.path.exists(p):
        return set()
    return set(x.strip() for x in open(p) if x.strip())


sfari_hc = load_set("sfari_hc.txt")
sfari_all = load_set("sfari_all.txt")
ddg2p_conf = load_set("ddg2p_confident.txt")
ddg2p_all = load_set("ddg2p_all.txt")
ndd = sfari_hc | sfari_all | ddg2p_conf | ddg2p_all
print("NDD/SFARI union gene set: %d genes" % len(ndd))

shet = {}
sp = os.path.join(T, "denovo_sv", "shet_by_symbol.tsv")
for ln in open(sp):
    f = ln.rstrip(NL).split(TAB)
    if len(f) >= 2:
        try:
            shet[f[0]] = float(f[1])
        except ValueError:
            pass
t1genes = set(g for g, v in shet.items() if v >= SHET_T1)
print("Tier 1 genes (s_het >= %.2f): %d" % (SHET_T1, len(t1genes)))


def gene_label(g):
    lab = []
    if g in sfari_hc:
        lab.append("SFARI_hc")
    elif g in sfari_all:
        lab.append("SFARI")
    if g in ddg2p_conf:
        lab.append("DDG2P_conf")
    elif g in ddg2p_all:
        lab.append("DDG2P")
    if g in t1genes:
        lab.append("TIER1_GENE")
    return ";".join(lab)


rows = []

# ---- SNV / indel
p = os.path.join(T, "denovo_tiered.tsv")
n_snv = 0
if os.path.exists(p):
    with open(p) as fh:
        hdr = fh.readline().rstrip(NL).split(TAB)
        ix = dict((h, i) for i, h in enumerate(hdr))
        for ln in fh:
            f = ln.rstrip(NL).split(TAB)
            g = f[ix["gene"]]
            if not g:
                continue
            n_snv += 1
            if g not in ndd and g not in t1genes:
                continue
            csq = f[ix["consequence"]]
            cls = "indel" if (len(f[ix["ref"]]) != len(f[ix["alt"]])) else "SNV"
            rows.append({"class": cls, "family": f[ix["family"]], "proband": f[ix["proband"]],
                         "chrom": f[ix["chrom"]], "start": f[ix["pos"]], "end": f[ix["pos"]],
                         "gene": g, "consequence": csq,
                         "s_het": ("%.4g" % shet[g]) if g in shet else "NA",
                         "lof_tier": f[ix["lof_tier"]], "miss_tier": f[ix["miss_tier"]],
                         "gene_flags": gene_label(g), "size": "", "flag": ""})

# ---- SV (already prioritised + masked + recurrence-flagged)
p = os.path.join(T, "denovo_sv", "denovo_sv_PRIORITIZED_v2.tsv")
if os.path.exists(p):
    with open(p) as fh:
        hdr = fh.readline().rstrip(NL).split(TAB)
        ix = dict((h, i) for i, h in enumerate(hdr))
        for ln in fh:
            f = ln.rstrip(NL).split(TAB)
            g = f[ix["gene"]]
            if g not in ndd and g not in t1genes:
                continue
            rows.append({"class": "SV_" + f[ix["svtype"]], "family": f[ix["family"]],
                         "proband": f[ix["proband"]], "chrom": f[ix["chrom"]],
                         "start": f[ix["start"]], "end": f[ix["end"]], "gene": g,
                         "consequence": "coding_" + f[ix["svtype"]],
                         "s_het": f[ix["s_het"]], "lof_tier": f[ix["lof_tier"]],
                         "miss_tier": "", "gene_flags": gene_label(g),
                         "size": f[ix["svlen"]], "flag": f[ix["flag"]]})

# ---- TR de novo
n_tr = 0
for tsv in sorted(glob.glob(os.path.join(R, "run_*", "analysis", "denovo_tr", "*.denovo_tr.tsv"))):
    if os.environ.get("EXCLUDE_FAMILY", "__none__") in tsv:
        continue
    b = os.path.basename(tsv).split(".")
    fam, pro = b[0], b[1]
    with open(tsv) as fh:
        hdr = fh.readline().rstrip(NL).split(TAB)
        ix = dict((h, i) for i, h in enumerate(hdr))
        for ln in fh:
            f = ln.rstrip(NL).split(TAB)
            n_tr += 1
            g = ""
            for k in ("gene", "SYMBOL", "symbol", "trid"):
                if k in ix and ix[k] < len(f):
                    g = f[ix[k]]
                    break
            if not g:
                continue
            # trid often looks like GENE_chr_start_end; take a leading symbol
            cand = g.split("_")[0]
            if cand in ndd or cand in t1genes:
                rows.append({"class": "TR", "family": fam, "proband": pro,
                             "chrom": f[ix["chrom"]] if "chrom" in ix else "",
                             "start": f[ix["pos"]] if "pos" in ix else "",
                             "end": f[ix["end"]] if "end" in ix else "",
                             "gene": cand, "consequence": "TR_expansion",
                             "s_het": ("%.4g" % shet[cand]) if cand in shet else "NA",
                             "lof_tier": "", "miss_tier": "",
                             "gene_flags": gene_label(cand), "size": "", "flag": ""})

print("de novo SNV/indel scanned: %d ; TR rows scanned: %d" % (n_snv, n_tr))

out = os.path.join(T, "denovo_COMBINED_priority.tsv")
cols = ["class", "family", "proband", "chrom", "start", "end", "size", "gene",
        "consequence", "s_het", "lof_tier", "miss_tier", "gene_flags", "flag"]


def sortkey(r):
    t1 = 0 if "TIER1_GENE" in r["gene_flags"] else 1
    hc = 0 if "SFARI_hc" in r["gene_flags"] or "DDG2P_conf" in r["gene_flags"] else 1
    sv = 0 if r["class"].startswith("SV") else 1
    return (t1, hc, sv, r["family"])


with open(out, "w") as fh:
    fh.write(TAB.join(cols) + NL)
    for r in sorted(rows, key=sortkey):
        fh.write(TAB.join(str(r.get(c, "")) for c in cols) + NL)

print("")
print("wrote %s (%d rows)" % (out, len(rows)))
print("=" * 78)
c = collections.Counter(r["class"] for r in rows)
for k, v in sorted(c.items()):
    print("  %-10s %d" % (k, v))
print("  --")
print("  in a TIER 1 gene (s_het>=0.18): %d" % sum(1 for r in rows if "TIER1_GENE" in r["gene_flags"]))
print("  in SFARI high-confidence      : %d" % sum(1 for r in rows if "SFARI_hc" in r["gene_flags"]))
print("  in DDG2P confident            : %d" % sum(1 for r in rows if "DDG2P_conf" in r["gene_flags"]))
print("=" * 78)
for r in sorted(rows, key=sortkey):
    if r["flag"] and r["flag"] != "pass":
        continue
    print("  %-8s %-11s %-11s %-10s %-26s s_het=%-8s %s" %
          (r["class"], r["family"], r["proband"], r["gene"],
           r["consequence"][:26], r["s_het"], r["gene_flags"]))
