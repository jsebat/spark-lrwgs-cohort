"""Tier every de novo MISSENSE with James's missense tiers.

miss_t1 = n_flag 4, miss_t2 = 3, miss_t3 = 2, miss_t4 = 1, where n_flag counts
how many of these clear their v2 forward-selection thresholds:
    ClinPred_rankscore         >= 0.4298
    AlphaMissense_rankscore    >= 0.9603
    popEVE_converted_rankscore >= 0.9209
    MPC_rankscore              >= 0.8947
A missing score counts as NOT passing (COALESCE(...,FALSE) in the lab SQL);
popEVE coverage is 100% in this build so that case does not arise.

JOIN FIX: dbNSFP stores '#chr' without the 'chr' prefix ('2', not 'chr2').
Joining on the prefixed name silently matched almost nothing.
"""
import duckdb
import os
import collections

# Paths default to the cohort freeze-1 tiering directory; the three environment variables let the same tiers be
# applied to another de novo call set (the phase-aware module's, say) without forking this script.
T = os.environ.get("TIERING_DIR", "/expanse/lustre/projects/ddp195/jsebat/longread-autism/tiering")
IN_TSV = os.environ.get("DENOVO_TSV", "")     # default: <TIERING_DIR>/denovo_tiered.tsv
OUT_TSV = os.environ.get("MISS_OUT", "")      # default: <TIERING_DIR>/denovo_missense_tiered.tsv
D = "/expanse/projects/sebat1/s3/data/sebat/resources/dbNSFP/5.3.1a/parquet_expanded_mane_select"
GS = "/expanse/projects/sebat1/s3/data/sebat/nf_rare_spark_wes/resources/gene_sets"
TAB, NL = chr(9), chr(10)
T_STARS = [("ClinPred_rankscore", 0.4298),
           ("AlphaMissense_rankscore", 0.9603),
           ("popEVE_converted_rankscore", 0.9209),
           ("MPC_rankscore", 0.8947)]


def load_set(fn):
    p = os.path.join(GS, fn)
    return set(x.strip() for x in open(p) if x.strip()) if os.path.exists(p) else set()


sfari_hc, sfari_all = load_set("sfari_hc.txt"), load_set("sfari_all.txt")
ddg2p_conf, ddg2p_all = load_set("ddg2p_confident.txt"), load_set("ddg2p_all.txt")
ndd = sfari_hc | sfari_all | ddg2p_conf | ddg2p_all

shet = {}
sp = os.path.join(T, "denovo_sv", "shet_by_symbol.tsv")
for ln in open(sp):
    f = ln.rstrip(NL).split(TAB)
    if len(f) >= 2:
        try:
            shet[f[0]] = float(f[1])
        except ValueError:
            pass

# ---- de novo missense
rows = []
with open(IN_TSV or os.path.join(T, "denovo_tiered.tsv")) as fh:
    hdr = fh.readline().rstrip(NL).split(TAB)
    ix = dict((h, i) for i, h in enumerate(hdr))
    for ln in fh:
        f = ln.rstrip(NL).split(TAB)
        if "missense" not in f[ix["consequence"]]:
            continue
        rows.append({"family": f[ix["family"]], "proband": f[ix["proband"]],
                     "chrom": f[ix["chrom"]], "pos": int(f[ix["pos"]]),
                     "ref": f[ix["ref"]], "alt": f[ix["alt"]],
                     "gene": f[ix["gene"]], "csq": f[ix["consequence"]]})
print("de novo missense variants: %d" % len(rows))

# ---- join scores, chrom prefix stripped
con = duckdb.connect()
con.execute("PRAGMA threads=4")
by_chrom = collections.defaultdict(list)
for r in rows:
    by_chrom[r["chrom"]].append(r)

scols = ", ".join('"%s"' % c for c, _ in T_STARS)
found = 0
for chrom, rs in sorted(by_chrom.items()):
    pf = os.path.join(D, "%s.parquet" % chrom)
    if not os.path.exists(pf):
        print("  no parquet for %s" % chrom)
        continue
    bare = chrom.replace("chr", "")
    pos = sorted(set(r["pos"] for r in rs))
    q = ('SELECT "pos(1-based)" AS p, ref, alt, %s FROM read_parquet(?) '
         'WHERE "#chr" = ? AND "pos(1-based)" IN (%s)'
         % (scols, ",".join(str(x) for x in pos)))
    hit = {}
    for rec in con.execute(q, [pf, bare]).fetchall():
        p, rr, aa = rec[0], rec[1], rec[2]
        hit[(int(p), rr, aa)] = rec[3:]
    for r in rs:
        v = hit.get((r["pos"], r["ref"], r["alt"]))
        if v is None:
            v = next((val for (pp, rr, aa), val in hit.items() if pp == r["pos"]), None)
        r["scores"] = v
        if v:
            found += 1
print("missense with dbNSFP scores: %d of %d" % (found, len(rows)))

# ---- n_flag and tier
out = OUT_TSV or os.path.join(T, "denovo_missense_tiered.tsv")
cols = ["family", "proband", "chrom", "pos", "ref", "alt", "gene", "n_flag", "miss_tier",
        "ClinPred", "AlphaMissense", "popEVE", "MPC", "s_het", "gene_sets"]


def label(g):
    L = []
    if g in sfari_hc:
        L.append("SFARI_hc")
    elif g in sfari_all:
        L.append("SFARI")
    if g in ddg2p_conf:
        L.append("DDG2P_conf")
    elif g in ddg2p_all:
        L.append("DDG2P")
    if shet.get(g, 0) >= 0.18:
        L.append("TIER1_GENE")
    return ";".join(L) or "-"


tal = collections.Counter()
with open(out, "w") as fh:
    fh.write(TAB.join(cols) + NL)
    recs = []
    for r in rows:
        v = r.get("scores")
        nf = 0
        vals = []
        for i, (name, thr) in enumerate(T_STARS):
            x = v[i] if v else None
            try:
                fx = float(x)
                vals.append("%.4f" % fx)
                if fx >= thr:
                    nf += 1
            except (TypeError, ValueError):
                vals.append("NA")
        tier = {4: "miss_t1", 3: "miss_t2", 2: "miss_t3", 1: "miss_t4"}.get(nf, "none")
        tal[tier] += 1
        recs.append((nf, r, vals, tier))
    for nf, r, vals, tier in sorted(recs, key=lambda x: -x[0]):
        fh.write(TAB.join([r["family"], r["proband"], r["chrom"], str(r["pos"]),
                           r["ref"][:20], r["alt"][:20], r["gene"], str(nf), tier] +
                          vals + [("%.4g" % shet[r["gene"]]) if r["gene"] in shet else "NA",
                                  label(r["gene"])]) + NL)

print("")
print("wrote %s" % out)
print("=" * 76)
for t in ("miss_t1", "miss_t2", "miss_t3", "miss_t4", "none"):
    if tal[t]:
        print("  %-9s %d" % (t, tal[t]))
print("=" * 76)
print("  %-9s %-11s %-11s %-9s %-6s %-8s %-8s %-8s %-8s %-8s %s" %
      ("tier", "family", "proband", "gene", "nflag", "ClinPred", "AlphaMis", "popEVE", "MPC", "s_het", "sets"))
for nf, r, vals, tier in sorted(recs, key=lambda x: -x[0]):
    print("  %-9s %-11s %-11s %-9s %-6d %-8s %-8s %-8s %-8s %-8s %s" %
          (tier, r["family"], r["proband"], r["gene"][:9], nf,
           vals[0], vals[1], vals[2], vals[3],
           ("%.4g" % shet[r["gene"]]) if r["gene"] in shet else "NA", label(r["gene"])))
