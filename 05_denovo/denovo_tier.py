"""Tier every de novo SNV/indel call in the long-read cohort.

Main result 1 (JS): all Tier 1 de novo LoF and missense, combined across
variant classes. This script does the SNV/indel half; SV and TR follow.

Tiers (JS 2026-09-05):
  lof_t1   LOFTEE HC and s_het >= 0.18
  lof_t2   LOFTEE HC and 0.03 <= s_het < 0.18
  lof_t3   remaining HC LoF, INCLUDING genes with no GeneBayes score
  miss_t1..t4  by n_flag = count of 4 rankscore thresholds met
Filters: lab standard segdup/simpleRepeat/rmsk mask; gnomAD global AF < 0.001.
"""
import duckdb
import glob
import os
import sys
import gzip

T = "/expanse/lustre/projects/ddp195/jsebat/longread-autism/tiering"
RUNROOT = "/expanse/lustre/projects/ddp195/jsebat/longread-autism"
DBNSFP = "/expanse/projects/sebat1/s3/data/sebat/resources/dbNSFP/5.3.1a/parquet_expanded_mane_select"
GB = "/expanse/projects/sebat1/s3/data/sebat/g2mh/scripts/scripts_for_rare_pipeline/resources/GeneBayes/output/Supplementary_Table_1.tsv"
MASK = os.path.join(T, "mask.segdup_repeat.bed.gz")
AF_MAX = 0.001
T_STARS = {"ClinPred_rankscore": 0.4298,
           "AlphaMissense_rankscore": 0.9603,
           "popEVE_converted_rankscore": 0.9209,
           "MPC_rankscore": 0.8947}

# ---- 1. collect de novo calls from every family
rows = []
for tsv in sorted(glob.glob(os.path.join(RUNROOT, "run_*", "analysis", "denovo", "*.denovo.hiconf.tsv"))):
    if os.environ.get("EXCLUDE_FAMILY", "__none__") in tsv:
        continue
    base = os.path.basename(tsv)
    parts = base.split(".")
    fam, proband = parts[0], parts[1]
    with open(tsv) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = dict((h, i) for i, h in enumerate(hdr))
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            af = f[ix["gnomad_af"]]
            try:
                afv = float(af)
            except ValueError:
                afv = -1.0
            rows.append({"family": fam, "proband": proband,
                         "chrom": f[ix["chrom"]], "pos": int(f[ix["pos"]]),
                         "ref": f[ix["ref"]], "alt": f[ix["alt"]],
                         "gene": f[ix["gene"]], "csq": f[ix["consequence"]],
                         "gnomad_af": afv})
print("de novo calls collected: %d from %d families" %
      (len(rows), len(set(r["family"] for r in rows))))

# ---- 2. rare filter (absent from gnomAD counts as rare)
rows = [r for r in rows if r["gnomad_af"] < AF_MAX]
print("after gnomAD global AF < %s: %d" % (AF_MAX, len(rows)))

# ---- 3. segdup / repeat mask
mask = {}
with gzip.open(MASK, "rt") as fh:
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        mask.setdefault(f[0], []).append((int(f[1]), int(f[2])))
for c in mask:
    mask[c].sort()


def masked(chrom, pos):
    import bisect
    iv = mask.get(chrom)
    if not iv:
        return False
    starts = [a for a, b in iv]
    i = bisect.bisect_right(starts, pos - 1) - 1
    while i >= 0 and i < len(iv):
        a, b = iv[i]
        if a < pos <= b:
            return True
        if b < pos - 1:
            break
        i -= 1
    return False


before = len(rows)
rows = [r for r in rows if not masked(r["chrom"], r["pos"])]
print("after segdup/repeat mask: %d (dropped %d)" % (len(rows), before - len(rows)))

# ---- 4. LoF tier from the cohort LoF table
lof = {}
tf = os.path.join(T, "lof_tiered.tsv")
if os.path.exists(tf):
    with open(tf) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = dict((h, i) for i, h in enumerate(hdr))
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            lof[(f[ix["chrom"]], int(f[ix["pos"]]))] = (f[ix["lof_tier"]], f[ix["s_het"]], f[ix["symbol"]])
print("cohort LoF sites available for join: %d" % len(lof))

# ---- 5. missense scores from dbNSFP
con = duckdb.connect()
con.execute("PRAGMA threads=4")
by_chrom = {}
for r in rows:
    by_chrom.setdefault(r["chrom"], []).append(r)
scores = {}
for chrom, rs in sorted(by_chrom.items()):
    pf = os.path.join(DBNSFP, "%s.parquet" % chrom)
    if not os.path.exists(pf):
        continue
    pos = sorted(set(r["pos"] for r in rs))
    cols = [r[0] for r in con.execute("DESCRIBE SELECT * FROM read_parquet('%s')" % pf).fetchall()]
    poscol = "pos(1-based)" if "pos(1-based)" in cols else ("pos" if "pos" in cols else None)
    refc = "ref" if "ref" in cols else None
    altc = "alt" if "alt" in cols else None
    if poscol is None:
        print("  %s: no position column found; skipping missense scores" % chrom)
        continue
    q = ('SELECT "%s" AS p, %s, %s, "ClinPred_rankscore", "AlphaMissense_rankscore", '
         '"popEVE_converted_rankscore", "MPC_rankscore" FROM read_parquet(?) '
         'WHERE "%s" IN (%s)') % (poscol,
                                  ('"%s" AS r' % refc) if refc else "NULL AS r",
                                  ('"%s" AS a' % altc) if altc else "NULL AS a",
                                  poscol, ",".join(str(p) for p in pos))
    try:
        for rec in con.execute(q, [pf]).fetchall():
            p, rr, aa, cp, am, pe, mpc = rec
            scores[(chrom, int(p), rr, aa)] = (cp, am, pe, mpc)
    except Exception as e:
        print("  %s: score query failed: %s" % (chrom, str(e)[:120]))
print("missense score records matched: %d" % len(scores))


def nflag(chrom, pos, ref, alt):
    v = scores.get((chrom, pos, ref, alt))
    if v is None:
        v = None
        for k in ((chrom, pos, None, None),):
            v = scores.get(k)
            if v:
                break
    if v is None:
        return None
    cp, am, pe, mpc = v
    n = 0
    for val, key in ((cp, "ClinPred_rankscore"), (am, "AlphaMissense_rankscore"),
                     (pe, "popEVE_converted_rankscore"), (mpc, "MPC_rankscore")):
        try:
            if val is not None and float(val) >= T_STARS[key]:
                n += 1
        except (TypeError, ValueError):
            pass
    return n


# ---- 6. assign tiers and write
out = os.path.join(T, "denovo_tiered.tsv")
cols = ["family", "proband", "chrom", "pos", "ref", "alt", "gene", "consequence",
        "gnomad_af", "lof_tier", "s_het", "miss_nflag", "miss_tier", "tier"]
n_t1_lof = n_t1_miss = 0
with open(out, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        lt, sh, sym = lof.get((r["chrom"], r["pos"]), ("", "", ""))
        mt = ""
        nf = ""
        if "missense" in r["csq"]:
            v = nflag(r["chrom"], r["pos"], r["ref"], r["alt"])
            if v is not None:
                nf = v
                mt = {4: "miss_t1", 3: "miss_t2", 2: "miss_t3", 1: "miss_t4"}.get(v, "")
        tier = lt if lt else mt
        if lt == "lof_t1":
            n_t1_lof += 1
        if mt == "miss_t1":
            n_t1_miss += 1
        fh.write("\t".join([r["family"], r["proband"], r["chrom"], str(r["pos"]),
                            r["ref"][:30], r["alt"][:30], r["gene"], r["csq"],
                            "%.6g" % r["gnomad_af"], lt, sh, str(nf), mt, tier]) + "\n")

print("")
print("wrote %s  (%d rows)" % (out, len(rows)))
print("=" * 66)
print("TIER 1 de novo LoF     : %d" % n_t1_lof)
print("TIER 1 de novo missense: %d" % n_t1_miss)
print("=" * 66)
