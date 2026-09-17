"""Quantify the ascertainment: what did WES already find in OUR probands?

JS: samples with a high-priority de novo SNV/indel were mostly NOT included in
the lrWGS set. If true, our probands should be depleted of high-tier de novo in
the SPARK WES SynthDNM callset relative to other SPARK probands.

Compares, in the SAME callset and with the SAME tiering:
   our 33 lrWGS probands   vs   all other SPARK WES probands
"""
import duckdb
import glob
import os
import collections

# paths from config/cohort.env (the repository README promises no cohort paths in code; these were one person's, X13)
S3 = os.environ.get("S3_DATA_ROOT", "/expanse/projects/sebat1/s3/data/sebat")
W = os.environ.get("WES_DNM_FILTERED", os.path.join(S3, "dnm_callsets_v1/outputs/spark_wes/filtered"))
DB = os.environ.get("DBNSFP_PARQUET", os.path.join(S3, "resources/dbNSFP/5.3.1a/parquet_expanded_mane_select"))
T = os.path.join(os.environ.get("DATA_ROOT", "/expanse/lustre/projects/ddp195/jsebat/longread-autism"), "tiering")
PED = os.environ.get("PED_SPARK_WES", os.path.join(S3, "SPARK_iWES_v3/SPARK_iWES_v3.ped"))
TAB, NL = chr(9), chr(10)
T_STARS = [("ClinPred_rankscore", 0.4298), ("AlphaMissense_rankscore", 0.9603),
           ("popEVE_converted_rankscore", 0.9209), ("MPC_rankscore", 0.8947)]

ours = set(x.strip() for x in open(os.path.join(T, "our_probands.txt")) if x.strip())
print("our lrWGS probands: %d" % len(ours))

# SPARK affected offspring = the comparison pool
spark_aff = set()
for ln in open(PED):
    f = ln.split()
    if len(f) >= 6 and f[2] != "0" and f[5] == "2":
        spark_aff.add(f[1])
print("SPARK WES affected offspring (comparison pool): %d" % len(spark_aff))

con = duckdb.connect()
con.execute("PRAGMA threads=4")
files = sorted(glob.glob(os.path.join(W, "chr*.parquet")))
print("de novo callset files: %d" % len(files))

flags = " + ".join(['CAST(COALESCE(TRY_CAST(d."%s" AS DOUBLE) >= %s, FALSE) AS INTEGER)'
                    % (c, t) for c, t in T_STARS])

q = """
SELECT v.SAMPLE AS s, ({flags}) AS n_flag
FROM read_parquet(?) v
JOIN read_parquet(?) d
  ON d."#chr" = replace(v."#CHROM", 'chr', '')
 AND d."pos(1-based)" = v.POS
 AND d.ref = v.REF AND d.alt = v.ALT
""".format(flags=flags)

per_sample = collections.defaultdict(lambda: collections.Counter())
tot = 0
for vf in files:
    ch = os.path.basename(vf).replace(".parquet", "")
    df = os.path.join(DB, ch + ".parquet")
    if not os.path.exists(df):
        continue
    try:
        for s, nf in con.execute(q, [vf, df]).fetchall():
            tot += 1
            tier = {4: "t1", 3: "t2", 2: "t3", 1: "t4"}.get(int(nf or 0), "none")
            per_sample[s][tier] += 1
    except Exception as e:
        print("  %s: %s" % (ch, str(e)[:110]))
print("de novo missense records tiered: %d across %d samples" % (tot, len(per_sample)))


def summarise(label, samples):
    n = len(samples)
    c = collections.Counter()
    carriers = collections.Counter()
    for s in samples:
        d = per_sample.get(s)
        if not d:
            continue
        for k, v in d.items():
            c[k] += v
        for k in ("t1", "t2"):
            if d.get(k):
                carriers[k] += 1
    print("")
    print("  %s  (n=%d probands)" % (label, n))
    for k in ("t1", "t2", "t3", "t4", "none"):
        print("    miss_%-5s total=%-7d  per-proband=%.4f" % (k, c[k], (1.0 * c[k] / n) if n else 0))
    for k in ("t1", "t2"):
        print("    probands carrying >=1 miss_%s : %d (%.1f%%)" % (k, carriers[k], 100.0 * carriers[k] / n if n else 0))


others = spark_aff - ours
summarise("OUR lrWGS PROBANDS", ours & spark_aff)
summarise("ALL OTHER SPARK WES PROBANDS", others)
