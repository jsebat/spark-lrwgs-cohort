"""Genome-wide missense tier distribution in dbNSFP 5.3.1a.

JS: what fraction of missense calls are tier 1?
This is the calibration for whether 0 miss_t1 among our 17 de novo missense is
meaningful or simply expected.
"""
import duckdb
import glob
import os

D = "/expanse/projects/sebat1/s3/data/sebat/resources/dbNSFP/5.3.1a/parquet_expanded_mane_select"
T = [("ClinPred_rankscore", 0.4298),
     ("AlphaMissense_rankscore", 0.9603),
     ("popEVE_converted_rankscore", 0.9209),
     ("MPC_rankscore", 0.8947)]

con = duckdb.connect()
con.execute("PRAGMA threads=4")
files = sorted(glob.glob(os.path.join(D, "chr*.parquet")))
print("parquet files: %d" % len(files))

flags = " + ".join(
    ['CAST(COALESCE(TRY_CAST("%s" AS DOUBLE) >= %s, FALSE) AS INTEGER)' % (c, t)
     for c, t in T])

q = """
WITH s AS (
  SELECT ({flags}) AS n_flag
  FROM read_parquet(?)
  WHERE "AlphaMissense_rankscore" IS NOT NULL
     OR "ClinPred_rankscore" IS NOT NULL
     OR "MPC_rankscore" IS NOT NULL
)
SELECT count(*) AS n,
       count(*) FILTER (WHERE n_flag = 4) AS t1,
       count(*) FILTER (WHERE n_flag = 3) AS t2,
       count(*) FILTER (WHERE n_flag = 2) AS t3,
       count(*) FILTER (WHERE n_flag = 1) AS t4,
       count(*) FILTER (WHERE n_flag = 0) AS none
FROM s
""".format(flags=flags)

n, t1, t2, t3, t4, none = con.execute(q, [files]).fetchone()
print("")
print("=" * 62)
print("  GENOME-WIDE MISSENSE TIER DISTRIBUTION (dbNSFP 5.3.1a)")
print("=" * 62)
print("  scored missense variants : {:,}".format(n))
for lab, v in (("miss_t1 (n_flag 4)", t1), ("miss_t2 (n_flag 3)", t2),
               ("miss_t3 (n_flag 2)", t3), ("miss_t4 (n_flag 1)", t4),
               ("no tier (n_flag 0)", none)):
    print("  {:<22}: {:>12,}  ({:.3f}%)".format(lab, v, 100.0 * v / n if n else 0))
print("=" * 62)
p1 = (1.0 * t1 / n) if n else 0
for k in (17, 36, 100, 1000):
    print("  expected miss_t1 among %4d missense calls : %.2f" % (k, k * p1))
print("=" * 62)
