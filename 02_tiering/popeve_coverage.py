"""popEVE coverage among missense variants, dbNSFP 5.3.1a.

JS: what fraction of missense are missing popEVE?
Determines whether a missing popEVE counts as a failed flag (capping the variant
at miss_t2, which is what the lab SQL does via COALESCE(...,FALSE)) or whether
n_flag should be scored over the metrics actually present.

No pandas, no line continuations.
"""
import duckdb
import glob
import os

D = "/expanse/projects/sebat1/s3/data/sebat/resources/dbNSFP/5.3.1a/parquet_expanded_mane_select"
WANT = ["popEVE_converted_rankscore", "ClinPred_rankscore",
        "AlphaMissense_rankscore", "MPC_rankscore"]
T_STARS = {"ClinPred_rankscore": 0.4298,
           "AlphaMissense_rankscore": 0.9603,
           "popEVE_converted_rankscore": 0.9209,
           "MPC_rankscore": 0.8947}

con = duckdb.connect()
con.execute("PRAGMA threads=4")

files = sorted(glob.glob(os.path.join(D, "chr*.parquet")))
print("parquet files: %d" % len(files))
cols = [r[0] for r in con.execute("DESCRIBE SELECT * FROM read_parquet('%s')" % files[0]).fetchall()]
present = [c for c in WANT if c in cols]
print("score columns present: %d of 4 -> %s" % (len(present), ", ".join(present)))

MP = [c for c in ("ClinPred_rankscore", "AlphaMissense_rankscore", "MPC_rankscore") if c in cols]
anyexpr = " OR ".join(['"%s" IS NOT NULL' % c for c in MP])
PE = '"popEVE_converted_rankscore"'

parts = ["count(*) AS n_rows"]
for c in WANT:
    parts.append('count(*) FILTER (WHERE "%s" IS NOT NULL) AS n_%s' % (c, c.split("_")[0]))
parts.append("count(*) FILTER (WHERE %s) AS any_miss" % anyexpr)
parts.append("count(*) FILTER (WHERE (%s) AND %s IS NULL) AS miss_no_pe" % (anyexpr, PE))
parts.append("count(*) FILTER (WHERE (%s) AND %s IS NOT NULL) AS miss_yes_pe" % (anyexpr, PE))
# how often would a missing popEVE actually change the tier?
# i.e. rows that clear the other three thresholds but have no popEVE
three = []
for c in ("ClinPred_rankscore", "AlphaMissense_rankscore", "MPC_rankscore"):
    three.append('TRY_CAST("%s" AS DOUBLE) >= %s' % (c, T_STARS[c]))
three_expr = " AND ".join(three)
parts.append("count(*) FILTER (WHERE %s) AS pass3" % three_expr)
parts.append("count(*) FILTER (WHERE (%s) AND %s IS NULL) AS pass3_no_pe" % (three_expr, PE))
parts.append("count(*) FILTER (WHERE (%s) AND TRY_CAST(%s AS DOUBLE) >= %s) AS pass4" % (three_expr, PE, T_STARS["popEVE_converted_rankscore"]))

Q = "SELECT " + ", ".join(parts) + " FROM read_parquet(?)"
KEYS = (["n_rows"] + ["n_" + c.split("_")[0] for c in WANT] +
        ["any_miss", "miss_no_pe", "miss_yes_pe", "pass3", "pass3_no_pe", "pass4"])


def pct(a, b):
    return (100.0 * a / b) if b else 0.0


def report(label, fl):
    row = con.execute(Q, [fl]).fetchone()
    r = dict(zip(KEYS, [int(x) if x is not None else 0 for x in row]))
    n = r["n_rows"]
    print("")
    print("=== %s ===" % label)
    print("  dbNSFP nsSNV rows                    : {:,}".format(n))
    for c in WANT:
        k = "n_" + c.split("_")[0]
        print("  with {:<32}: {:,} ({:.1f}%)".format(c.split("_")[0], r[k], pct(r[k], n)))
    am = r["any_miss"]
    print("  --")
    print("  rows with a missense predictor       : {:,}".format(am))
    print("  ... WITH popEVE                      : {:,} ({:.1f}%)".format(r["miss_yes_pe"], pct(r["miss_yes_pe"], am)))
    print("  ... MISSING popEVE                   : {:,} ({:.1f}%)".format(r["miss_no_pe"], pct(r["miss_no_pe"], am)))
    print("  --  does it actually change a tier?")
    print("  clear the other 3 thresholds         : {:,}".format(r["pass3"]))
    print("  ... and popEVE MISSING (capped t2)   : {:,} ({:.1f}% of those)".format(r["pass3_no_pe"], pct(r["pass3_no_pe"], r["pass3"])))
    print("  ... and popEVE also passes (= miss_t1): {:,}".format(r["pass4"]))


sub = [f for f in files if os.path.basename(f) in ("chr1.parquet", "chr21.parquet", "chrX.parquet")]
report("SUBSET chr1 + chr21 + chrX", sub)
report("ALL CHROMOSOMES", files)
