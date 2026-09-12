#!/usr/bin/env python3
"""Which X-linked genes are subject to X-inactivation in this cohort, and which escape it.

Whether the direction of skew matters for a gene depends on that gene being subject to inactivation. A gene
that escapes is expressed from both X chromosomes, so silencing the X carrying a damaging allele buys nothing.

Rather than import an escape list, this derives the call from the same methylation the skew estimate uses. A
promoter island of a gene subject to inactivation is unmethylated in males, who have one X and it is active,
and near half-methylated in females, who have one active and one inactive copy. An escaping gene is
unmethylated in both sexes. The separation is large and does not need a model.

Run this on any new cohort as a positive control before trusting the skew direction: recovering established
escapees such as KDM6A and DDX3X, while calling a gene like MECP2 subject, shows the pileups and the sex
assignments are right.
"""
import argparse
import csv
import glob
import gzip
import os
import sys

import numpy as np

GZIP_MAGIC = bytes([0x1F, 0x8B])


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ped", required=True, help="PLINK-format pedigree")
    p.add_argument("--pileup-dir", required=True, help="per-sample <sample>.combined.bed.gz pileups")
    p.add_argument("--islands", required=True, help="chrX CpG island BED")
    p.add_argument("--genes", required=True,
                   help="BED of X-linked genes: chrom, start, end, symbol (strand optional in column 6)")
    p.add_argument("--out", required=True, help="output TSV")
    p.add_argument("--min-cov", type=int, default=5)
    p.add_argument("--min-cpg", type=int, default=5)
    p.add_argument("--min-samples", type=int, default=10,
                   help="minimum males and females covering an island (default 10)")
    p.add_argument("--promoter-window", type=int, default=3000,
                   help="island must start within this distance of the TSS, or lie inside the gene body")
    p.add_argument("--male-max", type=float, default=0.20,
                   help="males below this counts as unmethylated (default 0.20)")
    p.add_argument("--female-min", type=float, default=0.25,
                   help="females at or above this counts as subject to inactivation (default 0.25)")
    return p.parse_args()


def read_ped(path):
    ped = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) >= 6:
                ped[f[1]] = f[4]
    return ped


def load_islands(path, intervals, min_cov, min_cpg):
    with open(path, "rb") as fh:
        gz = fh.read(2) == GZIP_MAGIC
    rows = []
    with (gzip.open if gz else open)(path, "rt") as fh:
        for line in fh:
            if line.startswith(("#", "track")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            try:
                if int(f[5]) >= min_cov:
                    rows.append((int(f[1]), float(f[3]) / 100.0))
            except ValueError:
                continue
    rows.sort()
    pos = np.array([r[0] for r in rows])
    beta = np.array([r[1] for r in rows])
    out = {}
    for start, end in intervals:
        i, j = np.searchsorted(pos, start), np.searchsorted(pos, end)
        if j - i >= min_cpg:
            out[(start, end)] = float(beta[i:j].mean())
    return out


def main():
    a = parse_args()
    ped = read_ped(a.ped)

    intervals = []
    with open(a.islands) as fh:
        for line in fh:
            if line.startswith(("#", "track")):
                continue
            f = line.split()
            if len(f) >= 3 and f[0] == "chrX":
                intervals.append((int(f[1]), int(f[2])))
    intervals.sort()

    genes = []
    with open(a.genes) as fh:
        for line in fh:
            if line.startswith(("#", "track")):
                continue
            f = line.split()
            if len(f) >= 4 and f[0] == "chrX":
                strand = f[5] if len(f) >= 6 and f[5] in ("+", "-") else "+"
                genes.append((f[3], int(f[1]), int(f[2]), strand))

    males, females = [], []
    for path in sorted(glob.glob(os.path.join(a.pileup_dir, "*.combined.bed.gz"))):
        sample = os.path.basename(path).split(".")[0]
        sex = ped.get(sample)
        if sex not in ("1", "2"):
            continue
        d = load_islands(path, intervals, a.min_cov, a.min_cpg)
        (males if sex == "1" else females).append(d)
    if not males or not females:
        sys.exit("need both sexes in the pileup set")
    sys.stderr.write("males %d, females %d, islands %d\n" % (len(males), len(females), len(intervals)))

    rows = []
    for symbol, start, end, strand in genes:
        tss = start if strand == "+" else end
        near = [iv for iv in intervals
                if abs(iv[0] - tss) <= a.promoter_window or (iv[0] >= start and iv[1] <= end)]
        best = None
        for iv in near:
            mv = [d[iv] for d in males if iv in d]
            fv = [d[iv] for d in females if iv in d]
            if len(mv) < a.min_samples or len(fv) < a.min_samples:
                continue
            mm, fm = float(np.median(mv)), float(np.median(fv))
            # Prefer the island with the largest female-minus-male separation: for a gene subject to
            # inactivation that is the promoter island carrying the inactive-X methylation.
            if best is None or (fm - mm) > (best[2] - best[1]):
                best = (iv, mm, fm, len(mv), len(fv))
        if best is None:
            rows.append(dict(gene=symbol, island="", male_median="", female_median="",
                             call="no island with coverage", n_male="", n_female=""))
            continue
        iv, mm, fm, nm, nf = best
        if mm < a.male_max and fm >= a.female_min:
            call = "subject to XCI"
        elif mm < a.male_max and fm < a.female_min:
            call = "escapes XCI"
        else:
            call = "uninformative"
        rows.append(dict(gene=symbol, island="chrX:%d-%d" % iv,
                         male_median=round(mm, 3), female_median=round(fm, 3),
                         call=call, n_male=nm, n_female=nf))

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    from collections import Counter
    for call, n in Counter(r["call"] for r in rows).most_common():
        sys.stderr.write("  %-24s %d\n" % (call, n))
    sys.stderr.write("wrote %s\n" % a.out)


if __name__ == "__main__":
    main()
