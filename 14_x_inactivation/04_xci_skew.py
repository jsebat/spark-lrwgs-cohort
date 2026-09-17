#!/usr/bin/env python3
"""X-inactivation skew for every female, by transmission phasing where possible and folded otherwise.

Model. At an informative X-linked CpG island the active X sits at methylation level u and the inactive X at m.
If a fraction f of cells carry haplotype 1 as the inactive X then

    hap1 = f*m + (1-f)*u        hap2 = (1-f)*m + f*u
    hap1 - hap2 = (2f - 1)(m - u)        combined = (m + u)/2

so m is recovered per sample as 2*combined - u, with u taken from males, who have one active X and no inactive
one. The skew then follows from the haplotype difference.

Informative islands are found empirically rather than from a published escape list: an island is informative
when males sit low (a single active X) and females sit intermediate (one active, one inactive). Islands that
escape inactivation are low in both sexes and drop out on their own.

Two modes.

  trio    Uses the per-block orientation from 03_trio_phase_x.py to put every island on chromosome-wide
          paternal/maternal labels, then takes the signed mean. This is directional: it names which parental X
          is silenced.

  folded  For samples with no orientation available. Uses median |hap1 - hap2|, which is invariant to the
          per-block label flips, then removes the noise bias. Because |x| is a folded quantity its expectation
          is above zero even when the true difference is zero, so the observed value is deconvolved in
          quadrature against a floor measured on autosomal CpG islands: those are diploid, not subject to
          inactivation, and measured in the same reads, so their haplotype asymmetry is pure assay noise.
          Controls are matched on methylation level, because sampling noise peaks near 0.5 and an unmethylated
          control would understate the floor.

Both are reported where both are computable; agreement between them is the check that licenses using the
folded estimate on samples that cannot be trio-phased.
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
    p.add_argument("--pileup-dir", required=True,
                   help="per-sample CpG pileups named <sample>.{hap1,hap2,combined}.bed.gz")
    p.add_argument("--islands", required=True,
                   help="directory holding chrX_cgi.bed and autosomal_cgi.bed from 01_make_island_sets.sh")
    p.add_argument("--orientation", default=None,
                   help="output of 03_trio_phase_x.py; without it every sample is scored folded only")
    p.add_argument("--out", required=True, help="output TSV")
    p.add_argument("--min-cov", type=int, default=5, help="minimum CpG coverage (default 5)")
    p.add_argument("--min-cpg", type=int, default=5, help="minimum covered CpGs per island (default 5)")
    p.add_argument("--min-islands", type=int, default=10,
                   help="minimum usable islands before a sample is scored (default 10)")
    p.add_argument("--male-max", type=float, default=0.20,
                   help="island is informative if the male median is below this (default 0.20)")
    p.add_argument("--female-range", type=float, nargs=2, default=[0.25, 0.75],
                   help="and the female median lies in this range (default 0.25 0.75)")
    p.add_argument("--control-tol", type=float, default=0.12,
                   help="autosomal control islands must match the female X level within this (default 0.12)")
    return p.parse_args()


def read_ped(path):
    ped = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) >= 6:
                ped[f[1]] = dict(fam=f[0], pat=f[2], mat=f[3], sex=f[4], aff=(f[5] == "2"))
    return ped


def load_pileup(path, min_cov):
    """(chrom, start, beta) for CpGs at or above min_cov; pb-CpG-tools bed layout"""
    with open(path, "rb") as fh:
        gz = fh.read(2) == GZIP_MAGIC
    opener = gzip.open if gz else open
    rows = []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(("#", "track")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            try:
                if int(f[5]) >= min_cov:
                    rows.append((f[0], int(f[1]), float(f[3]) / 100.0))
            except ValueError:
                continue
    return rows


def island_means(rows, intervals, min_cpg):
    by_chrom = {}
    for chrom, pos, beta in rows:
        by_chrom.setdefault(chrom, []).append((pos, beta))
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    arrays = {c: (np.array([x[0] for x in v]), np.array([x[1] for x in v]))
              for c, v in by_chrom.items()}
    out = {}
    for chrom, start, end in intervals:
        if chrom not in arrays:
            continue
        pos, beta = arrays[chrom]
        i, j = np.searchsorted(pos, start), np.searchsorted(pos, end)
        if j - i >= min_cpg:
            out[(chrom, start, end)] = float(beta[i:j].mean())
    return out


def read_bed(path):
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(("#", "track")):
                continue
            f = line.split()
            if len(f) >= 3:
                out.append((f[0], int(f[1]), int(f[2])))
    return sorted(out)


def main():
    a = parse_args()
    ped = read_ped(a.ped)
    x_iv = read_bed(os.path.join(a.islands, "chrX_cgi.bed"))
    auto_path = os.path.join(a.islands, "autosomal_cgi.bed")
    a_iv = read_bed(auto_path) if os.path.exists(auto_path) else []

    orient = {}
    blocks = {}
    if a.orientation:
        with open(a.orientation) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                sample = r["sample"]
                ps = int(r["phase_set"])
                orient.setdefault(sample, {})[ps] = 1 if r["paternal_hap"] == "1" else -1
                if r.get("block_start") not in (None, "",):
                    blocks.setdefault(sample, []).append(
                        (int(r["block_start"]), int(r["block_end"]), ps))
        for s in blocks:
            blocks[s].sort()
        sys.stderr.write("orientation loaded for %d samples\n" % len(orient))
        if orient and not blocks:
            sys.stderr.write("WARNING: the orientation table has no block_start/block_end (03_trio_phase_x.py was run without "
                             "--blocks-dir), so NO female can use the trio-directional estimator and every one falls to the "
                             "folded estimator (X4). Re-run 03 with --blocks-dir for the directional result.\n")

    XH1, XH2, XCB, AH1, AH2, ACB = {}, {}, {}, {}, {}, {}
    for path in sorted(glob.glob(os.path.join(a.pileup_dir, "*.hap1.bed.gz"))):
        sample = os.path.basename(path).split(".")[0]
        if sample not in ped:
            continue
        h2 = path.replace(".hap1.", ".hap2.")
        cb = path.replace(".hap1.", ".combined.")
        if not (os.path.exists(h2) and os.path.exists(cb)):
            continue
        r1, r2, rc = (load_pileup(p, a.min_cov) for p in (path, h2, cb))
        XH1[sample] = island_means(r1, x_iv, a.min_cpg)
        XH2[sample] = island_means(r2, x_iv, a.min_cpg)
        XCB[sample] = island_means(rc, x_iv, a.min_cpg)
        if a_iv:
            AH1[sample] = island_means(r1, a_iv, a.min_cpg)
            AH2[sample] = island_means(r2, a_iv, a.min_cpg)
            ACB[sample] = island_means(rc, a_iv, a.min_cpg)
    males = [s for s in XCB if ped[s]["sex"] == "1"]
    females = [s for s in XCB if ped[s]["sex"] == "2"]
    if not males:
        sys.exit("no males in the pileup set: the unmethylated reference level cannot be calibrated")
    sys.stderr.write("loaded %d samples (%d male, %d female)\n" % (len(XCB), len(males), len(females)))

    lo, hi = a.female_range
    informative = []
    for iv in x_iv:
        mv = [XCB[s][iv] for s in males if iv in XCB[s]]
        fv = [XCB[s][iv] for s in females if iv in XCB[s]]
        if len(mv) >= 15 and len(fv) >= 15 and np.median(mv) < a.male_max and lo <= np.median(fv) <= hi:
            informative.append(iv)
    if not informative:
        sys.exit("no informative X islands: check sex assignment and pileup coverage")
    u = float(np.median([np.median([XCB[s][iv] for iv in informative if iv in XCB[s]]) for s in males]))
    sys.stderr.write("informative X islands: %d; male unmethylated level u = %.3f\n" % (len(informative), u))

    controls = []
    if a_iv:
        fx = float(np.median([np.median([XCB[s][iv] for iv in informative if iv in XCB[s]])
                              for s in females]))
        for iv in a_iv:
            v = [ACB[s][iv] for s in ACB if iv in ACB[s]]
            if len(v) >= 30 and abs(float(np.median(v)) - fx) <= a.control_tol:
                controls.append(iv)
        sys.stderr.write("female X level %.3f; matched autosomal control islands: %d\n" % (fx, len(controls)))

    def block_of(sample, chrom, start, end):
        bl = blocks.get(sample)
        if not bl:
            return None
        starts = [b[0] for b in bl]
        i = int(np.searchsorted(starts, start, side="right")) - 1
        if 0 <= i < len(bl) and start >= bl[i][0] and end <= bl[i][1]:
            return bl[i][2]
        return None

    rows = []
    for s in females:
        signed, n_signed = None, 0
        if s in orient and s in blocks:
            d = []
            for iv in informative:
                if iv not in XH1[s] or iv not in XH2[s]:
                    continue
                ps = block_of(s, *iv)
                if ps is None or ps not in orient[s]:
                    continue
                d.append(orient[s][ps] * (XH1[s][iv] - XH2[s][iv]))   # paternal minus maternal
            if len(d) >= a.min_islands:
                d = np.array(d)
                signed = float(d.mean())
                n_signed = len(d)
                se = float(d.std(ddof=1) / np.sqrt(len(d)))

        dx = [abs(XH1[s][iv] - XH2[s][iv]) for iv in informative
              if iv in XH1[s] and iv in XH2[s]]
        da = [abs(AH1[s][iv] - AH2[s][iv]) for iv in controls
              if s in AH1 and iv in AH1[s] and iv in AH2[s]]
        if len(dx) < a.min_islands:
            continue
        combined = float(np.median([XCB[s][iv] for iv in informative if iv in XCB[s]]))
        span = (2 * combined - u) - u                 # m - u, the dynamic range available in this sample
        if span <= 0.05:
            sys.stderr.write("  %s: methylation span too small (%.3f), not scored\n" % (s, span))
            continue

        rec = dict(sample=s,
                   role=("child" if ped[s]["pat"] != "0" else "parent"),
                   affected="yes" if ped[s]["aff"] else "",
                   n_islands=len(dx))
        if signed is not None:
            f_pat = min(max(0.5 + signed / (2 * span), 0.0), 1.0)
            rec.update(method="trio-phased",
                       skew=round(max(f_pat, 1 - f_pat), 4),
                       silenced_X=("paternal" if signed > 0 else "maternal"),
                       delta_pat_minus_mat=round(signed, 4),
                       z=round(signed / se, 2) if se > 0 else "",
                       n_oriented_islands=n_signed)
        else:
            floor = float(np.median(da)) if len(da) >= a.min_islands else 0.0
            observed = float(np.median(dx))
            true = float(np.sqrt(max(observed ** 2 - floor ** 2, 0.0)))
            f = min(0.5 + true / (2 * span), 1.0)
            rec.update(method="folded",
                       skew=round(f, 4),
                       silenced_X="unknown",
                       observed_absdiff=round(observed, 4),
                       noise_floor=round(floor, 4) if da else "",
                       n_control_islands=len(da))
        rec["deviation"] = round(rec["skew"] - 0.5, 4)
        rec["band"] = ("extreme" if rec["skew"] >= 0.90 else
                       "marked" if rec["skew"] >= 0.80 else
                       "skewed" if rec["skew"] >= 0.70 else "")
        rows.append(rec)

    if not rows:
        sys.exit("no female was scorable")
    rows.sort(key=lambda r: -r["skew"])
    fields = []
    for r in rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(rows)

    for role in ("parent", "child"):
        grp = [r for r in rows if r["role"] == role]
        if grp:
            n = sum(1 for r in grp if r["skew"] >= 0.70)
            sys.stderr.write("%-7s n=%d, skewed (>=0.70): %d (%.0f%%)\n"
                             % (role, len(grp), n, 100.0 * n / len(grp)))
    sys.stderr.write("wrote %s (%d females)\n" % (a.out, len(rows)))


if __name__ == "__main__":
    main()
