#!/usr/bin/env python3
"""Assign an unphased structural variant to a parental haplotype, using breakpoint-spanning reads.

The problem. A hemizygous deletion removes one haplotype's sequence, so no heterozygous site survives inside
it and read-based phasers leave the SV unphased. Haplotype-tagged coverage inside the interval is also
uninformative for the same reason: with nothing to phase against, reads there carry no tag. So the deleted
haplotype is not directly observable, and without it a de novo deletion cannot be connected to the direction
of X-inactivation skew.

The solution. A read that crosses a breakpoint carries the deletion junction and also extends into flanking
sequence that is phased normally, where it picks up a haplotype tag. Those reads identify the haplotype
carrying the SV. Junction reads are found by their supplementary alignment: the primary alignment stops at one
breakpoint and the supplementary resumes at the other.

Each junction read is oriented in its own phase block, since a read near one breakpoint and a read near the
other may well sit in different blocks whose haplotype labels are unrelated. The script requires the reads to
agree, and reports NOT RESOLVED rather than guessing when they do not or when only untagged reads are found.

Works for deletions and duplications with a detectable junction. Balanced events and mobile-element insertions
need a different approach.
"""
import argparse
import csv
import glob
import os
import re
import subprocess
import sys
from collections import defaultdict

PAR_GRCH38 = [(10001, 2781479), (155701383, 156030895)]
PAR_GRCH37 = [(60001, 2699520), (154931044, 155260560)]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ped", required=True, help="PLINK-format pedigree")
    p.add_argument("--sample", required=True, help="the carrier")
    p.add_argument("--bam", required=True, help="haplotagged BAM for the carrier (HP and PS tags required)")
    p.add_argument("--vcf-dir", required=True, help="per-sample phased VCFs named <sample>.*.vcf.gz")
    p.add_argument("--region", required=True, help="the SV, as chrom:start-end")
    p.add_argument("--out", default=None, help="optional output TSV")
    p.add_argument("--build", choices=["GRCh38", "GRCh37"], default="GRCh38")
    p.add_argument("--flank", type=int, default=8000,
                   help="how far either side of a breakpoint to search for junction reads (default 8000)")
    p.add_argument("--window", type=int, default=2000,
                   help="a supplementary alignment within this distance of the other breakpoint counts "
                        "as a junction (default 2000)")
    p.add_argument("--orient-span", type=int, default=200000,
                   help="how far to search for informative sites when orienting a phase block (default 200kb)")
    p.add_argument("--min-gq", type=int, default=20)
    p.add_argument("--samtools", default="samtools")
    p.add_argument("--bcftools", default="bcftools")
    return p.parse_args()


def read_ped(path):
    ped = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) >= 6:
                ped[f[1]] = dict(fam=f[0], pat=f[2], mat=f[3], sex=f[4])
    return ped


def vcf_for(vcf_dir, sample):
    hits = sorted(glob.glob(os.path.join(vcf_dir, sample + ".*.vcf.gz")))
    return hits[0] if hits else None


def query(bcftools, vcf, fmt, region):
    r = subprocess.run([bcftools, "query", "-r", region, "-f", fmt, vcf],
                       capture_output=True, text=True)
    return r.stdout.splitlines() if r.returncode == 0 else []


def orient_block(a, ped, chrom, phase_set, kid_vcf, dad_vcf, mom_vcf, par):
    """+1 if hap1 is paternal in this block, -1 if hap2 is, 0 if undetermined. Also returns the evidence."""
    lo = max(phase_set - 5000, 1)
    hi = phase_set + a.orient_span
    region = "%s:%d-%d" % (chrom, lo, hi)

    paternal = {}
    for ln in query(a.bcftools, dad_vcf, "%POS\t[%GT]\t[%GQ]\n", region):
        pos, gt, gq = ln.split("\t")
        if gq in (".", "") or int(gq) < a.min_gq:
            continue
        al = [x for x in gt.replace("|", "/").split("/") if x != "."]
        if al and len(set(al)) == 1:
            paternal[int(pos)] = al[0]
    maternal = {}
    for ln in query(a.bcftools, mom_vcf, "%POS\t[%GT]\t[%GQ]\n", region):
        pos, gt, gq = ln.split("\t")
        if gq in (".", "") or int(gq) < a.min_gq:
            continue
        al = {x for x in gt.replace("|", "/").split("/") if x != "."}
        if al:
            maternal[int(pos)] = al

    votes, detail = [], []
    for ln in query(a.bcftools, kid_vcf, "%POS\t%REF\t%ALT\t[%GT]\t[%PS]\t[%GQ]\n", region):
        pos, ref, alt, gt, ps, gq = ln.split("\t")
        pos = int(pos)
        if ps in (".", "") or int(ps) != phase_set or "|" not in gt:
            continue
        if gq in (".", "") or int(gq) < a.min_gq:
            continue
        if any(x <= pos <= y for x, y in par):
            continue
        first, second = gt.split("|")
        if first == second or pos not in paternal or pos not in maternal:
            continue
        pat = paternal[pos]
        if pat not in (first, second):
            detail.append("%s:%d %s>%s skipped, Mendelian-inconsistent" % (chrom, pos, ref, alt))
            continue
        other = second if first == pat else first
        if other not in maternal[pos]:
            detail.append("%s:%d %s>%s skipped, allele in neither parent" % (chrom, pos, ref, alt))
            continue
        v = 1 if first == pat else -1
        votes.append(v)
        detail.append("%s:%d %s>%s child %s, father hemizygous %s -> hap%d paternal"
                      % (chrom, pos, ref, alt, gt, pat, 1 if v == 1 else 2))
    total = sum(votes)
    return (1 if total > 0 else -1 if total < 0 else 0), votes, detail


def main():
    a = parse_args()
    par = PAR_GRCH38 if a.build == "GRCh38" else PAR_GRCH37
    ped = read_ped(a.ped)
    if a.sample not in ped:
        sys.exit("%s is not in the pedigree" % a.sample)
    dad, mom = ped[a.sample]["pat"], ped[a.sample]["mat"]
    if dad == "0" or mom == "0" or dad not in ped or mom not in ped:
        sys.exit("%s is not a trio child; parental haplotypes cannot be named" % a.sample)

    m = re.match(r"^([\w.]+):(\d+)-(\d+)$", a.region)
    if not m:
        sys.exit("--region must look like chrom:start-end")
    chrom, left, right = m.group(1), int(m.group(2)), int(m.group(3))

    kid_vcf, dad_vcf, mom_vcf = (vcf_for(a.vcf_dir, s) for s in (a.sample, dad, mom))
    if not (kid_vcf and dad_vcf and mom_vcf):
        sys.exit("missing a phased VCF for the trio")

    # Junction reads: primary alignment near one breakpoint, supplementary near the other.
    search = "%s:%d-%d" % (chrom, max(left - a.flank, 1), right + a.flank)
    sam = subprocess.run([a.samtools, "view", a.bam, search], capture_output=True, text=True)
    if sam.returncode != 0:
        sys.exit("samtools view failed: %s" % sam.stderr.strip()[:300])

    junction = {}
    untagged = 0
    for line in sam.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        sa = hp = ps = None
        for tag in f[11:]:
            if tag.startswith("SA:Z:"):
                sa = tag[5:]
            elif tag.startswith("HP:i:"):
                hp = int(tag[5:])
            elif tag.startswith("PS:i:"):
                ps = int(tag[5:])
        if not sa:
            continue
        linked = None
        for seg in sa.split(";"):
            if not seg:
                continue
            g = seg.split(",")
            if len(g) >= 2 and g[0] == chrom:
                pos = int(g[1])
                if abs(pos - left) <= a.window or abs(pos - right) <= a.window:
                    linked = pos
        if linked is None:
            continue
        if hp is None or ps is None:
            untagged += 1
            continue
        junction[f[0]] = dict(pos=int(f[3]), hp=hp, ps=ps, mate=linked)

    print("=== junction reads ===")
    if not junction:
        print("  none with a haplotype tag (%d junction reads were untagged)" % untagged)
    for name, r in junction.items():
        print("  %s  aligned %d, supplementary %d, HP=%d PS=%d"
              % (name, r["pos"], r["mate"], r["hp"], r["ps"]))

    print("\n=== orienting the phase blocks those reads fall in ===")
    calls, rows = [], []
    for name, r in junction.items():
        o, votes, detail = orient_block(a, ped, chrom, r["ps"], kid_vcf, dad_vcf, mom_vcf, par)
        print("  PS=%d (%d informative sites)" % (r["ps"], len(votes)))
        for d in detail[:6]:
            print("      " + d)
        if o == 0:
            print("    -> could not orient this block")
            continue
        paternal_hap = 1 if o == 1 else 2
        parent = "paternal" if r["hp"] == paternal_hap else "maternal"
        print("    -> hap%d is paternal; read is HP=%d, so the SV is on the %s haplotype"
              % (paternal_hap, r["hp"], parent.upper()))
        calls.append(parent)
        rows.append(dict(sample=a.sample, region=a.region, read=name, phase_set=r["ps"],
                         read_hap=r["hp"], paternal_hap=paternal_hap,
                         n_informative=len(votes), call=parent))

    print("\n=== verdict ===")
    if calls and len(set(calls)) == 1:
        print("  the SV lies on the %s haplotype (%d of %d tagged junction reads, unanimous)"
              % (calls[0].upper(), len(calls), len(junction)))
    elif calls:
        print("  NOT RESOLVED: junction reads disagree (%s)" % ", ".join(calls))
    else:
        print("  NOT RESOLVED: no junction read could be both tagged and oriented")

    if a.out and rows:
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        print("\nwrote %s" % a.out)


if __name__ == "__main__":
    main()
