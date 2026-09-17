#!/usr/bin/env python3
"""Orient read-based X phase blocks onto chromosome-wide parental haplotypes, using transmission.

Read-based phasers emit hundreds to thousands of independent phase blocks per chromosome. Their haplotype
labels are arbitrary per block, so any signed per-haplotype statistic cancels when blocks are combined. In a
trio this is recoverable on the X at no cost: the father is hemizygous across the non-pseudoautosomal region,
so at every heterozygous site in a daughter the paternal allele is whatever single allele the father carries.
Comparing that with the child's phased genotype says, for each block, whether haplotype 1 or haplotype 2 is the
paternal X.

A site is used only when it is Mendelian-consistent. The father's allele identifies the paternal X only if the
child's other allele actually came from the mother; where the child carries an allele neither parent has, the
variant is de novo or miscalled and could sit on either haplotype, so the site is skipped rather than guessed.

The per-block vote consistency is the diagnostic. Votes within a block should be unanimous or nearly so; a
consistency near 0.5 means the pedigree, the sex assignment or the phasing is wrong, not that the sample is
unusual.

Outputs one row per (sample, phase_set) with the orientation and the vote count behind it.
"""
import argparse
import csv
import glob
import os
import subprocess
import sys
from collections import defaultdict

# GRCh38 pseudoautosomal regions: the father is diploid here, so his allele is not informative.
PAR_GRCH38 = [(10001, 2781479), (155701383, 156030895)]
PAR_GRCH37 = [(60001, 2699520), (154931044, 155260560)]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ped", required=True,
                   help="PLINK-format pedigree: family, individual, father, mother, sex (1=M,2=F), affected")
    p.add_argument("--vcf-dir", required=True,
                   help="directory of per-sample phased VCFs, each indexed and named <sample>.*.vcf.gz")
    p.add_argument("--blocks-dir", default=None,
                   help="optional directory of phase-block TSVs named <sample>.*blocks.tsv "
                        "(columns include chrom, start, end, phase_block_id); only used to report block spans")
    p.add_argument("--out", required=True, help="output TSV")
    p.add_argument("--chrom", default="chrX", help="chromosome to phase (default chrX)")
    p.add_argument("--build", choices=["GRCh38", "GRCh37"], default="GRCh38")
    p.add_argument("--min-gq", type=int, default=20, help="minimum genotype quality (default 20)")
    p.add_argument("--bcftools", default="bcftools", help="path to bcftools")
    p.add_argument("--min-consistency", type=float, default=0.90,
                   help="warn if a sample's mean within-block vote consistency falls below this")
    return p.parse_args()


def read_ped(path):
    ped = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) >= 6:
                ped[f[1]] = dict(fam=f[0], pat=f[2], mat=f[3], sex=f[4], aff=f[5])
    return ped


def vcf_for(vcf_dir, sample):
    hits = sorted(glob.glob(os.path.join(vcf_dir, sample + ".*.vcf.gz")))
    return hits[0] if hits else None


def query(bcftools, vcf, fmt, region):
    r = subprocess.run([bcftools, "query", "-r", region, "-f", fmt, vcf],
                       capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write("bcftools query failed on %s: %s\n" % (vcf, r.stderr.strip()[:300]))
        return []
    return r.stdout.splitlines()


def hom_allele(gt, gq, min_gq):
    """the single allele of a hemizygous or homozygous call, else None"""
    if gq in (".", "") or int(gq) < min_gq:
        return None
    alleles = [a for a in gt.replace("|", "/").split("/") if a != "."]
    if alleles and len(set(alleles)) == 1:
        return alleles[0]
    return None


def carried_alleles(gt, gq, min_gq):
    if gq in (".", "") or int(gq) < min_gq:
        return None
    alleles = {a for a in gt.replace("|", "/").split("/") if a != "."}
    return alleles or None


def block_spans(blocks_dir, sample, chrom):
    if not blocks_dir:
        return {}
    hits = sorted(glob.glob(os.path.join(blocks_dir, sample + ".*blocks.tsv")))
    if not hits:
        return {}
    spans = {}
    with open(hits[0]) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            ci = header.index("chrom")
            si = header.index("start")
            ei = header.index("end")
            bi = header.index("phase_block_id")
        except ValueError:
            return {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > max(ci, si, ei, bi) and f[ci] == chrom:
                spans[int(f[bi])] = (int(f[si]), int(f[ei]))
    return spans


def main():
    a = parse_args()
    par = PAR_GRCH38 if a.build == "GRCh38" else PAR_GRCH37

    def in_par(pos):
        return any(lo <= pos <= hi for lo, hi in par)

    ped = read_ped(a.ped)
    # Daughters with both parents present. Sons are hemizygous and have no X to skew.
    kids = sorted(s for s, d in ped.items()
                  if d["sex"] == "2" and d["pat"] in ped and d["mat"] in ped)
    if not kids:
        sys.exit("no female trio children found in %s" % a.ped)
    sys.stderr.write("female trio children: %d\n" % len(kids))

    rows = []
    for kid in kids:
        kv = vcf_for(a.vcf_dir, kid)
        dv = vcf_for(a.vcf_dir, ped[kid]["pat"])
        mv = vcf_for(a.vcf_dir, ped[kid]["mat"])
        if not (kv and dv and mv):
            sys.stderr.write("  %s: missing a trio VCF, skipped\n" % kid)
            continue

        paternal = {}
        for ln in query(a.bcftools, dv, "%POS\t[%GT]\t[%GQ]\n", a.chrom):
            pos, gt, gq = ln.split("\t")
            al = hom_allele(gt, gq, a.min_gq)
            if al is not None:
                paternal[int(pos)] = al

        maternal = {}
        for ln in query(a.bcftools, mv, "%POS\t[%GT]\t[%GQ]\n", a.chrom):
            pos, gt, gq = ln.split("\t")
            al = carried_alleles(gt, gq, a.min_gq)
            if al:
                maternal[int(pos)] = al

        votes = defaultdict(list)
        skipped_mendel = 0
        for ln in query(a.bcftools, kv, "%POS\t[%GT]\t[%PS]\t[%GQ]\n", a.chrom):
            pos, gt, ps, gq = ln.split("\t")
            pos = int(pos)
            if ps in (".", "") or "|" not in gt:
                continue
            if gq in (".", "") or int(gq) < a.min_gq or in_par(pos):
                continue
            first, second = gt.split("|")
            if first == second:
                continue
            if pos not in paternal or pos not in maternal:
                continue
            pat = paternal[pos]
            if pat not in (first, second):
                skipped_mendel += 1
                continue
            other = second if first == pat else first
            if other not in maternal[pos]:
                skipped_mendel += 1
                continue
            votes[int(ps)].append(1 if first == pat else -1)

        spans = block_spans(a.blocks_dir, kid, a.chrom)
        cons = []
        for ps, v in sorted(votes.items()):
            total = sum(v)
            if total == 0:
                continue                      # tied block: no orientation rather than a coin flip
            orient = 1 if total > 0 else -1
            agree = max(v.count(1), v.count(-1))
            cons.append(agree / len(v))
            lo, hi = spans.get(ps, ("", ""))
            rows.append(dict(sample=kid, chrom=a.chrom, phase_set=ps,
                             paternal_hap=(1 if orient == 1 else 2),
                             n_informative=len(v), n_agree=agree,
                             consistency=round(agree / len(v), 3),
                             block_start=lo, block_end=hi))
        mean_cons = sum(cons) / len(cons) if cons else float("nan")
        flag = "" if (cons and mean_cons >= a.min_consistency) else "   <-- LOW, CHECK PEDIGREE/SEX"
        sys.stderr.write("  %s: %d blocks oriented, mean consistency %.3f, %d sites skipped as "
                         "Mendelian-inconsistent%s\n"
                         % (kid, len(cons), mean_cons, skipped_mendel, flag))

    if not rows:
        sys.exit("no blocks could be oriented")
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    sys.stderr.write("wrote %s (%d oriented blocks across %d samples)\n"
                     % (a.out, len(rows), len({r["sample"] for r in rows})))


if __name__ == "__main__":
    main()
