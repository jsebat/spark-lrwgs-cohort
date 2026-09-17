"""VCF-level minitrio: a simulated trio with known haplotypes, crossovers, phase blocks, switch errors and
planted de novo sites, written in the SAME shape as the WDL outputs (one phased VCF per sample split from
a joint call, PS = block start, plus the joint unphased VCF and a manifest). Truth tables let the tests
check orientation (M1a), transmission and crossovers (M1b) and the VCF-level candidate generator (M2).

Pure Python and deterministic. Read-level fixtures (BAMs) are produced separately, in the container.
Identifiers are synthetic and chosen NOT to match the repository's PHI denylist patterns.
"""
from __future__ import annotations

import argparse
import gzip
import os
import random
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

FAMILY = "minifam"
CHILD, FATHER, MOTHER = "mini_child", "mini_father", "mini_mother"
CONTIGS = [("chr21", 2_000_000), ("chr22", 1_500_000)]


@dataclass
class Sim:
    seed: int = 7
    site_spacing: int = 800          # mean bp between polymorphic sites
    block_mean_bp: int = 200_000     # read-based phase block length (exponential)
    block_min_bp: int = 20_000
    p_switch_per_block: float = 0.10  # within-block switch error (deliberately high so tests see some)
    p_unphased_het: float = 0.03
    p_gt_error: float = 0.003
    p_low_gq: float = 0.02
    n_dnm_per_chrom: int = 3
    child_sex: str = "2"             # female by default: chrX-free fixture anyway (autosomes only)


@dataclass
class Truth:
    crossovers: List[Tuple[str, str, int, int, int]] = field(default_factory=list)   # chrom, parent, pos, from_hap, to_hap
    dnms: List[Tuple[str, int, str, str, str]] = field(default_factory=list)          # chrom, pos, ref, alt, hap (pat|mat)
    child_blocks: List[Tuple[str, int, int, int, str, int]] = field(default_factory=list)  # chrom, ps, start, end, HAP1_PAT|HAP1_MAT, switch_pos(0=none)
    parent_blocks: List[Tuple[str, str, int, int, int, int, int]] = field(default_factory=list)  # sample, chrom, ps, start, end, hap1_is_true_hap(1|2), switch_pos


BASES = "ACGT"


def _alt(ref: str, rng: random.Random) -> str:
    return rng.choice([b for b in BASES if b != ref])


def simulate(sim: Sim, out_dir: str) -> Truth:
    rng = random.Random(sim.seed)
    truth = Truth()
    os.makedirs(out_dir, exist_ok=True)
    # per sample: list of records (chrom, pos, ref, alt, gt_string, gq, dp, ad, ps)
    recs: Dict[str, List[tuple]] = {CHILD: [], FATHER: [], MOTHER: []}
    joint: List[tuple] = []

    for chrom, length in CONTIGS:
        # ---- polymorphic sites and founder haplotypes F1 F2 M1 M2
        pos = 1000
        sites = []
        while pos < length - 1000:
            pos += max(50, int(rng.expovariate(1.0 / sim.site_spacing)))
            if pos >= length - 1000:
                break
            af = rng.uniform(0.05, 0.95)
            ref = rng.choice(BASES)
            haps = tuple(1 if rng.random() < af else 0 for _ in range(4))   # F1 F2 M1 M2
            if sum(haps) in (0, 4):
                continue                                                   # monomorphic in the family
            sites.append((pos, ref, _alt(ref, rng), haps))
        # ---- one crossover per parent per chromosome, in the middle 60 %
        xo = {}
        for parent in ("F", "M"):
            xpos = rng.randint(int(0.2 * length), int(0.8 * length))
            first = rng.choice((1, 2))
            xo[parent] = (xpos, first, 3 - first)
            truth.crossovers.append((chrom, parent, xpos, first, 3 - first))

        def transmitted(parent: str, p: int, haps: tuple) -> int:
            xpos, first, second = xo[parent]
            h = first if p < xpos else second
            base = 0 if parent == "F" else 2
            return haps[base + h - 1]

        # ---- child true haplotypes (pat, mat) per site; then planted DNMs
        child_true: List[Tuple[int, str, str, int, int]] = []   # pos, ref, alt, pat_allele, mat_allele
        for p, ref, alt, haps in sites:
            child_true.append((p, ref, alt, transmitted("F", p, haps), transmitted("M", p, haps)))
        dnm_positions = set()
        for _ in range(sim.n_dnm_per_chrom):
            while True:
                p = rng.randint(2000, length - 2000)
                if all(abs(p - s[0]) > 5 for s in sites) and p not in dnm_positions:
                    break
            dnm_positions.add(p)
            ref = rng.choice(BASES)
            hap = rng.choice(("pat", "mat"))
            truth.dnms.append((chrom, p, ref, _alt(ref, rng), hap))
        # merge DNMs into the site list as child-only alt, parents hom-ref
        all_sites = [(p, ref, alt, haps, None) for p, ref, alt, haps in sites]
        for (c, p, ref, alt, hap) in truth.dnms:
            if c == chrom:
                all_sites.append((p, ref, alt, (0, 0, 0, 0), hap))
        all_sites.sort(key=lambda s: s[0])

        # ---- genotypes per sample at every site: (allele_hapA, allele_hapB) in TRUE haplotype order
        true_gt: Dict[str, List[Tuple[int, int]]] = {CHILD: [], FATHER: [], MOTHER: []}
        for p, ref, alt, haps, dnm_hap in all_sites:
            true_gt[FATHER].append((haps[0], haps[1]))
            true_gt[MOTHER].append((haps[2], haps[3]))
            if dnm_hap is None:
                true_gt[CHILD].append((transmitted("F", p, haps), transmitted("M", p, haps)))
            else:
                true_gt[CHILD].append((1, 0) if dnm_hap == "pat" else (0, 1))

        # ---- simulate read-based phasing per sample: blocks, orientation flips, switch errors
        for sample in (CHILD, FATHER, MOTHER):
            gts = true_gt[sample]
            het_idx = [i for i, g in enumerate(gts) if g[0] != g[1]]
            # block boundaries along het sites
            blocks: List[List[int]] = []
            i = 0
            while i < len(het_idx):
                blen = max(sim.block_min_bp, int(rng.expovariate(1.0 / sim.block_mean_bp)))
                start_pos = all_sites[het_idx[i]][0]
                blk = []
                while i < len(het_idx) and all_sites[het_idx[i]][0] < start_pos + blen:
                    blk.append(het_idx[i])
                    i += 1
                blocks.append(blk)
            phased_gt: Dict[int, Tuple[str, int]] = {}   # site index -> (gt string, ps)
            for blk in blocks:
                if len(blk) < 2:
                    continue                                   # singleton: left unphased
                flip = rng.random() < 0.5                      # hap1 = true hapB when flipped
                ps = all_sites[blk[0]][0]
                switch_at = 0
                if rng.random() < sim.p_switch_per_block and len(blk) >= 6:
                    k = rng.randint(2, len(blk) - 2)
                    switch_at = all_sites[blk[k]][0]
                cur_flip = flip
                for idx in blk:
                    p = all_sites[idx][0]
                    if switch_at and p >= switch_at:
                        cur_flip = not flip
                    a, b = gts[idx]
                    h1, h2 = (b, a) if cur_flip else (a, b)
                    phased_gt[idx] = ("%d|%d" % (h1, h2), ps)
                start, end = all_sites[blk[0]][0], all_sites[blk[-1]][0]
                if sample == CHILD:
                    truth.child_blocks.append((chrom, ps, start, end, "HAP1_MAT" if flip else "HAP1_PAT", switch_at))
                else:
                    truth.parent_blocks.append((sample, chrom, ps, start, end, 2 if flip else 1, switch_at))
            # ---- emit records with noise
            for idx, (p, ref, alt, haps, dnm_hap) in enumerate(all_sites):
                a, b = gts[idx]
                dp = max(4, int(rng.gauss(22, 5)))
                gq = rng.randint(3, 19) if rng.random() < sim.p_low_gq else rng.randint(25, 60)
                if idx in phased_gt and rng.random() >= sim.p_unphased_het:
                    gt, ps = phased_gt[idx]
                else:
                    gt, ps = ("%d/%d" % (min(a, b), max(a, b))), None
                if rng.random() < sim.p_gt_error:
                    gt, ps = rng.choice(("0/0", "0/1", "1/1")), None
                n_alt = gt.count("1")
                ad = (dp - (dp * n_alt) // 2, (dp * n_alt) // 2)
                recs[sample].append((chrom, p, ref, alt, gt, gq, dp, ad, ps))
        for idx, (p, ref, alt, haps, dnm_hap) in enumerate(all_sites):
            joint.append((chrom, p, ref, alt))

    # ---- write per-sample phased VCFs (WDL shape) and the joint unphased VCF + manifest + truth
    header = ["##fileformat=VCFv4.2", "##source=make_minitrio_vcf.py (simulated)"]
    header += ["##contig=<ID=%s,length=%d>" % (c, l) for c, l in CONTIGS]
    header += ['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
               '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
               '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
               '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
               '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">']
    paths = {}
    for sample in (CHILD, FATHER, MOTHER):
        path = os.path.join(out_dir, "%s.%s.joint.GRCh38.small_variants.phased.vcf.gz" % (sample, FAMILY))
        paths[sample] = path
        with gzip.open(path, "wt") as fh:
            fh.write("\n".join(header) + "\n")
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % sample)
            for chrom, p, ref, alt, gt, gq, dp, ad, ps in recs[sample]:
                fmt, val = "GT:GQ:DP:AD", "%s:%d:%d:%d,%d" % (gt, gq, dp, ad[0], ad[1])
                if ps is not None:
                    fmt, val = fmt + ":PS", val + ":%d" % ps
                fh.write("\t".join((chrom, str(p), ".", ref, alt, "50", "PASS", ".", fmt, val)) + "\n")
    # joint unphased VCF (GLnexus-like: GT:GQ:DP:AD for all three) — M2 candidate source
    jpath = os.path.join(out_dir, "%s.joint.GRCh38.small_variants.vcf.gz" % FAMILY)
    by_key = {s: {(r[0], r[1]): r for r in recs[s]} for s in recs}
    with gzip.open(jpath, "wt") as fh:
        fh.write("\n".join(header) + "\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\t%s\t%s\n" % (CHILD, FATHER, MOTHER))
        for chrom, p, ref, alt in joint:
            cols = []
            for s in (CHILD, FATHER, MOTHER):
                r = by_key[s][(chrom, p)]
                cols.append("%s:%d:%d:%d,%d" % (r[4].replace("|", "/"), r[5], r[6], r[7][0], r[7][1]))
            fh.write("\t".join((chrom, str(p), ".", ref, alt, "50", "PASS", ".", "GT:GQ:DP:AD", *cols)) + "\n")
    paths["joint"] = jpath
    with open(os.path.join(out_dir, "manifest.tsv"), "w") as fh:
        fh.write("sample_id\tfamily_id\tfather_id\tmother_id\tsex\taffected\trole\tlr_haplotagged_bam\t"
                 "lr_family_sv_vcf\tlr_family_smallvar_vcf\tlr_family_trgt_vcf\n")
        fh.write("\t".join((CHILD, FAMILY, FATHER, MOTHER, sim.child_sex, "affected", "offspring", "", "", paths[CHILD], "")) + "\n")
        fh.write("\t".join((FATHER, FAMILY, "0", "0", "1", "unaffected", "parent", "", "", paths[FATHER], "")) + "\n")
        fh.write("\t".join((MOTHER, FAMILY, "0", "0", "2", "unaffected", "parent", "", "", paths[MOTHER], "")) + "\n")
    with open(os.path.join(out_dir, "truth_crossovers.tsv"), "w") as fh:
        fh.write("chrom\tparent\tpos\tfrom_hap\tto_hap\n")
        for r in truth.crossovers:
            fh.write("\t".join(map(str, r)) + "\n")
    with open(os.path.join(out_dir, "truth_dnm.tsv"), "w") as fh:
        fh.write("chrom\tpos\tref\talt\thap\n")
        for r in truth.dnms:
            fh.write("\t".join(map(str, r)) + "\n")
    with open(os.path.join(out_dir, "truth_child_blocks.tsv"), "w") as fh:
        fh.write("chrom\tphase_block_id\tstart\tend\ttrue_orientation\tswitch_pos\n")
        for r in truth.child_blocks:
            fh.write("\t".join(map(str, r)) + "\n")
    with open(os.path.join(out_dir, "truth_parent_blocks.tsv"), "w") as fh:
        fh.write("sample\tchrom\tphase_block_id\tstart\tend\thap1_is_true_hap\tswitch_pos\n")
        for r in truth.parent_blocks:
            fh.write("\t".join(map(str, r)) + "\n")
    return truth


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("out_dir")
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    t = simulate(Sim(seed=a.seed), a.out_dir)
    print("minitrio written to %s: %d child blocks, %d crossovers, %d planted DNMs"
          % (a.out_dir, len(t.child_blocks), len(t.crossovers), len(t.dnms)))


if __name__ == "__main__":
    main()
