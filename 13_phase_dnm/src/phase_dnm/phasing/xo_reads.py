"""M1b2 — classify transmission change points as CROSSOVER or SWITCH_ERROR from the parent's reads (DESIGN P3).

A change of the transmitted haplotype inside a parent's phase block is a crossover if the parent's phasing is
intact across the change point and only the child's inheritance changes; it is a phase-switch error if HiPhase
had no read evidence bridging two consecutive heterozygous sites there. So, for each candidate (left_pos = last
informative site before the change, right_pos = first after), take the parent's PHASED heterozygous sites in the
block between the two, and for every consecutive pair count the parent's haplotagged reads (HP present, PS equal to
the block, MAPQ >= min_mapq, primary) whose reference span covers both sites. The WEAKEST LINK across the interval
decides: >= min_spanning reads at every gap -> CROSSOVER; a gap with 0 spanning tagged reads -> SWITCH_ERROR;
otherwise AMBIGUOUS. The weakest gap is reported, which is also the switch error's location when it is one.

pysam is imported lazily so the rest of the package (and its tests) stays pysam-free.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

COLUMNS = ["chrom", "parent", "parent_phase_block_id", "left_pos", "right_pos", "resolution_bp", "n_left", "n_right",
           "left_hap", "right_hap", "status", "n_parent_hets_in_interval", "n_gaps", "weakest_gap_start",
           "weakest_gap_end", "weakest_gap_spanning_reads", "min_spanning_reads_required"]


@dataclass
class Resolved:
    row: dict
    status: str
    n_hets: int
    n_gaps: int
    weakest: Tuple[int, int, int]   # (gap_start, gap_end, n_spanning)


def parent_phased_hets(vcf, sample: str, chrom: str, start: int, end: int, ps: int) -> List[int]:
    """Positions of the parent's phased heterozygous sites with phase set `ps` in [start, end]."""
    out = []
    for rec in vcf.fetch(chrom, max(0, start - 1), end):
        s = rec.samples[sample]
        gt = s.get("GT")
        if gt is None or len(gt) != 2 or gt[0] is None or gt[1] is None or gt[0] == gt[1] or not s.phased:
            continue
        if s.get("PS") != ps:
            continue
        if start <= rec.pos <= end:
            out.append(rec.pos)
    return out


def spanning_tagged_reads(bam, chrom: str, a: int, b: int, ps: int, min_mapq: int) -> int:
    n = 0
    for read in bam.fetch(chrom, a - 1, b):
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate:
            continue
        if read.mapping_quality < min_mapq or not read.has_tag("HP") or not read.has_tag("PS"):
            continue
        if read.get_tag("PS") != ps:
            continue
        if read.reference_start <= a - 1 and read.reference_end is not None and read.reference_end >= b:
            n += 1
    return n


def classify(changepoints_tsv: str, out_tsv: str, parent_bam: Dict[str, Tuple[str, Optional[str]]],
             parent_vcf: Dict[str, Tuple[str, Optional[str], str]], min_spanning: int = 3, min_mapq: int = 20,
             max_hets_per_interval: int = 400) -> Dict[str, int]:
    """parent_bam: {'F': (bam, bai), 'M': (bam, bai)}; parent_vcf: {'F': (vcf, index, sample_name), ...}."""
    import pysam
    bams = {k: (pysam.AlignmentFile(p, "rb", index_filename=i) if i else pysam.AlignmentFile(p, "rb")) for k, (p, i) in parent_bam.items()}
    vcfs = {k: (pysam.VariantFile(p, index_filename=i) if i else pysam.VariantFile(p), s) for k, (p, i, s) in parent_vcf.items()}
    counts: Dict[str, int] = {"CROSSOVER": 0, "SWITCH_ERROR": 0, "AMBIGUOUS": 0, "UNTESTABLE": 0}
    with open(changepoints_tsv, newline="") as fh, open(out_tsv, "w") as out:
        rd = csv.DictReader(fh, delimiter="\t")
        out.write("\t".join(COLUMNS) + "\n")
        for r in rd:
            parent, chrom, ps = r["parent"], r["chrom"], int(r["parent_phase_block_id"])
            left, right = int(r["left_pos"]), int(r["right_pos"])
            if parent not in bams or parent not in vcfs:
                status, n_hets, n_gaps, weakest = "UNTESTABLE", 0, 0, (left, right, -1)
            else:
                vcf, sample = vcfs[parent]
                hets = parent_phased_hets(vcf, sample, chrom, left, right, ps)
                if len(hets) > max_hets_per_interval:           # a huge interval: thin to the flanks + evenly spaced sites
                    step = len(hets) // max_hets_per_interval + 1
                    hets = hets[::step] + [hets[-1]]
                pts = sorted(set([left] + hets + [right]))
                gaps = list(zip(pts[:-1], pts[1:]))
                weakest = (left, right, 10 ** 9)
                for a, b in gaps:
                    n = spanning_tagged_reads(bams[parent], chrom, a, b, ps, min_mapq)
                    if n < weakest[2]:
                        weakest = (a, b, n)
                    if n == 0:
                        break
                n_hets, n_gaps = len(pts) - 2, len(gaps)
                if weakest[2] >= min_spanning:
                    status = "CROSSOVER"
                elif weakest[2] == 0:
                    status = "SWITCH_ERROR"
                else:
                    status = "AMBIGUOUS"
            counts[status] += 1
            out.write("\t".join(str(x) for x in (chrom, parent, ps, left, right, right - left, r["n_left"], r["n_right"],
                                                 r["left_hap"], r["right_hap"], status, n_hets, n_gaps, weakest[0], weakest[1],
                                                 weakest[2], min_spanning)) + "\n")
    for b in bams.values():
        b.close()
    return counts
