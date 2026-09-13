"""M1b2 — classify transmission change points as CROSSOVER or SWITCH_ERROR from the parent's reads (DESIGN P3).

A change of the transmitted haplotype inside a parent's phase block is a crossover if the parent's phasing is
intact across the change point and only the child's inheritance changes; it is a phase-switch error if the
parent's tagged phase is wrong across some gap between consecutive heterozygous sites. The test is ALLELE
CONCORDANCE, not read count: for each consecutive pair of the parent's phased SNV hets inside the candidate
interval, take the parent's primary reads (MAPQ >= min_mapq, HP and PS present, PS = block) whose reference span
covers both sites and read their base at each. A read is CONCORDANT when both bases support the same tagged
haplotype (the phase across the gap is what the read says it is) and DISCORDANT when they support different
haplotypes. Intact phasing gives discordance near 0 (sequencing error only); a switch error between the two sites
makes every correct read discordant relative to the tags. Verdict over the interval: any gap with >= min_spanning
informative reads and discordance >= switch_min_disc -> SWITCH_ERROR (located at that gap); all gaps with
>= min_spanning informative reads and discordance <= crossover_max_disc -> CROSSOVER; otherwise AMBIGUOUS (too
few informative reads at some gap, or mixed).

The first version of this step counted spanning tagged reads only; with 15 kb reads and ~11 reads per haplotype
almost every gap has >= 3, so it passed 137-179 "crossovers" per meiosis against an expected ~26-42 (2026-09-12).
Change points caused by a CHILD-side switch error (which flips the child's parental allele) cannot be seen in the
parent's reads at all - the parent's phasing is intact there; they are flagged by `child_boundary_bp` (distance to
the nearest boundary of the child's oriented segments) computed downstream, and by M2's dist_block_edge.

pysam is imported lazily so the rest of the package (and its tests) stays pysam-free.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

COLUMNS = ["chrom", "parent", "parent_phase_block_id", "left_pos", "right_pos", "resolution_bp", "n_left", "n_right",
           "left_hap", "right_hap", "status", "n_parent_hets_in_interval", "n_gaps", "weakest_gap_start",
           "weakest_gap_end", "weakest_gap_informative_reads", "weakest_gap_discordant", "max_gap_disc_frac",
           "min_spanning_reads_required", "child_switch_in_interval"]


@dataclass
class Resolved:
    row: dict
    status: str
    n_hets: int
    n_gaps: int
    weakest: Tuple[int, int, int]   # (gap_start, gap_end, n_spanning)


def parent_phased_snv_hets(vcf, sample: str, chrom: str, start: int, end: int, ps: int) -> List[Tuple[int, str, str]]:
    """(pos, hap1_base, hap2_base) of the parent's phased heterozygous SNVs with phase set `ps` in [start, end]."""
    out = []
    for rec in vcf.fetch(chrom, max(0, start - 1), end):
        if not (start <= rec.pos <= end):
            continue
        s = rec.samples[sample]
        gt = s.get("GT")
        if gt is None or len(gt) != 2 or gt[0] is None or gt[1] is None or gt[0] == gt[1] or not s.phased:
            continue
        if s.get("PS") != ps:
            continue
        alleles = rec.alleles
        a1, a2 = alleles[gt[0]], alleles[gt[1]]
        if len(a1) != 1 or len(a2) != 1:
            continue                                            # indel het: base lookup is not well defined
        out.append((rec.pos, a1, a2))
    return out


def _hap_of_base(base: Optional[str], h1: str, h2: str) -> Optional[int]:
    if base == h1:
        return 1
    if base == h2:
        return 2
    return None


def interval_gap_concordance(bam, chrom: str, hets: List[Tuple[int, str, str]], ps: int, min_mapq: int
                             ) -> List[Tuple[int, int]]:
    """For consecutive pairs of SNV hets, (concordant, discordant) counts over reads spanning both — ONE pass per read.

    Each read of the block (HP and PS present, PS = block, MAPQ >= min_mapq, primary) is walked once over its
    aligned pairs, collecting the base at every het position it covers (set lookup); a read then contributes to
    the gap between two consecutive hets it spans when both bases support a haplotype: concordant if the same
    haplotype, discordant otherwise. This replaces a per-gap re-fetch that walked each 15 kb read once per gap.
    """
    if len(hets) < 2:
        return []
    pos_index = {p - 1: i for i, (p, _, _) in enumerate(hets)}          # 0-based reference position -> het index
    n_gaps = len(hets) - 1
    conc = [0] * n_gaps
    disc = [0] * n_gaps
    first, last = hets[0][0], hets[-1][0]
    for read in bam.fetch(chrom, first - 1, last):
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate:
            continue
        if read.mapping_quality < min_mapq or not read.has_tag("HP") or not read.has_tag("PS") or read.get_tag("PS") != ps:
            continue
        seq = read.query_sequence
        if seq is None:
            continue
        support: Dict[int, int] = {}
        for qpos, rpos in read.get_aligned_pairs(matches_only=True):
            i = pos_index.get(rpos)
            if i is not None:
                h = _hap_of_base(seq[qpos], hets[i][1], hets[i][2])
                if h is not None:
                    support[i] = h
        for i in range(n_gaps):
            if i in support and i + 1 in support:
                if support[i] == support[i + 1]:
                    conc[i] += 1
                else:
                    disc[i] += 1
    return list(zip(conc, disc))


def classify(changepoints_tsv: str, out_tsv: str, parent_bam: Dict[str, Tuple[str, Optional[str]]],
             parent_vcf: Dict[str, Tuple[str, Optional[str], str]], min_spanning: int = 3, min_mapq: int = 20,
             crossover_max_disc: float = 0.2, switch_min_disc: float = 0.8, max_hets_per_interval: int = 400,
             flank_pad: int = 20000, child_orientation_tsv: Optional[str] = None) -> Dict[str, int]:
    """parent_bam: {'F': (bam, bai), 'M': (bam, bai)}; parent_vcf: {'F': (vcf, index, sample_name), ...}.

    child_orientation_tsv: the child's <child>.orientation.tsv; when given, each row records whether the interval
    contains a located CHILD switch position (`child_switch_in_interval`) - a child-side switch flips the child's
    parental allele and fakes a change in the parent's block that the parent's reads cannot reveal (2026-09-12:
    10.9% of CROSSOVER calls vs 4.4% of AMBIGUOUS ones overlapped a child switch).
    """
    import pysam
    child_switches: Dict[str, List[int]] = {}
    if child_orientation_tsv:
        with open(child_orientation_tsv, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if int(r.get("switch_pos", 0) or 0) > 0:
                    child_switches.setdefault(r["chrom"], []).append(int(r["switch_pos"]))
    bams = {k: (pysam.AlignmentFile(p, "rb", index_filename=i) if i else pysam.AlignmentFile(p, "rb")) for k, (p, i) in parent_bam.items()}
    vcfs = {k: (pysam.VariantFile(p, index_filename=i) if i else pysam.VariantFile(p), s) for k, (p, i, s) in parent_vcf.items()}
    counts: Dict[str, int] = {"CROSSOVER": 0, "SWITCH_ERROR": 0, "AMBIGUOUS": 0, "UNTESTABLE": 0}
    with open(changepoints_tsv, newline="") as fh, open(out_tsv, "w") as out:
        rd = csv.DictReader(fh, delimiter="\t")
        out.write("\t".join(COLUMNS) + "\n")
        for r in rd:
            parent, chrom, ps = r["parent"], r["chrom"], int(r["parent_phase_block_id"])
            left, right = int(r["left_pos"]), int(r["right_pos"])
            n_hets = n_gaps = 0
            weakest = (left, right, -1, -1)                    # (gap_start, gap_end, informative reads, discordant)
            max_disc = -1.0
            child_hit = "Y" if child_switches and any(left <= p <= right for p in child_switches.get(chrom, ())) else ("N" if child_switches else ".")
            if parent not in bams or parent not in vcfs:
                status = "UNTESTABLE"
            else:
                vcf, sample = vcfs[parent]
                # The change lies somewhere in [left, right]. The tested span must COVER that whole interval, so it
                # always runs from the nearest phased SNV het of the block at or before `left` to the nearest at or
                # after `right`; an endpoint that is an indel het (not testable by base lookup) is thereby bridged
                # rather than leaving an untested stretch that a parental switch could hide in (first cohort run
                # over-called paternal crossovers by ~10 per meiosis with the span starting at the first SNV inside).
                wider = parent_phased_snv_hets(vcf, sample, chrom, max(1, left - flank_pad), right + flank_pad, ps)
                lefts = [h for h in wider if h[0] <= left]
                rights = [h for h in wider if h[0] >= right]
                inside = [h for h in wider if left < h[0] < right]
                hets = ([lefts[-1]] if lefts else []) + inside + ([rights[0]] if rights else [])
                child_hit = "Y" if child_switches and any(left <= p <= right for p in child_switches.get(chrom, ())) else "N"
                if len(hets) > max_hets_per_interval:           # a huge interval: thin to evenly spaced sites + last
                    step = len(hets) // max_hets_per_interval + 1
                    hets = hets[::step] + [hets[-1]]
                n_hets, n_gaps = len(hets), max(0, len(hets) - 1)
                if n_gaps == 0:
                    status = "AMBIGUOUS"                        # fewer than two SNV hets to bridge: nothing to test
                else:
                    all_supported = True
                    status = None
                    weakest = (hets[0][0], hets[-1][0], 10 ** 9, 0)
                    per_gap = interval_gap_concordance(bams[parent], chrom, hets, ps, min_mapq)
                    for (a, b), (conc, disc) in zip(zip(hets[:-1], hets[1:]), per_gap):
                        n = conc + disc
                        frac = disc / n if n else 0.0
                        if n < weakest[2]:
                            weakest = (a[0], b[0], n, disc)
                        if n >= min_spanning:
                            max_disc = max(max_disc, frac)
                            if frac >= switch_min_disc:
                                status, weakest = "SWITCH_ERROR", (a[0], b[0], n, disc)
                                break
                            if frac > crossover_max_disc:
                                all_supported = False
                        else:
                            all_supported = False
                    if status is None:
                        status = "CROSSOVER" if all_supported else "AMBIGUOUS"
            counts[status] += 1
            out.write("\t".join(str(x) for x in (chrom, parent, ps, left, right, right - left, r["n_left"], r["n_right"],
                                                 r["left_hap"], r["right_hap"], status, n_hets, n_gaps, weakest[0], weakest[1],
                                                 weakest[2], weakest[3], round(max_disc, 3), min_spanning, child_hit)) + "\n")
    for b in bams.values():
        b.close()
    return counts
