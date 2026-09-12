"""M1b — transmission map: for each parent, inside each of the PARENT's phase blocks, which haplotype the child
received, and where that changes (DESIGN P3).

Inputs are the merged trio sites (io.vcf.iter_trio) and the child's oriented segments (orient.py). At a site the
child's paternal allele is known when the child is homozygous (both alleles equal - no phasing needed) or when the
child is a phased heterozygote inside an oriented segment (the allele on the paternal haplotype). At a site where
the father is a phased heterozygote, the transmitted haplotype is the one carrying the child's paternal allele:
a vote (+1 = father's hap1 transmitted, -1 = hap2) for the father's block. Blocks are segmented with the same
change-point dynamic programme as orient.segment_votes.

A within-block change of transmitted haplotype is EITHER a crossover OR a phase-switch error in the parent's
read-based phasing, and at the VCF level the two are indistinguishable. Switch errors outnumber crossovers roughly
ten to one (300-470 located per genome vs ~30-45 crossovers), so change points are emitted as CANDIDATES and
resolved by parental read support across the change point in the read-level step (M1b2). Crossovers falling
between two parent blocks are invisible: the crossover count is a lower bound and the fraction of the genome
inside resolved parent segments bounds it.
"""
from __future__ import annotations

import json
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple

from ..io.vcf import TrioSite, normalise_sex
from .orient import (AMBIGUOUS, HAP1_MAT, HAP1_PAT, PAR_GRCH38, SKIP_ALWAYS, OrientParams, _in_par,
                     segment_votes)

HAP1, HAP2, UNRESOLVED = "HAP1", "HAP2", "UNRESOLVED"


@dataclass
class TransmissionParams:
    xo_min_sites: int = 10          # informative sites on EACH side of a change point (isolated-run size for the DP)
    min_sites: int = 10             # informative sites needed to call a segment's transmitted haplotype
    min_frac: float = 0.95
    min_gq: int = 20
    max_splits: int = 6             # change points per parent block (a long block can hold a crossover AND switches)


class Orientation:
    """Child segments from <child>.orientation.tsv, queried by (chrom, phase_block_id, pos)."""

    def __init__(self, rows: Iterable[dict]):
        self.segs: Dict[Tuple[str, int], List[Tuple[int, int, str]]] = defaultdict(list)
        for r in rows:
            self.segs[(r["chrom"], int(r["phase_block_id"]))].append((int(r["start"]), int(r["end"]), r["orientation"]))
        for v in self.segs.values():
            v.sort()

    @classmethod
    def from_tsv(cls, path: str) -> "Orientation":
        import csv
        with open(path, newline="") as fh:
            return cls(csv.DictReader(fh, delimiter="\t"))

    def hap1_is(self, chrom: str, ps: int, pos: int) -> Optional[str]:
        """'P' or 'M' for the child's hap1 at this position, None if unoriented."""
        for start, end, orient in self.segs.get((chrom, ps), ()):
            if start <= pos <= end:
                return "P" if orient == HAP1_PAT else ("M" if orient == HAP1_MAT else None)
        return None


@dataclass
class ParentBlockVotes:
    parent: str                     # 'F' | 'M'
    chrom: str
    ps: int
    start: int = 1 << 62
    end: int = 0
    het_positions: List[int] = field(default_factory=list)
    votes: List[Tuple[int, int]] = field(default_factory=list)   # (pos, +1 hap1 transmitted / -1 hap2)

    def add_het(self, pos: int):
        self.het_positions.append(pos)
        self.start = min(self.start, pos)
        self.end = max(self.end, pos)


@dataclass
class Segment:
    parent: str
    chrom: str
    ps: int
    start: int
    end: int
    segment: int
    n_segments: int
    n_het_phased: int
    n_hap1: int
    n_hap2: int
    transmitted: str                # HAP1 | HAP2 | UNRESOLVED
    vote_frac: float
    reason: str
    left_boundary: str              # BLOCK_EDGE | CHANGE_POINT
    right_boundary: str
    change_pos: int                 # first informative position of this segment when left_boundary is CHANGE_POINT


@dataclass
class ChangePoint:
    parent: str
    chrom: str
    ps: int
    left_pos: int                   # last informative position before the change
    right_pos: int                  # first informative position after the change
    n_left: int
    n_right: int
    left_hap: str
    right_hap: str
    status: str = "CANDIDATE"       # CROSSOVER | SWITCH_ERROR after the read-level step (M1b2)


def child_parental_alleles(t: TrioSite, orientation: Orientation) -> Tuple[Optional[int], Optional[int]]:
    """(paternal allele, maternal allele) of the child at this site, or (None, None) when unknown."""
    c = t.child
    if c.alleles is None or len(c.alleles) != 2:
        return None, None
    a, b = c.alleles
    if a == b:
        return a, a
    if not c.phased or c.ps is None:
        return None, None
    h1 = orientation.hap1_is(t.chrom, c.ps, t.pos)
    if h1 == "P":
        return a, b
    if h1 == "M":
        return b, a
    return None, None


def build_transmission(trio_sites: Iterable[TrioSite], orientation: Orientation, child_sex: str,
                       p: TransmissionParams) -> Tuple[List[Segment], List[ChangePoint], Dict[str, Dict[str, int]]]:
    sex = normalise_sex(child_sex)
    blocks: Dict[Tuple[str, str, int], ParentBlockVotes] = {}
    stats: Dict[str, Dict[str, int]] = {"F": defaultdict(int), "M": defaultdict(int)}
    for t in trio_sites:
        if t.chrom in SKIP_ALWAYS:
            continue
        pat_allele, mat_allele = child_parental_alleles(t, orientation)
        for parent, ps_site, child_allele in (("F", t.father, pat_allele), ("M", t.mother, mat_allele)):
            if ps_site is None or not ps_site.is_het or not ps_site.phased or ps_site.ps is None:
                continue
            if parent == "F" and t.chrom in ("chrX", "X") and not _in_par(t.pos):
                continue                                          # father hemizygous: no paternal X phase blocks
            if parent == "M" and sex == "M" and t.chrom in ("chrX", "X") and not _in_par(t.pos):
                child_allele = t.child.alleles[0] if (t.child.alleles and len(set(t.child.alleles)) == 1) else None
            st = stats[parent]
            key = (parent, t.chrom, ps_site.ps)
            bv = blocks.get(key)
            if bv is None:
                bv = blocks[key] = ParentBlockVotes(parent=parent, chrom=t.chrom, ps=ps_site.ps)
            bv.add_het(t.pos)
            st["n_parent_het_phased"] += 1
            if ps_site.gq is not None and ps_site.gq < p.min_gq or (t.child.gq is not None and t.child.gq < p.min_gq):
                st["n_low_gq"] += 1
                continue
            if child_allele is None:
                st["n_child_allele_unknown"] += 1
                continue
            f1, f2 = ps_site.alleles
            if child_allele == f1 and child_allele != f2:
                bv.votes.append((t.pos, 1)); st["n_informative"] += 1
            elif child_allele == f2 and child_allele != f1:
                bv.votes.append((t.pos, -1)); st["n_informative"] += 1
            else:
                st["n_mendel_inconsistent"] += 1
    segments: List[Segment] = []
    changes: List[ChangePoint] = []
    dp_params = OrientParams(min_sites=p.min_sites, min_frac=p.min_frac, max_splits=p.max_splits,
                             split_min_sites=p.xo_min_sites, keep_dissent_positions=False)
    for _, bv in sorted(blocks.items(), key=lambda kv: (kv[0][0], kv[1].chrom, kv[1].start)):
        votes = sorted(bv.votes)
        cuts = segment_votes(votes, dp_params) if votes else []
        edges = [0] + cuts + [len(votes)]
        n_seg = len(edges) - 1
        st = stats[bv.parent]
        st["n_blocks"] += 1
        if n_seg > 1:
            st["n_blocks_with_change_points"] += 1
            st["n_change_points"] += n_seg - 1
        for i in range(n_seg):
            seg = votes[edges[i]:edges[i + 1]]
            n1 = sum(1 for _, v in seg if v > 0); n2 = len(seg) - n1; n = n1 + n2
            frac = max(n1, n2) / n if n else 0.0
            if n == 0:
                tr, reason = UNRESOLVED, "NO_INFORMATIVE_SITES"
            elif n < p.min_sites:
                tr, reason = UNRESOLVED, "LOW_SITES"
            elif frac < p.min_frac:
                tr, reason = UNRESOLVED, "MIXED_VOTES"
            else:
                tr, reason = (HAP1 if n1 > n2 else HAP2), "OK"
            start = bv.start if i == 0 else votes[edges[i]][0]
            end = bv.end if i == n_seg - 1 else votes[edges[i + 1]][0] - 1
            n_het = sum(1 for h in bv.het_positions if start <= h <= end)
            segments.append(Segment(parent=bv.parent, chrom=bv.chrom, ps=bv.ps, start=start, end=end, segment=i + 1,
                                    n_segments=n_seg, n_het_phased=n_het, n_hap1=n1, n_hap2=n2, transmitted=tr,
                                    vote_frac=round(frac, 4), reason=reason,
                                    left_boundary="BLOCK_EDGE" if i == 0 else "CHANGE_POINT",
                                    right_boundary="BLOCK_EDGE" if i == n_seg - 1 else "CHANGE_POINT",
                                    change_pos=0 if i == 0 else votes[edges[i]][0]))
            if tr == UNRESOLVED:
                st["n_segments_unresolved"] += 1; st["het_unresolved"] += n_het; st["bp_unresolved"] += max(0, end - start)
            else:
                st["n_segments_resolved"] += 1; st["het_resolved"] += n_het; st["bp_resolved"] += max(0, end - start)
        for i in range(1, n_seg):
            left = votes[edges[i - 1]:edges[i]]; right = votes[edges[i]:edges[i + 1]]
            lh = HAP1 if sum(v for _, v in left) > 0 else HAP2
            rh = HAP1 if sum(v for _, v in right) > 0 else HAP2
            changes.append(ChangePoint(parent=bv.parent, chrom=bv.chrom, ps=bv.ps, left_pos=left[-1][0],
                                       right_pos=right[0][0], n_left=len(left), n_right=len(right), left_hap=lh, right_hap=rh))
    return segments, changes, {k: dict(v) for k, v in stats.items()}


SEGMENT_COLUMNS = ["chrom", "start", "end", "parent", "transmitted", "parent_phase_block_id", "segment", "n_segments",
                   "n_het_phased", "n_hap1", "n_hap2", "vote_frac", "reason", "left_boundary", "right_boundary", "change_pos"]
CHANGE_COLUMNS = ["chrom", "parent", "parent_phase_block_id", "left_pos", "right_pos", "resolution_bp", "n_left", "n_right",
                  "left_hap", "right_hap", "status"]


def write_segments(segments: List[Segment], path: str):
    with open(path, "w") as fh:
        fh.write("\t".join(SEGMENT_COLUMNS) + "\n")
        for s in segments:
            fh.write("\t".join(str(x) for x in (s.chrom, s.start, s.end, s.parent, s.transmitted, s.ps, s.segment, s.n_segments,
                                                 s.n_het_phased, s.n_hap1, s.n_hap2, s.vote_frac, s.reason, s.left_boundary,
                                                 s.right_boundary, s.change_pos)) + "\n")


def write_changes(changes: List[ChangePoint], path: str):
    with open(path, "w") as fh:
        fh.write("\t".join(CHANGE_COLUMNS) + "\n")
        for c in changes:
            fh.write("\t".join(str(x) for x in (c.chrom, c.parent, c.ps, c.left_pos, c.right_pos, c.right_pos - c.left_pos,
                                                 c.n_left, c.n_right, c.left_hap, c.right_hap, c.status)) + "\n")


def summarise(stats: Dict[str, Dict[str, int]], p: TransmissionParams) -> dict:
    out = {"params": p.__dict__, "per_parent": {}}
    for parent, st in stats.items():
        het_all = st.get("het_resolved", 0) + st.get("het_unresolved", 0)
        bp_all = st.get("bp_resolved", 0) + st.get("bp_unresolved", 0)
        out["per_parent"][parent] = dict(st, frac_het_resolved=round(st.get("het_resolved", 0) / het_all, 4) if het_all else None,
                                         frac_bp_resolved=round(st.get("bp_resolved", 0) / bp_all, 4) if bp_all else None,
                                         mendel_inconsistent_per_informative=round(st.get("n_mendel_inconsistent", 0) / st["n_informative"], 5)
                                         if st.get("n_informative") else None)
    return out


def write_summary(summary: dict, path: str):
    with open(path, "w") as fh:
        json.dump(summary, fh, indent=1, sort_keys=True)
