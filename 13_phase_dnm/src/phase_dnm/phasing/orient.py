"""M1a — orient the child's read-based phase blocks to parent of origin by Mendelian vote (DESIGN P2).

Generalises 12_x_inactivation/03_trio_phase_x.py (father hemizygous on X) and the four-case rule in
09_methylation/pofo_local.py to any site: with the child heterozygous a|b, the assignment
"hap1 = paternal" is allowed iff a is among the father's alleles and b among the mother's; the
reverse assignment likewise. A site is informative iff exactly one assignment is allowed. A site
where neither is allowed is Mendelian-inconsistent (a de novo candidate or a genotype error) and is
counted, never voted. Sex chromosomes need no special rule for a female child (the father's haploid
call is an allele set of size one); for a male child, chrX outside the PAR and chrY are skipped.

The per-block vote consistency is the diagnostic, exactly as in the X script: a block whose votes
split near 0.5 is a phasing switch, a pedigree error or a sample swap, not an unusual sample.
"""
from __future__ import annotations

import json
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple

from ..io.vcf import TrioSite, normalise_sex

# GRCh38 pseudoautosomal regions (same values as 12_x_inactivation/03_trio_phase_x.py)
PAR_GRCH38 = ((10001, 2781479), (155701383, 156030895))
SKIP_ALWAYS = {"chrY", "Y", "chrM", "MT", "chrEBV"}

HAP1_PAT, HAP1_MAT, AMBIGUOUS = "HAP1_PAT", "HAP1_MAT", "AMBIGUOUS"


def assign(child_alleles: Tuple[int, int], father: set, mother: set) -> Optional[int]:
    """+1: hap1 paternal; -1: hap1 maternal; 0: uninformative; None: Mendelian-inconsistent."""
    a, b = child_alleles
    hap1_pat = (a in father) and (b in mother)
    hap1_mat = (b in father) and (a in mother)
    if hap1_pat and hap1_mat:
        return 0
    if hap1_pat:
        return 1
    if hap1_mat:
        return -1
    return None


@dataclass
class OrientParams:
    min_sites: int = 20        # informative sites needed to orient a block (or a segment after a split)
    min_frac: float = 0.95     # majority fraction needed
    small_min_sites: int = 10  # second tier: a segment with small_min_sites..min_sites-1 votes orients only if ...
    small_min_frac: float = 1.0  # ... its majority fraction reaches this (default: unanimous)
    min_gq: int = 20           # applied to child and both parents
    max_splits: int = 3        # phase-switch errors located per block; 0 disables splitting
    split_min_sites: int = 2   # minimum isolated minority run to accept a change point (a lone genotype error cannot split)
    keep_dissent_positions: bool = True


@dataclass
class BlockVotes:
    chrom: str
    ps: int
    start: int = 1 << 62
    end: int = 0
    het_positions: List[int] = field(default_factory=list)      # every phased het in the block
    votes: List[Tuple[int, int]] = field(default_factory=list)   # (pos, +1/-1), position-ordered

    def add_site(self, pos: int):
        self.het_positions.append(pos)
        self.start = min(self.start, pos)
        self.end = max(self.end, pos)

    @property
    def n_het_phased(self) -> int:
        return len(self.het_positions)


@dataclass
class BlockOrientation:
    chrom: str
    ps: int
    start: int
    end: int
    n_het_phased: int
    n_inf_pat: int
    n_inf_mat: int
    orientation: str
    vote_frac: float
    reason: str
    n_dissent: int
    dissent_positions: List[int]
    segment: int = 1           # 1-based index of this segment within the HiPhase block
    n_segments: int = 1        # number of segments the block was split into (1 = no switch found)
    switch_pos: int = 0        # first informative position AFTER the switch that starts this segment (0 = none)

    @property
    def n_informative(self) -> int:
        return self.n_inf_pat + self.n_inf_mat


def _majority(votes: List[Tuple[int, int]]) -> Tuple[int, int, float]:
    n_pat = sum(1 for _, v in votes if v > 0)
    n_mat = len(votes) - n_pat
    frac = max(n_pat, n_mat) / len(votes) if votes else 0.0
    return n_pat, n_mat, frac


def segment_votes(votes: List[Tuple[int, int]], p: OrientParams) -> List[int]:
    """Change points (indices where a new segment starts) that best explain the block's position-ordered votes
    as alternating runs of opposite sign, i.e. as phase-switch errors.

    Exact dynamic programme over RUN BOUNDARIES (an optimal cut never falls inside a run of equal votes), with
    at most max_splits cuts and a penalty of (split_min_sites - 0.5) per cut, so every accepted cut must raise
    the total majority count by at least split_min_sites: a lone dissenting vote cannot split a block, a run of
    two can. Adjacent segments always end up with opposite majority signs (merging same-sign neighbours costs
    nothing and saves a penalty). Purity and length of each segment are judged afterwards by min_sites /
    min_frac, so a switch close to a block end yields one oriented segment plus a short AMBIGUOUS:LOW_SITES
    tail rather than a discarded block, and a switch-and-back pattern becomes three segments.
    """
    n = len(votes)
    m = max(1, p.split_min_sites)
    if p.max_splits <= 0 or n < 2 * m:
        return []
    n_pat = sum(1 for _, v in votes if v > 0)
    if min(n_pat, n - n_pat) < m:
        return []                                     # fewer dissenters than a splittable run
    bounds = [0] + [i for i in range(1, n) if votes[i][1] != votes[i - 1][1]] + [n]
    cum = [0]
    for _, v in votes:
        cum.append(cum[-1] + (1 if v > 0 else 0))

    def maj(i: int, j: int) -> int:
        pat = cum[j] - cum[i]
        return max(pat, (j - i) - pat)

    B, S, penalty, NEG = len(bounds), p.max_splits + 1, m - 0.5, float("-inf")
    dp = [[NEG] * B for _ in range(S + 1)]
    back = [[-1] * B for _ in range(S + 1)]
    dp[0][0] = 0.0
    for s in range(1, S + 1):
        for b in range(1, B):
            best, arg = NEG, -1
            for a in range(0, b):
                if dp[s - 1][a] == NEG or bounds[b] - bounds[a] < m:
                    continue
                sc = dp[s - 1][a] + maj(bounds[a], bounds[b]) - (penalty if s > 1 else 0.0)
                if sc > best:
                    best, arg = sc, a
            dp[s][b], back[s][b] = best, arg
    s_best = max(range(1, S + 1), key=lambda s: dp[s][B - 1])
    cuts: List[int] = []
    b, s = B - 1, s_best
    while s > 1:
        a = back[s][b]
        cuts.append(bounds[a])
        b, s = a, s - 1
    return sorted(cuts)


def _orientation_row(bv: BlockVotes, votes: List[Tuple[int, int]], start: int, end: int, p: OrientParams,
                     segment: int, n_segments: int, switch_pos: int) -> BlockOrientation:
    n_pat, n_mat, frac = _majority(votes)
    n = n_pat + n_mat
    if n == 0:
        orient, reason, minority = AMBIGUOUS, "NO_INFORMATIVE_SITES", 0
    elif n < p.min_sites and (n < p.small_min_sites or frac < p.small_min_frac):
        orient, reason, minority = AMBIGUOUS, "LOW_SITES", 0
    elif n >= p.min_sites and frac < p.min_frac:
        orient, reason, minority = AMBIGUOUS, "MIXED_VOTES", 0
    else:
        orient = HAP1_PAT if n_pat > n_mat else HAP1_MAT
        reason = ("OK" if n >= p.min_sites else "OK_SMALL_UNANIMOUS") if n_segments == 1 else "SPLIT_AT_SWITCH"
        minority = -1 if orient == HAP1_PAT else 1
    dissent = [pos for pos, v in votes if v == minority] if minority else []
    n_het = sum(1 for h in bv.het_positions if start <= h <= end)
    return BlockOrientation(chrom=bv.chrom, ps=bv.ps, start=start, end=end, n_het_phased=n_het,
                            n_inf_pat=n_pat, n_inf_mat=n_mat, orientation=orient, vote_frac=round(frac, 4),
                            reason=reason, n_dissent=len(dissent),
                            dissent_positions=dissent if p.keep_dissent_positions else [],
                            segment=segment, n_segments=n_segments, switch_pos=switch_pos)


def decide(bv: BlockVotes, p: OrientParams) -> List[BlockOrientation]:
    """Orient one HiPhase block; returns one row, or several when phase-switch errors are located inside it."""
    votes = sorted(bv.votes)
    cuts = segment_votes(votes, p)
    if not cuts:
        return [_orientation_row(bv, votes, bv.start, bv.end, p, 1, 1, 0)]
    edges = [0] + cuts + [len(votes)]
    rows = []
    for i in range(len(edges) - 1):
        seg = votes[edges[i]:edges[i + 1]]
        sw = 0 if i == 0 else votes[edges[i]][0]                 # first informative position after the switch
        start = bv.start if i == 0 else sw
        end = bv.end if i == len(edges) - 2 else votes[edges[i + 1]][0] - 1
        rows.append(_orientation_row(bv, seg, start, end, p, i + 1, len(edges) - 1, sw))
    if all(r.orientation == AMBIGUOUS for r in rows):
        # splitting explained nothing: report the block once as mixed, do not count phantom switches
        return [_orientation_row(bv, votes, bv.start, bv.end, p, 1, 1, 0)]
    return rows


def _in_par(pos: int) -> bool:
    return any(lo <= pos <= hi for lo, hi in PAR_GRCH38)


def orient_child(trio_sites: Iterable[TrioSite], child_sex: str, params: OrientParams
                 ) -> Tuple[List[BlockOrientation], Dict[str, Dict[str, int]]]:
    """One pass over merged trio sites → oriented blocks + per-chromosome counters."""
    sex = normalise_sex(child_sex)
    blocks: Dict[Tuple[str, int], BlockVotes] = {}
    stats: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for t in trio_sites:
        st = stats[t.chrom]
        st["n_sites"] += 1
        c = t.child
        if not c.is_het:
            continue
        st["n_het"] += 1
        if t.chrom in SKIP_ALWAYS or (sex == "M" and t.chrom in ("chrX", "X") and not _in_par(t.pos)):
            st["n_skipped_sex_chrom"] += 1
            continue
        if not c.phased or c.ps is None:
            st["n_het_unphased"] += 1
            continue
        st["n_het_phased"] += 1
        key = (t.chrom, c.ps)
        bv = blocks.get(key)
        if bv is None:
            bv = blocks[key] = BlockVotes(chrom=t.chrom, ps=c.ps)
        bv.add_site(t.pos)
        if c.gq is not None and c.gq < params.min_gq:
            st["n_low_gq"] += 1
            continue
        if t.father is None or t.mother is None or t.father.alleles is None or t.mother.alleles is None:
            st["n_parent_missing"] += 1
            continue
        if (t.father.gq is not None and t.father.gq < params.min_gq) or \
           (t.mother.gq is not None and t.mother.gq < params.min_gq):
            st["n_low_gq"] += 1
            continue
        v = assign((c.alleles[0], c.alleles[1]), t.father.carried, t.mother.carried)
        if v is None:
            st["n_mendel_inconsistent"] += 1
        elif v == 0:
            st["n_uninformative"] += 1
        else:
            st["n_informative"] += 1
            bv.votes.append((t.pos, v))
    result: List[BlockOrientation] = []
    for _, bv in sorted(blocks.items(), key=lambda kv: (kv[1].chrom, kv[1].start)):
        rows = decide(bv, params)
        result.extend(rows)
        st = stats[bv.chrom]
        st["n_blocks"] += 1
        if len(rows) > 1:
            st["n_blocks_split"] += 1
            st["n_switches_located"] += len(rows) - 1
        for b in rows:
            st["n_segments"] += 1
            if b.orientation == AMBIGUOUS:
                st["n_segments_ambiguous"] += 1
                st["bp_ambiguous"] += max(0, b.end - b.start)
                st["het_ambiguous"] += b.n_het_phased
                if b.reason == "MIXED_VOTES":
                    st["bp_mixed_votes"] += max(0, b.end - b.start)
                    st["het_mixed_votes"] += b.n_het_phased
            else:
                st["n_segments_oriented"] += 1
                st["bp_oriented"] += max(0, b.end - b.start)
                st["het_oriented"] += b.n_het_phased
                st["n_dissent"] += b.n_dissent
    return result, {k: dict(v) for k, v in stats.items()}


ORIENTATION_COLUMNS = ["chrom", "phase_block_id", "segment", "n_segments", "start", "end", "switch_pos",
                       "n_het_phased", "n_inf_pat", "n_inf_mat", "n_informative", "vote_frac", "orientation",
                       "reason", "n_dissent"]


def write_orientation(blocks: List[BlockOrientation], path: str):
    with open(path, "w") as fh:
        fh.write("\t".join(ORIENTATION_COLUMNS) + "\n")
        for b in blocks:
            fh.write("\t".join(str(x) for x in (
                b.chrom, b.ps, b.segment, b.n_segments, b.start, b.end, b.switch_pos, b.n_het_phased,
                b.n_inf_pat, b.n_inf_mat, b.n_informative, b.vote_frac, b.orientation, b.reason,
                b.n_dissent)) + "\n")


def write_dissent(blocks: List[BlockOrientation], path: str):
    """Positions voting against their block's orientation: switch-error / genotype-error candidates."""
    with open(path, "w") as fh:
        fh.write("chrom\tpos\tphase_block_id\tblock_orientation\n")
        for b in blocks:
            for pos in b.dissent_positions:
                fh.write("%s\t%d\t%d\t%s\n" % (b.chrom, pos, b.ps, b.orientation))


def summarise(stats: Dict[str, Dict[str, int]], blocks: List[BlockOrientation], params: OrientParams) -> dict:
    tot: Dict[str, int] = defaultdict(int)
    for st in stats.values():
        for k, v in st.items():
            tot[k] += v
    n_segments = len(blocks)
    n_amb = sum(1 for b in blocks if b.orientation == AMBIGUOUS)
    bp_all = tot["bp_oriented"] + tot["bp_ambiguous"]
    het_all = tot["het_oriented"] + tot["het_ambiguous"]
    # Two weightings, deliberately: bp-weighted ambiguity is dominated by a handful of giant blocks spanning
    # pericentromeric het deserts (an 18 Mb chr1q12 block with 384 hets was 0.8% of one genome's bp);
    # het-weighted ambiguity measures what phasing quality actually is - the fraction of phased hets whose
    # parent of origin is unknown - and is the gated quantity (thresholds.yaml).
    return {
        "frac_het_ambiguous": round(tot["het_ambiguous"] / het_all, 4) if het_all else None,
        "frac_het_mixed_votes": round(tot["het_mixed_votes"] / het_all, 5) if het_all else None,
        "frac_bp_mixed_votes": round(tot["bp_mixed_votes"] / bp_all, 5) if bp_all else None,
        "params": params.__dict__,
        "total": dict(tot),
        "per_chrom": stats,
        "n_blocks": tot["n_blocks"],
        "n_blocks_split": tot["n_blocks_split"],
        "n_segments": n_segments,
        "n_segments_oriented": n_segments - n_amb,
        "frac_segments_ambiguous": round(n_amb / n_segments, 4) if n_segments else None,
        "frac_bp_ambiguous": round(tot["bp_ambiguous"] / bp_all, 4) if bp_all else None,
        "mendel_inconsistent_per_informative": round(tot["n_mendel_inconsistent"] / tot["n_informative"], 5)
        if tot["n_informative"] else None,
    }


def write_summary(summary: dict, path: str):
    with open(path, "w") as fh:
        json.dump(summary, fh, indent=1, sort_keys=True)
