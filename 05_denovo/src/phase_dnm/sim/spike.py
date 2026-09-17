"""Spike-in harness (DESIGN R3, §5; README `sim/spike.py`): plant de novo variants of every class into the REAL
haplotagged reads of a trio and check that the six-haplotype review recovers class, parent of origin and the
transmitted-haplotype evidence. Truth is known by construction; nothing here touches folds or thresholds.

Scenarios (per class SNV, INDEL, SV, TR):
  G    germline DNM on child haplotype h          child tagged reads of h: fraction 1.0, untagged 0.5
  CM   child postzygotic mosaic                    child tagged reads of h: fraction c, untagged c/2
  PM   parental mosaic, transmitted                as G, plus the ORIGIN parent's TRANSMITTED haplotype at fraction m
  IM   inherited, missed in the parent             as G, plus the origin parent's transmitted haplotype at 1.0 (untagged 0.5)
Parent of origin of h and the parent's transmitted haplotype come from the M1 label tables at the site, exactly as
the review reads them; sites without labels are not planted (spike-ins test the review, not the phasing).

Three steps, so the expensive one (apply) can be inspected before the review runs:
  plan      choose sites that pass a read-level QC (unanimous readable base on every haplotype, CIGAR-clean window
            for indels/SVs, consistent repeat lengths for TR) and write plan.tsv + candidates.tsv + truth.tsv
  apply     write per-sample SLICE BAMs (all reads in +-pad of every site, edited where the plan says) + .bai
  evaluate  join truth with the reviewed (and likelihood-scored) evidence table -> recovery per class x scenario

The edits themselves are `sim/edit.py` (pure functions, unit-tested); pysam is imported lazily here.
Identifiers never appear in code; read names are left alone (the review hashes them).
"""
from __future__ import annotations

import csv
import json
import os
import random
import statistics as st
from collections import defaultdict
from dataclasses import asdict, dataclass, fields
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from ..records import CandidateRecord, write_candidates
from . import edit as E

SCENARIOS = ("G", "CM", "PM", "IM")
EXPECTED_CLASS = {"G": "germline_DNM_phased", "CM": "child_postzygotic_mosaic",
                  "PM": "parental_mosaic_transmitted", "IM": "inherited_missed_in_parent"}
SUBTYPES = {"SNV": ("SNV",), "INDEL": ("DEL", "INS"), "SV": ("DEL", "INS", "BIGDEL"), "TR": ("EXP",)}
DROP = "DROP"   # a read of the deleted haplotype that lies inside a large deletion: it does not exist in that genome
LENGTHS = {("INDEL", "DEL"): (1, 2, 4, 8), ("INDEL", "INS"): (1, 3, 6, 12),
           ("SV", "DEL"): (200, 500, 1500), ("SV", "INS"): (100, 250, 400), ("TR", "EXP"): (3, 5, 8),
           # BIGDEL spans the size range the pedigree swap cannot supply: after the rarity gate the synthetic positives
           # hold 43 deletions between 10 and 50 kb and 13 above, so the classifier has no way to learn this class from
           # them. These are planted by removing the deleted haplotype's reads and clipping the ones that cross a
           # breakpoint, which is what the genome actually looks like (no read spans the event).
           ("SV", "BIGDEL"): (5000, 20000, 50000)}
FRACTIONS = {"CM": (0.15, 0.3), "PM": (0.1, 0.25)}
BASES = "ACGT"


@dataclass
class PlanRow:
    variant_id: str
    chrom: str
    pos: int                      # 1-based anchor (VCF POS); SNV: the base itself
    variant_class: str            # SNV | INDEL | SV | TR
    subtype: str                  # SNV | DEL | INS | EXP
    length: int                   # bp (TR: bp = units * motif_unit_bp)
    ref: str                      # SNV/INDEL: VCF-style REF; SV/TR: "."
    alt: str                      # SNV/INDEL: VCF-style ALT; INS/SV-INS: inserted sequence; TR: motif; SV-DEL: "."
    scenario: str
    child_hap: int
    child_frac: float
    parent: str                   # F | M | .
    parent_hap: int               # transmitted haplotype of the origin parent (0 when unused)
    parent_frac: float
    expected_poo: str             # paternal | maternal
    expected_class: str
    child_ps: int
    parent_ps: int
    trid: str = "."
    motif_unit_bp: int = 0
    child_base_al: str = "."      # TR: per-haplotype repeat lengths before planting, "h1,h2"
    father_al: str = "."
    mother_al: str = "."
    seed: int = 0
    locus_end: int = 0            # TR: 1-based locus end (insertion goes to the locus midpoint)


PLAN_COLS = [f.name for f in fields(PlanRow)]


# ----------------------------------------------------------------------------------------------
# helpers on reads (pysam objects in, plain data out)
# ----------------------------------------------------------------------------------------------
def aln_of(read) -> E.Aln:
    q = read.query_qualities
    return E.Aln(read.reference_start, list(read.cigartuples or []), read.query_sequence or "", list(q) if q is not None else None)


def write_back(read, a: E.Aln) -> None:
    # reference_start too: apply_breakpoint(keep="right") moves the alignment start to the breakpoint, and until
    # 2026-09-16 that new start was never written, so every right-hand junction read of a planted deletion kept its
    # ORIGINAL start with a leading soft clip -- its aligned block sat one clip-length left of the breakpoint, no
    # junction matcher could place it there, and the far breakpoint of every 20/50 kb planted deletion showed ~0
    # junction reads. That, not the six-haplotype matrix, was the deficit P31 recorded.
    read.reference_start = a.reference_start
    read.cigartuples = a.cigar
    read.query_sequence = a.seq
    read.query_qualities = a.qual
    for tag in ("MD",):
        if read.has_tag(tag):
            read.set_tag(tag, None)


def _hp(read) -> Optional[int]:
    v = read.get_tag("HP") if read.has_tag("HP") else None
    return v if v in (1, 2) else None


def _ps(read) -> Optional[int]:
    return read.get_tag("PS") if read.has_tag("PS") else None


def _primary(read, min_mapq: int) -> bool:
    return not (read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate) and read.mapping_quality >= min_mapq


def clean_window(a: E.Aln, s0: int, e0: int) -> bool:
    """Aligned through [s0, e0) by M ops only (no indel, no clip inside)."""
    ref = a.reference_start
    if ref > s0 or a.reference_end < e0:
        return False
    for op, L in a.cigar:
        if op in (E.M, E.EQ, E.X):
            nxt = ref + L
            if ref <= s0 and nxt >= e0:
                return True
            ref = nxt
        elif op in (E.D, E.N):
            if ref < e0 and ref + L > s0:
                return False
            ref += L
        elif op == E.I:
            if s0 < ref < e0:
                return False
        if ref >= e0:
            break
    return False


def consensus_base(alns: Sequence[E.Aln], pos0: int) -> Tuple[Optional[str], int, float]:
    """(majority base, n readable, majority fraction) at a 0-based position."""
    c: Dict[str, int] = defaultdict(int)
    for a in alns:
        qi = E.query_index_at(a, pos0)
        if qi is not None:
            c[a.seq[qi].upper()] += 1
    n = sum(c.values())
    if not n:
        return None, 0, 0.0
    b = max(c, key=c.get)
    return b, n, c[b] / n


def consensus_seq(alns: Sequence[E.Aln], s0: int, e0: int) -> Optional[str]:
    out = []
    for p in range(s0, e0):
        b, n, f = consensus_base(alns, p)
        if b is None or f < 0.9 or b == "N":
            return None
        out.append(b)
    return "".join(out)


# ----------------------------------------------------------------------------------------------
# planning
# ----------------------------------------------------------------------------------------------
@dataclass
class SiteQC:
    k: int = 5                    # readable reads per haplotype per sample
    unanimity: float = 0.97
    clean_frac: float = 0.9       # fraction of reads (all samples) with an indel-free window
    margin: int = 30
    min_mapq: int = 20
    pad: int = 2500               # slice half-width around a site (>= review windows + read overhang)
    # An SV deletion is the only class whose evidence is measured OUTSIDE the event: sv_interval_evidence compares
    # depth inside the interval against a flank probe at start - 2300 .. start - 300 (bp_window + probe_bp; note that
    # its `flank_bp=20000` argument is declared and never used, so the probe is near, not 22 kb out). A read covering
    # that probe may begin a full read length earlier, and reads outside the slice were never written -- so a narrow
    # slice truncates the flank depth at its own edge and biases every ratio low. The pad is therefore a read length.
    big_pad: int = 25000          # slice half-width for SV deletions: one HiFi read length beyond the flank probe
    big_mask_frac: float = 0.5    # BIGDEL: breakpoints must be mask-free; this much INTERIOR mask overlap is allowed
    spacing: int = 60000          # min distance between planted sites (no read spans two)
    max_tries_per_site: int = 400


def _slice_window(pos1: int, cls: str, sub: str, length: int, qc: SiteQC) -> Tuple[int, int]:
    """The 0-based BAM slice apply_plan writes for a site. The planner and apply_plan MUST agree on it: the planner
    keeps these windows `spacing` apart, which is what makes "no read is written twice" true."""
    pad = qc.big_pad if (cls == "SV" and sub in ("DEL", "BIGDEL")) else qc.pad
    return pos1 - 1 - pad, pos1 - 1 + max(0, length) + pad


def _clear(used_chrom: List[Tuple[int, int]], lo: int, hi: int, spacing: int) -> bool:
    """Is a slice [lo, hi) at least `spacing` away from every slice already reserved on this contig?

    apply_plan writes one fetch window per site in coordinate order, and fetch returns every read OVERLAPPING a
    window -- including reads that begin far to its left. Two windows closer together than a read length therefore
    make apply_plan write the same read twice AND emit it out of order, leaving a BAM that samtools cannot index.
    Anchor-to-anchor distance does not capture this: a 50 kb deletion's slice reaches 75 kb past its anchor, more
    than the whole spacing budget, which is how a site 60,490 bp away still landed 13 kb inside its neighbour."""
    return all(hi + spacing <= u_lo or lo >= u_hi + spacing for u_lo, u_hi in used_chrom)


def _reads_at(bam, chrom: str, s0: int, e0: int, qc: SiteQC):
    return [r for r in bam.fetch(chrom, max(0, s0), e0) if _primary(r, qc.min_mapq) and r.reference_start <= s0 and (r.reference_end or 0) >= e0]


def _labels_at(labels, chrom: str, pos1: int, child_reads, parent_reads: Dict[str, list]):
    """(child_hap1_is, child_ps, {F: (transmitted, ps), M: ...}) from the majority PS of the spanning reads."""
    def maj_ps(reads):
        c: Dict[int, int] = defaultdict(int)
        for r in reads:
            p = _ps(r)
            if p is not None and _hp(r) is not None:
                c[p] += 1
        return max(c, key=c.get) if c else None
    cps = maj_ps(child_reads)
    h1 = labels.child_hap1_is(chrom, cps, pos1)
    par = {}
    for role in ("F", "M"):
        pps = maj_ps(parent_reads[role])
        par[role] = (labels.parent_transmitted(role, chrom, pps, pos1), pps)
    return h1, cps, par


def _scenario_grid(rng: random.Random, n_per: int) -> List[Tuple[str, str, str, int, float, float]]:
    """(class, subtype, scenario, length, child_frac, parent_frac) x n_per, shuffled."""
    grid = []
    for cls, subs in SUBTYPES.items():
        for sub in subs:
            lengths = LENGTHS.get((cls, sub), (1,))
            for sc in SCENARIOS:
                for i in range(n_per):
                    L = lengths[i % len(lengths)]
                    cf = 1.0 if sc in ("G", "PM", "IM") else FRACTIONS["CM"][i % 2]
                    pf = {"G": 0.0, "CM": 0.0, "PM": FRACTIONS["PM"][i % 2], "IM": 1.0}[sc]
                    grid.append((cls, sub, sc, L, cf, pf))
    rng.shuffle(grid)
    return grid


def plan_sites(bams: Dict[str, object], tr_bams: Dict[str, object], labels, regions: List[Tuple[str, int, int]],
               n_per: int, seed: int, family: str, child: str, qc: SiteQC, trgt_vcf: Optional[str] = None,
               mask=None, log=None, trgt_vcf_index: Optional[str] = None) -> Tuple[List[PlanRow], List[CandidateRecord], Dict[str, int]]:
    rng = random.Random(seed)
    grid = _scenario_grid(rng, n_per)
    used: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    plan: List[PlanRow] = []
    cands: List[CandidateRecord] = []
    stats: Dict[str, int] = defaultdict(int)
    tr_loci = _tr_loci(trgt_vcf, regions, child, trgt_vcf_index) if trgt_vcf else []
    rng.shuffle(tr_loci)
    tr_iter = iter(tr_loci)
    stats["planned"] = len(grid)
    for n_item, (cls, sub, sc, L, cf, pf) in enumerate(grid):
        placed = False
        for _ in range(qc.max_tries_per_site):
            stats["tries"] += 1
            if cls == "TR":
                locus = next(tr_iter, None)
                if locus is None:
                    break
                lo, hi = (locus["start"] - qc.pad, int(locus.get("end") or locus["start"]) + qc.pad)
                if not _clear(used[locus["chrom"]], lo, hi, qc.spacing):
                    stats["skip_spacing"] += 1
                    continue
                site = _try_tr_site(locus, tr_bams, labels, rng, qc, sc, L, cf, pf, family, child, n_item, seed)
            else:
                chrom, rs, re_ = regions[rng.randrange(len(regions))]
                # the anchor must leave room for the event AND its whole slice at both ends
                lo_pad = qc.big_pad if (cls == "SV" and sub in ("DEL", "BIGDEL")) else qc.pad
                hi = re_ - L - lo_pad + 1
                if hi <= rs + lo_pad:
                    stats["skip_window"] += 1
                    continue
                pos1 = rng.randrange(rs + lo_pad, hi)
                w_lo, w_hi = _slice_window(pos1, cls, sub, L, qc)
                if not _clear(used[chrom], w_lo, w_hi, qc.spacing):
                    stats["skip_spacing"] += 1
                    continue
                if mask is not None:
                    if sub == "BIGDEL":
                        # Requiring zero mask overlap across a 50 kb interval is effectively impossible in this
                        # genome, which is why the first run that placed BIGDEL at all placed 24 at 5 kb and none at
                        # 20 or 50 kb: skip_mask was the largest bucket. What the evidence layer actually needs is
                        # clean BREAKPOINTS; interior segdup is tolerable up to a fraction, and is flagged anyway.
                        bad = (mask.overlap_bp(chrom, pos1 - qc.margin, pos1 + qc.margin) > 0
                               or mask.overlap_bp(chrom, pos1 + L - qc.margin, pos1 + L + qc.margin) > 0
                               or mask.overlap_bp(chrom, pos1, pos1 + L) > qc.big_mask_frac * L)
                    else:
                        bad = mask.overlap_bp(chrom, pos1, pos1 + max(1, L)) > 0
                    if bad:
                        stats["skip_mask"] += 1
                        continue
                site = _try_site(chrom, pos1, cls, sub, L, bams, labels, rng, qc, sc, cf, pf, family, child, n_item, seed, stats)
            if site is None:
                continue
            row, rec = site
            used[row.chrom].append(_slice_window(row.pos, row.variant_class, row.subtype, row.length, qc))
            plan.append(row); cands.append(rec)
            placed = True
            break
        stats["placed" if placed else "unplaced"] += 1
        if log and n_item % 25 == 0:
            log("spike plan: %d/%d items, %d placed" % (n_item + 1, len(grid), len(plan)))
    return plan, cands, dict(stats)


def _window_qc(chrom, w0, w1, anchor0, bams, qc, stats, check_clean):
    """QC one window: every haplotype of every sample readable and unanimous at `anchor0`, and optionally a clean
    indel-free alignment across [w0, w1). Returns (reads, alns, allr) or None."""
    reads = {role: _reads_at(b, chrom, w0, w1, qc) for role, b in bams.items()}
    alns = {role: [aln_of(r) for r in rs] for role, rs in reads.items()}
    for role in ("C", "F", "M"):
        by_hp: Dict[int, List[E.Aln]] = defaultdict(list)
        for r, a in zip(reads[role], alns[role]):
            h = _hp(r)
            if h:
                by_hp[h].append(a)
        for h in (1, 2):
            b, n, f = consensus_base(by_hp[h], anchor0)
            if n < qc.k or f < qc.unanimity or b is None or b == "N":
                stats["skip_hap_qc"] += 1
                return None
    allr = [a for role in alns for a in alns[role]]
    if check_clean:
        clean = sum(1 for a in allr if clean_window(a, w0, w1))
        if not allr or clean / len(allr) < qc.clean_frac:
            stats["skip_window"] += 1
            return None
    return reads, alns, allr


def _bp_homology(left_alns, right_alns, l0: int, r0: int, limit: int) -> int:
    """Microhomology at a deletion's two breakpoints: how many bases just inside the left breakpoint are identical to
    the bases just inside the right one, read off the same consensus the site was planted from.

    That count is the number of positions the identical deletion could equally be placed at, which is what a caller
    reports as HOMLEN. The planter hardcoded it to None, so bp_homology_len was empty for every synthetic SV while
    real sawfish calls carry it -- the classifier could only ever see that column on one side of the comparison."""
    n = 0
    while n < limit:
        a = consensus_base(left_alns, l0 + n)[0]
        b = consensus_base(right_alns, r0 + n)[0]
        if a is None or b is None or a == "N" or a != b:
            break
        n += 1
    return n


def _try_site(chrom, pos1, cls, sub, L, bams, labels, rng, qc, sc, cf, pf, family, child, n_item, seed, stats):
    pos0 = pos1 - 1
    if cls == "SV":
        # An SV is validated at its BREAKPOINTS, never across its span. _reads_at keeps only reads covering the
        # window end to end and clean_window then demands that whole span be indel-free, so the gate gets harder the
        # longer the event is and is outright unsatisfiable past a read length. e449b8c granted BIGDEL the breakpoint
        # treatment; leaving the rest of the SV class on the span gate hid the same defect one size down. Of the SV
        # deletions that did place, EVERY one was the shortest length in the grid (23 of 23 at 200 bp, none at 500 or
        # 1500) -- and 200 bp is below readers.min_interval_bp, so no SV row in the whole spike ever reached the
        # interval evidence at all. The span was never the thing being planted: apply_plan edits the breakpoints.
        got = _window_qc(chrom, pos0 - qc.margin, pos0 + qc.margin + 1, pos0, bams, qc, stats, True)
        if got is None:
            return None
        right = None
        if sub in ("DEL", "BIGDEL"):
            right = _window_qc(chrom, pos0 + L - qc.margin, pos0 + L + qc.margin + 1, pos0 + L, bams, qc, stats, True)
            if right is None:
                return None
        reads, alns, allr = got
        homlen = _bp_homology(allr, right[2], pos0, pos0 + L, qc.margin) if right is not None else 0
    else:
        s0, e0 = pos0 - qc.margin, pos0 + max(1, L) + qc.margin + 1
        got = _window_qc(chrom, s0, e0, pos0, bams, qc, stats, cls != "SNV")
        if got is None:
            return None
        reads, alns, allr = got
    h1, cps, par = _labels_at(labels, chrom, pos1, reads["C"], {"F": reads["F"], "M": reads["M"]})
    if h1 is None or cps is None:
        stats["skip_unoriented"] += 1
        return None
    child_hap = rng.choice((1, 2))
    origin = ("F" if h1 == "P" else "M") if child_hap == 1 else ("M" if h1 == "P" else "F")
    t_hap, pps = par[origin]
    if sc in ("PM", "IM") and (t_hap is None or pps is None):
        stats["skip_untransmitted"] += 1
        return None
    # the variant itself, from the reads' consensus (reference-free)
    # Only an INDEL deletion writes its deleted bases into REF. An SV deletion is symbolic (<DEL> plus END), so it
    # needs the anchor base alone -- asking for a unanimous L+1 base consensus was a second span-wide gate on exactly
    # the sites the first one had already thinned.
    ref_seq = consensus_seq(allr, pos0, pos0 + (L + 1 if (cls, sub) == ("INDEL", "DEL") else 1))
    if ref_seq is None:
        stats["skip_consensus"] += 1
        return None
    if cls == "SNV":
        alt = rng.choice([b for b in BASES if b != ref_seq[0]])
        ref, vid, end = ref_seq[0], "spike:%s:%d:%s:%s" % (chrom, pos1, ref_seq[0], alt), pos1
        payload = {"spike": True}
        rec_ref, rec_alt = ref, alt
    elif cls == "INDEL":
        if sub == "DEL":
            ref, alt = ref_seq, ref_seq[0]
        else:
            ins = "".join(rng.choice(BASES) for _ in range(L))
            ref, alt = ref_seq[0], ref_seq[0] + ins
        vid, end = "spike:%s:%d:%s:%s" % (chrom, pos1, ref, alt), pos1 + (L if sub == "DEL" else 0)
        payload = {"spike": True}
        rec_ref, rec_alt = ref, alt
    else:  # SV
        if sub in ("DEL", "BIGDEL"):
            # BIGDEL is this planner's size label, not a callable type: what a caller emits for a 5-50 kb event is a
            # DEL with an END. Recording it as anything else forfeits the five depth/junction features outright,
            # because sv_interval_evidence returns {} unless svtype is an interval type AND end > start. Every
            # BIGDEL planted so far fell through to the insertion branch and was written as <INS:5000bp> with
            # end == start, so the one instrument built for large deletions never saw a single one of them.
            svtype, svlen, end, alt = "DEL", -L, pos1 + L, "<DEL>"
        else:
            svtype, svlen, end, alt = sub, L, pos1, "".join(rng.choice(BASES) for _ in range(L))
        vid = "spike:%s:%d:%s:%d" % (chrom, pos1, sub, L)     # the id keeps BIGDEL: recall is reported by event size
        payload = {"svtype": svtype, "svlen": svlen, "end": end, "imprecise": False, "mateid": None, "svclaim": None,
                   "caller_id": vid, "child_cn": None, "father_cn": None, "mother_cn": None, "insseq_len": L if sub == "INS" else None,
                   "homlen": homlen, "spike": True}
        rec_ref, rec_alt = ref_seq[0], (alt if sub in ("DEL", "BIGDEL") else "<INS:%dbp>" % L)
    row = PlanRow(variant_id=vid, chrom=chrom, pos=pos1, variant_class=cls, subtype=sub, length=L, ref=rec_ref if cls != "SV" else ".",
                  alt=(alt if cls != "SV" or sub == "INS" else "."), scenario=sc, child_hap=child_hap, child_frac=cf,
                  parent=origin if sc in ("PM", "IM") else ".", parent_hap=t_hap or 0, parent_frac=pf,
                  expected_poo="paternal" if origin == "F" else "maternal", expected_class=EXPECTED_CLASS[sc],
                  child_ps=cps, parent_ps=pps or 0, seed=seed)
    rec = CandidateRecord(family_id=family, sample_id=child, variant_id=vid, chrom=chrom, start=pos1, end=end, ref=rec_ref, alt=rec_alt,
                          variant_class=cls, caller="spike", class_payload=payload)
    return row, rec


def _tr_loci(trgt_vcf: str, regions, child: str, index: Optional[str] = None) -> List[dict]:
    import pysam
    out = []
    vf = pysam.VariantFile(trgt_vcf, index_filename=index) if index else pysam.VariantFile(trgt_vcf)
    if child not in list(vf.header.samples):
        return out
    for chrom, s, e in regions:
        try:
            it = vf.fetch(chrom, s, e)
        except ValueError:
            continue
        for rec in it:
            motifs = rec.info.get("MOTIFS")
            motifs = ",".join(motifs) if isinstance(motifs, tuple) else (motifs or "")
            unit = max(1, min((len(x) for x in motifs.split(",") if x), default=1))
            out.append({"trid": rec.info.get("TRID"), "chrom": chrom, "start": rec.pos, "end": rec.stop, "motifs": motifs, "unit": unit})
    vf.close()
    return out


def _try_tr_site(locus, tr_bams, labels, rng, qc, sc, n_units, cf, pf, family, child, n_item, seed):
    chrom, start1, end1, trid, unit = locus["chrom"], locus["start"], locus["end"], locus["trid"], locus["unit"]
    from ..evidence.readers import tr_read_length
    lengths: Dict[str, Dict[int, List[int]]] = {r: defaultdict(list) for r in ("C", "F", "M")}
    reads: Dict[str, list] = {}
    for role, bam in tr_bams.items():
        rs = []
        for r in bam.fetch(chrom, max(0, start1 - 1), end1 + 1):
            if not _primary(r, 0) or not r.has_tag("TR") or r.get_tag("TR") != trid:
                continue
            L = tr_read_length(r)
            if L is None:
                continue
            rs.append(r)
            h = _hp(r)
            if h:
                lengths[role][h].append(L)
        reads[role] = rs
    for role in ("C", "F", "M"):
        for h in (1, 2):
            v = lengths[role][h]
            if len(v) < qc.k or (len(v) > 1 and st.pstdev(v) > 2.0):
                return None
    al = {role: (int(st.median(lengths[role][1])), int(st.median(lengths[role][2]))) for role in ("C", "F", "M")}
    h1, cps, par = _labels_at(labels, chrom, start1, reads["C"], {"F": reads["F"], "M": reads["M"]})
    if h1 is None or cps is None:
        return None
    child_hap = rng.choice((1, 2))
    origin = ("F" if h1 == "P" else "M") if child_hap == 1 else ("M" if h1 == "P" else "F")
    t_hap, pps = par[origin]
    if sc in ("PM", "IM") and (t_hap is None or pps is None):
        return None
    delta = n_units * unit
    base = al["C"][child_hap - 1]
    planted = base + delta
    other = al["C"][2 - child_hap]
    # separable from every competing allele by more than the tolerance (support_tr caps it at half the gap)
    if min(abs(planted - o) for o in (other, *al["F"], *al["M"])) < 2 * unit + 2:
        return None
    motif = locus["motifs"].split(",")[0] or "A"
    vid = "spike:%s:EXP:%d" % (trid, delta)
    payload = {"trid": trid, "motifs": locus["motifs"], "struc": None, "motif_unit_bp": unit,
               "child_AL": [other, planted], "father_AL": list(al["F"]), "mother_AL": list(al["M"]),
               "outlier_allele_idx": 1, "direction": "expansion", "delta_bp": planted - max(al["F"] + al["M"]),
               "delta_units": round((planted - max(al["F"] + al["M"])) / unit, 2), "spike": True}
    row = PlanRow(variant_id=vid, chrom=chrom, pos=start1, variant_class="TR", subtype="EXP", length=delta, ref=".", alt=motif,
                  scenario=sc, child_hap=child_hap, child_frac=cf, parent=origin if sc in ("PM", "IM") else ".",
                  parent_hap=t_hap or 0, parent_frac=pf, expected_poo="paternal" if origin == "F" else "maternal",
                  expected_class=EXPECTED_CLASS[sc], child_ps=cps, parent_ps=pps or 0, trid=trid, motif_unit_bp=unit,
                  child_base_al="%d,%d" % al["C"], father_al="%d,%d" % al["F"], mother_al="%d,%d" % al["M"], seed=seed, locus_end=end1)
    rec = CandidateRecord(family_id=family, sample_id=child, variant_id=vid, chrom=chrom, start=start1, end=end1, ref="<TR>",
                          alt="<AL=%d>" % planted, variant_class="TR", caller="spike", class_payload=payload)
    return row, rec


# ----------------------------------------------------------------------------------------------
# applying the plan to reads
# ----------------------------------------------------------------------------------------------
def _edit_one(a: E.Aln, row: PlanRow, rng: random.Random) -> Optional[E.Aln]:
    pos0 = row.pos - 1
    if row.variant_class == "SNV":
        return E.apply_snv(a, pos0, row.alt)
    if row.variant_class == "TR":
        motif = row.alt if row.alt not in (".", "") else "A"
        ins = (motif * (row.length // len(motif) + 1))[:row.length]
        mid = (row.pos + max(row.pos, row.locus_end)) // 2
        return E.apply_insertion(a, mid - 1, ins, margin=10)
    if row.subtype == "BIGDEL":
        start, end = row.pos, row.pos + row.length
        r_s, r_e = a.reference_start, a.reference_end
        if r_s >= start and r_e <= end:
            return DROP                                             # wholly inside the deletion: this read does not exist
        if r_s < start and r_e > end:
            return E.apply_deletion(a, start, row.length)           # a read long enough to span it keeps a D operation
        if r_s < start < r_e:
            return E.apply_breakpoint(a, start, "left")             # crosses the left breakpoint -> clipped junction read
        if r_s < end < r_e:
            return E.apply_breakpoint(a, end, "right")
        return None
    if row.subtype == "DEL":
        return E.apply_deletion(a, row.pos, row.length)             # deleted bases are 0-based [pos, pos+L): right after the anchor
    ins = row.alt[1:] if row.variant_class == "INDEL" else row.alt
    return E.apply_insertion(a, pos0, ins)


def _decide(read, role: str, row: PlanRow, rng: random.Random) -> bool:
    """Does this read get the edit?"""
    h = _hp(read)
    if role == "C":
        if h == row.child_hap:
            return rng.random() < row.child_frac
        if h is None:
            return rng.random() < row.child_frac * 0.5
        return False
    if role == row.parent and row.parent_frac > 0:
        if h == row.parent_hap:
            return rng.random() < row.parent_frac
        if h is None:
            return rng.random() < row.parent_frac * 0.5
    return False


def apply_plan(plan: List[PlanRow], bams: Dict[str, object], out_dir: str, tag: str, qc: SiteQC, seed: int,
               tr: bool = False, log=None, ledger_path: Optional[str] = None) -> Dict[str, int]:
    """Write <out_dir>/<role>.<tag>.bam with every primary/secondary/supplementary read in +-pad of each planted site of
    the matching class group (tr=False: SNV/INDEL/SV from the genome BAMs; tr=True: TR from the TRGT spanning BAMs),
    edited where the plan says. The planner keeps these windows >= spacing apart (see _clear), which is what makes
    "no read is written twice" true and the output coordinate-sorted.

    `ledger_path` additionally records, per site and per role, how many reads carry the planted allele (edited, plus
    the ones removed because they fell inside a large deletion) and how many carry the reference. Those are the exact
    allelic depths of the planted genotype, and sim/genotype.py turns them into the caller block for the SV and TR
    classes, which no caller can be re-run on for a slice (P29)."""
    import pysam
    rows = [r for r in plan if (r.variant_class == "TR") == tr]
    stats: Dict[str, int] = defaultdict(int)
    ledger: Dict[Tuple[str, str], Dict[str, int]] = defaultdict(lambda: {"alt": 0, "ref": 0, "dropped": 0})
    rng = random.Random(seed + (1 if tr else 0))
    os.makedirs(out_dir, exist_ok=True)
    for role, bam in bams.items():
        out_path = os.path.join(out_dir, "%s.%s.bam" % (role, tag))
        rows.sort(key=lambda r: (bam.get_tid(r.chrom), r.pos))
        # Fail here, naming the two sites, rather than 6 minutes of planning later inside samtools index with
        # "Unsorted positions on sequence #20". Overlapping windows do not merely break the sort order: they write
        # the reads in the overlap TWICE, which inflates exactly the depth the SV features are read from.
        prev_chrom, prev_lo, prev_hi, prev_id = None, 0, 0, ""
        for row in rows:
            lo, hi = _slice_window(row.pos, row.variant_class, row.subtype, row.length, qc)
            if row.chrom == prev_chrom and lo < prev_hi:
                raise ValueError("spike apply: slice windows overlap on %s -- %s [%d,%d) then %s [%d,%d); the "
                                 "planner must keep them %d bp apart"
                                 % (row.chrom, prev_id, prev_lo, prev_hi, row.variant_id, lo, hi, qc.spacing))
            prev_chrom, prev_lo, prev_hi, prev_id = row.chrom, lo, hi, row.variant_id
        # Reads are written in fetch order, i.e. sorted by their ORIGINAL start. write_back moves a right-clipped junction
        # read to the edited start, so the stream is no longer coordinate-sorted: sort before indexing (the slices are small).
        unsorted = out_path + ".unsorted.bam"
        with pysam.AlignmentFile(unsorted, "wb", template=bam) as out:
            for row in rows:
                s0, e0 = _slice_window(row.pos, row.variant_class, row.subtype, row.length, qc)
                for read in bam.fetch(row.chrom, max(0, s0), e0):
                    if read.is_unmapped:
                        continue
                    stats["reads_%s" % role] += 1
                    if tr and (not read.has_tag("TR") or read.get_tag("TR") != row.trid):
                        out.write(read); continue
                    spanning = read.reference_start <= row.pos - 1 and (read.reference_end or 0) >= row.pos
                    if _primary(read, 0) and _decide(read, role, row, rng):
                        a = _edit_one(aln_of(read), row, rng)
                        if a is DROP:
                            stats["dropped_%s" % role] += 1
                            ledger[(row.variant_id, role)]["alt"] += 1
                            ledger[(row.variant_id, role)]["dropped"] += 1
                            continue                                # the read is not written: that haplotype is deleted here
                        if a is not None:
                            ledger[(row.variant_id, role)]["alt"] += 1
                            write_back(read, a)
                            if row.subtype == "BIGDEL":
                                # a clipped junction read carries a supplementary alignment at the other breakpoint;
                                # the SA tag is what the SV adapter and sawfish both read as junction evidence
                                other = row.pos + row.length if read.reference_start < row.pos else row.pos
                                read.set_tag("SA", "%s,%d,%s,%dM,60,0;" % (row.chrom, other + 1, "-" if read.is_reverse else "+",
                                                                           max(1, read.query_length // 2)))
                            stats["edited_%s" % role] += 1
                        else:
                            stats["edit_refused_%s" % role] += 1
                            if spanning:
                                ledger[(row.variant_id, role)]["ref"] += 1
                    elif spanning and _primary(read, 0):
                        ledger[(row.variant_id, role)]["ref"] += 1
                    out.write(read)
        pysam.sort("-o", out_path, unsorted)
        os.remove(unsorted)
        pysam.index(out_path)
        if log:
            log("spike apply: %s -> %s (%d reads, %d edited)" % (role, os.path.basename(out_path), stats["reads_%s" % role], stats["edited_%s" % role]))
    if ledger_path:
        new = not os.path.exists(ledger_path)
        with open(ledger_path, "a", newline="\n") as fh:
            if new:
                fh.write("variant_id\trole\talt\tref\tdropped\n")
            for (vid, role), d in sorted(ledger.items()):
                fh.write("%s\t%s\t%d\t%d\t%d\n" % (vid, role, d["alt"], d["ref"], d["dropped"]))
    return dict(stats)


# ----------------------------------------------------------------------------------------------
# tables
# ----------------------------------------------------------------------------------------------
def write_plan(rows: List[PlanRow], path: str) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=PLAN_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(asdict(r))


def read_plan(path: str) -> List[PlanRow]:
    out = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            kw = {}
            for f in fields(PlanRow):
                v = r[f.name]
                kw[f.name] = (int(v) if f.type == "int" else float(v) if f.type == "float" else v)
            out.append(PlanRow(**kw))
    return out


# ----------------------------------------------------------------------------------------------
# evaluation
# ----------------------------------------------------------------------------------------------
def evaluate(plan: List[PlanRow], evidence_rows: Iterable[dict], k: int = 5) -> Tuple[List[dict], List[dict]]:
    """Per planted site: observed? class recovered? parent of origin recovered? Then a summary per class x scenario.
    'observable' = all six haplotypes observed at k (hap_obs_k<k>); recovery is reported both overall and among
    observable sites, because a spike-in on an unobservable haplotype tests the phasing coverage, not the review."""
    ev = {r["variant_id"]: r for r in evidence_rows}
    present = {r.get("variant_class") for r in ev.values()}
    if present:
        plan = [p for p in plan if p.variant_class in present]       # one evidence table per class group
    post_key = {"G": "lik_post_germline", "CM": "lik_post_child_mosaic", "PM": "lik_post_parental_mosaic", "IM": "lik_post_inherited"}
    per_site = []
    for p in plan:
        r = ev.get(p.variant_id)
        d = dict(variant_id=p.variant_id, variant_class=p.variant_class, subtype=p.subtype, length=p.length, scenario=p.scenario,
                 child_frac=p.child_frac, parent_frac=p.parent_frac, expected_class=p.expected_class, expected_poo=p.expected_poo,
                 reviewed=int(r is not None))
        if r is None:
            per_site.append(d); continue
        obs = r.get("hap_obs_k%d" % k, "")
        cls = r.get("phase_class", "")
        d.update(phase_class=cls, parent_of_origin=r.get("parent_of_origin", ""), rule_score=r.get("rule_score", ""),
                 phase_score=r.get("phase_score", ""), flags=r.get("flags", ""), hap_obs=obs,
                 observable=int(obs not in ("", None) and int(float(obs)) == 6),
                 class_ok=int(cls == p.expected_class),
                 class_ok_lenient=int(cls == p.expected_class or (p.scenario == "G" and cls == "germline_DNM_unphased")),
                 poo_ok=int(r.get("parent_of_origin") == p.expected_poo),
                 post_expected=r.get(post_key[p.scenario], ""),
                 c_alt_hapA=r.get("c_alt_hapA", ""), c_dp_hapA=r.get("c_dp_hapA", ""), t_alt_reads=r.get("t_alt_reads", ""), t_dp=r.get("t_dp", ""))
        per_site.append(d)
    summary = []
    groups: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    for d in per_site:
        groups[(d["variant_class"], d["scenario"])].append(d)
    for (cls, sc), rows in sorted(groups.items()):
        rev = [d for d in rows if d["reviewed"]]
        obs = [d for d in rev if d.get("observable")]
        def frac(xs, key):
            return round(sum(d[key] for d in xs) / len(xs), 3) if xs else None
        def mean(xs, key):
            v = [float(d[key]) for d in xs if d.get(key) not in ("", None)]
            return round(st.mean(v), 3) if v else None
        called: Dict[str, int] = defaultdict(int)
        for d in rev:
            called[d.get("phase_class") or "."] += 1
        summary.append(dict(variant_class=cls, scenario=sc, expected_class=EXPECTED_CLASS[sc], n_planted=len(rows), n_reviewed=len(rev),
                            n_observable=len(obs), class_ok_all=frac(rev, "class_ok"), class_ok_lenient_all=frac(rev, "class_ok_lenient"),
                            class_ok_observable=frac(obs, "class_ok"), poo_ok_all=frac(rev, "poo_ok"), poo_ok_observable=frac(obs, "poo_ok"),
                            mean_phase_score=mean(rev, "phase_score"), mean_phase_score_observable=mean(obs, "phase_score"),
                            mean_post_expected=mean(rev, "post_expected"), mean_rule_score=mean(rev, "rule_score"),
                            called_classes=";".join("%s:%d" % kv for kv in sorted(called.items(), key=lambda kv: -kv[1]))))
    return per_site, summary


def write_rows(rows: List[dict], path: str) -> None:
    if not rows:
        open(path, "w").close(); return
    cols = list(rows[0].keys())
    for r in rows[1:]:
        for c in r:
            if c not in cols:
                cols.append(c)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: ("" if r.get(c) is None else r.get(c)) for c in cols})


def parse_regions(spec: str) -> List[Tuple[str, int, int]]:
    out = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        chrom, rng = part.split(":")
        s, e = rng.replace(",", "").split("-")
        out.append((chrom, int(s), int(e)))
    return out
