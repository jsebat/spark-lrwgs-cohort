"""M2 read-level layer (pysam): turn one CandidateRecord into ReadObs for the three samples. The ONLY class-specific
code in Module 2 lives here, in the three `support_*` hooks that answer "does this read support the alt allele?"
(README §2.2, DESIGN P17). Everything downstream (hapmatrix.py) is class-agnostic.

  SNV      base at the position from the read's aligned pairs; ALT if it equals the alt base, REF if the ref base,
           AMB otherwise (sequencing error, deletion at the site).
  INDEL    the read's CIGAR at the site: an insertion of the alt length starting within `indel_tol` bp of the
           position, or a deletion of the ref length; REF when the read is aligned through the site with no
           indel within the tolerance; AMB otherwise (a different indel length, or a soft-clip at the site).
  SV       ALT when the read name is in sawfish's supporting-read list for this SV id (any sample), or when it
           carries a junction signature at a breakpoint (SA tag, or a CIGAR deletion/insertion >= min_junction_bp
           within `bp_window` of a breakpoint - the rule of 06_clinical/sv_read_review.sh); REF when it spans a
           breakpoint window without a signature; AMB when it neither spans nor signals. For DEL/DUP the per-
           haplotype depth inside vs. flanks comes from hapdepth bins, not from this hook.
  TR       the per-read `AL:i` tag from TRGT's spanning-read BAM: ALT when within `tr_al_tol` bp (or one motif
           unit) of the outlier allele length, REF when within tolerance of another called allele, AMB otherwise.

Read names are hashed (blake2b, 8 bytes, salted) before they leave; identifiers are values, never code.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
from typing import Dict, Iterable, List, Optional, Set, Tuple

from ..records import CandidateRecord
from .hapmatrix import ReadObs

_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def rid_hash(name: str, salt: str) -> str:
    return hashlib.blake2b((salt + name).encode(), digest_size=8).hexdigest()


def _nm_rate(read) -> Optional[float]:
    try:
        nm = read.get_tag("NM")
    except KeyError:
        return None
    L = read.query_alignment_length or 0
    return nm / L if L else None


def _base_at(read, pos1: int) -> Optional[str]:
    """Base of the read aligned to 1-based reference position pos1, None if deleted/unaligned there.

    Walks the CIGAR (a few hundred ops on a HiFi read) instead of materialising all ~15k aligned pairs; the
    first version did the latter and cost ~0.45 s per candidate."""
    seq = read.query_sequence
    if seq is None or read.cigartuples is None:
        return None
    target = pos1 - 1
    ref, q = read.reference_start, 0
    for op, L in read.cigartuples:
        if op in (0, 7, 8):                          # M, =, X: consume both
            if ref <= target < ref + L:
                return seq[q + (target - ref)]
            ref += L; q += L
        elif op == 2 or op == 3:                    # D, N: consume reference only
            if ref <= target < ref + L:
                return None
            ref += L
        elif op == 1 or op == 4:                    # I, S: consume query only
            q += L
        # H (5) and P (6): consume nothing
        if ref > target:
            break
    return None


def _indel_at(read, pos1: int, tol: int) -> Tuple[Optional[str], int]:
    """('I'|'D', length) for an indel whose reference anchor lies within tol bp of pos1 (VCF anchor base), else (None, 0)."""
    if read.cigartuples is None:
        return None, 0
    ref = read.reference_start                      # 0-based
    for op, L in read.cigartuples:
        if op in (0, 7, 8):
            ref += L
        elif op == 2:                               # D: starts right after the anchor base
            if abs(ref - pos1) <= tol:
                return "D", L
            ref += L
        elif op == 1:                               # I
            if abs(ref - pos1) <= tol:
                return "I", L
        elif op == 3:
            ref += L
        if ref > pos1 + tol:
            break
    return None, 0


def _spans(read, a: int, b: int) -> bool:
    return read.reference_start <= a - 1 and read.reference_end is not None and read.reference_end >= b


# ----------------------------------------------------------------------------------------------
# support hooks
# ----------------------------------------------------------------------------------------------
def support_snv(read, rec: CandidateRecord) -> str:
    b = _base_at(read, rec.start)
    if b is None:
        return "AMB"
    return "ALT" if b == rec.alt else ("REF" if b == rec.ref else "AMB")


def support_indel(read, rec: CandidateRecord, tol: int = 3) -> str:
    op, L = _indel_at(read, rec.start, tol)
    if len(rec.alt) > len(rec.ref):                 # insertion of len(alt)-len(ref)
        want = ("I", len(rec.alt) - len(rec.ref))
    else:                                           # deletion of len(ref)-len(alt)
        want = ("D", len(rec.ref) - len(rec.alt))
    if op is None:
        return "REF" if _spans(read, rec.start - tol, rec.end + tol) else "AMB"
    return "ALT" if (op, L) == want else "AMB"


def support_sv(read, rec: CandidateRecord, supporting: Set[str], bp_window: int = 300, len_tol: float = 0.3,
               sa_tol: int = 2000) -> str:
    """ALT: named by sawfish as supporting, or carries a junction signature MATCHING THIS EVENT (an indel of the
    right length starting within bp_window of the breakpoint, or a supplementary alignment landing near the other
    breakpoint). REF: spans a breakpoint window with no matching signature. AMB: otherwise.
    The first version accepted ANY SA tag or large indel near a breakpoint and called 448 of 627 SVs on the first
    family 'phase_conflict_artifact' - parents' unrelated split reads in repeats look like alt on two haplotypes."""
    if read.query_name in supporting:
        return "ALT"
    pl = rec.class_payload
    svtype, svlen = pl.get("svtype"), pl.get("svlen")
    bps = [rec.start] + ([rec.end] if rec.end and rec.end != rec.start else [])
    near_bp = any(_spans(read, bp - bp_window, bp + bp_window) for bp in bps)
    touches = any(read.reference_start <= bp <= (read.reference_end or 0) for bp in bps)
    if not (near_bp or touches):
        return "AMB"
    if _sv_signature_matches(read, rec, svtype, svlen, bps, bp_window, len_tol, sa_tol):
        return "ALT"
    return "REF" if near_bp else "AMB"


def _sv_signature_matches(read, rec: CandidateRecord, svtype, svlen, bps, bp_window: int, len_tol: float, sa_tol: int) -> bool:
    return sv_signature(read.reference_start, read.reference_end, read.cigartuples,
                        read.get_tag("SA") if read.has_tag("SA") else None, rec.chrom, svtype, svlen, bps,
                        bp_window, len_tol, sa_tol)


def _cigar_ref_span(cigar: str) -> int:
    """Reference bases consumed by a CIGAR string (M/D/N/=/X), for the span of a supplementary alignment."""
    n = 0
    for L, op in _CIGAR.findall(cigar):
        if op in "MDN=X":
            n += int(L)
    return n


def sv_signature(ref_start: int, ref_end: Optional[int], cigartuples, sa_tag: Optional[str], chrom: str, svtype, svlen,
                 bps, bp_window: int, len_tol: float, sa_tol: int) -> bool:
    """Does this alignment carry a junction signature for THIS event? Pure function of the alignment's coordinates,
    CIGAR and SA tag, so it can be tested without a BAM.

    Two signatures. (1) A CIGAR indel of the event's length whose reference anchor lies within bp_window of a
    breakpoint. (2) A split read: this alignment ENDS (either end) within 3 x bp_window of one breakpoint and its
    supplementary alignment lands (either end) within sa_tol of the OTHER breakpoint.

    The first version of (2) compared only reference_start with the breakpoints, so a junction read whose primary
    segment lies to the LEFT of a breakpoint -- clipped at its 3' end, reference_end == breakpoint, start up to a
    read length away -- could never match; and an inner `any(... for bp in bps)` re-bound `bp`, so the "other
    breakpoint" constraint was a no-op. Every junction read not named in sawfish's supporting-read list went through
    this path: all of them in the spike arm and in the surrogate parents of synthetic trios, and any real read sawfish
    omitted. The planted-deletion junction deficit above 20 kb recorded in DESIGN P31 is this function's signature
    (review 2026-09-16, D4)."""
    L_want = abs(int(svlen)) if svlen not in (None, "", ".") else None
    # 1. CIGAR indel at a breakpoint with the event's length (DEL -> D, INS/DUP -> I)
    if L_want and svtype in ("DEL", "INS", "DUP") and cigartuples:
        want_op = 2 if svtype == "DEL" else 1
        ref = ref_start
        for op, L in cigartuples:
            if op in (0, 7, 8):
                ref += L
            elif op == 2:
                if op == want_op and abs(L - L_want) <= max(20, len_tol * L_want) and any(abs(ref - bp) <= bp_window for bp in bps):
                    return True
                ref += L
            elif op == 1:
                if op == want_op and abs(L - L_want) <= max(20, len_tol * L_want) and any(abs(ref - bp) <= bp_window for bp in bps):
                    return True
            elif op == 3:
                ref += L
            if ref > max(bps) + bp_window:
                break
    # 2. split read: primary ends at one breakpoint, supplementary lands at the other
    if sa_tag:
        ends_here = [ref_start] + ([ref_end] if ref_end is not None else [])
        for part in str(sa_tag).split(";"):
            f = part.split(",")
            if len(f) < 2:
                continue
            sa_chrom, sa_pos = f[0], int(f[1])
            if sa_chrom != chrom:
                continue
            sa_ends = [sa_pos] + ([sa_pos + _cigar_ref_span(f[3])] if len(f) > 3 and f[3] else [])
            for bp in bps:
                if min(abs(e - bp) for e in ends_here) > bp_window * 3:
                    continue                                  # this alignment does not end at this breakpoint
                others = [b for b in bps if b != bp] or [bp]
                if any(abs(se - o) <= sa_tol for se in sa_ends for o in others):
                    return True
    return False


def tr_read_length(read) -> Optional[int]:
    """Repeat length carried by one TRGT spanning read.

    TRGT trims each spanning read to the repeat plus two flanks whose lengths are in `FL:B:I`, so the repeat
    length is query_length - FL[0] - FL[1]. Verified on 4,385 reads at 300 loci (2026-09-12): it equals the
    sample's called AL of the read's assigned allele within +-2 bp in 98.4% of reads (median difference 0).
    Semantics of the other tags, measured the same way: `AL:i` is the ALLELE INDEX the read was assigned to
    within its own sample's call (0/1), not a length; `SO`/`EO` are the read's offsets relative to the repeat
    (SO is negative), so EO-SO is the read span, not the repeat; sum(MC_i * motif_i) is off by hundreds of bp
    on multi-motif loci."""
    if read.has_tag("FL") and read.query_length:
        fl = list(read.get_tag("FL"))
        L = read.query_length - sum(int(x) for x in fl)
        return L if L >= 0 else None
    return None


def support_tr(read, rec: CandidateRecord, tol_bp: int, role: str = "C") -> str:
    """Child: ALT if the read's repeat length is within tol of the outlier allele (the read's own AL index must
    agree when present); REF if within tol of another child allele. Parents: ALT if within tol of the CHILD's
    outlier allele (a parental read carrying the child's allele = missed inheritance or parental mosaic), REF if
    within tol of one of the parent's own called alleles; AMB otherwise."""
    L = tr_read_length(read)
    if L is None:
        return "AMB"
    pl = rec.class_payload
    idx = pl["outlier_allele_idx"]
    target = pl["child_AL"][idx]
    others = [a for i, a in enumerate(pl["child_AL"]) if i != idx]
    own = pl.get("father_AL" if role == "F" else "mother_AL") or []
    # The tolerance can never exceed half the distance to the nearest competing allele, otherwise stutter reads of
    # the child's OTHER allele count as alt (smoke run #2: 56% of TR rows had 'alt on both child haplotypes' for
    # the 1-2-unit candidates that dominate the raw set). Competing alleles: the child's other allele(s) for the
    # child, the parent's own alleles for a parent.
    compet = others if role == "C" else own
    if compet:
        gap = min(abs(target - o) for o in compet)
        tol_bp = min(tol_bp, max(0, (gap - 1) // 2))
    near_target = abs(L - target) <= tol_bp
    near_compet = any(abs(L - o) <= tol_bp for o in compet)
    if near_target and near_compet:
        return "AMB"
    if near_target:
        return "ALT"
    if near_compet:
        return "REF"
    return "AMB"


# ----------------------------------------------------------------------------------------------
# per-candidate extraction
# ----------------------------------------------------------------------------------------------
class TrioBams:
    """Open BAMs per role; for TR the TRGT spanning-read BAMs are used instead of the genome BAMs."""

    def __init__(self, bams: Dict[str, Tuple[str, Optional[str]]], tr_bams: Optional[Dict[str, Tuple[str, Optional[str]]]] = None):
        import pysam
        op = lambda p, i: pysam.AlignmentFile(p, "rb", index_filename=i) if i else pysam.AlignmentFile(p, "rb")
        self.bam = {role: op(*pi) for role, pi in bams.items()}
        self.tr = {role: op(*pi) for role, pi in (tr_bams or {}).items()}

    def close(self):
        for b in list(self.bam.values()) + list(self.tr.values()):
            b.close()


def load_supporting_reads(path: Optional[str], sv_id: str) -> Set[str]:
    """sawfish supporting_reads.json.gz: {sv_id: {sample: [read names]}} -> all names for this id."""
    if not path or not os.path.exists(path):
        return set()
    if not hasattr(load_supporting_reads, "_cache") or load_supporting_reads._cache[0] != path:   # type: ignore[attr-defined]
        with gzip.open(path, "rt") as fh:
            load_supporting_reads._cache = (path, json.load(fh))                                  # type: ignore[attr-defined]
    d = load_supporting_reads._cache[1].get(sv_id, {})                                           # type: ignore[attr-defined]
    return {n for names in d.values() for n in names}


def read_obs(rec: CandidateRecord, bams: TrioBams, salt: str, supporting_json: Optional[str] = None,
             pad: int = 0, min_mapq_keep: int = 0, tr_tol_bp: Optional[int] = None) -> List[ReadObs]:
    """ReadObs for every primary read of the three samples spanning the candidate (SV: touching a breakpoint window)."""
    out: List[ReadObs] = []
    cls = rec.variant_class
    if cls == "TR":
        pl = rec.class_payload
        tol = tr_tol_bp if tr_tol_bp is not None else max(2, pl.get("motif_unit_bp", 1))
        trid = pl.get("trid")
        for role, bam in bams.tr.items():
            for read in bam.fetch(rec.chrom, max(0, rec.start - 1), rec.end + 1):
                if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate:
                    continue
                if trid and read.has_tag("TR") and read.get_tag("TR") != trid:
                    continue                                    # a read of an overlapping/neighbouring locus
                out.append(_obs(read, role, support_tr(read, rec, tol, role), salt, al=tr_read_length(read)))
        return out
    supporting = load_supporting_reads(supporting_json, rec.class_payload.get("caller_id", rec.variant_id)) if cls == "SV" else set()
    if cls == "SV":
        windows = [(rec.start - 300, rec.start + 300)] + ([(rec.end - 300, rec.end + 300)] if rec.end != rec.start else [])
    else:
        windows = [(rec.start - pad, rec.end + pad)]
    for role, bam in bams.bam.items():
        seen: Set[str] = set()
        for a, b in windows:
            for read in bam.fetch(rec.chrom, max(0, a - 1), b):
                if read.is_unmapped or read.is_secondary or read.is_duplicate or read.query_name in seen:
                    continue
                if cls != "SV" and read.is_supplementary:
                    continue
                seen.add(read.query_name)
                if cls == "SNV":
                    if not _spans(read, rec.start, rec.start):
                        continue
                    sup = support_snv(read, rec)
                elif cls == "INDEL":
                    if not _spans(read, rec.start - 3, rec.end + 3):
                        continue
                    sup = support_indel(read, rec)
                else:
                    sup = support_sv(read, rec, supporting)
                out.append(_obs(read, role, sup, salt))
    return out


# ----------------------------------------------------------------------------------------------
# SV interval evidence (P17, implemented 2026-09-14 — it had been specified, declared in the registry and left unbuilt).
#
# A joint-caller genotype on a handful of junction reads does not establish a constitutional heterozygous SV, and the
# six-haplotype matrix at a breakpoint is a point-variant instrument: a 35 kb heterozygous deletion does not present as
# alt reads at a position, it presents as DEPTH LOST ACROSS THE INTERVAL ON ONE HAPLOTYPE. This is the read-level review
# of 06_clinical/sv_read_review.sh, done per haplotype for all six haplotypes of the trio:
#   (1) reads spanning each breakpoint that carry a junction signature, by haplotype;
#   (2) depth inside the interval vs the flanks, by haplotype, for the child and both parents;
#   (3) heterozygous-SNV persistence inside vs the flanks (the child's phased small-variant VCF), which separates a
#       constitutional deletion (hets vanish on the deleted haplotype) from mosaicism or chimeric reads (hets persist).
# Parental depth across the interval is what the original single-sample script could not see and is the point of doing
# it here: a parent that also loses depth carries the event, and the candidate is inherited, not de novo.
# ----------------------------------------------------------------------------------------------
SV_INTERVAL_TYPES = ("DEL", "DUP", "CNV", "INV")


def _hp_of(read) -> Optional[int]:
    hp = read.get_tag("HP") if read.has_tag("HP") else None
    return hp if hp in (1, 2) else None


def _depth_by_hap(bam, chrom: str, a: int, b: int, min_mapq: int = 5, n_points: int = 5) -> Dict[Optional[int], float]:
    """MEAN PER-BASE DEPTH over [a, b), by haplotype tag (1, 2, None = untagged).

    Depth must be measured at points, not as reads-per-window: a 15 kb HiFi read overlaps every short window it touches,
    so counting overlapping reads and dividing by the window width inflates short windows (it returned ~26x for a 302 bp
    event). Depth is therefore averaged over `n_points` evenly spaced positions.

    A point counts only where the read is ALIGNED there, which is not the same as lying between its start and end: a
    read carrying the deletion as a single D operation still spans the interval, and counting the span reported full
    depth inside a deletion that the reads plainly showed. Measured on a planted 5 kb heterozygous deletion, this
    function returned 0.968 and 0.938 for the two haplotypes where samtools depth over the same interval gives 0.397 --
    so the depth rule was blind to exactly the size class an aligner represents as one D in a spanning read, roughly
    300 bp up to the read length. Larger events were unaffected only because their reads go missing altogether, which
    is why the 35 kb MECP2 deletion scored correctly while mid-size ones could not."""
    out: Dict[Optional[int], float] = {1: 0.0, 2: 0.0, None: 0.0}
    if b <= a:
        return out
    pts = [a + int((i + 0.5) * (b - a) / n_points) for i in range(max(1, n_points))]
    for pt in pts:
        try:
            it = bam.fetch(chrom, max(0, pt), pt + 1)
        except (ValueError, KeyError):
            continue
        for read in it:
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate:
                continue
            if read.mapping_quality < min_mapq:
                continue
            if read.reference_start > pt or (read.reference_end or 0) <= pt:
                continue
            if not any(s0 <= pt < e0 for s0, e0 in read.get_blocks()):
                continue                      # spans the point inside a D/N: not aligned there, so not depth
            out[_hp_of(read)] += 1.0
    n = float(len(pts))
    return {k: v / n for k, v in out.items()}


def _junction_by_hap(bam, rec: CandidateRecord, pos: int, bp_window: int, supporting: Set[str]) -> Dict[Optional[int], int]:
    """Reads at one breakpoint carrying a junction signature for THIS event, tallied by haplotype."""
    out: Dict[Optional[int], int] = {1: 0, 2: 0, None: 0}
    try:
        it = bam.fetch(rec.chrom, max(0, pos - bp_window - 1), pos + bp_window)
    except (ValueError, KeyError):
        return out
    seen: Set[str] = set()
    for read in it:
        if read.is_unmapped or read.is_secondary or read.is_duplicate or read.query_name in seen:
            continue
        seen.add(read.query_name)
        if read.query_name in supporting or support_sv(read, rec, supporting) == "ALT":
            out[_hp_of(read)] += 1
    return out


def sv_interval_evidence(rec: CandidateRecord, bams: TrioBams, supporting_json: Optional[str] = None,
                         flank_bp: int = 20000, bp_window: int = 300, probe_bp: int = 2000,
                         min_interval_bp: int = 300, min_flank_reads: int = 5) -> Dict[str, object]:
    """Per-haplotype depth inside the interval vs the flanks, and junction reads by haplotype, for all six haplotypes.

    Depth is measured in probe windows (up to three inside, one either side outside) and expressed per base, so the
    ratio is independent of how long the event is. Returns {} for events with no interval (insertions, breakends) and
    for intervals below `min_interval_bp`, where the breakpoint matrix already carries the evidence."""
    pl = rec.class_payload or {}
    svtype = str(pl.get("svtype") or "")
    start, end = int(rec.start), int(rec.end or rec.start)
    length = end - start
    if svtype not in SV_INTERVAL_TYPES or length < min_interval_bp:
        return {}
    supporting = load_supporting_reads(supporting_json, pl.get("caller_id", rec.variant_id))
    inner_lo, inner_hi = start + bp_window, end - bp_window
    if inner_hi - inner_lo < 200:
        inner_lo, inner_hi = start + length // 4, end - length // 4
    span = max(inner_hi - inner_lo, 1)
    n_probe = 1 if span <= probe_bp else min(3, max(1, span // probe_bp))
    inside = []
    for i in range(n_probe):
        c = inner_lo + int((i + 0.5) * span / n_probe)
        w = min(probe_bp, span) // 2
        inside.append((c - w, c + w))
    flanks = [(max(0, start - bp_window - probe_bp), max(0, start - bp_window)),
              (end + bp_window, end + bp_window + probe_bp)]
    out: Dict[str, object] = {"sv_interval_len": length, "sv_interval_probes": n_probe}
    for role, bam in bams.bam.items():
        ins: Dict[Optional[int], float] = {1: 0.0, 2: 0.0, None: 0.0}
        for a, b in inside:
            c = _depth_by_hap(bam, rec.chrom, a, b)
            for k in ins:
                ins[k] += c[k] / len(inside)
        fl: Dict[Optional[int], float] = {1: 0.0, 2: 0.0, None: 0.0}
        for a, b in flanks:
            c = _depth_by_hap(bam, rec.chrom, a, b)
            for k in fl:
                fl[k] += c[k] / len(flanks)
        for hp in (1, 2, None):
            tag = "u" if hp is None else str(hp)
            out["%s_sv_in_hap%s" % (role, tag)] = round(ins[hp], 2)
            out["%s_sv_fl_hap%s" % (role, tag)] = round(fl[hp], 2)
        tot_in, tot_fl = sum(ins.values()), sum(fl.values())
        out["%s_sv_dp_in" % role] = round(tot_in, 2)
        out["%s_sv_dp_fl" % role] = round(tot_fl, 2)
        for hp in (1, 2):
            out["%s_sv_ratio_hap%d" % (role, hp)] = round(ins[hp] / fl[hp], 4) if fl[hp] >= min_flank_reads else None
        out["%s_sv_ratio_all" % role] = round(tot_in / tot_fl, 4) if tot_fl >= min_flank_reads else None
        # Inside a heterozygous deletion one haplotype is absent and the remaining one has no heterozygous site to be
        # phased against, so HiPhase cannot tag reads there at all: haplotype-resolved depth is structurally undefined
        # exactly where the event is, and the LOSS of haplotype resolution is itself the evidence (JS 2026-09-14).
        tag_in = (ins[1] + ins[2]) / tot_in if tot_in else None
        tag_fl = (fl[1] + fl[2]) / tot_fl if tot_fl else None
        out["%s_sv_tagged_frac_in" % role] = round(tag_in, 4) if tag_in is not None else None
        out["%s_sv_tagged_frac_fl" % role] = round(tag_fl, 4) if tag_fl is not None else None
        out["%s_sv_tagged_loss" % role] = round(tag_fl - tag_in, 4) if (tag_in is not None and tag_fl is not None) else None
        j1 = _junction_by_hap(bam, rec, start, bp_window, supporting)
        j2 = _junction_by_hap(bam, rec, end, bp_window, supporting) if end != start else {1: 0, 2: 0, None: 0}
        for hp in (1, 2, None):
            tag = "u" if hp is None else str(hp)
            out["%s_sv_junc_hap%s" % (role, tag)] = j1[hp] + j2[hp]
        out["%s_sv_junc_both_ends" % role] = int(sum(j1.values()) >= 2 and sum(j2.values()) >= 2)
    return out


def vcf_index_for(vcf_path: str):
    """The tabix/CSI index for a VCF, including the WDL layout where the VCF is a symlink and its index lives in a
    SIBLING <name>_index/ directory rather than next to the file. Returns None when there is none, so a caller can say
    so instead of silently fetching nothing."""
    import os as _os
    for cand in (vcf_path + ".tbi", vcf_path + ".csi",
                 vcf_path.replace("/joint_small_variants_vcf/", "/joint_small_variants_vcf_index/") + ".tbi",
                 vcf_path.replace("/joint_small_variants_vcf/", "/joint_small_variants_vcf_index/") + ".csi"):
        if _os.path.exists(cand):
            return cand
    return None


def sv_het_persistence(rec: CandidateRecord, smallvar_vcf: Optional[str], child: str,
                       flank_bp: int = 20000, bp_window: int = 300, min_interval_bp: int = 5000,
                       min_sites: int = 20) -> Dict[str, object]:
    """Heterozygous-SNV persistence inside the interval vs the flanks, from the child's phased small-variant VCF.

    A constitutional heterozygous deletion removes one haplotype, so heterozygous sites inside it collapse to
    hemizygous calls and the het fraction falls towards zero. Het sites persisting at the flanking rate mean both
    haplotypes are present inside the interval, which excludes a constitutional deletion and leaves mosaicism or a
    chimeric read as the explanation (the rule of 06_clinical/sv_read_review.sh)."""
    pl = rec.class_payload or {}
    if str(pl.get("svtype") or "") not in ("DEL", "CNV") or not smallvar_vcf or not os.path.exists(smallvar_vcf):
        return {}
    start, end = int(rec.start), int(rec.end or rec.start)
    if end - start < min_interval_bp:
        return {}
    import pysam
    idx = vcf_index_for(smallvar_vcf)
    if idx is None:
        return {"sv_het_index_missing": 1}
    def het_frac(a: int, b: int):
        n = h = 0
        try:
            vf = pysam.VariantFile(smallvar_vcf, index_filename=idx)
        except (OSError, ValueError):
            return None, 0
        try:
            for v in vf.fetch(rec.chrom, max(0, a), b):
                s = v.samples.get(child)
                if s is None:
                    continue
                al = [x for x in (s.get("GT") or ()) if x is not None]
                if len(al) < 2:
                    continue
                n += 1
                h += int(al[0] != al[1])
        except (ValueError, KeyError):
            return None, 0
        finally:
            vf.close()
        return (h / n if n else None), n
    f_in, n_in = het_frac(start + bp_window, end - bp_window)
    f_l, n_l = het_frac(max(0, start - flank_bp), max(0, start - bp_window))
    f_r, n_r = het_frac(end + bp_window, end + flank_bp)
    fl = [x for x in (f_l, f_r) if x is not None]
    if not fl or (n_l + n_r) < min_sites:
        return {"sv_het_sites_inside": n_in, "sv_het_sites_flank": n_l + n_r}
    flank = sum(fl) / len(fl)
    if f_in is None or n_in == 0:
        # no called site inside at all: for a deletion that is loss of heterozygosity, i.e. the positive signal, and it
        # must not be reported as missing data (it was, until 2026-09-14). Only meaningful when the flanks are informative.
        return {"sv_het_sites_inside": n_in, "sv_het_sites_flank": n_l + n_r, "sv_het_frac_inside": 0.0,
                "sv_het_frac_flank": round(flank, 4), "sv_het_snv_persistence": 0.0 if flank > 0 else None,
                "sv_het_no_sites_inside": 1}
    return {"sv_het_sites_inside": n_in, "sv_het_sites_flank": n_l + n_r,
            "sv_het_frac_inside": round(f_in, 4), "sv_het_frac_flank": round(flank, 4),
            "sv_het_snv_persistence": round(f_in / flank, 4) if flank > 0 else None}


def _obs(read, role: str, support: str, salt: str, al: Optional[int] = None) -> ReadObs:
    hp = read.get_tag("HP") if read.has_tag("HP") else None
    ps = read.get_tag("PS") if read.has_tag("PS") else None
    clipped = bool(read.cigartuples) and (read.cigartuples[0][0] in (4, 5) or read.cigartuples[-1][0] in (4, 5))
    rq = read.get_tag("rq") if read.has_tag("rq") else None
    return ReadObs(role=role, hp=hp if hp in (1, 2) else None, ps=ps, support=support, mapq=read.mapping_quality,
                   nm_rate=_nm_rate(read), clipped=clipped, supplementary=read.is_supplementary,
                   read_len=read.query_length or 0, al=al, rid=rid_hash(read.query_name, salt), rq=float(rq) if rq is not None else None)
