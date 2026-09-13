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
           within `bp_window` of a breakpoint - the rule of 05_denovo/sv_read_review.sh); REF when it spans a
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
    L_want = abs(int(svlen)) if svlen not in (None, "", ".") else None
    # 1. CIGAR indel at a breakpoint with the event's length (DEL -> D, INS/DUP -> I)
    if L_want and svtype in ("DEL", "INS", "DUP") and read.cigartuples:
        want_op = 2 if svtype == "DEL" else 1
        ref = read.reference_start
        for op, L in read.cigartuples:
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
    # 2. supplementary alignment landing near the OTHER breakpoint (or, for BND/INV, near any breakpoint on the same contig)
    if read.has_tag("SA"):
        for part in read.get_tag("SA").split(";"):
            f = part.split(",")
            if len(f) < 2:
                continue
            chrom, pos = f[0], int(f[1])
            if chrom != rec.chrom:
                continue
            here = read.reference_start
            for bp in bps:
                other = [b for b in bps if b != bp] or [bp]
                if any(abs(here - bp) <= bp_window * 3 for bp in bps) and any(abs(pos - o) <= sa_tol for o in other):
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
    if role == "C":
        if abs(L - target) <= tol_bp:
            return "ALT"
        others = [a for i, a in enumerate(pl["child_AL"]) if i != idx]
        return "REF" if any(abs(L - o) <= tol_bp for o in others) else "AMB"
    own = pl.get("father_AL" if role == "F" else "mother_AL") or []
    if abs(L - target) <= tol_bp and not any(abs(L - o) <= tol_bp for o in own):
        return "ALT"
    if any(abs(L - o) <= tol_bp for o in own):
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


def _obs(read, role: str, support: str, salt: str, al: Optional[int] = None) -> ReadObs:
    hp = read.get_tag("HP") if read.has_tag("HP") else None
    ps = read.get_tag("PS") if read.has_tag("PS") else None
    clipped = bool(read.cigartuples) and (read.cigartuples[0][0] in (4, 5) or read.cigartuples[-1][0] in (4, 5))
    return ReadObs(role=role, hp=hp if hp in (1, 2) else None, ps=ps, support=support, mapq=read.mapping_quality,
                   nm_rate=_nm_rate(read), clipped=clipped, supplementary=read.is_supplementary,
                   read_len=read.query_length or 0, al=al, rid=rid_hash(read.query_name, salt))
