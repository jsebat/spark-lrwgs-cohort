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
    """Base of the read aligned to 1-based reference position pos1, None if deleted/unaligned there."""
    seq = read.query_sequence
    if seq is None:
        return None
    for q, r in read.get_aligned_pairs(matches_only=True):
        if r == pos1 - 1:
            return seq[q]
        if r is not None and r > pos1 - 1:
            break
    return None


def _indel_at(read, pos1: int, tol: int) -> Tuple[Optional[str], int]:
    """('I'|'D', length) for an indel whose reference anchor lies within tol bp of pos1 (VCF anchor base), else (None, 0)."""
    ref = read.reference_start                      # 0-based
    for L, op in _CIGAR.findall(read.cigarstring or ""):
        L = int(L)
        if op in "M=X":
            ref += L
        elif op == "D":
            if abs(ref - (pos1)) <= tol:            # deletion starts right after the anchor base
                return "D", L
            ref += L
        elif op == "I":
            if abs(ref - (pos1)) <= tol:
                return "I", L
        elif op == "N":
            ref += L
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


def support_sv(read, rec: CandidateRecord, supporting: Set[str], min_junction_bp: int = 200, bp_window: int = 300) -> str:
    if read.query_name in supporting:
        return "ALT"
    bps = [rec.start] + ([rec.end] if rec.end and rec.end != rec.start else [])
    near_bp = any(_spans(read, bp - bp_window, bp + bp_window) for bp in bps)
    if not near_bp and not any(read.reference_start <= bp <= (read.reference_end or 0) for bp in bps):
        return "AMB"
    sig = read.has_tag("SA")
    if not sig:
        for L, op in _CIGAR.findall(read.cigarstring or ""):
            if op in "DI" and int(L) >= min_junction_bp:
                sig = True
                break
    if sig:
        return "ALT"
    return "REF" if near_bp else "AMB"


def support_tr(read, rec: CandidateRecord, tol_bp: int) -> str:
    if not read.has_tag("AL"):
        return "AMB"
    al = int(read.get_tag("AL"))
    pl = rec.class_payload
    target = pl["child_AL"][pl["outlier_allele_idx"]]
    others = [a for i, a in enumerate(pl["child_AL"]) if i != pl["outlier_allele_idx"]]
    if abs(al - target) <= tol_bp:
        return "ALT"
    if any(abs(al - o) <= tol_bp for o in others):
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
        for role, bam in bams.tr.items():
            for read in bam.fetch(rec.chrom, max(0, rec.start - 1), rec.end + 1):
                if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate:
                    continue
                out.append(_obs(read, role, support_tr(read, rec, tol), salt, al=int(read.get_tag("AL")) if read.has_tag("AL") else None))
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
