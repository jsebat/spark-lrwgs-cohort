"""Alignment surgery for the spike-in harness (DESIGN R3 / §5): plant an SNV, a deletion or an insertion into ONE
aligned read, as pure functions on the read's (reference_start, CIGAR, sequence, qualities). No pysam here, so the
edits are unit-tested on the Windows side; `sim/spike.py` wraps them around pysam.AlignedSegment.

Coordinates are 0-based reference positions. Edits split CIGAR operations at the edit boundaries instead of
expanding the alignment column by column (a HiFi read is ~15 kb; a few hundred ops is what we walk).

  apply_snv        substitute the base aligned to `pos` (read must be aligned through the position, else None)
  apply_deletion   delete reference bases [start, start+length): aligned query bases are removed, inserted bases
                   inside the interval are removed too, deletions inside merge; one D op results
  apply_insertion  insert `seq` into the read between the bases aligned to `pos` and `pos + 1`

Every function returns a new Aln or None when the read does not cover the edit with `margin` aligned reference
bases on both sides (the read is then left as it is — it becomes a REF/AMB observation naturally). `=`/`X` ops are
normalised to `M`; MD/NM tags are the wrapper's problem (it drops MD and recomputes nothing: NM stays approximate
and is documented as such in features.yaml).
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import List, Optional, Tuple

M, I, D, N, S, H, P, EQ, X = 0, 1, 2, 3, 4, 5, 6, 7, 8
_REF_OPS = {M, D, N, EQ, X}
_QRY_OPS = {M, I, S, EQ, X}

Cigar = List[Tuple[int, int]]


@dataclass
class Aln:
    reference_start: int
    cigar: Cigar
    seq: str
    qual: Optional[List[int]] = None

    @property
    def reference_end(self) -> int:
        return self.reference_start + sum(L for op, L in self.cigar if op in _REF_OPS)

    @property
    def query_length(self) -> int:
        return sum(L for op, L in self.cigar if op in _QRY_OPS)


def normalise(cig: Cigar) -> Cigar:
    """=/X -> M, then merge adjacent identical ops and drop zero-length ops."""
    out: Cigar = []
    for op, L in cig:
        if L <= 0:
            continue
        if op in (EQ, X):
            op = M
        if out and out[-1][0] == op:
            out[-1] = (op, out[-1][1] + L)
        else:
            out.append((op, L))
    return out


def query_index_at(a: Aln, pos: int) -> Optional[int]:
    """Query index of the base aligned to reference position `pos`, None if deleted / clipped / outside."""
    ref, q = a.reference_start, 0
    for op, L in a.cigar:
        if op in (M, EQ, X):
            if ref <= pos < ref + L:
                return q + (pos - ref)
            ref += L; q += L
        elif op in (D, N):
            if ref <= pos < ref + L:
                return None
            ref += L
        elif op in (I, S):
            q += L
        if ref > pos:
            return None
    return None


def covers(a: Aln, start: int, end: int, margin: int) -> bool:
    """Aligned (M) through [start - margin, start) and [end, end + margin)."""
    if margin <= 0:
        return a.reference_start <= start and end <= a.reference_end
    for p in (start - margin, start - 1, end, end + margin - 1):
        if query_index_at(a, p) is None:
            return False
    return True


def apply_snv(a: Aln, pos: int, alt: str, margin: int = 1) -> Optional[Aln]:
    qi = query_index_at(a, pos)
    if qi is None or not covers(a, pos, pos + 1, margin):
        return None
    if a.seq[qi].upper() == alt.upper():
        return None                               # already alt: nothing to plant (the caller counts it)
    seq = a.seq[:qi] + alt.upper() + a.seq[qi + 1:]
    return replace(a, cigar=normalise(a.cigar), seq=seq)


def apply_deletion(a: Aln, start: int, length: int, margin: int = 20) -> Optional[Aln]:
    end = start + length
    if not covers(a, start, end, margin):
        return None
    cig: Cigar = []
    seq_parts: List[str] = []
    qual_parts: List[List[int]] = []
    ref, q = a.reference_start, 0
    for op, L in normalise(a.cigar):
        if op == M:
            o_s, o_e = max(ref, start), min(ref + L, end)          # overlap with the deletion in ref space
            if o_s < o_e:
                left = o_s - ref
                right = (ref + L) - o_e
                if left:
                    cig.append((M, left)); seq_parts.append(a.seq[q:q + left])
                    if a.qual is not None: qual_parts.append(a.qual[q:q + left])
                cig.append((D, o_e - o_s))
                if right:
                    cig.append((M, right)); seq_parts.append(a.seq[q + L - right:q + L])
                    if a.qual is not None: qual_parts.append(a.qual[q + L - right:q + L])
            else:
                cig.append((op, L)); seq_parts.append(a.seq[q:q + L])
                if a.qual is not None: qual_parts.append(a.qual[q:q + L])
            ref += L; q += L
        elif op in (D, N):
            cig.append((D if (start <= ref < end or start < ref + L <= end) else op, L))
            ref += L
        elif op == I:
            if start < ref <= end:                                    # inserted bases inside the deleted interval vanish
                cig.append((D, 0))
            else:
                cig.append((op, L)); seq_parts.append(a.seq[q:q + L])
                if a.qual is not None: qual_parts.append(a.qual[q:q + L])
            q += L
        elif op == S:
            cig.append((op, L)); seq_parts.append(a.seq[q:q + L])
            if a.qual is not None: qual_parts.append(a.qual[q:q + L])
            q += L
        else:                                                          # H, P: consume nothing
            cig.append((op, L))
    return Aln(a.reference_start, normalise(cig), "".join(seq_parts),
               [x for part in qual_parts for x in part] if a.qual is not None else None)


def apply_insertion(a: Aln, pos: int, ins: str, margin: int = 20, ins_qual: int = 40) -> Optional[Aln]:
    """Insert `ins` between the base aligned to `pos` and the one aligned to `pos + 1`."""
    if not covers(a, pos + 1, pos + 1, margin):
        return None
    cig: Cigar = []
    ref, q = a.reference_start, 0
    seq, qual = a.seq, a.qual
    done = False
    for op, L in normalise(a.cigar):
        if not done and op == M and ref <= pos < ref + L:
            left = pos - ref + 1                                       # bases up to and including pos
            cig.append((M, left)); cig.append((I, len(ins))); cig.append((M, L - left))
            cut = q + left
            seq = seq[:cut] + ins.upper() + seq[cut:]
            if qual is not None:
                qual = qual[:cut] + [ins_qual] * len(ins) + qual[cut:]
            done = True
        else:
            cig.append((op, L))
        if op in _REF_OPS:
            ref += L
        if op in _QRY_OPS:
            q += L
    if not done:
        return None
    return Aln(a.reference_start, normalise(cig), seq, qual)


def cigar_string(cig: Cigar) -> str:
    return "".join("%d%s" % (L, "MIDNSHP=X"[op]) for op, L in cig)
