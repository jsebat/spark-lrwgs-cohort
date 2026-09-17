"""Alignment surgery for spike-ins: each edit must (i) keep the read internally consistent (query length == bases
consumed by the CIGAR, reference span moves by exactly the planted length) and (ii) read back through the SAME
functions the M2 adapters use (a CIGAR walk), so a planted variant is seen the way a real one is."""
import random

from phase_dnm.sim import edit as E

REF = "ACGTTGCAACGTAGCTAGCTTAGGCTAAGCTTGCA" * 20          # 700 bp of "reference"


def read_from_ref(start, length, clip=0, indel=None):
    """A perfectly aligned read over REF[start:start+length]; optional soft clip and an existing indel (op, at, L)."""
    seq = REF[start:start + length]
    cig = [(E.M, length)]
    if indel:
        op, at, L = indel
        if op == E.D:
            seq = seq[:at] + seq[at + L:]
            cig = [(E.M, at), (E.D, L), (E.M, length - at - L)]
        else:
            seq = seq[:at] + "T" * L + seq[at:]
            cig = [(E.M, at), (E.I, L), (E.M, length - at)]
    if clip:
        seq = "N" * clip + seq
        cig = [(E.S, clip)] + cig
    return E.Aln(start, cig, seq, [30] * len(seq))


def consistent(a):
    assert a.query_length == len(a.seq), (a.query_length, len(a.seq), E.cigar_string(a.cigar))
    if a.qual is not None:
        assert len(a.qual) == len(a.seq)
    return True


def base_at(a, pos):
    qi = E.query_index_at(a, pos)
    return None if qi is None else a.seq[qi]


def test_snv_substitution_reads_back_and_leaves_neighbours():
    a = read_from_ref(100, 300, clip=7)
    pos = 250
    ref_base = REF[pos]
    alt = "A" if ref_base != "A" else "C"
    b = E.apply_snv(a, pos, alt)
    assert consistent(b) and base_at(b, pos) == alt and base_at(a, pos) == ref_base
    assert base_at(b, pos - 1) == REF[pos - 1] and base_at(b, pos + 1) == REF[pos + 1]
    assert b.reference_end == a.reference_end
    assert E.apply_snv(a, 100 - 1, "A") is None            # outside the read
    assert E.apply_snv(a, pos, ref_base) is None           # already the requested base


def test_snv_at_a_deleted_position_is_not_planted():
    a = read_from_ref(100, 300, indel=(E.D, 50, 4))        # ref 150..153 deleted in this read
    assert E.apply_snv(a, 151, "A") is None
    assert E.apply_snv(a, 160, "A" if REF[160] != "A" else "C") is not None


def test_deletion_inside_a_match_block():
    a = read_from_ref(100, 300)
    b = E.apply_deletion(a, 200, 25)
    assert consistent(b)
    assert b.cigar == [(E.M, 100), (E.D, 25), (E.M, 175)]
    assert b.reference_end == a.reference_end and len(b.seq) == len(a.seq) - 25
    # bases flanking the deletion are the reference bases flanking it
    assert base_at(b, 199) == REF[199] and base_at(b, 225) == REF[225] and base_at(b, 210) is None
    # the adapters' view: a D op of the right length starting at the breakpoint
    assert (E.D, 25) in b.cigar


def test_deletion_merges_with_existing_indels_and_respects_margins():
    # existing 4 bp deletion at read offset 50 (ref 150..153) and a 5 bp insertion after ref 220
    a = read_from_ref(100, 300, indel=(E.D, 50, 4))
    b = E.apply_deletion(a, 140, 30)                       # spans the existing deletion -> one merged D of 30
    assert consistent(b) and b.cigar == [(E.M, 40), (E.D, 30), (E.M, 230)]
    c = read_from_ref(100, 300, indel=(E.I, 120, 5))       # insertion between ref 219 and 220
    d = E.apply_deletion(c, 210, 20)                       # inserted bases fall inside the deleted interval: gone
    assert consistent(d) and d.cigar == [(E.M, 110), (E.D, 20), (E.M, 170)] and len(d.seq) == 300 - 20 + 5 - 5
    # not enough aligned sequence on the right: refused, read untouched
    assert E.apply_deletion(a, 380, 30, margin=20) is None
    assert E.apply_deletion(a, 385, 15, margin=5) is None      # end + margin runs past the read
    assert E.apply_deletion(a, 370, 15, margin=5) is not None


def test_insertion_reads_back_as_an_I_op_at_the_breakpoint():
    a = read_from_ref(100, 300, clip=3)
    ins = "GATTACA" * 10
    b = E.apply_insertion(a, 249, ins)                      # between ref 249 and 250
    assert consistent(b)
    assert b.cigar == [(E.S, 3), (E.M, 150), (E.I, 70), (E.M, 150)]
    assert b.reference_end == a.reference_end and len(b.seq) == len(a.seq) + 70
    qi = E.query_index_at(b, 249)
    assert b.seq[qi + 1:qi + 71] == ins and base_at(b, 250) == REF[250]
    assert b.qual[qi + 1] == 40 and b.qual[qi] == 30
    assert E.apply_insertion(a, 395, ins) is None           # margin violated near the read end


def test_random_edits_stay_consistent():
    rng = random.Random(3)
    for _ in range(300):
        start = rng.randrange(0, 200)
        L = rng.randrange(120, 400)
        kind = rng.choice([None, (E.D, rng.randrange(20, 80), rng.randrange(1, 10)), (E.I, rng.randrange(20, 80), rng.randrange(1, 10))])
        a = read_from_ref(start, L, clip=rng.choice([0, 0, 5]), indel=kind)
        pos = rng.randrange(start + 25, start + L - 60)
        for out in (E.apply_snv(a, pos, "A" if REF[pos] != "A" else "G"),
                    E.apply_deletion(a, pos, rng.randrange(1, 30)),
                    E.apply_insertion(a, pos, "ACGT" * rng.randrange(1, 15))):
            if out is not None:
                consistent(out)
                assert out.reference_start == a.reference_start
