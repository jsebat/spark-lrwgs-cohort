"""Planting a LARGE deletion: no HiFi read spans a 35 kb event, so it cannot be planted by editing one CIGAR.
The deleted haplotype's reads disappear and the ones crossing a breakpoint survive as clipped junction reads."""
from phase_dnm.sim import edit as E
from phase_dnm.sim import spike as SP

M, I, D, S = 0, 1, 2, 4


def aln(start=1000, cigar=None, n=500):
    return E.Aln(reference_start=start, cigar=cigar or [(M, n)], seq="A" * n, qual=[30] * n)


def test_breakpoint_clip_keeps_the_read_and_moves_the_edge():
    a = aln()
    left = E.apply_breakpoint(a, 1300, "left")
    assert left.reference_start == 1000 and left.reference_end == 1300
    assert left.query_length == a.query_length          # the bases are still there, they are soft-clipped
    assert left.cigar[-1][0] == S
    right = E.apply_breakpoint(a, 1300, "right")
    assert right.reference_start == 1300 and right.reference_end == 1500
    assert right.query_length == a.query_length and right.cigar[0][0] == S


def test_breakpoint_clip_refuses_a_read_that_barely_touches():
    assert E.apply_breakpoint(aln(), 1005, "left") is None
    assert E.apply_breakpoint(aln(), 1495, "right") is None
    assert E.apply_breakpoint(aln(), 900, "left") is None            # outside the alignment entirely


def test_breakpoint_clip_survives_existing_clips_and_insertions():
    a = E.Aln(2000, [(S, 20), (M, 300), (I, 5), (M, 200), (S, 10)], "C" * 535, [30] * 535)
    left = E.apply_breakpoint(a, 2400, "left")
    assert left.reference_end == 2400 and left.query_length == a.query_length


def test_large_deletion_removes_inside_reads_and_clips_crossing_ones():
    row = SP.PlanRow(variant_id="bd1", chrom="chr1", pos=10000, variant_class="SV", subtype="BIGDEL", length=20000,
                     ref=".", alt=".", scenario="G", child_hap=1, child_frac=1.0, parent="F", parent_hap=1,
                     parent_frac=0.0, expected_poo="paternal", expected_class="germline_DNM_phased",
                     child_ps=0, parent_ps=0)
    import random
    rng = random.Random(0)
    inside = aln(start=15000, n=1000)                                 # wholly inside -> the read does not exist
    assert SP._edit_one(inside, row, rng) is SP.DROP
    crossing_left = aln(start=9000, n=2000)                           # 9000-11000 crosses the left breakpoint
    e = SP._edit_one(crossing_left, row, rng)
    assert e is not None and e is not SP.DROP and e.reference_end == 10000
    crossing_right = aln(start=29000, n=2000)                         # 29000-31000 crosses the right breakpoint
    e = SP._edit_one(crossing_right, row, rng)
    assert e is not None and e.reference_start == 30000
    outside = aln(start=40000, n=500)
    assert SP._edit_one(outside, row, rng) is None                    # untouched
    spanning = aln(start=5000, n=40000)                               # long enough to span: keeps a D operation
    e = SP._edit_one(spanning, row, rng)
    assert e is not None and any(op == D and L == 20000 for op, L in e.cigar)


def test_bigdel_sizes_cover_what_the_pedigree_swap_cannot_supply():
    assert "BIGDEL" in SP.SUBTYPES["SV"]
    assert min(SP.LENGTHS[("SV", "BIGDEL")]) >= 5000
    assert max(SP.LENGTHS[("SV", "BIGDEL")]) >= 50000


def test_write_back_moves_the_read_start_for_a_right_clip():
    """The BAM read must receive the new reference_start, not only the CIGAR: until 2026-09-16 a keep="right" junction
    read kept its original start with a leading soft clip and sat one clip-length left of the breakpoint."""

    class Read:  # the four attributes write_back touches, plus the tag interface it calls
        def __init__(self):
            self.reference_start, self.cigartuples, self.query_sequence, self.query_qualities = 1000, [(M, 500)], "A" * 500, [30] * 500
        def has_tag(self, t): return False
        def set_tag(self, t, v): pass

    r = Read()
    right = E.apply_breakpoint(aln(1000, n=500), 1300, "right")
    SP.write_back(r, right)
    assert r.reference_start == 1300, r.reference_start
    assert r.cigartuples[0][0] == S                                     # leading soft clip
    assert r.reference_start + sum(L for op, L in r.cigartuples if op in (0, 2, 3, 7, 8)) == 1500
