"""The split-read junction signature (readers.sv_signature), pinned on the geometry that was mis-handled.

A 35 kb deletion at [S, E). A junction read from the deleted haplotype aligns either
  (a) to the LEFT of S, clipped at its 3' end: reference_end == S, reference_start = S - ~10 kb, SA at E; or
  (b) to the RIGHT of E, clipped at its 5' end: reference_start == E, SA at S.
The first version compared only reference_start with the breakpoints, so (a) never matched (review 2026-09-16, D4).
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from phase_dnm.evidence.readers import sv_signature  # noqa: E402

S, E = 1_000_000, 1_035_454
BPS = [S, E]
MATCH = [(0, 12_000)]                       # a plain 12 kb match, no indel


def sig(ref_start, ref_end, sa):
    return sv_signature(ref_start, ref_end, MATCH, sa, "chrX", "DEL", E - S, BPS, bp_window=300, len_tol=0.3, sa_tol=2000)


def test_left_flank_read_clipped_at_its_3prime_end_matches():
    # primary ends AT the left breakpoint; supplementary starts at the right breakpoint
    assert sig(S - 12_000, S, "chrX,%d,+,12000S3000M,60,0;" % E)


def test_right_flank_read_clipped_at_its_5prime_end_matches():
    assert sig(E, E + 12_000, "chrX,%d,+,9000M12000S,60,0;" % (S - 9_000))


def test_supplementary_that_ENDS_at_the_other_breakpoint_matches():
    # SA record gives its leftmost position; its reference END is what lands on the breakpoint here
    assert sig(E, E + 12_000, "chrX,%d,+,9000M,60,0;" % (S - 9_000))


def test_supplementary_far_from_the_other_breakpoint_does_not_match():
    assert not sig(S - 12_000, S, "chrX,%d,+,3000M,60,0;" % (E + 50_000))


def test_primary_not_ending_at_any_breakpoint_does_not_match():
    # a read wholly inside the left flank, 20 kb from S, with an SA somewhere near E: not a junction read
    assert not sig(S - 40_000, S - 28_000, "chrX,%d,+,3000M,60,0;" % E)


def test_other_contig_supplementary_ignored():
    assert not sig(S - 12_000, S, "chr7,%d,+,3000M,60,0;" % E)


def test_same_breakpoint_both_sides_requires_the_OTHER_breakpoint():
    # primary ends at S and the SA also lands at S: not a deletion junction (the old code accepted this)
    assert not sig(S - 12_000, S, "chrX,%d,+,3000M,60,0;" % (S + 100))


def test_cigar_deletion_of_the_event_length_at_a_breakpoint_matches():
    cig = [(0, 5000), (2, 5000), (0, 5000)]          # a 5 kb D at S
    assert sv_signature(S - 5000, S + 10_000, cig, None, "chr1", "DEL", 5000, [S, S + 5000], 300, 0.3, 2000)
    assert not sv_signature(S - 5000, S + 10_000, cig, None, "chr1", "DEL", 20_000, [S, S + 20_000], 300, 0.3, 2000)


if __name__ == "__main__":
    import inspect
    bad = 0
    for n, f in sorted(globals().items()):
        if n.startswith("test_") and inspect.isfunction(f):
            try:
                f(); print("ok  ", n)
            except AssertionError:
                bad += 1; print("FAIL", n)
    sys.exit(1 if bad else 0)
