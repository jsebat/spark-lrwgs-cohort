"""apply_plan end to end on a tiny BAM: a planted large deletion clips reads at both breakpoints; the left-clipped ones
(crossing the RIGHT breakpoint) move to a later start (write_back, 8abfdca), so the written stream is no longer in
fetch order. The output must still be coordinate-sorted and indexable (v8 spike: 32/33 tasks died in pysam.index)."""
import os
import random

import pytest

pysam = pytest.importorskip("pysam")

from phase_dnm.sim import spike as SP

CHROM, CLEN, RL = "chr1", 400000, 3000
DEL_POS, DEL_LEN = 200001, 5000          # anchor; deleted bases 0-based [200001, 205001)


def _bam(path):
    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": CHROM, "LN": CLEN}]}
    rng = random.Random(1)
    reads = []
    # reads every 250 bp across the site on both haplotypes, so that plenty cross each breakpoint
    for i, start in enumerate(range(DEL_POS - 30000, DEL_POS + DEL_LEN + 30000, 250)):
        for hp in (1, 2):
            a = pysam.AlignedSegment()
            a.query_name = "r%d_h%d" % (i, hp)
            a.reference_id, a.reference_start, a.mapping_quality = 0, start, 60
            a.cigartuples = [(0, RL)]
            a.query_sequence = "".join(rng.choice("ACGT") for _ in range(RL))
            a.query_qualities = pysam.qualitystring_to_array("I" * RL)
            a.set_tag("HP", hp); a.set_tag("PS", 1)
            a.flag = 0 if rng.random() < 0.5 else 16
            reads.append(a)
    reads.sort(key=lambda r: r.reference_start)
    with pysam.AlignmentFile(path, "wb", header=hdr) as out:
        for r in reads:
            out.write(r)
    pysam.index(path)
    return path


def _row():
    return SP.PlanRow(variant_id="G_BIGDEL_1", chrom=CHROM, pos=DEL_POS, variant_class="SV", subtype="BIGDEL", length=DEL_LEN,
                      ref=".", alt=".", scenario="G", child_hap=1, child_frac=1.0, parent=".", parent_hap=0, parent_frac=0.0,
                      expected_poo="paternal", expected_class="child_dnm", child_ps=1, parent_ps=0)


def test_apply_plan_output_is_sorted_and_indexed_after_junction_reads_move(tmp_path):
    src = _bam(str(tmp_path / "C.bam"))
    with pysam.AlignmentFile(src) as bam:
        stats = SP.apply_plan([_row()], {"C": bam}, str(tmp_path), "spiked", SP.SiteQC(), seed=0)
    out = str(tmp_path / "C.spiked.bam")
    assert os.path.exists(out) and os.path.exists(out + ".bai")
    assert not os.path.exists(out + ".unsorted.bam")
    assert stats["edited_C"] > 0 and stats["dropped_C"] > 0
    prev, n, moved = -1, 0, 0
    with pysam.AlignmentFile(out) as bam:
        for r in bam.fetch(until_eof=True):
            assert r.reference_start >= prev, "not coordinate-sorted at %s" % r.query_name
            prev = r.reference_start; n += 1
            # a left-clipped junction read (crossing the right breakpoint) now starts AT the right breakpoint
            if r.has_tag("SA") and r.cigartuples[0][0] == 4:
                moved += 1
                assert r.reference_start == DEL_POS + DEL_LEN     # deleted bases are 0-based [pos, pos+L) (_edit_one)
    assert n > 0 and moved > 0
    # and the index is usable for a region query across the deletion
    with pysam.AlignmentFile(out) as bam:
        assert sum(1 for _ in bam.fetch(CHROM, DEL_POS - 100, DEL_POS + 100)) > 0
