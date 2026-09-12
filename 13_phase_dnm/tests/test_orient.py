import json
import os
import subprocess
import sys

from conftest import read_tsv
from phase_dnm.io.vcf import VcfReader, iter_trio
from phase_dnm.phasing import orient as O


def test_assign_rule():
    # child 0|1; father hom-ref, mother carries alt -> hap1 (0) paternal
    assert O.assign((0, 1), {0}, {0, 1}) == 1
    # child 1|0; father hom-ref, mother carries alt -> hap1 (1) maternal
    assert O.assign((1, 0), {0}, {0, 1}) == -1
    # both parents het -> uninformative
    assert O.assign((0, 1), {0, 1}, {0, 1}) == 0
    # father hom-alt, mother het: alt must be paternal
    assert O.assign((1, 0), {1}, {0, 1}) == 1
    # hemizygous father {1} (chrX daughter), mother hom-ref: hap2 paternal
    assert O.assign((0, 1), {1}, {0}) == -1
    # both parents hom-ref: de novo candidate -> Mendelian-inconsistent
    assert O.assign((0, 1), {0}, {0}) is None
    # multi-allelic: child 1|2, father {0,1}, mother {0,2}
    assert O.assign((1, 2), {0, 1}, {0, 2}) == 1


def _run_orient(minitrio, **kw):
    p = minitrio["paths"]
    trio = iter_trio(VcfReader(p["child"]), VcfReader(p["father"]), VcfReader(p["mother"]))
    params = O.OrientParams(**kw)
    return O.orient_child(trio, "2", params), params


def test_orientation_matches_truth_on_clean_blocks(minitrio):
    (blocks, stats), params = _run_orient(minitrio)
    truth = {(c, ps): (orient, sw) for c, ps, s, e, orient, sw in minitrio["truth"].child_blocks}
    assert blocks, "no blocks"
    n_checked = n_wrong = 0
    for b in blocks:
        t = truth.get((b.chrom, b.ps))
        assert t is not None, "block %s:%d not in truth" % (b.chrom, b.ps)
        true_orient, switch_pos = t
        if switch_pos == 0 and b.orientation != O.AMBIGUOUS:
            n_checked += 1
            if b.orientation != true_orient:
                n_wrong += 1
    assert n_checked >= 5
    assert n_wrong == 0, "%d of %d clean oriented blocks wrong" % (n_wrong, n_checked)


def test_switch_blocks_are_split_or_flagged_never_silently_oriented(minitrio):
    (blocks, stats), params = _run_orient(minitrio)
    truth = {(c, ps): (orient, sw) for c, ps, s, e, orient, sw in minitrio["truth"].child_blocks}
    flipped = {O.HAP1_PAT: O.HAP1_MAT, O.HAP1_MAT: O.HAP1_PAT}
    switched_ps = {k for k, (o, sw) in truth.items() if sw != 0}
    assert switched_ps, "fixture produced no switch-error blocks; raise p_switch_per_block"
    n_split = 0
    for key in switched_ps:
        rows = [b for b in blocks if (b.chrom, b.ps) == key]
        true_orient, true_switch = truth[key]
        if len(rows) == 1:
            b = rows[0]
            assert b.orientation == O.AMBIGUOUS or b.n_dissent >= 1, (key, b.vote_frac)
            continue
        n_split += 1
        assert len(rows) == 2 and [r.segment for r in rows] == [1, 2] and all(r.n_segments == 2 for r in rows)
        left, right = rows
        # each side is either oriented correctly for its side of the switch, or too short to orient
        for seg, expected in ((left, true_orient), (right, flipped[true_orient])):
            if seg.orientation == O.AMBIGUOUS:
                assert seg.reason == "LOW_SITES", (key, seg.segment, seg.reason)
            else:
                assert seg.reason == "SPLIT_AT_SWITCH" and seg.orientation == expected, (key, seg.segment)
        assert not (left.orientation == right.orientation == O.AMBIGUOUS), "split produced two unusable halves"
        # segments are contiguous (a read is never unlabelled): right starts at the first INFORMATIVE site after
        # the switch, which is at or after the first het site after it (the truth position)
        assert left.start < true_switch <= right.start == right.switch_pos == left.end + 1, \
            (key, left.start, true_switch, right.start)
        assert right.start <= right.end
    assert n_split >= 1, "no switched block was split; change-point search is not working"


def test_clean_blocks_never_get_a_wrongly_oriented_segment(minitrio):
    (blocks, stats), params = _run_orient(minitrio)
    truth = {(c, ps): (orient, sw) for c, ps, s, e, orient, sw in minitrio["truth"].child_blocks}
    for b in blocks:
        true_orient, sw = truth[(b.chrom, b.ps)]
        if sw == 0 and b.n_segments > 1:
            # two adjacent genotype errors at a block edge are indistinguishable from a switch; the only
            # acceptable outcome is a short AMBIGUOUS tail, never an oriented segment that contradicts truth
            assert b.orientation in (true_orient, O.AMBIGUOUS), ("spurious oriented segment", b.chrom, b.ps)
            if b.orientation == O.AMBIGUOUS:
                assert b.reason == "LOW_SITES"
    # and splitting can be switched off
    (blocks0, _), _ = _run_orient(minitrio, max_splits=0)
    assert all(b.n_segments == 1 for b in blocks0)


def test_small_block_tier():
    p = O.OrientParams(min_sites=20, min_frac=0.95, small_min_sites=10, small_min_frac=1.0)
    mk = lambda signs: O.BlockVotes(chrom="chr1", ps=1, start=1, end=1000 * len(signs),
                                    het_positions=[i * 1000 for i in range(len(signs))],
                                    votes=[(i * 1000, s) for i, s in enumerate(signs)])
    # 12 unanimous votes: oriented under the second tier
    r = O.decide(mk([-1] * 12), p)[0]
    assert r.orientation == O.HAP1_MAT and r.reason == "OK_SMALL_UNANIMOUS"
    # 12 votes with one dissenter: not unanimous -> LOW_SITES (never MIXED_VOTES below min_sites)
    r = O.decide(mk([-1] * 11 + [1]), p)[0]
    assert r.orientation == O.AMBIGUOUS and r.reason == "LOW_SITES"
    # 9 unanimous votes: below the second tier
    r = O.decide(mk([-1] * 9), p)[0]
    assert r.orientation == O.AMBIGUOUS and r.reason == "LOW_SITES"
    # 25 votes with 2 dissenters interior: full tier, frac 0.92 < 0.95 -> MIXED_VOTES (no clean split: interior pair)
    r = O.decide(mk([-1] * 12 + [1, 1] + [-1] * 11), p)
    assert len(r) == 1 and r[0].reason == "MIXED_VOTES"
    # tier disabled by setting small_min_sites == min_sites
    p0 = O.OrientParams(min_sites=20, small_min_sites=20)
    assert O.decide(mk([-1] * 12), p0)[0].reason == "LOW_SITES"


def test_segment_votes_rules():
    p = O.OrientParams(min_sites=5, split_min_sites=2, max_splits=3)
    P, M = 1, -1
    mk = lambda signs: [(i * 1000, s) for i, s in enumerate(signs)]
    # clean block: no cut
    assert O.segment_votes(mk([M] * 30), p) == []
    # a lone dissenter anywhere: no cut
    assert O.segment_votes(mk([M] * 10 + [P] + [M] * 10), p) == []
    assert O.segment_votes(mk([P] + [M] * 20), p) == []
    # a run of two at the head: one cut at index 2 (the case that motivated the rule)
    assert O.segment_votes(mk([P, P] + [M] * 58), p) == [2]
    # a switch in the middle
    assert O.segment_votes(mk([P] * 20 + [M] * 20), p) == [20]
    # switch and back: two cuts
    assert O.segment_votes(mk([P] * 20 + [M] * 6 + [P] * 20), p) == [20, 26]
    # max_splits respected: an interior run needs two cuts; with one allowed, no single cut beats the
    # unsplit majority (20 + 20 - penalty < 40), so the block is left whole with 6 dissenters
    p1 = O.OrientParams(min_sites=5, split_min_sites=2, max_splits=1)
    assert O.segment_votes(mk([P] * 20 + [M] * 6 + [P] * 20), p1) == []
    assert O.segment_votes(mk([P] * 20 + [M] * 20 + [P] * 5), p1) == [20]   # one cut, tail stays as dissent
    # disabled
    assert O.segment_votes(mk([P] * 20 + [M] * 20), O.OrientParams(max_splits=0)) == []


def test_planted_dnms_are_counted_as_mendelian_inconsistent_not_voted(minitrio):
    (blocks, stats), params = _run_orient(minitrio)
    total = sum(st.get("n_mendel_inconsistent", 0) for st in stats.values())
    # planted DNMs that landed phased in a block with GQ >= 20 in all three samples must be here;
    # simulated genotype errors add a few more, so >= not ==
    assert total >= 1
    # and they never contribute votes: every vote position is a truth site with a parent carrying the allele
    dnm_pos = {(c, p) for c, p, r, a, h in minitrio["truth"].dnms}
    for b in blocks:
        for pos in b.dissent_positions:
            assert (b.chrom, pos) not in dnm_pos


def test_summary_and_ambiguity_gate(minitrio):
    (blocks, stats), params = _run_orient(minitrio)
    s = O.summarise(stats, blocks, params)
    assert s["n_segments"] == len(blocks)
    assert s["n_blocks"] == len({(b.chrom, b.ps) for b in blocks})
    assert s["n_blocks_split"] >= 1
    assert s["frac_bp_ambiguous"] is not None and s["frac_bp_ambiguous"] < 0.5
    assert s["frac_het_ambiguous"] is not None and 0 <= s["frac_het_ambiguous"] < 0.5
    assert s["frac_het_mixed_votes"] is not None and s["frac_bp_mixed_votes"] is not None
    het_all = s["total"]["het_oriented"] + s["total"]["het_ambiguous"]
    assert het_all == s["total"]["n_het_phased"], "every phased het must be in exactly one segment"
    assert s["total"]["n_informative"] > 100


def test_cli_end_to_end(minitrio, tmp_path):
    p = minitrio["paths"]
    c, f, m = minitrio["ids"]
    out = tmp_path / "phase"
    env = dict(os.environ)
    env["PYTHONPATH"] = os.path.join(os.path.dirname(__file__), "..", "src")
    cmd = [sys.executable, "-m", "phase_dnm.cli", "orient", "--child", c, "--manifest", p["manifest"],
           "--child-vcf", p["child"], "--father-vcf", p["father"], "--mother-vcf", p["mother"],
           "--out-dir", str(out)]
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    assert r.returncode == 0, r.stderr
    rows = read_tsv(out / ("%s.orientation.tsv" % c))
    assert rows and set(O.ORIENTATION_COLUMNS) <= set(rows[0])
    summ = json.load(open(out / ("%s.orientation.summary.json" % c)))
    assert summ["child"] == c and summ["father"] == f and summ["mother"] == m
    assert os.path.exists(out / ("%s.orientation.dissent.tsv" % c))
    # planned sub-commands refuse cleanly
    r2 = subprocess.run([sys.executable, "-m", "phase_dnm.cli", "haplotag"], capture_output=True, text=True, env=env)
    assert r2.returncode != 0 and "not implemented" in (r2.stderr + r2.stdout)
