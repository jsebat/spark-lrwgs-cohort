import os
import subprocess
import sys

from conftest import read_tsv
from trio_phase.io.vcf import VcfReader, iter_trio
from trio_phase.phasing import orient as O
from trio_phase.phasing import transmission as T


def _orient_rows(minitrio):
    p = minitrio["paths"]
    trio = iter_trio(VcfReader(p["child"]), VcfReader(p["father"]), VcfReader(p["mother"]))
    blocks, _ = O.orient_child(trio, "2", O.OrientParams())
    return [dict(chrom=b.chrom, phase_block_id=b.ps, start=b.start, end=b.end, orientation=b.orientation) for b in blocks]


def _run(minitrio, **kw):
    p = minitrio["paths"]
    orientation = T.Orientation(_orient_rows(minitrio))
    trio = iter_trio(VcfReader(p["child"]), VcfReader(p["father"]), VcfReader(p["mother"]))
    params = T.TransmissionParams(**kw)
    return T.build_transmission(trio, orientation, "2", params), params


def _expected_label(truth_block, xo, chrom, pos):
    """Block-coordinate label (HAP1/HAP2) of the truly transmitted haplotype at pos, in a parent's block."""
    sample, c, ps, start, end, hap1_is_true, switch_pos = truth_block
    xpos, first, second = xo[(chrom, "F" if "father" in sample else "M")]
    true_hap = first if pos < xpos else second               # 1 or 2 in TRUE haplotype numbering
    label = true_hap if hap1_is_true == 1 else 3 - true_hap  # block hap1 == true hap1 when hap1_is_true == 1
    if switch_pos and pos >= switch_pos:
        label = 3 - label
    return "HAP1" if label == 1 else "HAP2"


def test_resolved_segments_match_truth(minitrio):
    (segments, changes, stats), params = _run(minitrio)
    truth = minitrio["truth"]
    xo = {(c, par): (pos, f, t) for c, par, pos, f, t in truth.crossovers}
    pblocks = {("F" if "father" in s else "M", c, ps): (s, c, ps, st, en, h, sw) for s, c, ps, st, en, h, sw in truth.parent_blocks}
    n_checked = n_wrong = n_straddle = 0
    for s in segments:
        if s.transmitted == T.UNRESOLVED:
            continue
        tb = pblocks[(s.parent, s.chrom, s.ps)]
        xpos = xo[(s.chrom, s.parent)][0]
        switch_pos = tb[6]
        # a segment that contains a TRUE change point (crossover or parental switch) can only be right on one side;
        # whether the change was located is the change-point test's job. Segments without one must match exactly.
        if s.start < xpos <= s.end or (switch_pos and s.start < switch_pos <= s.end):
            n_straddle += 1
            continue
        for pos in (s.start, s.end):
            n_checked += 1
            if _expected_label(tb, xo, s.chrom, pos) != s.transmitted:
                n_wrong += 1
    assert n_checked >= 10
    assert n_wrong == 0, "%d of %d segment ends disagree with truth (segments without an internal change)" % (n_wrong, n_checked)
    # straddling segments exist only when a change had < xo_min_sites votes on one side; they must be few
    assert n_straddle <= len(segments) // 4, "too many resolved segments straddle an undetected change: %d" % n_straddle


def test_change_points_are_crossovers_or_parent_switches(minitrio):
    (segments, changes, stats), params = _run(minitrio)
    truth = minitrio["truth"]
    xo_pos = {(c, par): pos for c, par, pos, f, t in truth.crossovers}
    switches = {("F" if "father" in s else "M", c, ps): sw for s, c, ps, st, en, h, sw in truth.parent_blocks if sw}
    assert changes, "no change points found; the fixture plants one crossover per parent per chromosome"
    unexplained = []
    for cp in changes:
        ok_xo = cp.left_pos < xo_pos[(cp.chrom, cp.parent)] <= cp.right_pos
        sw = switches.get((cp.parent, cp.chrom, cp.ps), 0)
        ok_sw = bool(sw) and cp.left_pos < sw <= cp.right_pos
        if not (ok_xo or ok_sw):
            unexplained.append((cp.chrom, cp.parent, cp.ps, cp.left_pos, cp.right_pos))
    assert not unexplained, unexplained
    # every change point is a candidate until the read-level step resolves it
    assert all(cp.status == "CANDIDATE" for cp in changes)


def test_crossovers_inside_blocks_are_recovered(minitrio):
    (segments, changes, stats), params = _run(minitrio)
    truth = minitrio["truth"]
    recovered = 0
    for c, par, pos, f, t in truth.crossovers:
        if any(cp.chrom == c and cp.parent == par and cp.left_pos < pos <= cp.right_pos for cp in changes):
            recovered += 1
    # crossovers between parent blocks are invisible by design (P3); require at least half recovered on this fixture
    assert recovered >= len(truth.crossovers) // 2, "recovered %d of %d planted crossovers" % (recovered, len(truth.crossovers))


def test_summary_and_cli(minitrio, tmp_path):
    (segments, changes, stats), params = _run(minitrio)
    s = T.summarise(stats, params)
    for parent in ("F", "M"):
        pp = s["per_parent"][parent]
        assert pp["n_informative"] > 100 and pp["frac_het_resolved"] is not None
        assert pp["mendel_inconsistent_per_informative"] < 0.05
    p = minitrio["paths"]; c, f, m = minitrio["ids"]
    out = tmp_path / "phase"
    env = dict(os.environ); env["PYTHONPATH"] = os.path.join(os.path.dirname(__file__), "..", "src")
    base = [sys.executable, "-m", "trio_phase.cli"]
    r = subprocess.run(base + ["orient", "--child", c, "--manifest", p["manifest"], "--child-vcf", p["child"],
                               "--father-vcf", p["father"], "--mother-vcf", p["mother"], "--out-dir", str(out)],
                       capture_output=True, text=True, env=env)
    assert r.returncode == 0, r.stderr
    r = subprocess.run(base + ["transmission", "--child", c, "--manifest", p["manifest"], "--child-vcf", p["child"],
                               "--father-vcf", p["father"], "--mother-vcf", p["mother"], "--out-dir", str(out)],
                       capture_output=True, text=True, env=env)
    assert r.returncode == 0, r.stderr
    seg = read_tsv(out / ("%s.transmission.tsv" % c))
    cps = read_tsv(out / ("%s.changepoints.tsv" % c))
    assert seg and set(T.SEGMENT_COLUMNS) <= set(seg[0])
    assert cps and set(T.CHANGE_COLUMNS) <= set(cps[0])
    assert os.path.exists(out / ("%s.transmission.summary.json" % c))
