"""The audit must pass on a well-formed $PHASE_DIR and fail, non-zero, on each invariant it declares.

A checker that has never been seen to fail is not evidence of anything, so every check below is
exercised in both directions on the same fixture."""
import os
import shutil
import subprocess
import sys

import pytest

from conftest import read_tsv
from make_minitrio_vcf import CHILD, FAMILY
from trio_phase import audit as A
from trio_phase.phasing import xo_reads as X


def _messages(phase_dir, manifest, **kw):
    lines = []
    rc = A.run(phase_dir, manifest, log=lines.append, **kw)
    return rc, "\n".join(lines)


@pytest.fixture(scope="module")
def phase_dir(minitrio, tmp_path_factory):
    """A $PHASE_DIR with every artefact a complete trio must have, produced by the real commands."""
    root = tmp_path_factory.mktemp("phase")
    out = root / FAMILY
    p = minitrio["paths"]
    env = dict(os.environ)
    env["PYTHONPATH"] = os.path.join(os.path.dirname(__file__), "..", "src")
    for cmd in ("orient", "transmission"):
        r = subprocess.run([sys.executable, "-m", "trio_phase.cli", cmd, "--child", CHILD,
                            "--manifest", p["manifest"], "--child-vcf", p["child"],
                            "--father-vcf", p["father"], "--mother-vcf", p["mother"],
                            "--out-dir", str(out)], capture_output=True, text=True, env=env)
        assert r.returncode == 0, r.stderr
    # xo-reads needs BAMs, which the minitrio does not carry; the resolved table is built here with the
    # writer's own column list so the schema check is still checking the real declaration.
    rows = read_tsv(str(out / (CHILD + ".changepoints.tsv")))
    assert rows, "the minitrio must produce change points for this test to mean anything"
    with open(out / (CHILD + ".changepoints.resolved.tsv"), "w") as fh:
        fh.write("\t".join(X.COLUMNS) + "\n")
        for i, r in enumerate(rows):
            r = dict(r)
            r["status"] = "CROSSOVER" if i == 0 else "AMBIGUOUS"
            r.update(n_parent_hets_in_interval=10, n_gaps=9, weakest_gap_start=r["left_pos"],
                     weakest_gap_end=r["right_pos"], weakest_gap_informative_reads=5,
                     weakest_gap_discordant=0, max_gap_disc_frac=0.0, min_spanning_reads_required=3,
                     child_switch_in_interval="N")
            fh.write("\t".join(str(r[c]) for c in X.COLUMNS) + "\n")
    return str(root), p["manifest"]


def test_clean_phase_dir_passes(phase_dir):
    root, manifest = phase_dir
    rc, msg = _messages(root, manifest)
    assert rc == 0, msg
    for check in ("trios", "nonempty", "schema", "vocabulary", "overlap", "summaries"):
        assert "OK    %-22s" % check in msg, (check, msg)
    assert not [l for l in msg.splitlines() if l.startswith("FAIL")], msg
    assert "audit: 0 FAIL" in msg


def test_missing_table_fails(phase_dir, tmp_path):
    root, manifest = phase_dir
    d = tmp_path / "missing"
    shutil.copytree(root, d)
    os.remove(d / FAMILY / (CHILD + ".transmission.tsv"))
    rc, msg = _messages(str(d), manifest)
    assert rc == 1 and "no .transmission.tsv" in msg


def test_empty_table_fails(phase_dir, tmp_path):
    root, manifest = phase_dir
    d = tmp_path / "empty"
    shutil.copytree(root, d)
    p = d / FAMILY / (CHILD + ".orientation.tsv")
    header = open(p).readline()
    open(p, "w").write(header)
    rc, msg = _messages(str(d), manifest)
    assert rc == 1 and "header and no rows" in msg


def test_unclassified_change_points_fail(phase_dir, tmp_path):
    """Every change point `transmission` emits must come back from the read-level step."""
    root, manifest = phase_dir
    d = tmp_path / "short"
    shutil.copytree(root, d)
    p = d / FAMILY / (CHILD + ".changepoints.resolved.tsv")
    lines = open(p).read().splitlines(True)
    open(p, "w").writelines(lines[:-1])
    rc, msg = _messages(str(d), manifest)
    assert rc == 1 and "did not classify every change point" in msg


def test_no_crossover_fails(phase_dir, tmp_path):
    root, manifest = phase_dir
    d = tmp_path / "nox"
    shutil.copytree(root, d)
    p = d / FAMILY / (CHILD + ".changepoints.resolved.tsv")
    text = open(p).read()
    open(p, "w").write(text.replace("CROSSOVER", "AMBIGUOUS"))
    rc, msg = _messages(str(d), manifest)
    assert rc == 1 and "no CROSSOVER at all" in msg


def test_unknown_status_fails(phase_dir, tmp_path):
    root, manifest = phase_dir
    d = tmp_path / "vocab"
    shutil.copytree(root, d)
    p = d / FAMILY / (CHILD + ".changepoints.resolved.tsv")
    text = open(p).read()
    open(p, "w").write(text.replace("AMBIGUOUS", "MAYBE_A_CROSSOVER", 1))
    rc, msg = _messages(str(d), manifest)
    assert rc == 1 and "OUTSIDE the declared vocabulary" in msg


def test_changed_header_fails(phase_dir, tmp_path):
    root, manifest = phase_dir
    d = tmp_path / "schema"
    shutil.copytree(root, d)
    p = d / FAMILY / (CHILD + ".transmission.tsv")
    lines = open(p).read().splitlines(True)
    lines[0] = lines[0].replace("transmitted", "transmitted_hap")
    open(p, "w").writelines(lines)
    rc, msg = _messages(str(d), manifest)
    assert rc == 1 and "differs from the writer's declaration" in msg


def test_overlapping_segments_within_a_block_fail(phase_dir, tmp_path):
    """Two rows of one phase block covering the same base would give that base two orientations."""
    root, manifest = phase_dir
    d = tmp_path / "overlap"
    shutil.copytree(root, d)
    p = d / FAMILY / (CHILD + ".orientation.tsv")
    lines = open(p).read().splitlines(True)
    open(p, "a").write(lines[1])                      # the same segment of the same block, twice
    rc, msg = _messages(str(d), manifest)
    assert rc == 1 and "overlap inside one phase block" in msg


def test_summary_without_thresholds_version_fails(phase_dir, tmp_path):
    root, manifest = phase_dir
    d = tmp_path / "vers"
    shutil.copytree(root, d)
    import json
    p = d / FAMILY / (CHILD + ".orientation.summary.json")
    j = json.load(open(p))
    j.pop("thresholds_version", None)
    json.dump(j, open(p, "w"))
    rc, msg = _messages(str(d), manifest)
    assert rc == 1 and "no thresholds_version" in msg


def test_missing_hapdepth_fails(phase_dir, tmp_path):
    root, manifest = phase_dir
    empty = tmp_path / "hapdepth_empty"
    empty.mkdir()
    rc, msg = _messages(root, manifest, hapdepth_dir=str(empty))
    assert rc == 1 and "no depth table" in msg


def test_cli_audit_exit_code(phase_dir):
    root, manifest = phase_dir
    env = dict(os.environ)
    env["PYTHONPATH"] = os.path.join(os.path.dirname(__file__), "..", "src")
    r = subprocess.run([sys.executable, "-m", "trio_phase.cli", "audit", "--phase-dir", root,
                        "--manifest", manifest], capture_output=True, text=True, env=env)
    assert r.returncode == 0, r.stderr
    assert "audit: 0 FAIL" in r.stderr
