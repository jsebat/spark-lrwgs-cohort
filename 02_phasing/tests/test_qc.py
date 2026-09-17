import json
import os
import subprocess
import sys

from conftest import read_tsv
from trio_phase.phasing import qc as Q


def _run_m1_vcf_level(minitrio, out):
    p = minitrio["paths"]; c, f, m = minitrio["ids"]
    env = dict(os.environ); env["PYTHONPATH"] = os.path.join(os.path.dirname(__file__), "..", "src")
    base = [sys.executable, "-m", "trio_phase.cli"]
    for cmd in ("orient", "transmission"):
        r = subprocess.run(base + [cmd, "--child", c, "--manifest", p["manifest"], "--child-vcf", p["child"],
                                   "--father-vcf", p["father"], "--mother-vcf", p["mother"], "--out-dir", str(out)],
                           capture_output=True, text=True, env=env)
        assert r.returncode == 0, r.stderr
    return env, base


def test_child_qc_tolerates_missing_read_level_inputs(minitrio, tmp_path):
    c, f, m = minitrio["ids"]
    fam_dir = tmp_path / "minifam"; fam_dir.mkdir()
    _run_m1_vcf_level(minitrio, fam_dir)
    gates = dict(max_frac_het_ambiguous=0.10, max_frac_het_mixed_votes=0.005, max_mendel_inconsistent_per_informative=0.05,
                 min_frac_het_transmission_resolved=0.5, crossovers_per_meiosis=[0, 999], min_total_depth=12)
    q = Q.child_qc(str(fam_dir), c, f, m, None, gates, "test")
    assert q["orientation"]["n_blocks"] > 0
    assert set(q["transmission"]) == {"F", "M"}
    assert "changepoints_resolved" in q["missing"] and "hapdepth" in q["missing"]
    assert "depth" not in q or not q["depth"]
    # gates only judge what is present: no crossover or depth flags can fire without those inputs
    assert not any(fl.startswith(("CROSSOVERS", "LOW_DEPTH", "LOW_TAGGED")) for fl in q["flags"])


def test_gates_fire_on_bad_values():
    q = {"orientation": dict(sex="M", frac_het_ambiguous=0.2, frac_het_mixed_votes=0.02, mendel_inconsistent_per_informative=0.05),
         "transmission": {"F": dict(frac_het_resolved=0.5)},
         "change_points": {"M": dict(candidates=100, CROSSOVER=120, crossover_excl_child_switch=110)},
         "depth": {"child": dict(dp_total=8.0, tagged_frac=0.6, hap_balance=0.5, chrx_autosome_ratio=1.0)}}
    g = dict(max_frac_het_ambiguous=0.10, max_frac_het_mixed_votes=0.005, max_mendel_inconsistent_per_informative=0.01,
             min_frac_het_transmission_resolved=0.90, crossovers_per_meiosis=[15, 65], min_total_depth=12,
             min_tagged_frac=0.75, min_hap_balance=0.9)
    flags = Q.apply_gates(q, g)
    for expected in ("ORIENT_AMBIGUOUS_HET", "ORIENT_MIXED_VOTES", "MENDEL", "TRANSMISSION_UNRESOLVED_F", "CROSSOVERS_M=110",
                     "LOW_DEPTH_child=8.0", "LOW_TAGGED_child", "HAP_IMBALANCE_child", "SEX_DEPTH_MISMATCH"):
        assert expected in flags, (expected, flags)
    # a female with a female-like ratio and everything else fine passes
    ok = {"orientation": dict(sex="F", frac_het_ambiguous=0.03, frac_het_mixed_votes=0.001, mendel_inconsistent_per_informative=0.001),
          "depth": {"child": dict(dp_total=23.0, tagged_frac=0.88, hap_balance=1.0, chrx_autosome_ratio=1.05)}}
    assert Q.apply_gates(ok, g) == []


def test_phase_qc_cli_cohort_table(minitrio, tmp_path):
    c, f, m = minitrio["ids"]
    phase_dir = tmp_path / "phase"; fam_dir = phase_dir / "minifam"; fam_dir.mkdir(parents=True)
    env, base = _run_m1_vcf_level(minitrio, fam_dir)
    r = subprocess.run(base + ["phase-qc", "--manifest", minitrio["paths"]["manifest"], "--phase-dir", str(phase_dir)],
                       capture_output=True, text=True, env=env)
    assert r.returncode == 0, r.stderr
    rows = read_tsv(phase_dir / "cohort_phase_qc.tsv")
    assert len(rows) == 1 and rows[0]["child"] == c and rows[0]["family"] == "minifam"
    assert set(Q.COHORT_COLUMNS) <= set(rows[0])
    q = json.load(open(fam_dir / ("%s.phase_qc.json" % c)))
    assert q["qc"] in ("PASS",) or all(not fl.startswith(("LOW_DEPTH", "CROSSOVERS")) for fl in q["flags"])
    assert os.path.exists(phase_dir / "cohort_phase_qc.summary.json")
