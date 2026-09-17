"""Spike-in harness, pysam-free parts: the scenario grid, the per-read edit decision, plan round-trip, window QC
and the evaluator on a synthetic evidence table."""
import random

from phase_dnm.sim import edit as E
from phase_dnm.sim import spike as SP


def test_scenario_grid_covers_every_class_subtype_scenario():
    grid = SP._scenario_grid(random.Random(1), 3)
    combos = {(c, s, sc) for c, s, sc, *_ in grid}
    assert ("SNV", "SNV", "G") in combos and ("SV", "DEL", "PM") in combos and ("TR", "EXP", "IM") in combos and ("INDEL", "INS", "CM") in combos
    assert len(grid) == 3 * sum(len(v) for v in SP.SUBTYPES.values()) * len(SP.SCENARIOS)
    for c, s, sc, L, cf, pf in grid:
        assert (sc != "CM") == (cf == 1.0)
        assert (pf > 0) == (sc in ("PM", "IM"))
        assert pf == 1.0 if sc == "IM" else True


class _R:
    def __init__(self, hp):
        self._hp = hp
    def has_tag(self, t):
        return t == "HP" and self._hp is not None
    def get_tag(self, t):
        return self._hp


def _row(sc="G", child_hap=1, parent=".", parent_hap=0, pf=0.0, cf=1.0):
    return SP.PlanRow("v", "chr1", 100, "SNV", "SNV", 1, "A", "C", sc, child_hap, cf, parent, parent_hap, pf, "paternal",
                      SP.EXPECTED_CLASS[sc], 1, 1)


def test_edit_decision_follows_haplotype_and_fraction():
    rng = random.Random(0)
    g = _row("G", child_hap=2)
    assert all(SP._decide(_R(2), "C", g, rng) for _ in range(50))          # target haplotype: always
    assert not any(SP._decide(_R(1), "C", g, rng) for _ in range(50))      # other haplotype: never
    assert not any(SP._decide(_R(2), "F", g, rng) for _ in range(50))      # parents untouched in G
    unt = sum(SP._decide(_R(None), "C", g, rng) for _ in range(2000)) / 2000
    assert 0.42 < unt < 0.58                                                # untagged child reads: half
    pm = _row("PM", child_hap=1, parent="F", parent_hap=2, pf=0.2)
    f2 = sum(SP._decide(_R(2), "F", pm, rng) for _ in range(4000)) / 4000
    assert 0.16 < f2 < 0.24 and not any(SP._decide(_R(1), "F", pm, rng) for _ in range(50))
    assert not any(SP._decide(_R(2), "M", pm, rng) for _ in range(50))     # the other parent never
    cm = _row("CM", child_hap=1, cf=0.3)
    c1 = sum(SP._decide(_R(1), "C", cm, rng) for _ in range(4000)) / 4000
    assert 0.25 < c1 < 0.35


def test_plan_round_trip(tmp_path):
    rows = [_row("G"), _row("IM", parent="M", parent_hap=1, pf=1.0)]
    p = str(tmp_path / "plan.tsv")
    SP.write_plan(rows, p)
    back = SP.read_plan(p)
    assert back == rows


def test_clean_window_and_consensus():
    ref = "ACGT" * 100
    a = E.Aln(0, [(E.M, 400)], ref, None)
    b = E.Aln(0, [(E.M, 100), (E.D, 5), (E.M, 295)], ref[:100] + ref[105:], None)
    c = E.Aln(0, [(E.M, 200), (E.I, 3), (E.M, 200)], ref[:200] + "TTT" + ref[200:], None)
    assert SP.clean_window(a, 50, 150) and not SP.clean_window(b, 50, 150) and SP.clean_window(b, 150, 250)
    assert not SP.clean_window(c, 150, 250) and SP.clean_window(c, 0, 100)
    assert not SP.clean_window(a, -5, 10) and not SP.clean_window(a, 390, 405)
    base, n, f = SP.consensus_base([a, b, c], 120)
    assert base == ref[120] and n == 3 and f == 1.0
    base, n, f = SP.consensus_base([a, b, c], 102)     # deleted in b
    assert n == 2
    assert SP.consensus_seq([a, b, c], 300, 310) == ref[300:310]


def test_evaluate_scores_class_and_parent_of_origin():
    plan = [SP.PlanRow("s1", "chr1", 100, "SNV", "SNV", 1, "A", "C", "G", 1, 1.0, ".", 0, 0.0, "paternal", "germline_DNM_phased", 1, 1),
            SP.PlanRow("s2", "chr1", 200, "SNV", "SNV", 1, "A", "C", "G", 2, 1.0, ".", 0, 0.0, "maternal", "germline_DNM_phased", 1, 1),
            SP.PlanRow("s3", "chr1", 300, "SV", "DEL", 500, ".", ".", "PM", 1, 1.0, "F", 2, 0.2, "paternal", "parental_mosaic_transmitted", 1, 1),
            SP.PlanRow("s4", "chr1", 400, "SV", "DEL", 500, ".", ".", "PM", 1, 1.0, "F", 2, 0.2, "paternal", "parental_mosaic_transmitted", 1, 1)]
    ev = [dict(variant_id="s1", variant_class="SNV", phase_class="germline_DNM_phased", parent_of_origin="paternal", hap_obs_k5="6", phase_score="0.98", rule_score="6", flags="", lik_post_germline="0.98"),
          dict(variant_id="s2", variant_class="SNV", phase_class="germline_DNM_unphased", parent_of_origin="maternal", hap_obs_k5="5", phase_score="0.7", rule_score="4", flags="LOW_HAP_DEPTH", lik_post_germline="0.7"),
          dict(variant_id="s3", variant_class="SV", phase_class="parental_mosaic_transmitted", parent_of_origin="paternal", hap_obs_k5="6", phase_score="0.01", rule_score="3", flags="", lik_post_parental_mosaic="0.9")]
    per_site, summary = SP.evaluate(plan, ev)
    by = {d["variant_id"]: d for d in per_site}
    assert by["s1"]["class_ok"] == 1 and by["s1"]["observable"] == 1 and by["s1"]["poo_ok"] == 1
    assert by["s2"]["class_ok"] == 0 and by["s2"]["class_ok_lenient"] == 1 and by["s2"]["observable"] == 0
    assert by["s4"]["reviewed"] == 0
    s = {(r["variant_class"], r["scenario"]): r for r in summary}
    assert s[("SNV", "G")]["n_planted"] == 2 and s[("SNV", "G")]["n_observable"] == 1 and s[("SNV", "G")]["class_ok_observable"] == 1.0
    assert s[("SNV", "G")]["class_ok_all"] == 0.5 and s[("SNV", "G")]["class_ok_lenient_all"] == 1.0
    assert s[("SV", "PM")]["n_reviewed"] == 1 and s[("SV", "PM")]["poo_ok_all"] == 1.0 and s[("SV", "PM")]["mean_phase_score"] == 0.01
    assert s[("SV", "PM")]["mean_post_expected"] == 0.9 and s[("SNV", "G")]["mean_post_expected"] == 0.84
    assert s[("SNV", "G")]["called_classes"] == "germline_DNM_phased:1;germline_DNM_unphased:1"
    # a table restricted to one class group evaluates only that group's planted sites
    _, summary_sv = SP.evaluate(plan, [ev[2]])
    assert {r["variant_class"] for r in summary_sv} == {"SV"}


def test_parse_regions():
    assert SP.parse_regions("chr20:1000-2000, chr21:5-9") == [("chr20", 1000, 2000), ("chr21", 5, 9)]
