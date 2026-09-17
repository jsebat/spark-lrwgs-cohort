"""The QUALITY verdict: is this variant real and de novo? Decided from the evidence 05_denovo already measured.

Recorded independently of clinical relevance (DESIGN P32): a variant can be quality-PASS and clinically uninteresting,
and that is a complete result. The verdict is a technical judgement about depth, junctions, heterozygosity, parental
coverage and phase, so unlike a clinical call it is reusable -- arm A's verdicts are an unbiased sample and may serve
as validation truth for the evidence (P33); arm B's are conditioned on the variant being interesting and may not.

Three verdicts. PASS: the evidence supports a constitutional de novo event. REVIEW: consistent but incomplete, worth
an analyst's eyes (this is where DNMT3A lands: depth halved and phased, but no junction read at both ends of a 302 bp
event). FAIL: the evidence says inherited, artefact, common, or absent. Every verdict carries its reasons, so a reader
can disagree with a rule rather than with a label.

The constants are the P8 deletion-path constants (`del_max_ratio` 0.7, `del_parent_min_ratio` 0.85,
`del_max_het_persistence` 0.35). P31 records that they were chosen with MECP2 and DNMT3A in view, which is why these
two events are an ACCEPTANCE test for the module and not a validation of it.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

PASS, REVIEW, FAIL = "PASS", "REVIEW", "FAIL"

# phase classes that say "not a constitutional de novo event"
NOT_DE_NOVO = {"phase_conflict_artifact": "PHASE_CONFLICT_ARTIFACT",
               "inherited_missed_in_parent": "INHERITED_MISSED_IN_PARENT",
               "parental_mosaic_transmitted": "PARENTAL_MOSAIC_TRANSMITTED"}


def _f(x) -> Optional[float]:
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _i(x) -> int:
    v = _f(x)
    return int(v) if v is not None else 0


@dataclass
class Verdict:
    verdict: str
    reasons: List[str] = field(default_factory=list)
    met: int = 0
    of: int = 0

    def as_row(self, prefix: str = "quality") -> Dict[str, object]:
        return {prefix + "_verdict": self.verdict, prefix + "_reasons": ";".join(self.reasons),
                prefix + "_criteria": "%d/%d" % (self.met, self.of) if self.of else ""}


def _common_fails(r: Dict[str, object]) -> List[str]:
    """Rejections that apply to every class: not de novo by phase, recurrent in unrelated people, masked."""
    out = []
    pc = str(r.get("phase_class") or "")
    if pc in NOT_DE_NOVO:
        out.append(NOT_DE_NOVO[pc])
    if _i(r.get("cohort_AC_loo")) > 0:
        out.append("SEEN_IN_OTHER_FAMILY(cohort_AC_loo=%d)" % _i(r.get("cohort_AC_loo")))
    if _i(r.get("pon_founder_recurrence_loo")) > 0:
        out.append("IN_FOUNDER_PANEL(n=%d)" % _i(r.get("pon_founder_recurrence_loo")))
    if str(r.get("site_filter_fail") or "") not in ("", "0", "0.0", "None"):
        out.append("SITE_FILTER_FAIL")
    return out


def _mask_note(r: Dict[str, object], frac_fail: float = 0.5) -> Tuple[bool, Optional[str]]:
    """Mask overlap is a FLAG in the de novo module (P5); here an event that is mostly repeat fails, a touch is noted.
    For SVs `mask_overlap` is fractional; for point variants it is 0/1."""
    m = _f(r.get("mask_overlap"))
    if m is None or m <= 0:
        return False, None
    if m >= frac_fail:
        return True, "MASKED(%.2f)" % m
    return False, "TOUCHES_MASK(%.2f)" % m


# ----------------------------------------------------------------------------------------------
# structural variants
# ----------------------------------------------------------------------------------------------
@dataclass
class SvParams:
    del_max_ratio: float = 0.7
    dup_min_ratio: float = 1.3
    del_parent_min_ratio: float = 0.85
    del_max_het_persistence: float = 0.35
    min_interval_for_het: int = 5000
    min_child_junction_reads: int = 2
    max_parent_junction_reads: int = 1
    min_child_support: int = 3
    min_rule_score_pass: int = 5
    min_rule_score_review: int = 3
    min_rf_review: float = 0.3


def sv_quality(r: Dict[str, object], p: SvParams = SvParams()) -> Verdict:
    """r: a final-table SV row merged with its evidence-table row (raw six-haplotype SV columns)."""
    reasons = list(_common_fails(r))
    svtype = str(r.get("svtype") or "")
    length = abs(_i(r.get("svlen"))) or max(0, _i(r.get("end")) - _i(r.get("start")))
    masked, note = _mask_note(r)
    if note:
        reasons.append(note)
    if masked:
        reasons.append("FAIL_MOSTLY_REPEAT")
    # parents carrying the event: junction reads named by sawfish, or depth lost across the interval
    pj = max(_i(r.get("F_sv_junc_hap1")) + _i(r.get("F_sv_junc_hap2")) + _i(r.get("F_sv_junc_hapu")),
             _i(r.get("M_sv_junc_hap1")) + _i(r.get("M_sv_junc_hap2")) + _i(r.get("M_sv_junc_hapu")),
             _i(r.get("max_parent_sv_support")))
    if pj > p.max_parent_junction_reads:
        reasons.append("PARENT_JUNCTION_READS(%d)" % pj)
    pr = [_f(r.get("F_sv_ratio_all")), _f(r.get("M_sv_ratio_all"))]
    pr = [x for x in pr if x is not None]
    if svtype in ("DEL", "CNV") and pr and any(x <= p.del_max_ratio for x in pr):
        reasons.append("PARENT_DEPTH_DEPLETED(%.2f)" % min(pr))
    if svtype == "DUP" and pr and any(x >= p.dup_min_ratio for x in pr):
        reasons.append("PARENT_DEPTH_GAINED(%.2f)" % max(pr))
    if any(x.startswith(("PHASE_CONFLICT", "INHERITED", "PARENTAL_MOSAIC", "SEEN_IN_OTHER", "IN_FOUNDER", "PARENT_", "FAIL_", "SITE_FILTER")) for x in reasons):
        return Verdict(FAIL, reasons)

    ratio = _f(r.get("sv_depth_ratio_inside_flank"))
    het = _f(r.get("sv_het_snv_persistence"))
    both = _i(r.get("C_sv_junc_both_ends") if r.get("C_sv_junc_both_ends") not in (None, "") else r.get("c_sv_junc_both_ends"))
    cj = _i(r.get("C_sv_junc_hap1")) + _i(r.get("C_sv_junc_hap2")) + _i(r.get("C_sv_junc_hapu"))
    csup = max(cj, _i(r.get("child_sv_support")))
    rule = _i(r.get("rule_score"))
    rf = _f(r.get("rf_prob"))

    if svtype in ("DEL", "CNV", "DUP") and ratio is not None:
        # the interval path (P8 amendment): depth, parents, junctions, heterozygosity
        crit: List[Tuple[str, bool]] = []
        if svtype == "DUP":
            crit.append(("CHILD_DEPTH_GAINED", ratio >= p.dup_min_ratio))
        else:
            crit.append(("CHILD_DEPTH_HALVED", ratio <= p.del_max_ratio))
        crit.append(("PARENTS_NOT_DEPLETED", bool(pr) and all(p.del_max_ratio < x for x in pr) and all(x >= p.del_parent_min_ratio for x in pr) if svtype != "DUP" else bool(pr)))
        crit.append(("BOTH_PARENTS_MEASURED", len(pr) == 2))
        crit.append(("JUNCTION_BOTH_ENDS", both > 0))
        crit.append(("CHILD_JUNCTION_READS>=%d" % p.min_child_junction_reads, csup >= p.min_child_junction_reads))
        if svtype != "DUP" and length >= p.min_interval_for_het:
            crit.append(("LOSS_OF_HETEROZYGOSITY", het is not None and het <= p.del_max_het_persistence))
        met = sum(1 for _, ok in crit if ok)
        for name, ok in crit:
            reasons.append(("+" if ok else "-") + name)
        depth_ok = crit[0][1]
        parents_ok = crit[1][1]
        hets_persist = svtype != "DUP" and het is not None and het > p.del_max_het_persistence and length >= p.min_interval_for_het
        if depth_ok and parents_ok and (both > 0 or (svtype != "DUP" and het is not None and het <= p.del_max_het_persistence)) and not hets_persist:
            return Verdict(PASS, reasons, met, len(crit))
        if not depth_ok and hets_persist:
            reasons.append("NO_DEPTH_LOSS_HETS_PERSIST")
            return Verdict(FAIL, reasons, met, len(crit))
        if hets_persist and depth_ok:
            reasons.append("DEPTH_LOST_BUT_HETS_PERSIST(mosaic_or_chimeric)")
            return Verdict(REVIEW, reasons, met, len(crit))
        if depth_ok or (both > 0 and csup >= p.min_child_junction_reads):
            return Verdict(REVIEW, reasons, met, len(crit))
        return Verdict(FAIL, reasons + ["NO_DEPTH_CHANGE_NO_JUNCTIONS"], met, len(crit))

    # insertions, inversions, breakends and short intervals: the breakpoint matrix and the rule score carry it
    crit = [("CHILD_SUPPORT>=%d" % p.min_child_support, csup >= p.min_child_support),
            ("PARENTS_NO_SUPPORT", pj == 0),
            ("ALT_CONFINED_TO_ONE_HAPLOTYPE", str(r.get("c_alt_confined") or "") in ("1", "1.0", "True")),
            ("RULE_SCORE>=%d" % p.min_rule_score_pass, rule >= p.min_rule_score_pass),
            ("GERMLINE_PHASE_CLASS", str(r.get("phase_class") or "").startswith("germline_DNM"))]
    met = sum(1 for _, ok in crit if ok)
    for name, ok in crit:
        reasons.append(("+" if ok else "-") + name)
    if crit[0][1] and crit[1][1] and crit[3][1] and crit[4][1]:
        return Verdict(PASS, reasons, met, len(crit))
    if (rule >= p.min_rule_score_review and crit[1][1]) or (rf is not None and rf >= p.min_rf_review):
        return Verdict(REVIEW, reasons, met, len(crit))
    return Verdict(FAIL, reasons, met, len(crit))


# ----------------------------------------------------------------------------------------------
# SNV / indel
# ----------------------------------------------------------------------------------------------
@dataclass
class SmallParams:
    min_alt_reads: int = 3
    min_rule_score_pass: int = 5
    min_rule_score_review: int = 3
    min_rf_review: float = 0.3
    max_gnomad_af: float = 0.001
    min_child_dp: int = 8


def smallvar_quality(r: Dict[str, object], p: SmallParams = SmallParams()) -> Verdict:
    reasons = list(_common_fails(r))
    masked, note = _mask_note(r, frac_fail=1.0)
    if note:
        reasons.append(note)
    af = _f(r.get("gnomad_af"))
    if af is not None and af >= p.max_gnomad_af:
        reasons.append("COMMON_IN_GNOMAD(af=%.4g)" % af)
    palt = _i(r.get("p_max_alt_any_hap"))
    if palt >= 2 or _i(r.get("t_alt_reads")) >= 2:
        reasons.append("PARENT_ALT_READS(max_hap=%d,transmitted=%d)" % (palt, _i(r.get("t_alt_reads"))))
    if any(x.startswith(("PHASE_CONFLICT", "INHERITED", "PARENTAL_MOSAIC", "SEEN_IN_OTHER", "IN_FOUNDER", "COMMON_", "PARENT_", "SITE_FILTER")) for x in reasons):
        return Verdict(FAIL, reasons)
    pc = str(r.get("phase_class") or "")
    rule = _i(r.get("rule_score"))
    rf = _f(r.get("rf_prob"))
    alt = _i(r.get("c_alt_hapA")) + _i(r.get("c_untagged_alt"))
    dp = _f(r.get("child_DP"))
    crit = [("ALT_READS>=%d" % p.min_alt_reads, alt >= p.min_alt_reads),
            ("ALT_CONFINED_TO_ONE_HAPLOTYPE", str(r.get("c_alt_confined") or "") in ("1", "1.0", "True")),
            ("PARENTS_ALT_FREE", palt <= 1),
            ("RULE_SCORE>=%d" % p.min_rule_score_pass, rule >= p.min_rule_score_pass),
            ("GERMLINE_PHASE_CLASS", pc.startswith("germline_DNM")),
            ("CHILD_DP>=%d" % p.min_child_dp, dp is not None and dp >= p.min_child_dp)]
    met = sum(1 for _, ok in crit if ok)
    for name, ok in crit:
        reasons.append(("+" if ok else "-") + name)
    if pc == "child_postzygotic_mosaic":
        reasons.append("POSTZYGOTIC_MOSAIC")
        return Verdict(REVIEW, reasons, met, len(crit))
    if all(ok for _, ok in crit):
        return Verdict(PASS, reasons, met, len(crit))
    if (pc.startswith("germline_DNM") and rule >= p.min_rule_score_review and alt >= p.min_alt_reads) or (rf is not None and rf >= p.min_rf_review):
        return Verdict(REVIEW, reasons, met, len(crit))
    return Verdict(FAIL, reasons, met, len(crit))


# ----------------------------------------------------------------------------------------------
# tandem repeats
# ----------------------------------------------------------------------------------------------
def tr_quality(r: Dict[str, object]) -> Verdict:
    reasons = list(_common_fails(r))
    if any(x.startswith(("PHASE_CONFLICT", "INHERITED", "PARENTAL_MOSAIC", "SEEN_IN_OTHER", "IN_FOUNDER")) for x in reasons):
        return Verdict(FAIL, reasons)
    rule = _i(r.get("rule_score"))
    rf = _f(r.get("rf_prob"))
    call = str(r.get("dnm_call") or "")
    crit = [("CALLED_TIER1", call == "YES"),
            ("RULE_SCORE>=5", rule >= 5),
            ("GERMLINE_PHASE_CLASS", str(r.get("phase_class") or "").startswith("germline_DNM")),
            ("ALT_CONFINED_TO_ONE_HAPLOTYPE", str(r.get("c_alt_confined") or "") in ("1", "1.0", "True"))]
    met = sum(1 for _, ok in crit if ok)
    for name, ok in crit:
        reasons.append(("+" if ok else "-") + name)
    if call == "YES" and crit[2][1]:
        return Verdict(PASS, reasons, met, len(crit))
    if call == "CANDIDATE" or rule >= 3 or (rf is not None and rf >= 0.3):
        return Verdict(REVIEW, reasons, met, len(crit))
    return Verdict(FAIL, reasons, met, len(crit))
