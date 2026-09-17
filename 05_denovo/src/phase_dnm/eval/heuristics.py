"""The original pipeline's heuristic de novo filters as scorers (DESIGN P22, P25). Each arm returns

    pass_        the pipeline's actual decision with its thresholds verbatim
    sweep_score  the rule's natural knob with every other criterion held fixed; rows that fail a fixed criterion score
                 below every passing row (a large negative offset minus the knob), so the ROC sweep's operating point
                 is the pipeline's decision and the curve is labelled a sweep, never "the pipeline"

Inputs are the module's own rows: a CandidateRecord-shaped dict (caller fields) and, for the cohort arms, the
annotation columns of the feature table. Sources, verbatim:
  SNV/indel H1  01_qc/family_downstream.sb, slivar:  kid.het && mom.hom_ref && dad.hom_ref && kid.GQ>=20 && mom.GQ>=20 &&
                dad.GQ>=20 && kid.DP>=10 && mom.DP>=10 && dad.DP>=10 && kid.AB>=0.25 && kid.AB<=0.75 && mom.AD[1]==0 && dad.AD[1]==0
  SNV/indel H2  H1 and outside the segdup/repeat mask and gnomAD AF < 0.001  (hiconf)
  SNV/indel H3  H2 and PON leave-one-out AC < 3 and not shared by both siblings  (04_panel/apply_pon_loo.sb, 05_denovo)
  SV H1         child GT carries the allele, both parents 0/0, no quality filter; knob = child GQ
  SV H2         H1 and mask overlap < 50 % of the interval and founder-panel recurrence AC < 6 (reciprocal overlap 0.5)
  TR H1         child's longest allele >= longest parental + 50 bp and >= 1.5x it (family_downstream.sb); knob = gain_bp
  TR H2         03_tiering/tr_outliers.py: SD >= 5 on every allele of all three, child allele > max(parents) + max(6 bp,
                3 x shortest motif), and > founder maximum; knob = excess over max(parents) in bp
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

FAIL = -1e6          # offset for rows failing a fixed criterion: below any passing row, ordered by the knob within


def _f(x) -> Optional[float]:
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _ints(s) -> Optional[List[int]]:
    if s in (None, "", "."):
        return None
    try:
        return [int(float(x)) for x in str(s).split(",")]
    except ValueError:
        return None


def _gt_alleles(gt: str) -> Optional[Tuple[int, ...]]:
    if gt in (None, "", ".", "./.", ".|."):
        return None
    try:
        return tuple(int(a) for a in str(gt).replace("|", "/").split("/") if a != ".")
    except ValueError:
        return None


def _het(gt: str) -> bool:
    a = _gt_alleles(gt)
    return a is not None and len(a) == 2 and a[0] != a[1] and 0 in a


def _homref(gt: str) -> bool:
    a = _gt_alleles(gt)
    return a is not None and len(a) >= 1 and all(x == 0 for x in a)


def _carries(gt: str) -> bool:
    a = _gt_alleles(gt)
    return a is not None and any(x > 0 for x in a)


def allele_balance(ad: Optional[List[int]], allele_index: int = 1) -> Optional[float]:
    if not ad or sum(ad) == 0 or allele_index >= len(ad):
        return None
    return ad[allele_index] / sum(ad)


# ----------------------------------------------------------------------------------------------
# SNV / indel
# ----------------------------------------------------------------------------------------------
def slivar_h1(rec: Dict[str, object], min_gq: int = 20, min_dp: int = 10, ab: Tuple[float, float] = (0.25, 0.75)) -> Tuple[bool, float]:
    """H1: the slivar trio expression, knob = min GQ of the trio (higher = more confident)."""
    c_gt, f_gt, m_gt = rec.get("caller_gt"), rec.get("father_gt"), rec.get("mother_gt")
    c_gq, f_gq, m_gq = _f(rec.get("caller_gq")), _f(rec.get("father_gq")), _f(rec.get("mother_gq"))
    c_dp, f_dp, m_dp = _f(rec.get("caller_dp")), _f(rec.get("father_dp")), _f(rec.get("mother_dp"))
    c_ad, f_ad, m_ad = _ints(rec.get("child_ad")), _ints(rec.get("father_ad")), _ints(rec.get("mother_ad"))
    ai = 1
    try:
        pl = rec.get("class_payload") or {}
        ai = int(pl.get("allele_index", 1)) if isinstance(pl, dict) else 1
    except (TypeError, ValueError):
        ai = 1
    if c_dp is None and c_ad:
        c_dp = float(sum(c_ad))
    if f_dp is None and f_ad:
        f_dp = float(sum(f_ad))
    if m_dp is None and m_ad:
        m_dp = float(sum(m_ad))
    c_ab = allele_balance(c_ad, ai)
    fixed = (_het(c_gt) and _homref(f_gt) and _homref(m_gt)
             and c_dp is not None and f_dp is not None and m_dp is not None and c_dp >= min_dp and f_dp >= min_dp and m_dp >= min_dp
             and c_ab is not None and ab[0] <= c_ab <= ab[1]
             and f_ad is not None and m_ad is not None and len(f_ad) > ai and len(m_ad) > ai and f_ad[ai] == 0 and m_ad[ai] == 0)
    gqs = [g for g in (c_gq, f_gq, m_gq) if g is not None]
    knob = min(gqs) if len(gqs) == 3 else -1.0
    pass_ = bool(fixed and knob >= min_gq)
    return pass_, (knob if fixed else FAIL + knob)


def snv_h2(rec: Dict[str, object], feat: Dict[str, object], af_max: float = 0.001) -> Tuple[bool, float]:
    """H2 hiconf: H1 and outside the mask and gnomAD AF < 0.001 (absent from gnomAD counts as rare)."""
    p1, s1 = slivar_h1(rec)
    mask = str(feat.get("segdup_overlap", rec.get("mask_overlap", "0"))) == "1"
    af = _f(feat.get("gnomad_af"))
    rare = af is None or af < af_max
    ok = p1 and not mask and rare
    return ok, (s1 if (not mask and rare and s1 > FAIL / 2) else FAIL + min(s1, 0))


def snv_h3(rec: Dict[str, object], feat: Dict[str, object], pon_max_ac: int = 3) -> Tuple[bool, float]:
    """H3 cohort: H2 and PON leave-one-out AC < 3 and not shared by both siblings."""
    p2, s2 = snv_h2(rec, feat)
    pon = _f(feat.get("pon_founder_recurrence_loo"))
    sib = str(feat.get("sib_shared", "0")) == "1"
    ok = p2 and (pon is None or pon < pon_max_ac) and not sib
    return ok, (s2 if ((pon is None or pon < pon_max_ac) and not sib and s2 > FAIL / 2) else FAIL + min(s2, 0))


# ----------------------------------------------------------------------------------------------
# SV
# ----------------------------------------------------------------------------------------------
def sv_h1(rec: Dict[str, object]) -> Tuple[bool, float]:
    """H1: child carries the allele, both parents 0/0; knob = child GQ (no quality filter in the pipeline)."""
    fixed = _carries(rec.get("caller_gt")) and _homref(rec.get("father_gt")) and _homref(rec.get("mother_gt"))
    gq = _f(rec.get("caller_gq"))
    knob = gq if gq is not None else 0.0
    return bool(fixed), (knob if fixed else FAIL + knob)


def sv_h2(rec: Dict[str, object], feat: Dict[str, object], mask_frac_max: float = 0.5, panel_ac_max: int = 6) -> Tuple[bool, float]:
    """H2 cohort: H1 and mask overlap < 50 % of the interval and founder-panel recurrence AC < 6."""
    p1, s1 = sv_h1(rec)
    mf = _f(feat.get("mask_frac", feat.get("segdup_overlap")))
    ac = _f(feat.get("pon_founder_recurrence_loo"))
    ok_fixed = (mf is None or mf < mask_frac_max) and (ac is None or ac < panel_ac_max)
    return bool(p1 and ok_fixed), (s1 if (ok_fixed and s1 > FAIL / 2) else FAIL + min(s1, 0))


# ----------------------------------------------------------------------------------------------
# TR
# ----------------------------------------------------------------------------------------------
def _tr_alleles(rec: Dict[str, object]) -> Tuple[Optional[List[int]], Optional[List[int]], Optional[List[int]], Dict]:
    pl = rec.get("class_payload") or {}
    if not isinstance(pl, dict):
        return None, None, None, {}
    return pl.get("child_AL"), pl.get("father_AL"), pl.get("mother_AL"), pl


def tr_h1(rec: Dict[str, object], min_gain_bp: int = 50, min_ratio: float = 1.5) -> Tuple[bool, float]:
    """H1 per family: child's longest allele >= longest parental + 50 bp and >= 1.5x; knob = gain in bp."""
    c, f, m, _ = _tr_alleles(rec)
    if not c or not f or not m:
        return False, FAIL - 1
    cmax, pmax = max(c), max(f + m)
    gain = cmax - pmax
    fixed = pmax > 0 and cmax >= min_ratio * pmax
    return bool(fixed and gain >= min_gain_bp), (gain if fixed else FAIL + gain)


def tr_h2(rec: Dict[str, object], feat: Dict[str, object], min_sd: int = 5, min_excess_units: int = 3, min_excess_bp: int = 6) -> Tuple[bool, float]:
    """H2 cohort (tr_outliers.py): SD >= 5 on every allele of all three, child allele > max(parents) + max(6, 3 x unit),
    and > founder maximum (feature `trgt_founder_max_loo`; absent -> criterion not evaluable, treated as passed and
    flagged by the caller). Knob = excess over max(parents) in bp."""
    c, f, m, pl = _tr_alleles(rec)
    if not c or not f or not m:
        return False, FAIL - 1
    def sd_ok(key):
        v = _ints(pl.get(key))
        return v is not None and all(x >= min_sd for x in v)
    unit = int(pl.get("motif_unit_bp") or 1)
    excess_min = max(min_excess_bp, min_excess_units * unit)
    idx = int(pl.get("outlier_allele_idx", 0) or 0)
    allele = c[idx] if idx < len(c) else max(c)
    excess = allele - max(f + m)
    founder_max = _f(feat.get("trgt_founder_max_loo"))
    fixed = sd_ok("child_SD") and sd_ok("father_SD") and sd_ok("mother_SD") and (founder_max is None or allele > founder_max)
    return bool(fixed and excess >= excess_min), (excess if fixed else FAIL + excess)


ARMS = {
    "snv_indel": {"H1_slivar": lambda r, f: slivar_h1(r), "H2_hiconf": snv_h2, "H3_cohort": snv_h3},
    "sv": {"H1_genotype": lambda r, f: sv_h1(r), "H2_cohort": sv_h2},
    "tr": {"H1_family": lambda r, f: tr_h1(r), "H2_cohort": tr_h2},
}


def score_row(class_group: str, rec: Dict[str, object], feat: Dict[str, object]) -> Dict[str, object]:
    out: Dict[str, object] = {}
    for name, fn in ARMS[class_group].items():
        p, s = fn(rec, feat)
        out[name + "_pass"] = int(p)
        out[name + "_score"] = s
    return out
