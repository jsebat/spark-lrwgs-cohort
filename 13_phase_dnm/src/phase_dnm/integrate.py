"""Module 3 — integration: one UNFILTERED final table per class per family (README §2.3), the P15 decision, and the
class-specific columns. Inputs are the M2 tables that already exist for every candidate: the likelihood-scored
evidence table (rule class, parent of origin, six-haplotype counts, posterior) and the registered feature table.
`rf_prob` comes from Module 4; until a model exists the decision runs in a documented PROVISIONAL mode.

P15 (DESIGN):
  dnm_call = YES  iff rf_prob >= tau_class,
               or (rf_prob >= tau_rescue_class and phase_class == germline_DNM_phased and hap_obs_k5 == 6)   [rescue]
  dnm_call = NO   whenever phase_class in {phase_conflict_artifact, inherited_missed_in_parent}, whatever rf_prob    [demotion]
  mosaic classes are never YES; they are reported separately (mosaic_flag).
Provisional mode (call_mode = phase_only; rf_prob missing or tau unset): YES iff phase_class == germline_DNM_phased,
hap_obs_k5 == 6 and phase_score >= phase_only_min_score — the intersection that carried the parent-of-origin signal
on the cohort (paternal fraction 0.72 vs 0.50 for either layer alone, PLAN 2026-09-13). Every row keeps the reason.
"""
from __future__ import annotations

import csv
import json
import os
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple

FINAL_CORE = ["family_id", "sample_id", "chrom", "start", "end", "ref", "alt", "variant_class", "caller", "caller_gt", "caller_qual",
              "rf_prob", "dnm_call", "call_mode", "decision_reason", "mosaic_flag",
              "parent_of_origin", "poo_reason", "poo_confidence",
              "phase_class", "rule_score", "hap_obs_k3", "hap_obs_k5", "child_alt_hap_frac", "child_alt_other_hap",
              "transmitted_parent_alt_reads", "untransmitted_parent_alt_reads", "transmitted_parent_dp", "phase_score",
              "lik_best_alternative", "lik_log10lr_germline", "flags", "source_tier", "source_list", "mask_overlap",
              "variant_id", "thresholds_version"]
CLASS_COLS = {"SV": ["svtype", "svlen", "bp_precision"],
              "TR": ["trid", "motif", "child_AL_pat", "child_AL_mat", "father_AL_T", "father_AL_U", "mother_AL_T", "mother_AL_U"]}
DEMOTE = ("phase_conflict_artifact", "inherited_missed_in_parent")
MOSAIC = ("child_postzygotic_mosaic", "parental_mosaic_transmitted")


@dataclass
class FinalParams:
    tau: Dict[str, Optional[float]] = field(default_factory=dict)          # per variant_class, from M4 outer folds
    tau_rescue: Dict[str, Optional[float]] = field(default_factory=dict)
    phase_only_min_score: float = 0.9
    require_hap_obs: int = 6


def _f(x) -> Optional[float]:
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def decide(row: Dict[str, object], p: FinalParams, rf_prob: Optional[float]) -> Dict[str, object]:
    """P15. Returns dnm_call, call_mode, decision_reason, mosaic_flag."""
    cls = str(row.get("phase_class") or "")
    vclass = str(row.get("variant_class") or "")
    hap_obs = _f(row.get("hap_obs_k5"))
    six = hap_obs is not None and int(hap_obs) >= p.require_hap_obs
    mosaic = int(cls in MOSAIC)
    if cls in DEMOTE:
        return dict(dnm_call="NO", call_mode="rf+phase" if rf_prob is not None else "phase_only", decision_reason="DEMOTED:%s" % cls, mosaic_flag=mosaic)
    if mosaic:
        return dict(dnm_call="NO", call_mode="rf+phase" if rf_prob is not None else "phase_only", decision_reason="MOSAIC:%s" % cls, mosaic_flag=1)
    tau, tau_r = p.tau.get(vclass), p.tau_rescue.get(vclass)
    if rf_prob is not None and tau is not None:
        if rf_prob >= tau:
            return dict(dnm_call="YES", call_mode="rf+phase", decision_reason="RF", mosaic_flag=0)
        if tau_r is not None and rf_prob >= tau_r and cls == "germline_DNM_phased" and six:
            return dict(dnm_call="YES", call_mode="rf+phase", decision_reason="RESCUED", mosaic_flag=0)
        return dict(dnm_call="NO", call_mode="rf+phase", decision_reason="BELOW_TAU", mosaic_flag=0)
    ps = _f(row.get("phase_score"))
    if cls == "germline_DNM_phased" and six and ps is not None and ps >= p.phase_only_min_score:
        return dict(dnm_call="YES", call_mode="phase_only", decision_reason="PHASE_ONLY", mosaic_flag=0)
    why = "NOT_PHASED_GERMLINE" if cls != "germline_DNM_phased" else ("HAP_UNOBSERVED" if not six else "LOW_POSTERIOR")
    return dict(dnm_call="NO", call_mode="phase_only", decision_reason=why, mosaic_flag=0)


# ----------------------------------------------------------------------------------------------
# class-specific columns from the evidence row's class_payload
# ----------------------------------------------------------------------------------------------
def _payload(row: Dict[str, object]) -> dict:
    try:
        return json.loads(str(row.get("class_payload") or "{}"))
    except ValueError:
        return {}


def _phased_alleles(gt: str, al: Optional[list]) -> Optional[Tuple[int, int]]:
    """(hap1 allele length, hap2 allele length) from a PHASED TRGT genotype 'a|b' and the AL list; None otherwise."""
    if not al or not isinstance(gt, str) or "|" not in gt:
        return None
    try:
        a, b = (int(x) for x in gt.split("|"))
        return al[a], al[b]
    except (ValueError, IndexError, TypeError):
        return None


def class_columns(row: Dict[str, object]) -> Dict[str, object]:
    vclass = row.get("variant_class")
    pl = _payload(row)
    if vclass == "SV":
        return dict(svtype=pl.get("svtype"), svlen=pl.get("svlen"), bp_precision=None if pl.get("imprecise") is None else ("IMPRECISE" if pl.get("imprecise") else "PRECISE"))
    if vclass == "TR":
        out: Dict[str, object] = dict(trid=pl.get("trid"), motif=(pl.get("motifs") or "").split(",")[0] or None,
                                      child_AL_pat=None, child_AL_mat=None, father_AL_T=None, father_AL_U=None, mother_AL_T=None, mother_AL_U=None)
        h1 = row.get("child_hap1_is")
        c = _phased_alleles(str(row.get("caller_gt") or ""), pl.get("child_AL"))
        if c and h1 in ("P", "M"):
            out["child_AL_pat"], out["child_AL_mat"] = (c[0], c[1]) if h1 == "P" else (c[1], c[0])
        for par, key, gt_key in (("F", "father", "father_gt"), ("M", "mother", "mother_gt")):
            t = str(row.get(par + "_transmitted_hap") or "")
            a = _phased_alleles(str(pl.get(gt_key) or row.get(gt_key) or ""), pl.get(key + "_AL"))
            if a and t in ("1", "2"):
                ti = int(t) - 1
                out[key + "_AL_T"], out[key + "_AL_U"] = a[ti], a[1 - ti]
        return out
    return {}


# ----------------------------------------------------------------------------------------------
# assembly
# ----------------------------------------------------------------------------------------------
RENAME = {"child_alt_hap_frac": "c_alt_hap_frac", "child_alt_other_hap": "c_alt_hapO", "transmitted_parent_alt_reads": "t_alt_reads",
          "untransmitted_parent_alt_reads": "u_alt_reads", "transmitted_parent_dp": "t_dp"}


def final_columns(vclass_group: str, feature_names: List[str]) -> List[str]:
    cols = list(FINAL_CORE)
    for vc in (("SV",) if vclass_group == "sv" else ("TR",) if vclass_group == "tr" else ()):
        cols += CLASS_COLS.get(vc, [])
    cols += [f for f in feature_names if f not in cols]
    return cols


def load_rf_probs(path: Optional[str], column: str = "rf_prob") -> Dict[str, float]:
    """variant_id -> score. `column` is rf_prob (calibrated probability) or rf_q (fold-quantile score, train/rescore.py)."""
    out: Dict[str, float] = {}
    if not path or not os.path.exists(path):
        return out
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        col = column if column in (rd.fieldnames or []) else "rf_prob"
        for r in rd:
            v = _f(r.get(col))
            if v is not None:
                out[r["variant_id"]] = v
    return out


def integrate_table(evidence_lik_tsv: str, features_tsv: Optional[str], out_tsv: str, p: FinalParams, vclass_group: str,
                    rf_probs: Optional[Dict[str, float]] = None) -> Dict[str, object]:
    feats: Dict[str, dict] = {}
    feat_names: List[str] = []
    if features_tsv and os.path.exists(features_tsv):
        with open(features_tsv, newline="") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            id_cols = {"family_id", "sample_id", "variant_id", "chrom", "start", "end", "ref", "alt", "variant_class", "caller"}
            feat_names = [c for c in (rd.fieldnames or []) if c not in id_cols]
            for r in rd:
                feats[r["variant_id"]] = r
    cols = final_columns(vclass_group, feat_names)
    counts: Dict[str, int] = {"rows": 0, "YES": 0, "NO": 0, "mosaic": 0}
    reasons: Dict[str, int] = {}
    per_class_yes: Dict[str, int] = {}
    rf_probs = rf_probs or {}
    os.makedirs(os.path.dirname(os.path.abspath(out_tsv)), exist_ok=True)
    with open(evidence_lik_tsv, newline="") as fh, open(out_tsv, "w", newline="") as out:
        rd = csv.DictReader(fh, delimiter="\t")
        w = csv.DictWriter(out, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rd:
            row: Dict[str, object] = dict(r)
            for new, old in RENAME.items():
                row[new] = r.get(old)
            rf = rf_probs.get(r["variant_id"])
            row["rf_prob"] = "" if rf is None else rf
            d = decide(row, p, rf)
            row.update(d)
            row.update(class_columns(r))
            f = feats.get(r["variant_id"], {})
            for k in feat_names:
                if k not in row or row.get(k) in (None, ""):
                    row[k] = f.get(k, "")
            w.writerow({c: ("" if row.get(c) is None else row.get(c)) for c in cols})
            counts["rows"] += 1
            counts[d["dnm_call"]] += 1
            counts["mosaic"] += d["mosaic_flag"]
            reasons[d["decision_reason"]] = reasons.get(d["decision_reason"], 0) + 1
            if d["dnm_call"] == "YES":
                per_class_yes[r["variant_class"]] = per_class_yes.get(r["variant_class"], 0) + 1
    return {"counts": counts, "reasons": reasons, "yes_by_class": per_class_yes, "columns": len(cols), "features": len(feat_names)}


def write_parquet(tsv: str) -> Optional[str]:
    """Parquet twin when pyarrow is available; otherwise None (documented in the family summary)."""
    try:
        import pyarrow.csv as pc  # type: ignore
        import pyarrow.parquet as pq  # type: ignore
    except ImportError:
        return None
    out = tsv[:-4] + ".parquet" if tsv.endswith(".tsv") else tsv + ".parquet"
    t = pc.read_csv(tsv, parse_options=pc.ParseOptions(delimiter="\t"), convert_options=pc.ConvertOptions(strings_can_be_null=True))
    pq.write_table(t, out)
    return out
