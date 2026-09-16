#!/usr/bin/env python3
"""P28 evidence: recall against planted event size, for the classifier score path and the deterministic depth rule.

DESIGN P28 asserts that a classifier carries small and mid-size variants while large deletions are called from the
interval evidence by a deterministic rule, on the grounds that the pedigree-swap positives hold almost no large
deletions for the classifier to learn from. That was an argument, not a measurement. This tool measures it: planted
deletions of known size are scored by the fold model that held the child's family out, and the SAME rows are put
through integrate.decide() twice -- once with the depth rule enabled and once without -- so both paths are the
production decision function rather than a reimplementation of it.

  recall_classifier  YES with sv_depth_rule_tier1 = False  (score >= tau, germline review class, rule layer passed)
  recall_depth_rule  decision_reason == TIER1_SV_DEPTH with the rule enabled
  recall_combined    YES with the rule enabled (what the pipeline actually calls)

Only planted GERMLINE (scenario G) rows count as positives, matching the harness. Rows the review could not see at
all are reported as n_unobservable rather than silently lowering the denominator.

Usage:
  python tools/sv_size_recall.py --evidence-dir DIR --harness-dir DIR --folds-dir DIR --manifest TSV \
      --seed 0 --thresholds config/thresholds.json --out sv_size_recall.tsv
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from phase_dnm import integrate as I                     # noqa: E402
from phase_dnm.eval import external as EX                # noqa: E402
from phase_dnm.train import rescore as RS                # noqa: E402

# the columns integrate.rules_fail and integrate.decide read that do not travel in the harness EV_COLS
RULE_COLS = ["variant_id", "gnomad_af", "lr_sv_catalog_af", "pon_founder_recurrence_loo", "cohort_AC_loo",
             "segdup_overlap", "mask_overlap", "rule_score", "flags", "phase_class", "hap_obs_k5"]


def _read_some(path, cols):
    """Those of `cols` the file actually has, as strings; an empty frame when the file is absent."""
    if not os.path.exists(path):
        return pd.DataFrame(columns=cols)
    head = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
    use = [c for c in cols if c in head]
    if "variant_id" not in use:
        return pd.DataFrame(columns=cols)
    df = pd.read_csv(path, sep="\t", usecols=use, dtype=str, keep_default_na=False)
    for c in cols:
        if c not in df.columns:
            df[c] = ""
    return df


def _params(thresholds_path, rule_on, tau_value, tau_tier2):
    thr = json.load(open(thresholds_path)) if thresholds_path and os.path.exists(thresholds_path) else {}
    f = thr.get("final", {}) or {}
    rules = dict(I.DEFAULT_RULES)
    rules.update(f.get("rules") or {})
    return I.FinalParams(tau={"SV": tau_value}, tau_rescue={},
                         tau_tier2={"SV": tau_tier2} if tau_tier2 is not None else {},
                         phase_only_min_score=f.get("phase_only_min_score", 0.9),
                         require_hap_obs=f.get("require_hap_obs", 6),
                         tr_rescue_min_units=float(f.get("tr_rescue_min_units", 3)),
                         rules=rules, apply_rules=bool(f.get("apply_rules", True)),
                         sv_depth_rule_tier1=rule_on,
                         sv_depth_min_rule_score=float(f.get("sv_depth_min_rule_score", 6)))


def collect(a, fam_of, fold_of, fms):
    """Every planted SV row, scored by the fold model that held its family out, plus the real rows for the ecdf."""
    rows, real = [], []
    for sd in sorted(glob.glob(os.path.join(a.evidence_dir, "*", "spike", "*"))):
        child = os.path.basename(sd)
        fam = fam_of.get(child)
        if fam is None or fam not in fold_of or fold_of[fam] not in fms:
            continue
        fm = fms[fold_of[fam]]
        X, side = EX.spike_rows(sd, child, "sv")
        if X.empty:
            continue
        prob, _raw = fm.predict(X)
        extra = _read_some(os.path.join(sd, "sv.features.tsv"), RULE_COLS)
        side = side.merge(extra, on="variant_id", how="left", suffixes=("", "_f")).fillna("")
        for i, r in side.reset_index(drop=True).iterrows():
            vid = str(r["variant_id"])
            bits = vid.split(":")                        # spike:<chrom>:<pos>:<subtype>:<length>
            if len(bits) < 5 or not bits[4].isdigit():
                continue
            d = {k: r.get(k, "") for k in RULE_COLS}
            d.update(child=child, fold=fold_of[fam], variant_id=vid, subtype=bits[3], length=int(bits[4]),
                     scenario=str(r.get("scenario", "")), variant_class="SV", prob=float(prob[i]))
            rows.append(d)
        rf = os.path.join(a.evidence_dir, fam, "features", "%s.sv.features.rf.tsv" % child)
        if os.path.exists(rf):
            Xr = pd.read_csv(rf, sep="\t", dtype=str, keep_default_na=False)
            pr, _ = fm.predict(Xr)
            real.extend({"fold": fold_of[fam], "prob": float(x)} for x in pr)
    return rows, real


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--evidence-dir", required=True)
    ap.add_argument("--harness-dir", required=True)
    ap.add_argument("--folds-dir", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--thresholds")
    ap.add_argument("--out", required=True)
    ap.add_argument("--per-site-out", help="also write the per-row table behind the curve")
    a = ap.parse_args()

    man = {r["sample_id"]: r for r in csv.DictReader(open(a.manifest), delimiter="\t")}
    fam_of = {sid: r["family_id"] for sid, r in man.items()}
    fold_of = {r["family_id"]: int(r["outer_fold"])
               for r in csv.DictReader(open(os.path.join(a.folds_dir, "folds.seed%d.tsv" % a.seed)), delimiter="\t")}
    fms = EX.load_fold_models(a.harness_dir, "sv", a.seed)
    if not fms:
        sys.stderr.write("no sv fold models for seed %d in %s\n" % (a.seed, a.harness_dir))
        return 2

    rows, real = collect(a, fam_of, fold_of, fms)
    if not rows:
        sys.stderr.write("no spike sv rows found under %s\n" % a.evidence_dir)
        return 2
    df = pd.DataFrame(rows)
    rdf = pd.DataFrame(real) if real else pd.DataFrame(columns=["fold", "prob"])

    # fold-quantile score against the fold's REAL rows -- the construction tau_q was swept on
    df["rf_q"] = np.nan
    for k, g in df.groupby("fold"):
        ref = RS.ecdf_ref(rdf.loc[rdf["fold"] == k, "prob"].to_numpy())
        df.loc[g.index, "rf_q"] = RS.quantile(ref, g["prob"].to_numpy())

    tj = os.path.join(a.harness_dir, "tau.sv.json")
    t = json.load(open(tj)) if os.path.exists(tj) else {}
    use_q = t.get("score_column") == "rf_q"
    if use_q and t.get("tau_q") is not None:
        tau = float(t["tau_q"])
    elif t.get("tau") is not None:
        tau = float(t["tau"])
    else:
        sys.stderr.write("no tau in %s; refusing to invent one\n" % tj)
        return 2
    tau2 = t.get("tau_q_rescue") if use_q else t.get("tau_tier2")
    score = df["rf_q"].to_numpy() if use_q else df["prob"].to_numpy()
    sys.stderr.write("scoring with %s, tau=%.6g (%d spike rows, %d real rows for the ecdf)\n"
                     % ("rf_q" if use_q else "rf_prob", tau, len(df), len(rdf)))

    p_on = _params(a.thresholds, True, tau, float(tau2) if tau2 is not None else None)
    p_off = _params(a.thresholds, False, tau, float(tau2) if tau2 is not None else None)
    recs = df.to_dict("records")
    dec_on = [I.decide(r, p_on, float(score[i])) for i, r in enumerate(recs)]
    dec_off = [I.decide(r, p_off, float(score[i])) for i, r in enumerate(recs)]
    df["score"] = score
    df["call_rule_on"] = [d["dnm_call"] for d in dec_on]
    df["reason_rule_on"] = [d["decision_reason"] for d in dec_on]
    df["tier_rule_on"] = [d["dnm_tier"] for d in dec_on]
    df["call_rule_off"] = [d["dnm_call"] for d in dec_off]
    df["reason_rule_off"] = [d["decision_reason"] for d in dec_off]

    if a.per_site_out:
        df.to_csv(a.per_site_out, sep="\t", index=False)

    g = df[df["scenario"] == "G"]
    order = {"DEL": 0, "BIGDEL": 1}
    out = []
    for (sub, ln), grp in sorted(g.groupby(["subtype", "length"]), key=lambda kv: (order.get(kv[0][0], 9), kv[0][1])):
        n = len(grp)
        out.append(dict(subtype=sub, length_bp=ln, n_planted_G=n,
                        n_unobservable=int((grp["phase_class"] == "").sum()),
                        recall_classifier=round(int((grp["call_rule_off"] == "YES").sum()) / n, 4),
                        recall_depth_rule=round(int((grp["reason_rule_on"] == "TIER1_SV_DEPTH").sum()) / n, 4),
                        recall_combined=round(int((grp["call_rule_on"] == "YES").sum()) / n, 4),
                        n_classifier=int((grp["call_rule_off"] == "YES").sum()),
                        n_depth_rule=int((grp["reason_rule_on"] == "TIER1_SV_DEPTH").sum()),
                        n_combined=int((grp["call_rule_on"] == "YES").sum()),
                        median_score=round(float(np.median(grp["score"])), 4),
                        frac_with_depth_evidence=round(float(grp["flags"].str.contains("SV_DEPTH_EVIDENCE").mean()), 4)))
    od = pd.DataFrame(out)
    od.to_csv(a.out, sep="\t", index=False)
    sys.stderr.write(od.to_string(index=False) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
