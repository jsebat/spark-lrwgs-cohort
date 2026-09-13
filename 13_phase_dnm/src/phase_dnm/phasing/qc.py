"""M1d — phase QC per child and per cohort (README §3.1), from the M1 outputs already on disk.

Per child: orientation summary (M1a), transmission summary and resolved change points (M1b/M1b2), per-sample
haplotype depth for child and both parents (hapdepth), optional HiPhase per-chromosome stats. Every quantity is
reported; gates from thresholds.yaml turn a subset into PASS / flag strings. Missing inputs are tolerated and named
in `missing`, so the QC can run on a family before its read-level steps have finished.
"""
from __future__ import annotations

import csv
import glob
import json
import os
import statistics as st
from typing import Dict, List, Optional

AUTOSOMES = {"chr%d" % i for i in range(1, 23)}


def _load_json(path: str) -> Optional[dict]:
    return json.load(open(path)) if os.path.exists(path) else None


def _hapdepth_means(path: str) -> Optional[dict]:
    s = _load_json(path)
    if not s:
        return None
    pc = {c: v for c, v in s["per_chrom"].items() if c in AUTOSOMES}
    if not pc:
        return None
    W = sum(v["n_bins"] for v in pc.values())
    mean = lambda k: sum(v[k] * v["n_bins"] for v in pc.values()) / W
    h1, h2, un, lo = mean("mean_dp_hap1"), mean("mean_dp_hap2"), mean("mean_dp_untagged"), mean("mean_dp_lowmapq")
    x = s["per_chrom"].get("chrX")
    xt = (x["mean_dp_hap1"] + x["mean_dp_hap2"] + x["mean_dp_untagged"] + x["mean_dp_lowmapq"]) if x else None
    tot = h1 + h2 + un + lo
    return dict(dp_hap1=round(h1, 2), dp_hap2=round(h2, 2), dp_untagged=round(un, 2), dp_lowmapq=round(lo, 2),
                dp_total=round(tot, 2), tagged_frac=round((h1 + h2) / max(1e-9, h1 + h2 + un), 3),
                hap_balance=round(min(h1, h2) / max(h1, h2, 1e-9), 3),
                chrx_autosome_ratio=round(xt / tot, 3) if (xt is not None and tot) else None)


def _changepoints(path: str) -> Optional[dict]:
    if not os.path.exists(path):
        return None
    out: Dict[str, Dict[str, int]] = {"F": {}, "M": {}}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            d = out[r["parent"]]
            d["candidates"] = d.get("candidates", 0) + 1
            d[r["status"]] = d.get(r["status"], 0) + 1
            if r["status"] == "CROSSOVER" and r.get("child_switch_in_interval", "N") != "Y":
                d["crossover_excl_child_switch"] = d.get("crossover_excl_child_switch", 0) + 1
    return out


def child_qc(family_dir: str, child: str, father: str, mother: str, hapdepth_dir: Optional[str], gates: dict,
             thresholds_version: Optional[str] = None) -> dict:
    stem = os.path.join(family_dir, child)
    missing: List[str] = []
    q: dict = {"child": child, "father": father, "mother": mother, "thresholds_version": thresholds_version}
    o = _load_json(stem + ".orientation.summary.json")
    if o:
        t = o["total"]
        q["orientation"] = dict(sex=o.get("sex"), n_blocks=o["n_blocks"], n_blocks_split=o["n_blocks_split"],
                                n_switches_located=t.get("n_switches_located", 0), n_informative=t.get("n_informative", 0),
                                mendel_inconsistent_per_informative=o["mendel_inconsistent_per_informative"],
                                frac_het_ambiguous=o.get("frac_het_ambiguous"), frac_het_mixed_votes=o.get("frac_het_mixed_votes"),
                                frac_bp_ambiguous=o.get("frac_bp_ambiguous"), frac_bp_mixed_votes=o.get("frac_bp_mixed_votes"))
    else:
        missing.append("orientation")
    tr = _load_json(stem + ".transmission.summary.json")
    if tr:
        q["transmission"] = {p: dict(n_blocks=v.get("n_blocks", 0), n_change_points=v.get("n_change_points", 0),
                                     frac_het_resolved=v.get("frac_het_resolved"), frac_bp_resolved=v.get("frac_bp_resolved"),
                                     mendel_inconsistent_per_informative=v.get("mendel_inconsistent_per_informative"))
                             for p, v in tr["per_parent"].items()}
    else:
        missing.append("transmission")
    cp = _changepoints(stem + ".changepoints.resolved.tsv")
    if cp:
        q["change_points"] = cp
    else:
        missing.append("changepoints_resolved")
    if hapdepth_dir:
        q["depth"] = {}
        for role, sample in (("child", child), ("father", father), ("mother", mother)):
            d = _hapdepth_means(os.path.join(hapdepth_dir, sample + ".hapdepth.summary.json"))
            if d:
                q["depth"][role] = d
            else:
                missing.append("hapdepth_" + role)
    else:
        missing.append("hapdepth")
    q["missing"] = missing
    q["flags"] = apply_gates(q, gates)
    q["qc"] = "PASS" if not q["flags"] else ";".join(q["flags"])
    return q


def apply_gates(q: dict, g: dict) -> List[str]:
    flags: List[str] = []
    o = q.get("orientation")
    if o:
        if o.get("frac_het_ambiguous") is not None and o["frac_het_ambiguous"] > g.get("max_frac_het_ambiguous", 0.10):
            flags.append("ORIENT_AMBIGUOUS_HET")
        if o.get("frac_het_mixed_votes") is not None and o["frac_het_mixed_votes"] > g.get("max_frac_het_mixed_votes", 0.005):
            flags.append("ORIENT_MIXED_VOTES")
        if o.get("mendel_inconsistent_per_informative") is not None and \
                o["mendel_inconsistent_per_informative"] > g.get("max_mendel_inconsistent_per_informative", 0.01):
            flags.append("MENDEL")
    for p, v in (q.get("transmission") or {}).items():
        if v.get("frac_het_resolved") is not None and v["frac_het_resolved"] < g.get("min_frac_het_transmission_resolved", 0.90):
            flags.append("TRANSMISSION_UNRESOLVED_" + p)
    lo, hi = g.get("crossovers_per_meiosis", [15, 60])
    for p, v in (q.get("change_points") or {}).items():
        x = v.get("crossover_excl_child_switch", v.get("CROSSOVER", 0))
        if v.get("candidates") and not (lo <= x <= hi):
            flags.append("CROSSOVERS_%s=%d" % (p, x))
    for role, d in (q.get("depth") or {}).items():
        if d["dp_total"] < g.get("min_total_depth", 12):
            flags.append("LOW_DEPTH_%s=%.1f" % (role, d["dp_total"]))
        if d["tagged_frac"] < g.get("min_tagged_frac", 0.75):
            flags.append("LOW_TAGGED_%s" % role)
        if d["hap_balance"] < g.get("min_hap_balance", 0.9):
            flags.append("HAP_IMBALANCE_%s" % role)
    sex = (q.get("orientation") or {}).get("sex")
    ratio = ((q.get("depth") or {}).get("child") or {}).get("chrx_autosome_ratio")
    if sex in ("M", "F") and ratio is not None:
        if (sex == "M" and ratio > 0.75) or (sex == "F" and ratio < 0.85):
            flags.append("SEX_DEPTH_MISMATCH")
    return flags


COHORT_COLUMNS = ["family", "child", "sex", "qc", "flags", "n_blocks", "n_switches_located", "frac_het_ambiguous",
                  "frac_het_mixed_votes", "mendel_per_informative", "F_frac_het_resolved", "M_frac_het_resolved",
                  "F_crossovers", "M_crossovers", "F_crossovers_excl_child_switch", "M_crossovers_excl_child_switch",
                  "F_switch_errors", "M_switch_errors", "child_dp_total", "child_dp_hap1", "child_dp_hap2", "child_tagged_frac",
                  "father_dp_total", "mother_dp_total", "child_chrx_ratio", "missing"]


def cohort_row(family: str, q: dict) -> dict:
    o = q.get("orientation") or {}
    tr = q.get("transmission") or {}
    cp = q.get("change_points") or {}
    dp = q.get("depth") or {}
    g = lambda d, k: d.get(k) if d else None
    return {
        "family": family, "child": q["child"], "sex": o.get("sex"), "qc": q["qc"], "flags": ";".join(q["flags"]),
        "n_blocks": o.get("n_blocks"), "n_switches_located": o.get("n_switches_located"),
        "frac_het_ambiguous": o.get("frac_het_ambiguous"), "frac_het_mixed_votes": o.get("frac_het_mixed_votes"),
        "mendel_per_informative": o.get("mendel_inconsistent_per_informative"),
        "F_frac_het_resolved": g(tr.get("F"), "frac_het_resolved"), "M_frac_het_resolved": g(tr.get("M"), "frac_het_resolved"),
        "F_crossovers": g(cp.get("F"), "CROSSOVER"), "M_crossovers": g(cp.get("M"), "CROSSOVER"),
        "F_crossovers_excl_child_switch": g(cp.get("F"), "crossover_excl_child_switch"),
        "M_crossovers_excl_child_switch": g(cp.get("M"), "crossover_excl_child_switch"),
        "F_switch_errors": g(cp.get("F"), "SWITCH_ERROR"), "M_switch_errors": g(cp.get("M"), "SWITCH_ERROR"),
        "child_dp_total": g(dp.get("child"), "dp_total"), "child_dp_hap1": g(dp.get("child"), "dp_hap1"),
        "child_dp_hap2": g(dp.get("child"), "dp_hap2"), "child_tagged_frac": g(dp.get("child"), "tagged_frac"),
        "father_dp_total": g(dp.get("father"), "dp_total"), "mother_dp_total": g(dp.get("mother"), "dp_total"),
        "child_chrx_ratio": g(dp.get("child"), "chrx_autosome_ratio"), "missing": ";".join(q["missing"]),
    }


def write_cohort_table(rows: List[dict], path: str):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COHORT_COLUMNS, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else r[k]) for k in COHORT_COLUMNS})


def cohort_summary(rows: List[dict]) -> dict:
    def dist(k):
        v = [float(r[k]) for r in rows if r.get(k) not in (None, "")]
        return dict(n=len(v), median=round(st.median(v), 4), min=round(min(v), 4), max=round(max(v), 4)) if v else None
    return {"n_children": len(rows), "n_pass": sum(1 for r in rows if r["qc"] == "PASS"),
            "flags": sorted({f for r in rows for f in r["flags"].split(";") if f}),
            "distributions": {k: dist(k) for k in COHORT_COLUMNS if k not in ("family", "child", "sex", "qc", "flags", "missing")}}
