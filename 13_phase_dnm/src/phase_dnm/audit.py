"""Assert that what the configuration DECLARES is what the data actually CONTAINS (P30).

Every defect in this module so far has had one shape. Something was declared and never produced, and nothing checked:

    P17   SV interval evidence declared in the registry, never implemented          0 of 22434 rows
    TR    population columns declared, the annotate step never run for the class    0 of 406565 rows
    P24   site_qual admitted as a feature; positives never restricted to rare ones  69 % / 85 % common
    P29   planted candidates carried no caller block at all                         42 of 72 columns empty
    -     qd_child derived at load for training, not at scoring                     scored as NaN
    -     Geometry read `read_class`; the pipeline writes `status`                  crossover features empty
    -     m2_features_family.sb never passed --annot                                whole population block empty
    -     annotate never passed --strchive; the TR branch fell through              catalogue never consulted
    -     the spike site QC demanded a read spanning the event                      zero large deletions planted

None of these failed loudly. Several produced numbers that looked reasonable and were reported. The common cause is
not carelessness at any one site: it is that the pipeline had no place where a declaration is compared against the
artefact it describes. This module is that place. It is meant to be run before any number leaves the cluster, and to
exit non-zero when a claim is not backed by data.

Each check answers one question of the form "the config says X; does the data agree?" and reports the evidence either
way, so a green result is a measurement rather than an absence of complaints.
"""
from __future__ import annotations

import collections
import csv
import glob
import json
import os
from typing import Dict, List, Optional, Tuple

TAB = "\t"


class Finding:
    __slots__ = ("level", "check", "detail")

    def __init__(self, level: str, check: str, detail: str):
        self.level, self.check, self.detail = level, check, detail

    def __str__(self):
        return "%-5s %-26s %s" % (self.level, self.check, self.detail)


def _cols_and_counts(path: str, limit: Optional[int] = None) -> Tuple[List[str], Dict[str, collections.Counter], int]:
    """Column names, the value counter per column, and the row count. Values are capped so a wide table stays cheap."""
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter=TAB)
        cols = list(rd.fieldnames or [])
        counts = {c: collections.Counter() for c in cols}
        n = 0
        for r in rd:
            n += 1
            for c in cols:
                v = r.get(c)
                ctr = counts[c]
                if len(ctr) < 50 or v in ctr:
                    ctr[v] += 1
            if limit and n >= limit:
                break
    return cols, counts, n


EMPTYISH = (None, "", ".", "NA", "nan", "None")


def audit_features(evidence_dir: str, registry, class_groups=("snv_indel", "sv", "tr"),
                   max_children: int = 0) -> List[Finding]:
    """Every registry feature applicable to a class must be present, non-empty somewhere, and take more than one value.

    A constant column is reported separately from a missing one because they fail differently: a missing column is a
    wiring bug, while a constant column is usually a feature that cannot vary in this cohort and should be dropped or
    narrowed. Both end up discarded by the presence-leak guard, so neither is visible in a model's metrics."""
    out: List[Finding] = []
    for cls in class_groups:
        paths = sorted(glob.glob(os.path.join(evidence_dir, "*", "features", "*.%s.features.rf.tsv" % cls)))
        if max_children:
            paths = paths[:max_children]
        if not paths:
            out.append(Finding("FAIL", "features:%s" % cls, "no feature matrices found under %s" % evidence_dir))
            continue
        declared = sorted(f.name for f in registry.for_class_group(cls)) if hasattr(registry, "for_class_group") \
            else sorted({f.name for vc in _vclasses(cls) for f in registry.for_class(vc)})
        seen_cols: set = set()
        nonempty: collections.Counter = collections.Counter()
        values: Dict[str, set] = collections.defaultdict(set)
        rows = 0
        for p in paths:
            cols, counts, n = _cols_and_counts(p)
            rows += n
            seen_cols |= set(cols)
            for c in cols:
                for v, k in counts[c].items():
                    if v not in EMPTYISH:
                        nonempty[c] += k
                        if len(values[c]) < 8:
                            values[c].add(v)
        missing = [f for f in declared if f not in seen_cols]
        empty = [f for f in declared if f in seen_cols and nonempty[f] == 0]
        constant = [f for f in declared if f in seen_cols and nonempty[f] > 0 and len(values[f]) == 1]
        out.append(Finding("INFO", "features:%s" % cls,
                           "%d children, %d rows, %d declared features" % (len(paths), rows, len(declared))))
        if missing:
            out.append(Finding("FAIL", "features:%s" % cls, "declared but NOT A COLUMN: %s" % ", ".join(missing)))
        if empty:
            out.append(Finding("FAIL", "features:%s" % cls, "column present but empty in every row: %s" % ", ".join(empty)))
        if constant:
            out.append(Finding("WARN", "features:%s" % cls,
                               "constant in every row (the leak guard will discard these): %s"
                               % ", ".join("%s=%s" % (c, next(iter(values[c]))) for c in constant)))
        if not (missing or empty or constant):
            out.append(Finding("OK", "features:%s" % cls, "all %d declared features present and varying" % len(declared)))
    return out


def _vclasses(cls: str):
    return {"snv_indel": ("SNV", "INDEL"), "sv": ("SV",), "tr": ("TR",)}[cls]


def audit_model_vs_matrix(harness_dir: str, evidence_dir: str, class_groups=("snv_indel", "sv", "tr")) -> List[Finding]:
    """Every column a frozen model was trained on must be obtainable from the matrices it will score.

    This is the train/score skew check. qd_child is derived at load time, so it is legitimately absent from the file
    and present after the loader runs; a column that is absent BOTH ways is scored as NaN against a model that was
    fitted with it."""
    out: List[Finding] = []
    try:
        from .train import nested_cv as CV
    except ImportError as e:
        return [Finding("WARN", "model-vs-matrix", "cannot import the loader (%s)" % e)]
    for cls in class_groups:
        man = os.path.join(harness_dir, "models", "%s.training_manifest.json" % cls)
        if not os.path.exists(man):
            out.append(Finding("WARN", "model:%s" % cls, "no training manifest in %s" % harness_dir))
            continue
        cols = json.load(open(man)).get("cols") or json.load(open(man)).get("features") or []
        paths = sorted(glob.glob(os.path.join(evidence_dir, "*", "features", "*.%s.features.rf.tsv" % cls)))
        if not cols or not paths:
            out.append(Finding("WARN", "model:%s" % cls, "manifest has %d columns, %d matrices found" % (len(cols), len(paths))))
            continue
        df = CV._read_matrix(paths[0])
        missing = [c for c in cols if c not in df.columns]
        if missing:
            out.append(Finding("FAIL", "model:%s" % cls,
                               "trained on columns the scorer cannot supply, they become NaN: %s" % ", ".join(missing)))
        else:
            out.append(Finding("OK", "model:%s" % cls, "all %d trained columns obtainable after the loader" % len(cols)))
    return out


def audit_spike_plan(evidence_dir: str) -> List[Finding]:
    """Every subtype and length the planner was asked for must actually appear in the plan.

    The large-deletion experiment ran four times and planted no large deletion: the site QC rejected every candidate
    and the run still reported success, because nothing compared what was planned against what was placed."""
    out: List[Finding] = []
    plans = sorted(glob.glob(os.path.join(evidence_dir, "*", "spike", "*", "plan.tsv")))
    if not plans:
        return [Finding("WARN", "spike:plan", "no spike plans found")]
    placed: collections.Counter = collections.Counter()
    for p in plans:
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh, delimiter=TAB):
                placed[(r.get("variant_class"), r.get("subtype"), r.get("length"))] += 1
    try:
        from .sim.spike import SUBTYPES, LENGTHS
    except ImportError as e:
        return [Finding("WARN", "spike:plan", "cannot import the planner (%s)" % e)]
    out.append(Finding("INFO", "spike:plan", "%d plans, %d placed sites" % (len(plans), sum(placed.values()))))
    for cls, subs in SUBTYPES.items():
        for sub in subs:
            lengths = LENGTHS.get((cls, sub), (None,))
            for L in lengths:
                n = sum(v for (c, s, l), v in placed.items()
                        if c == cls and s == sub and (L is None or str(l) == str(L)))
                if n == 0:
                    out.append(Finding("FAIL", "spike:plan",
                                       "%s/%s length %s was requested and NEVER placed" % (cls, sub, L)))
    if not any(f.level == "FAIL" for f in out):
        out.append(Finding("OK", "spike:plan", "every requested subtype and length was placed"))
    return out


def audit_annot(evidence_dir: str, class_groups=("snv_indel", "sv", "tr")) -> List[Finding]:
    """An annot table that is full while the matching feature matrix column is empty means the join never happened."""
    out: List[Finding] = []
    for cls in class_groups:
        apaths = sorted(glob.glob(os.path.join(evidence_dir, "*", "annot", "*.%s.annot.tsv" % cls)))
        if not apaths:
            out.append(Finding("WARN", "annot:%s" % cls, "no annot tables found"))
            continue
        filled: collections.Counter = collections.Counter()
        rows = 0
        for p in apaths:
            _, counts, n = _cols_and_counts(p)
            rows += n
            for c, ctr in counts.items():
                filled[c] += sum(k for v, k in ctr.items() if v not in EMPTYISH)
        detail = ", ".join("%s %.0f%%" % (c, 100.0 * filled[c] / max(rows, 1))
                           for c in sorted(filled) if c != "variant_id")
        out.append(Finding("INFO", "annot:%s" % cls, "%d tables, %d rows: %s" % (len(apaths), rows, detail)))
    return out


def run(evidence_dir: str, registry, harness_dir: Optional[str] = None, spike: bool = True,
        class_groups=("snv_indel", "sv", "tr"), log=print) -> int:
    findings: List[Finding] = []
    findings += audit_features(evidence_dir, registry, class_groups)
    findings += audit_annot(evidence_dir, class_groups)
    if harness_dir:
        findings += audit_model_vs_matrix(harness_dir, evidence_dir, class_groups)
    if spike:
        findings += audit_spike_plan(evidence_dir)
    for f in findings:
        log(str(f))
    n_fail = sum(1 for f in findings if f.level == "FAIL")
    n_warn = sum(1 for f in findings if f.level == "WARN")
    log("audit: %d FAIL, %d WARN, %d checks" % (n_fail, n_warn, len(findings)))
    return 1 if n_fail else 0
