"""P9 — likelihood layer: a posterior over hypotheses from the per-haplotype read counts.

Rows (alt, ref): A = the child's alt-carrying haplotype, O = the child's other haplotype, T / U = transmitted /
untransmitted haplotype of the parent of origin, N1 / N2 = the two haplotypes of the non-transmitting parent,
UNT = the child's untagged reads (either haplotype). Each row is binomial with a hypothesis-specific expected alt
fraction; nuisance fractions are integrated on a grid over uniform priors; ε is the per-read false-alt rate
(estimated per class from the cohort, default in thresholds.yaml); δ is the fraction of reads that carry the true
allele but land on the wrong haplotype row (mis-tagging), so a haplotype that "should" be pure ref shows alt at
ε + δ and one that should be pure alt shows alt at 1 − ε − δ.

Hypotheses (expected alt fraction per row):
  germline    A: 1-ε-δ   O: ε+δ   T, U, N1, N2: ε      UNT: 1/2
  inherited   A: 1-ε-δ   O: ε+δ   T: 1-ε-δ   U: ε+δ   N1, N2: ε   UNT: 1/2
  parental_mosaic  as germline but T: m,  m ~ U(0.02, 0.5)
  child_mosaic     A: c (c ~ U(0.05, 0.8))   O: ε+δ   parents: ε   UNT: c/2
  artefact    every row: a,  a ~ U(0, 1)   (a shared systematic alt fraction on all haplotypes)
When T/U are unresolved the parent-of-origin rows are treated as an unordered pair (both assignments averaged);
when the parent of origin is undetermined the inherited / parental-mosaic hypotheses consider both parents.

Outputs: posterior per hypothesis, log10 likelihood ratio germline vs best alternative, and `phase_score` =
posterior of germline. The rule layer (P8) and this layer are always both reported (P9).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

HYPS = ("germline", "inherited", "parental_mosaic", "child_mosaic", "artefact")


@dataclass
class LikParams:
    eps: float = 0.01                # per-read false-alt rate on a haplotype that carries no alt
    delta: float = 0.03              # mis-tagging: true-allele reads landing on the wrong haplotype row
    grid: int = 40
    prior: Dict[str, float] = None   # type: ignore[assignment]

    def __post_init__(self):
        if self.prior is None:
            # per-candidate priors on the RAW set: ~70 true DNMs among ~32k small-variant candidates per child;
            # inherited-and-missed and artefacts dominate; transmitted parental mosaics and postzygotic mosaics are
            # each a minority (~5-10%) of apparent DNMs (order of magnitude; unverified figures). At ~10 reads per
            # haplotype a low-fraction parental mosaic on T is nearly indistinguishable from germline, so the prior
            # ratio matters there - report sensitivity (R7, P10).
            self.prior = {"germline": 0.003, "inherited": 0.10, "parental_mosaic": 0.0003, "child_mosaic": 0.0005, "artefact": 0.8962}


Counts = Tuple[int, int]   # (alt, ref)


def _logbinom(k: int, n: int, p: float) -> float:
    if n == 0:
        return 0.0
    p = min(max(p, 1e-9), 1 - 1e-9)
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
            + k * math.log(p) + (n - k) * math.log(1 - p))


def _lse(xs: List[float]) -> float:
    m = max(xs)
    return m + math.log(sum(math.exp(x - m) for x in xs)) if xs else float("-inf")


def _row(c: Optional[Counts], p: float) -> float:
    if c is None:
        return 0.0
    return _logbinom(c[0], c[0] + c[1], p)


def loglik(hyp: str, rows: Dict[str, Optional[Counts]], prm: LikParams, po_resolved: bool) -> float:
    """rows: A, O, T, U, N1, N2, UNT (each (alt, ref) or None). po_resolved: T/U labels are real (else T/U is an
    unordered pair of the origin parent's haplotypes)."""
    e, d = prm.eps, prm.delta
    hi, lo = 1 - e - d, e + d
    A, O, T, U, N1, N2, UNT = (rows.get(k) for k in ("A", "O", "T", "U", "N1", "N2", "UNT"))
    base_child = _row(A, hi) + _row(O, lo) + _row(UNT, 0.5)
    parents_ref = _row(N1, e) + _row(N2, e)

    def tu(pt: float, pu: float) -> float:
        if po_resolved:
            return _row(T, pt) + _row(U, pu)
        return _lse([_row(T, pt) + _row(U, pu), _row(T, pu) + _row(U, pt)]) - math.log(2)

    if hyp == "germline":
        return base_child + tu(e, e) + parents_ref
    if hyp == "inherited":
        return base_child + tu(hi, lo) + parents_ref
    if hyp == "parental_mosaic":
        ms = [0.02 + (0.5 - 0.02) * (i + 0.5) / prm.grid for i in range(prm.grid)]
        return base_child + parents_ref + _lse([tu(m, e) for m in ms]) - math.log(prm.grid)
    if hyp == "child_mosaic":
        cs = [0.05 + (0.8 - 0.05) * (i + 0.5) / prm.grid for i in range(prm.grid)]
        return _lse([_row(A, c) + _row(O, lo) + _row(UNT, c / 2) + tu(e, e) + parents_ref for c in cs]) - math.log(prm.grid)
    if hyp == "artefact":
        as_ = [(i + 0.5) / prm.grid for i in range(prm.grid)]
        return _lse([sum(_row(r, a) for r in (A, O, T, U, N1, N2, UNT)) for a in as_]) - math.log(prm.grid)
    raise ValueError(hyp)


def posterior(rows: Dict[str, Optional[Counts]], prm: LikParams, po_resolved: bool) -> Dict[str, object]:
    ll = {h: loglik(h, rows, prm, po_resolved) for h in HYPS}
    lp = {h: ll[h] + math.log(prm.prior[h]) for h in HYPS}
    z = _lse(list(lp.values()))
    post = {h: math.exp(lp[h] - z) for h in HYPS}
    best_alt = max((h for h in HYPS if h != "germline"), key=lambda h: ll[h])
    out: Dict[str, object] = {"lik_post_" + h: round(post[h], 6) for h in HYPS}
    out["lik_best_alternative"] = best_alt
    out["lik_log10lr_germline"] = round((ll["germline"] - ll[best_alt]) / math.log(10), 3)
    out["phase_score"] = round(post["germline"], 6)
    return out


# ----------------------------------------------------------------------------------------------
# from an evidence row (review.py columns) — so the layer can run as a post-step over existing tables
# ----------------------------------------------------------------------------------------------
def _c(row: dict, pre: str) -> Optional[Counts]:
    try:
        a, r = int(float(row[pre + "_alt"])), int(float(row[pre + "_ref"]))
    except (KeyError, ValueError, TypeError):
        return None
    return (a, r)


def rows_from_evidence(row: dict) -> Tuple[Dict[str, Optional[Counts]], bool]:
    c1, c2 = _c(row, "C1"), _c(row, "C2")
    if c1 is None or c2 is None:
        return {}, False
    A, O = (c1, c2) if c1[0] >= c2[0] else (c2, c1)
    hap1_is = row.get("child_hap1_is", ".")
    a_is_hap1 = c1[0] >= c2[0]
    # parent of origin of A from the child's orientation
    if hap1_is in ("P", "M"):
        origin = ("F" if hap1_is == "P" else "M") if a_is_hap1 else ("M" if hap1_is == "P" else "F")
    else:
        origin = None
    unt = None
    try:
        u_alt = int(float(row["C_untagged_alt"]))
        # readable REF among untagged reads when the table carries it (2026-09-16 on); dp - alt before that counted
        # unreadable AMB reads as REF, unlike the tagged rows (D13)
        u_ref = str(row.get("C_untagged_ref") or "")
        unt = (u_alt, int(float(u_ref)) if u_ref not in ("", ".") else int(float(row["C_untagged_dp"])) - u_alt)
    except (KeyError, ValueError, TypeError):
        pass
    rows: Dict[str, Optional[Counts]] = {"A": A, "O": O, "UNT": unt}
    if origin is None:
        # undetermined: treat the parent with more alt as "origin" for the inherited/mosaic hypotheses (unordered)
        f1, f2, m1, m2 = _c(row, "F1"), _c(row, "F2"), _c(row, "M1"), _c(row, "M2")
        fa = sum(x[0] for x in (f1, f2) if x); ma = sum(x[0] for x in (m1, m2) if x)
        origin = "F" if fa >= ma else "M"
    other = "M" if origin == "F" else "F"
    tr = row.get(origin + "_transmitted_hap", ".")
    po_resolved = tr in ("1", "2", 1, 2)
    if po_resolved:
        t, u = int(tr), 3 - int(tr)
        rows["T"], rows["U"] = _c(row, "%s%d" % (origin, t)), _c(row, "%s%d" % (origin, u))
    else:
        rows["T"], rows["U"] = _c(row, origin + "1"), _c(row, origin + "2")
    rows["N1"], rows["N2"] = _c(row, other + "1"), _c(row, other + "2")
    return rows, po_resolved


def score_evidence_row(row: dict, prm: LikParams) -> Dict[str, object]:
    rows, po = rows_from_evidence(row)
    if not rows:
        return {"lik_post_" + h: None for h in HYPS} | {"lik_best_alternative": None, "lik_log10lr_germline": None, "phase_score": None}
    return posterior(rows, prm, po)
