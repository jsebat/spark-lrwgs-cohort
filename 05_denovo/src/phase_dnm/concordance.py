"""P18 — concordance of the module's final calls with the original pipeline's de novo sets, per class.

Baselines (copied from the Lustre analysis directory to <baselines>/, 2026-09-13):
  denovo_tiered.tsv          small variants, the filtered set: family proband chrom pos ref alt gene ... tier   (1,332 rows)
  denovo_sv_all.bed          the pipeline's de novo SV list, no header: chrom start end family proband svtype svlen  (7,736 rows)
  denovo_sv_PRIORITIZED_v2.tsv   the prioritised coding SVs (reported alongside; 4 rows)
  tr_denovo_expansions.tsv   TR expansions: trid chrom pos end motifs cds_genes child family affected ...   (205 rows)
Matching: small variants by (proband, chrom, pos, ref, alt); SVs by (proband, svtype) with breakpoints within
`sv_bp_tol` bp or reciprocal overlap >= `sv_recip`; TRs by (child, trid).

For every class: concordant YES, original-only (with the module's phase_class / decision_reason so a demotion is
explained), module-only (by source_tier and mask_overlap so a rescue is qualified), and per-proband counts before /
after (the R9 guard: a rescue that moves ~38 per proband towards ~70 with the same parent-of-origin ratio is
credible; one that overshoots is not). Identifiers stay in the tables on the filer; the summary is count-only.
"""
from __future__ import annotations

import csv
import os
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Tuple


def norm_allele(pos: int, ref: str, alt: str) -> Tuple[int, str, str]:
    """Canonical key (position of the first changed base, deleted sequence, inserted sequence): trim the shared prefix
    (moving pos), then the shared suffix, allowing EMPTY cores - so `GAT>GT` at 100 and `AT>T` at 101 (padded vs
    anchor-after forms of the same 1-bp deletion) both become (101, "A", ""). The original pipeline matched after
    `bcftools norm`; the joint VCF's multi-allelic-split records keep their padding, which hid 120 indels on
    2026-09-13. Left-alignment through a repeat needs the reference and is covered by the position-window fallback."""
    ref, alt = ref.upper(), alt.upper()
    while ref and alt and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    return pos, ref, alt


def load_small(path: str) -> Dict[Tuple[str, str, int, str, str], dict]:
    """Keyed by the NORMALISED (proband, chrom, pos, ref, alt)."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                p, a, b = norm_allele(int(r["pos"]), r["ref"], r["alt"])
                out[(r["proband"], r["chrom"], p, a, b)] = r
            except (KeyError, ValueError):
                continue
    return out


def small_match(small: Dict, sid: str, chrom: str, pos: int, ref: str, alt: str, window: int = 20) -> Optional[Tuple]:
    """Exact normalised match, else (indels only) the same length change within `window` bp for the same proband."""
    p, a, b = norm_allele(pos, ref, alt)
    key = (sid, chrom, p, a, b)
    if key in small:
        return key
    if len(a) == len(b):
        return None
    delta = len(b) - len(a)
    for k in small:
        if k[0] != sid or k[1] != chrom:
            continue
        # the tiered table truncates alleles at 30 characters: a same-position indel whose baseline core is a long
        # prefix of ours is the same event (117 of the 1,332 originals are >= 30 bp insertions, 2026-09-13)
        if k[2] == p and ((len(k[4]) >= 20 and b.startswith(k[4]) and k[3] == a) or (len(k[3]) >= 20 and a.startswith(k[3]) and k[4] == b)):
            return k
        if abs(k[2] - p) <= window and (len(k[4]) - len(k[3])) == delta:
            return k
    return None


def load_sv_bed(path: str) -> Dict[str, List[dict]]:
    """proband -> [ {chrom, start, end, svtype, svlen} ] (bed: chrom start end family proband svtype svlen)."""
    out: Dict[str, List[dict]] = defaultdict(list)
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6 or f[0].startswith("#") or f[0] == "chrom":
                continue
            try:
                out[f[4]].append(dict(chrom=f[0], start=int(f[1]), end=int(f[2]), family=f[3], svtype=f[5], svlen=int(f[6]) if len(f) > 6 and f[6] not in ("", ".") else None))
            except ValueError:
                continue
    return out


def load_tr(path: str) -> Dict[Tuple[str, str], dict]:
    out = {}
    if not os.path.exists(path):
        return out
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if "trid" in r and "child" in r:
                out[(r["child"], r["trid"])] = r
    return out


def sv_match(row: dict, cands: List[dict], bp_tol: int = 500, recip: float = 0.5) -> Optional[dict]:
    chrom, s, e, t = row["chrom"], int(row["start"]), int(row["end"]), (row.get("svtype") or "").upper()
    best = None
    for c in cands:
        if c["chrom"] != chrom:
            continue
        ct = c["svtype"].upper()
        if t and ct and t != ct and not ({t, ct} <= {"DUP", "INS"}):
            continue
        if abs(c["start"] - s) <= bp_tol and abs(c["end"] - e) <= bp_tol:
            return c
        if t in ("DEL", "DUP", "INV") and e > s and c["end"] > c["start"]:
            ov = min(e, c["end"]) - max(s, c["start"])
            if ov > 0 and ov / (e - s) >= recip and ov / (c["end"] - c["start"]) >= recip:
                best = c
        elif abs(c["start"] - s) <= bp_tol:
            best = best or c
    return best


def concordance(final_rows: Iterable[dict], vclass_group: str, baselines_dir: str, sv_bp_tol: int = 500, sv_recip: float = 0.5
                ) -> Tuple[List[dict], List[dict], Dict[str, object]]:
    """Returns (per-row concordance table, per-proband table, summary)."""
    rows = list(final_rows)
    small = load_small(os.path.join(baselines_dir, "denovo_tiered.tsv")) if vclass_group == "snv_indel" else {}
    svb = load_sv_bed(os.path.join(baselines_dir, "denovo_sv_all.bed")) if vclass_group == "sv" else {}
    trb = load_tr(os.path.join(baselines_dir, "tr_denovo_expansions.tsv")) if vclass_group == "tr" else {}
    matched_baseline: set = set()
    per_row = []
    for r in rows:
        sid = r["sample_id"]
        hit = None
        if vclass_group == "snv_indel":
            key = small_match(small, sid, r["chrom"], int(r["start"]), r["ref"], r["alt"])
            hit = small.get(key) if key is not None else None
            if hit is not None:
                matched_baseline.add(key)
        elif vclass_group == "sv":
            hit = sv_match(r, svb.get(sid, []), sv_bp_tol, sv_recip)
            if hit is not None:
                matched_baseline.add((sid, hit["chrom"], hit["start"], hit["end"], hit["svtype"]))
        else:
            key = (sid, r.get("trid") or r.get("variant_id"))
            hit = trb.get(key)
            if hit is not None:
                matched_baseline.add(key)
        yes = r.get("dnm_call") == "YES"
        cand = r.get("dnm_call") == "CANDIDATE"          # tier 2 (0.3.0): an original call recovered at tier 2 is neither concordant nor lost
        status = ("concordant_YES" if yes and hit is not None else "original_tier2" if cand and hit is not None else
                  "original_only" if hit is not None else "module_only" if yes else "module_tier2" if cand else "neither")
        per_row.append(dict(sample_id=sid, variant_id=r.get("variant_id"), variant_class=r.get("variant_class"), status=status,
                            dnm_call=r.get("dnm_call"), dnm_tier=r.get("dnm_tier", ""), phase_class=r.get("phase_class"), decision_reason=r.get("decision_reason"),
                            parent_of_origin=r.get("parent_of_origin"), source_tier=r.get("source_tier"), mask_overlap=r.get("mask_overlap"),
                            baseline_tier=(hit or {}).get("tier", "") if vclass_group == "snv_indel" else ""))
    # baseline rows the module never saw as a candidate at all (not in the unfiltered set)
    if vclass_group == "snv_indel":
        unseen = [k for k in small if k not in matched_baseline and k[0] in {r["sample_id"] for r in rows}]
    elif vclass_group == "sv":
        unseen = [(sid, c["chrom"], c["start"], c["end"], c["svtype"]) for sid, cs in svb.items() if sid in {r["sample_id"] for r in rows}
                  for c in cs if (sid, c["chrom"], c["start"], c["end"], c["svtype"]) not in matched_baseline]
    else:
        unseen = [k for k in trb if k not in matched_baseline and k[0] in {r["sample_id"] for r in rows}]
    # per proband
    pp: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for d in per_row:
        pp[d["sample_id"]][d["status"]] += 1
        if d["status"] in ("concordant_YES", "original_only", "original_tier2"):
            pp[d["sample_id"]]["original_seen"] += 1
        if d["dnm_call"] == "YES":
            pp[d["sample_id"]]["module_YES"] += 1
        if d["dnm_call"] == "CANDIDATE":
            pp[d["sample_id"]]["module_CANDIDATE"] += 1   # counter, distinct from the status module_tier2
        if d["dnm_call"] == "YES" and d["parent_of_origin"] == "paternal":
            pp[d["sample_id"]]["module_YES_paternal"] += 1
    for k in unseen:
        pp[k[0]]["original_unseen"] += 1
    per_proband = [dict(sample_id=s, **{k: v for k, v in sorted(c.items())}) for s, c in sorted(pp.items())]
    n = lambda st: sum(1 for d in per_row if d["status"] == st)
    by_class_orig_only: Dict[str, int] = defaultdict(int)
    by_reason_orig_only: Dict[str, int] = defaultdict(int)
    by_tier_mod_only: Dict[str, int] = defaultdict(int)
    mask_mod_only = 0
    yes_pat = yes_mat = 0
    for d in per_row:
        if d["status"] == "original_only":
            by_class_orig_only[d["phase_class"] or "."] += 1
            by_reason_orig_only[d["decision_reason"] or "."] += 1
        if d["status"] == "module_only":
            by_tier_mod_only[d["source_tier"] or "."] += 1
            mask_mod_only += int(str(d["mask_overlap"]) == "1")
        if d["dnm_call"] == "YES":
            yes_pat += d["parent_of_origin"] == "paternal"; yes_mat += d["parent_of_origin"] == "maternal"
    summary = dict(vclass_group=vclass_group, n_rows=len(rows), n_probands=len(pp), concordant_YES=n("concordant_YES"),
                   original_only=n("original_only"), module_only=n("module_only"), original_unseen_as_candidate=len(unseen),
                   original_tier2=n("original_tier2"), module_tier2=n("module_tier2"),
                   per_proband_module_tier2_median=_median([c.get("module_CANDIDATE", 0) for c in pp.values()]),
                   original_only_by_phase_class=dict(by_class_orig_only), original_only_by_reason=dict(by_reason_orig_only),
                   module_only_by_tier=dict(by_tier_mod_only), module_only_in_mask=mask_mod_only,
                   module_YES_paternal_fraction=round(yes_pat / (yes_pat + yes_mat), 3) if (yes_pat + yes_mat) else None,
                   per_proband_original_median=_median([c.get("original_seen", 0) + c.get("original_unseen", 0) for c in pp.values()]),
                   per_proband_module_YES_median=_median([c.get("module_YES", 0) for c in pp.values()]))
    return per_row, per_proband, summary


def _median(xs: List[float]) -> Optional[float]:
    xs = sorted(xs)
    if not xs:
        return None
    m = len(xs) // 2
    return xs[m] if len(xs) % 2 else (xs[m - 1] + xs[m]) / 2
