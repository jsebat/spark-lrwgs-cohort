"""The two arms of the targeted clinical workflow (DESIGN P31-P33), built from 05_denovo's tables.

  arm A  rare + large SVs: selected by SIZE and RARITY, quality first, then clinical annotation.
         Its quality verdicts are an unbiased sample and may serve as validation truth for the evidence.
  arm B  clinically led: selected by GENE PANEL and predicted IMPACT (any dnm_call, including below threshold),
         clinical relevance first, then quality. Its verdicts describe the panel, not the caller.

Both arms read every candidate 05_denovo scored, not just the called ones: DNMT3A is a 302 bp coding deletion the
genome-wide caller leaves BELOW_TAU (rule 5/6, rf 0.72), and a workflow restricted to called variants would not see
it. Every row carries discovery_mode = targeted_clinical and never enters a recall, FDR or rate estimate (P31).

Two verdicts per row, recorded independently: quality (is it real and de novo) and clinical (does it matter).
"""
from __future__ import annotations

import csv
import glob
import os
from typing import Dict, Iterable, List, Optional, Set, Tuple

from . import quality as Q
from .genes import GeneModel, Panel, shet_band

TAB = chr(9)
DISCOVERY_MODE = "targeted_clinical"

# SV loss-of-function rule (JS-agreed, from denovo_sv_priority.sb):
#   DEL overlapping CDS -> coding LoF;  INS with a breakpoint inside CDS -> coding LoF;
#   DUP overlapping CDS -> coding, flagged (a whole-gene DUP is a gain, not LoF);
#   INV with a breakpoint inside the gene -> coding;  BND -> excluded (not validatable).
# UTR overlaps are kept and labelled: a deletion of a 5' UTR or of a last exon's 3' UTR can abolish expression
# without touching a coding base, and a CDS-only rule cannot see it (JS, 2026-09-16).


def sv_impact(svtype: str, kinds: Set[str]) -> str:
    coding = "CDS" in kinds
    utr = bool(kinds & {"UTR5", "UTR3"})
    if svtype == "BND":
        return "excluded_BND"
    if svtype in ("DEL", "CNV"):
        return "coding_LoF" if coding else ("UTR_loss" if utr else "intronic_or_exon_noncoding")
    if svtype == "INS":
        return "coding_LoF_insertion" if coding else ("UTR_insertion" if utr else "intronic_or_exon_noncoding")
    if svtype == "DUP":
        return "coding_duplication(gain,flag)" if coding else ("UTR_duplication" if utr else "intronic_or_exon_noncoding")
    if svtype == "INV":
        return "coding_inversion" if coding else ("UTR_inversion" if utr else "intronic_or_exon_noncoding")
    return "other"


def _f(x) -> Optional[float]:
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _i(x) -> int:
    v = _f(x)
    return int(v) if v is not None else 0


def read_tsv(path: str) -> List[Dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter=TAB))


def load_final(final_dir: str, cls: str) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for p in sorted(glob.glob(os.path.join(final_dir, "*.%s.dnm.tsv" % cls))):
        rows.extend(read_tsv(p))
    return rows


def load_sv_evidence(evidence_dir: str) -> Dict[str, Dict[str, str]]:
    """variant_id -> raw six-haplotype SV evidence row (the final table keeps only the derived features)."""
    out: Dict[str, Dict[str, str]] = {}
    for p in sorted(glob.glob(os.path.join(evidence_dir, "*", "evidence", "*.sv.evidence.tsv"))):
        for r in read_tsv(p):
            out[r["variant_id"]] = r
    return out


def merge_sv(final_row: Dict[str, str], ev: Dict[str, Dict[str, str]]) -> Dict[str, object]:
    r: Dict[str, object] = dict(ev.get(final_row.get("variant_id", ""), {}))
    r.update(final_row)          # the final table's decision columns win where both carry a value
    return r


# ----------------------------------------------------------------------------------------------
# shared row construction
# ----------------------------------------------------------------------------------------------
BASE_COLS = ["arm", "discovery_mode", "family_id", "sample_id", "variant_class", "chrom", "start", "end", "ref", "alt",
             "svtype", "svlen", "variant_id",
             "gene", "impact", "gene_sets", "panel_rank", "s_het", "s_het_band", "clinical_relevance",
             "quality_verdict", "quality_reasons", "quality_criteria",
             "dnm_call", "dnm_tier", "decision_reason", "rf_prob", "phase_class", "rule_score", "parent_of_origin",
             "cohort_AC_loo", "pon_founder_recurrence_loo", "gnomad_af", "mask_overlap",
             "sv_depth_ratio_inside_flank", "sv_het_snv_persistence", "junction_both_ends", "child_sv_support",
             "F_sv_ratio_all", "M_sv_ratio_all", "t_alt_reads", "p_max_alt_any_hap", "c_alt_hapA", "child_DP",
             "consequence", "lof", "lof_filter", "missense_tier", "missense_nflag", "notes"]


def clinical_relevance(gene: str, impact: str, panel: Panel, shet: Optional[float]) -> str:
    """A label for HOW the variant matters, never a gate on quality. Panel membership is the criterion (P31/P33)."""
    rank = panel.rank(gene) if gene else 9
    hi_impact = impact.startswith(("coding_LoF", "LoF_HC", "stop_gained", "frameshift", "splice_donor", "splice_acceptor",
                                   "start_lost", "coding_inversion", "coding_dup"))
    mid_impact = impact.startswith(("missense", "UTR_", "coding_duplication", "inframe", "splice_region", "TR_"))
    if rank <= 2 and hi_impact:
        return "HIGH: LoF-type variant in a high-confidence panel gene"
    if rank <= 2 and mid_impact:
        return "MODERATE: non-truncating variant in a high-confidence panel gene"
    if rank <= 4 and hi_impact:
        return "MODERATE: LoF-type variant in a panel gene"
    if rank <= 4:
        return "LOW: panel gene, uncertain impact"
    if shet is not None and shet >= 0.18 and hi_impact:
        return "MODERATE: LoF-type variant in a constrained gene off-panel"
    return "NONE: off-panel"


def base_row(arm: str, r: Dict[str, object], gene: str, impact: str, panel: Panel, shet: Dict[str, float],
             verdict: Q.Verdict, extra: Optional[Dict[str, object]] = None) -> Dict[str, object]:
    sh = shet.get(gene) if gene else None
    row: Dict[str, object] = {k: r.get(k, "") for k in BASE_COLS}
    row.update(arm=arm, discovery_mode=DISCOVERY_MODE, gene=gene, impact=impact,
               gene_sets=panel.label(gene) if gene else "-", panel_rank=panel.rank(gene) if gene else 9,
               s_het=("" if sh is None else "%.4g" % sh), s_het_band=shet_band(sh),
               clinical_relevance=clinical_relevance(gene, impact, panel, sh),
               junction_both_ends=r.get("C_sv_junc_both_ends", r.get("c_sv_junc_both_ends", "")))
    row.update(verdict.as_row("quality"))
    if extra:
        row.update(extra)
    return row


# ----------------------------------------------------------------------------------------------
# arm A: rare + large SVs
# ----------------------------------------------------------------------------------------------
def arm_a(sv_rows: Iterable[Dict[str, str]], ev: Dict[str, Dict[str, str]], gm: GeneModel, panel: Panel,
          shet: Dict[str, float], min_svlen: int = 5000) -> List[Dict[str, object]]:
    """Selected by size and rarity only. Quality first; clinical annotation is attached afterwards and does not
    affect membership. Inherited/artefact phase classes are kept and FAIL, so the set stays an unbiased sample of
    what the caller saw."""
    out = []
    for fr in sv_rows:
        svtype = fr.get("svtype", "")
        length = abs(_i(fr.get("svlen"))) or max(0, _i(fr.get("end")) - _i(fr.get("start")))
        if length < min_svlen or svtype == "BND":
            continue
        if _i(fr.get("cohort_AC_loo")) > 0 or _i(fr.get("pon_founder_recurrence_loo")) > 0:
            continue
        r = merge_sv(fr, ev)
        v = Q.sv_quality(r)
        s0, e0 = _i(r.get("start")) - 1, max(_i(r.get("end")), _i(r.get("start")) + 1)
        hits = gm.hits(r["chrom"], s0, e0) if svtype != "INS" else gm.hits(r["chrom"], s0, s0 + 1)
        # one row per variant; the gene column lists panel genes first, then the rest, each with its impact
        genes = sorted(hits, key=lambda g: (panel.rank(g), g))
        gene = genes[0] if genes else ""
        impact = sv_impact(svtype, hits[gene]) if gene else "intergenic"
        extra = {"notes": "; ".join("%s:%s" % (g, sv_impact(svtype, hits[g])) for g in genes[1:8])}
        out.append(base_row("A_rare_large_sv", r, gene, impact, panel, shet, v, extra))
    return out


# ----------------------------------------------------------------------------------------------
# arm B: clinically led
# ----------------------------------------------------------------------------------------------
def arm_b_sv(sv_rows: Iterable[Dict[str, str]], ev: Dict[str, Dict[str, str]], gm: GeneModel, panel: Panel,
             shet: Dict[str, float]) -> List[Dict[str, object]]:
    """Every SV candidate, any dnm_call, whose interval (or insertion point) touches an exon or UTR of a panel gene."""
    out = []
    union = panel.union
    for fr in sv_rows:
        svtype = fr.get("svtype", "")
        if svtype == "BND":
            continue
        s0, e0 = _i(fr.get("start")) - 1, max(_i(fr.get("end")), _i(fr.get("start")) + 1)
        hits = gm.hits(fr["chrom"], s0, e0) if svtype != "INS" else gm.hits(fr["chrom"], s0, s0 + 1)
        pg = [g for g in hits if g in union]
        if not pg:
            continue
        r = merge_sv(fr, ev)
        v = Q.sv_quality(r)
        for g in sorted(pg, key=lambda g: (panel.rank(g), g)):
            imp = sv_impact(svtype, hits[g])
            if imp == "intronic_or_exon_noncoding":
                continue
            out.append(base_row("B_clinically_led", r, g, imp, panel, shet, v))
    return out


LOF_CSQ = ("stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant", "start_lost", "stop_lost",
           "transcript_ablation")
MID_CSQ = ("missense_variant", "inframe_insertion", "inframe_deletion", "splice_region_variant", "protein_altering_variant")


def arm_b_smallvar(rows: Iterable[Dict[str, str]], annot: Dict[Tuple[str, int, str, str], Dict[str, str]], panel: Panel,
                   shet: Dict[str, float], gm: Optional[GeneModel] = None) -> List[Dict[str, object]]:
    """SNV/indel candidates (any dnm_call) with a panel-gene consequence. `annot` maps (chrom, pos, ref, alt) to a
    record with SYMBOL, Consequence, LoF, LoF_filter, gnomAD AF and, when available, missense_tier / missense_nflag
    (see smallvar_annot.py). Candidates the annotation source does not cover are reported separately by the caller."""
    out = []
    union = panel.union
    for r in rows:
        key = (r["chrom"], _i(r.get("start")), r.get("ref", ""), r.get("alt", ""))
        a = annot.get(key)
        if not a:
            continue
        gene = a.get("SYMBOL") or a.get("gene") or ""
        if gene not in union:
            continue
        csq = a.get("Consequence") or a.get("consequence") or ""
        lof = a.get("LoF") or ""
        if lof == "HC":
            impact = "LoF_HC"
        elif any(c in csq for c in LOF_CSQ):
            impact = csq.split("&")[0] + ("(LoF_%s)" % lof if lof else "(LOFTEE_unflagged)")
        elif any(c in csq for c in MID_CSQ):
            impact = csq.split("&")[0]
            if a.get("missense_tier"):
                impact += "(%s)" % a["missense_tier"]
        else:
            continue                         # synonymous, intronic, UTR point changes: below this arm's threshold
        rr: Dict[str, object] = dict(r)
        af = a.get("gnomAD4.1_joint_AF") or a.get("gnomad_af")
        if af not in (None, "", "."):
            rr["gnomad_af"] = af
        v = Q.smallvar_quality(rr)
        extra = {"consequence": csq, "lof": lof, "lof_filter": a.get("LoF_filter", ""),
                 "missense_tier": a.get("missense_tier", ""), "missense_nflag": a.get("missense_nflag", "")}
        out.append(base_row("B_clinically_led", rr, gene, impact, panel, shet, v, extra))
    return out


def arm_b_tr(rows: Iterable[Dict[str, str]], gm: GeneModel, panel: Panel, shet: Dict[str, float]) -> List[Dict[str, object]]:
    """TR candidates at a STRchive pathogenic locus, or expanding inside an exon/UTR of a panel gene."""
    out = []
    union = panel.union
    seen = set()
    for r in rows:
        s0, e0 = _i(r.get("start")) - 1, max(_i(r.get("end")), _i(r.get("start")) + 1)
        hits = gm.hits(r["chrom"], s0, e0)
        pg = [g for g in hits if g in union and hits[g] & {"CDS", "UTR5", "UTR3"}]
        strc = str(r.get("strchive_locus") or "")
        if not pg and strc in ("", ".", "0", "None"):
            continue
        # one row per child per locus: the candidate table can carry two rows for one locus (one per outlier allele),
        # and the first run listed the same TR twice for the same child
        key = (r.get("sample_id"), r["chrom"], s0, e0)
        if key in seen:
            continue
        seen.add(key)
        v = Q.tr_quality(r)
        gene = sorted(pg, key=lambda g: (panel.rank(g), g))[0] if pg else strc.split(":")[0]
        impact = "TR_expansion_" + ("+".join(sorted(hits[gene])) if gene in hits else "STRchive_locus")
        extra = {"notes": ("STRchive:%s" % strc) if strc not in ("", ".", "0", "None") else ""}
        out.append(base_row("B_clinically_led", r, gene, impact, panel, shet, v, extra))
    return out


def write_rows(rows: List[Dict[str, object]], path: str):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=BASE_COLS, delimiter=TAB, lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in BASE_COLS})
