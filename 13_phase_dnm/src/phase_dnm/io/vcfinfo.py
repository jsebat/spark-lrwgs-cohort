"""VCF INFO equivalents of the final table (README §2.3): one header block, one INFO string per row, and a
sites-only VCF per class so the calls can travel with any VCF-based downstream step. Values are copied from the
final table verbatim; nothing is recomputed here."""
from __future__ import annotations

import csv
from typing import Dict, Iterable, List

INFO_FIELDS = [
    ("PDNM_PROB", "1", "Float", "Module 4 classifier probability (rf_prob); missing until a model is frozen"),
    ("PDNM_CALL", "1", "String", "Final de novo call after the P15 phase layer: YES or NO"),
    ("PDNM_MODE", "1", "String", "Decision mode: rf+phase or phase_only (provisional, no classifier)"),
    ("PDNM_WHY", "1", "String", "Decision reason: RF, RESCUED, BELOW_TAU, RF_UNSUPPORTED:<class>, PHASE_ONLY, DEMOTED:<class>, MOSAIC:<class>, NOT_PHASED_GERMLINE, HAP_UNOBSERVED, LOW_POSTERIOR"),
    ("PDNM_POO", "1", "String", "Parent of origin of the alt-carrying child haplotype: paternal, maternal, undetermined"),
    ("PDNM_POOR", "1", "String", "Parent-of-origin reason code"),
    ("PDNM_POOC", "1", "Float", "Parent-of-origin confidence (fraction of tagged alt reads on the origin haplotype)"),
    ("PDNM_CLASS", "1", "String", "Phase class (P8 rule layer)"),
    ("PDNM_RULE", "1", "Integer", "Rule score 0-6 (germline criteria met)"),
    ("PDNM_HAPOBS", "1", "Integer", "Six-haplotype observability at k=5 (0-6, readable reads)"),
    ("PDNM_CHF", "1", "Float", "Child alt fraction on the alt-carrying haplotype"),
    ("PDNM_CAO", "1", "Integer", "Child alt reads on the other haplotype"),
    ("PDNM_TALT", "1", "Integer", "Alt reads on the origin parent's TRANSMITTED haplotype"),
    ("PDNM_UALT", "1", "Integer", "Alt reads on the origin parent's UNTRANSMITTED haplotype"),
    ("PDNM_TDP", "1", "Integer", "Depth on the origin parent's transmitted haplotype"),
    ("PDNM_SCORE", "1", "Float", "P9 posterior of the germline hypothesis (phase_score)"),
    ("PDNM_ALT", "1", "String", "P9 best alternative hypothesis"),
    ("PDNM_LR", "1", "Float", "log10 likelihood ratio germline vs best alternative"),
    ("PDNM_MOSAIC", "0", "Flag", "Row is a mosaic class (never YES; reported separately)"),
    ("PDNM_FLAGS", ".", "String", "Review flags"),
    ("PDNM_TIER", "1", "String", "Source tier of the candidate: HIGH, LOW, UNFILTERED (P5)"),
    ("PDNM_MASK", "1", "Integer", "Overlaps the pipeline's region mask (flag, never a filter)"),
]
COLUMN_OF = {"PDNM_PROB": "rf_prob", "PDNM_CALL": "dnm_call", "PDNM_MODE": "call_mode", "PDNM_WHY": "decision_reason",
             "PDNM_POO": "parent_of_origin", "PDNM_POOR": "poo_reason", "PDNM_POOC": "poo_confidence", "PDNM_CLASS": "phase_class",
             "PDNM_RULE": "rule_score", "PDNM_HAPOBS": "hap_obs_k5", "PDNM_CHF": "child_alt_hap_frac", "PDNM_CAO": "child_alt_other_hap",
             "PDNM_TALT": "transmitted_parent_alt_reads", "PDNM_UALT": "untransmitted_parent_alt_reads", "PDNM_TDP": "transmitted_parent_dp",
             "PDNM_SCORE": "phase_score", "PDNM_ALT": "lik_best_alternative", "PDNM_LR": "lik_log10lr_germline",
             "PDNM_FLAGS": "flags", "PDNM_TIER": "source_tier", "PDNM_MASK": "mask_overlap"}


def header_lines(source: str = "phase_dnm") -> List[str]:
    out = ["##fileformat=VCFv4.3", "##source=%s" % source]
    for k, n, t, d in INFO_FIELDS:
        out.append('##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">' % (k, n, t, d))
    return out


def _clean(v: object) -> str:
    s = str(v).replace(";", "|").replace(" ", "_").replace("=", ":").replace(",", "|")
    return s


def info_string(row: Dict[str, object]) -> str:
    parts = []
    for k, n, t, d in INFO_FIELDS:
        if k == "PDNM_MOSAIC":
            if str(row.get("mosaic_flag", "0")) in ("1", "True", "true"):
                parts.append("PDNM_MOSAIC")
            continue
        v = row.get(COLUMN_OF[k])
        if v in (None, "", "."):
            continue
        parts.append("%s=%s" % (k, _clean(v)))
    return ";".join(parts) if parts else "."


def write_sites_vcf(rows: Iterable[Dict[str, object]], path: str, contigs: Iterable[str] = ()) -> int:
    """Sites-only VCF (no samples): CHROM POS ID REF ALT QUAL FILTER INFO. Symbolic ALT for SV/TR rows."""
    n = 0
    with open(path, "w") as fh:
        for line in header_lines():
            fh.write(line + "\n")
        for c in contigs:
            fh.write("##contig=<ID=%s>\n" % c)
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for r in rows:
            vclass = r.get("variant_class")
            ref, alt = str(r.get("ref") or "N"), str(r.get("alt") or ".")
            if vclass == "SV":
                alt = "<%s>" % (r.get("svtype") or "SV") if not alt.startswith("<") and len(alt) > 50 else alt
                ref = ref[:1] or "N"
            elif vclass == "TR":
                ref, alt = "N", "<TR>" if not alt.startswith("<") else alt
            qual = r.get("caller_qual")
            fh.write("\t".join([str(r.get("chrom")), str(r.get("start")), str(r.get("variant_id") or "."), ref, alt,
                                "." if qual in (None, "") else str(qual), ".", info_string(r)]) + "\n")
            n += 1
    return n


def rows_from_final(path: str) -> List[Dict[str, object]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))
