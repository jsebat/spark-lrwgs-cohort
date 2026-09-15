"""M2 orchestration: candidates -> six-haplotype evidence rows for one child (README §2.2 outputs).

One pass per child: open the trio's haplotagged BAMs (and TRGT spanning-read BAMs for TR, sawfish supporting reads
for SV) once, load the child's M1 label tables once, then for every CandidateRecord: ReadObs -> matrix ->
features -> transmission features -> P8 class. Output columns are stable and documented here so M3 can rely on them.
"""
from __future__ import annotations

import csv
import gzip
import json
import os
from collections import Counter
from dataclasses import asdict
from typing import Dict, Iterable, List, Optional

from ..records import CandidateRecord, read_candidates
from . import hapmatrix as H
from .readers import TrioBams, read_obs, sv_het_persistence, sv_interval_evidence

ROW_PREFIX = {("C", 1): "C1", ("C", 2): "C2", ("F", 1): "F1", ("F", 2): "F2", ("M", 1): "M1", ("M", 2): "M2"}
SV_INTERVAL_COLS = (["sv_interval_len", "sv_interval_probes"]
                    + ["%s_sv_%s_hap%s" % (r, w, h) for r in ("C", "F", "M") for w in ("in", "fl", "junc") for h in ("1", "2", "u")]
                    + ["%s_sv_%s" % (r, s) for r in ("C", "F", "M") for s in ("dp_in", "dp_fl", "ratio_hap1", "ratio_hap2", "ratio_all", "tagged_frac_in", "tagged_frac_fl", "tagged_loss", "junc_both_ends")]
                    + ["sv_depth_ratio_inside_flank", "c_sv_hap_depth_change_A", "c_sv_hap_depth_change_O",
                       "p_sv_max_hap_depth_change", "c_sv_junction_hap_concentration", "c_sv_junc_untagged_frac",
                       "c_sv_junc_both_ends", "sv_het_sites_inside", "sv_het_sites_flank", "sv_het_frac_inside",
                       "sv_het_frac_flank", "sv_het_snv_persistence", "sv_het_no_sites_inside", "sv_het_index_missing"])
ROW_FIELDS = ("dp", "alt", "ref", "amb", "mapq_mean", "mapq0_frac", "nm_alt_mean", "nm_ref_mean", "clip_alt_frac", "al_mean", "al_sd", "al_n")
CORE = ["family_id", "sample_id", "variant_id", "chrom", "start", "end", "ref", "alt", "variant_class", "caller", "caller_gt",
        "caller_gq", "caller_dp", "caller_qual", "caller_filter", "source_tier", "source_list", "mask_overlap",
        "phase_class", "rule_score", "flags", "parent_of_origin", "poo_reason", "poo_confidence",
        "child_hap1_is", "F_transmitted_hap", "M_transmitted_hap", "n_reads_used", "thresholds_version"]


def evidence_columns(k: Iterable[int]) -> List[str]:
    cols = list(CORE)
    for key in ROW_PREFIX.values():
        cols += ["%s_%s" % (key, f) for f in ROW_FIELDS]
    cols += ["C_untagged_dp", "C_untagged_alt", "F_untagged_dp", "F_untagged_alt", "M_untagged_dp", "M_untagged_alt"]
    cols += ["t_alt_reads", "t_dp", "u_alt_reads", "u_dp", "nt_parent_alt_reads", "t_hap_resolved"]
    # D-block features in a fixed order (registry names)
    cols += ["c_alt_hapA", "c_alt_hapO", "c_dp_hapA", "c_dp_hapO", "c_alt_hap_frac", "c_alt_confined", "c_alt_tagged_frac",
             "c_untagged_dp", "c_untagged_alt", "c_alt_mapq_mean", "c_alt_nm_rate", "c_ref_nm_rate", "c_alt_clip_frac",
             "c_alt_mapq0_frac", "c_alt_supp_frac", "c_alt_readlen_median", "c_alt_rq_mean",
             "c_tr_al_hapA_mean", "c_tr_al_hapA_sd", "c_tr_al_hapO_mean",
             "p_min_hap_dp", "p_max_alt_any_hap", "p_sum_alt_all_haps", "p_n_haps_with_alt", "p_max_alt_hap_frac",
             "p_untagged_dp_max", "p_untagged_alt_max", "c_amb_frac_hapA", "p_amb_frac_max"]
    for kk in k:
        cols += ["c_both_haps_obs_k%d" % kk, "p_n_haps_obs_k%d" % kk, "hap_obs_k%d" % kk]
    # SV interval evidence (P17, built 2026-09-14): per-haplotype depth inside vs flanks for all six haplotypes,
    # junction reads by haplotype, and heterozygous-SNV persistence. Empty for events without an interval.
    cols += SV_INTERVAL_COLS
    cols += ["class_payload"]
    return cols


def review_child(candidates_tsv: str, out_tsv: str, bams: TrioBams, labels: H.LabelTables, hp: H.HapParams, cp: H.ClassParams,
                 salt: str, supporting_json: Optional[str] = None, reads_jsonl_gz: Optional[str] = None,
                 max_per_class: Optional[int] = None, thresholds_version: Optional[str] = None,
                 class_filter: Optional[Iterable[str]] = None, smallvar_vcf: Optional[str] = None,
                 child_sample: Optional[str] = None) -> Dict[str, int]:
    cols = evidence_columns(hp.k)
    counts: Counter = Counter()
    per_class: Counter = Counter()
    keep = set(class_filter) if class_filter else None
    rj = gzip.open(reads_jsonl_gz, "wt") if reads_jsonl_gz else None
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for rec in read_candidates(candidates_tsv):
            if keep and rec.variant_class not in keep:
                continue
            if max_per_class is not None and per_class[rec.variant_class] >= max_per_class:
                continue
            per_class[rec.variant_class] += 1
            obs = read_obs(rec, bams, salt, supporting_json=supporting_json)
            sv_raw: Dict[str, object] = {}
            if rec.variant_class == "SV":
                # P17: alt-read confinement is a point-variant instrument. A read crossing a large deletion is split
                # between a primary and a supplementary alignment and has no heterozygous sites inside the deleted
                # segment to phase against, so its HP tag is often absent or inconsistent (measured on the cohort's
                # 35 kb MECP2 deletion: 5 junction reads per breakpoint, 3 untagged, the other 2 on OPPOSITE
                # haplotypes). Depth across the interval, per haplotype and for all six haplotypes, is the instrument
                # that works, and comparing the parents' haplotype depth is what says inherited rather than de novo.
                sv_raw = sv_interval_evidence(rec, bams, supporting_json=supporting_json)
                if sv_raw:
                    sv_raw.update(sv_het_persistence(rec, smallvar_vcf, child_sample or rec.sample_id))
            m = H.build_matrix(obs, labels, rec.chrom, rec.start, hp)
            f = H.features(m, hp)
            t = H.transmission_features(m, hp)
            c = H.classify(m, f, t, hp, cp, labels, rec.chrom, rec.start, sv=H.sv_depth_features(sv_raw) | sv_raw if sv_raw else None)
            row = {k: v for k, v in asdict(rec).items() if k in cols}
            row["class_payload"] = json.dumps(rec.class_payload, separators=(",", ":"), sort_keys=True)
            row.update(c); row.update(t); row.update(f)
            if sv_raw:
                row.update(sv_raw); row.update(H.sv_depth_features(sv_raw))
            row.update(child_hap1_is=m.child_hap1_is or ".", F_transmitted_hap=m.transmitted.get("F") or ".",
                       M_transmitted_hap=m.transmitted.get("M") or ".", n_reads_used=m.n_reads, thresholds_version=thresholds_version)
            for key, pre in ROW_PREFIX.items():
                s = m.row(*key).summary()
                for fld in ROW_FIELDS:
                    row["%s_%s" % (pre, fld)] = s.get(fld)
            for role in ("C", "F", "M"):
                row["%s_untagged_dp" % role] = m.untagged[role].dp
                row["%s_untagged_alt" % role] = m.untagged[role].alt
            w.writerow({k: ("" if row.get(k) is None else row.get(k)) for k in cols})
            counts[c["phase_class"]] += 1
            counts["rows"] += 1
            if rj is not None:
                rj.write(json.dumps({"variant_id": rec.variant_id, "class": rec.variant_class,
                                     "reads": [dict(rid=o.rid, role=o.role, hp=o.hp, ps=o.ps, support=o.support, mapq=o.mapq,
                                                    nm=o.nm_rate, clip=o.clipped, al=o.al) for o in obs]}, separators=(",", ":")) + "\n")
    if rj is not None:
        rj.close()
    return dict(counts)


def reclassify_table(evidence_in: str, out_tsv: str, hp: H.HapParams, cp: H.ClassParams,
                     thresholds_version: Optional[str] = None) -> Dict[str, int]:
    """Re-run features/transmission/classify from the count columns of an existing evidence table (no BAMs).
    Columns added since the table was written (registry D block) are appended; everything else keeps its place."""
    counts: Counter = Counter()
    with open(evidence_in, newline="") as fh, open(out_tsv, "w", newline="") as out:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = list(rd.fieldnames or [])
        for c in evidence_columns(hp.k):
            if c not in cols:
                cols.insert(cols.index("class_payload") if "class_payload" in cols else len(cols), c)
        w = csv.DictWriter(out, fieldnames=cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rd:
            row = H.reclassify_row(r, hp, cp, thresholds_version)
            w.writerow({k: ("" if row.get(k) is None else row.get(k)) for k in cols})
            counts[row["phase_class"]] += 1
            counts["rows"] += 1
    return dict(counts)
