"""Feature extraction: evidence row + candidate record -> one registered feature vector (DESIGN P11, features.yaml).

Blocks and where each comes from:
  A caller   the candidate record (GT/GQ/DP/AD/PL of child and parents; SynthDNM's universal 21 derived exactly as
             synthdnm/scripts/preprocess_features.py does: AR = ref/(alt+1), min/max over parents, PL0/1/2 for the
             allele, AB, indel_flag, haploid_flag from sex + PAR), class fields from the payload.
  B context  mask overlap from a BED (a FLAG, P5), sequence context from the reference FASTA when given
             (homopolymer run, 21-bp entropy, CpG); cohort/population columns are produced by a separate
             annotation step and left empty here (reported in the coverage summary).
  C reads    child read-quality summaries already in the evidence row.
  D phase    the six-haplotype block already in the evidence row; rf_safe by construction (hapmatrix.features).
  transmission (rf_safe False): copied through so M2/M3 tables are complete; never enters an RF matrix
             (registry.assert_rf_safe).

Two outputs per child and class: `<child>.<class>.features.tsv` (every registered feature that applies, id columns,
labels) and, on request, the classifier matrix restricted to rf_safe columns (`.rf.tsv`).
"""
from __future__ import annotations

import bisect
import csv
import gzip
import json
import math
import os
from collections import Counter
from typing import Dict, Iterable, List, Optional, Tuple

from ..io.vcf import normalise_sex
from ..records import CandidateRecord, read_candidates
from .registry import Registry

# GRCh38 pseudoautosomal regions, for haploid_flag. Same values as 12_x_inactivation/03_trio_phase_x.py
# and as ../../01_phasing/src/trio_phase/phasing/orient.py, which is where this used to be imported
# from; it is restated here rather than imported so that this module has no import dependency on the
# phasing module, whose tables it consumes by path.
PAR_GRCH38 = ((10001, 2781479), (155701383, 156030895))

ID_COLUMNS = ["family_id", "sample_id", "variant_id", "chrom", "start", "end", "ref", "alt", "variant_class", "caller",
              "source_tier", "source_list", "phase_class", "rule_score", "flags", "parent_of_origin", "poo_reason", "poo_confidence"]


# ----------------------------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------------------------
def _ints(s: Optional[str]) -> Optional[List[int]]:
    if s in (None, "", "."):
        return None
    out = []
    for p in str(s).split(","):
        if p in (".", ""):
            return None
        try:
            out.append(int(float(p)))
        except ValueError:
            return None
    return out


def _f(x) -> Optional[float]:
    if x in (None, "", "."):
        return None
    try:
        return float(x)
    except ValueError:
        return None


def allele_ratio(ad: Optional[List[int]], a: int) -> Optional[float]:
    """SynthDNM AR = ref / (alt + 1) for allele index a."""
    if not ad or a >= len(ad):
        return None
    return ad[0] / (ad[a] + 1)


def allele_balance(ad: Optional[List[int]], a: int) -> Optional[float]:
    if not ad or a >= len(ad) or (ad[0] + ad[a]) == 0:
        return None
    return ad[a] / (ad[0] + ad[a])


def pl_triplet(pl: Optional[List[int]], a: int) -> Tuple[Optional[int], Optional[int], Optional[int]]:
    """PL for genotypes 0/0, 0/a, a/a (VCF ordering: index of j/k = k(k+1)/2 + j)."""
    if not pl:
        return None, None, None
    i0, i1, i2 = 0, a * (a + 1) // 2, a * (a + 1) // 2 + a
    try:
        return pl[i0], pl[i1], pl[i2]
    except IndexError:
        return None, None, None


def _mm(vals):
    v = [x for x in vals if x is not None]
    return (min(v), max(v)) if v else (None, None)


class Geometry:
    """Child phase-segment geometry (orientation.tsv: chrom start end per segment; segments of one child do not overlap)
    and resolved crossover positions (changepoints.resolved.tsv, class CROSSOVER). Child-only quantities: the segment
    length and the distance to its edge are rf_safe; the crossover distances are transmission-derived and rf_safe: false
    (registry) - computed here for the rule/likelihood layer and the final table, never for the classifier matrix."""

    def __init__(self, orientation_tsv: Optional[str] = None, changepoints_tsv: Optional[str] = None):
        import bisect as _b
        self._b = _b
        self.seg: Dict[str, List[Tuple[int, int]]] = {}
        self.xo: Dict[str, List[int]] = {}
        if orientation_tsv and os.path.exists(orientation_tsv):
            with open(orientation_tsv, newline="") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    try:
                        self.seg.setdefault(r["chrom"], []).append((int(r["start"]), int(r["end"])))
                    except (KeyError, ValueError):
                        continue
            for c in self.seg:
                self.seg[c].sort()
        if changepoints_tsv and os.path.exists(changepoints_tsv):
            with open(changepoints_tsv, newline="") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    if str(r.get("status") or r.get("read_class") or r.get("class") or "").upper().startswith("CROSSOVER"):
                        try:
                            self.xo.setdefault(r["chrom"], []).append((int(r["left_pos"]) + int(r["right_pos"])) // 2)
                        except (KeyError, ValueError):
                            continue
            for c in self.xo:
                self.xo[c].sort()

    def at(self, chrom: str, pos1: int) -> Dict[str, object]:
        out: Dict[str, object] = {}
        segs = self.seg.get(chrom)
        if segs is not None:
            i = self._b.bisect_right([s for s, _ in segs], pos1) - 1
            if 0 <= i < len(segs) and segs[i][0] <= pos1 <= segs[i][1]:
                s, e = segs[i]
                out["c_block_len_log10"] = round(math.log10(max(1, e - s + 1)), 3)
                out["c_dist_block_edge_log10"] = round(math.log10(max(1, min(pos1 - s, e - pos1) + 1)), 3)
            else:
                out["c_block_len_log10"] = None
                out["c_dist_block_edge_log10"] = -1.0
        xs = self.xo.get(chrom)
        if xs:
            i = self._b.bisect_left(xs, pos1)
            d = min(abs(xs[j] - pos1) for j in (i - 1, i) if 0 <= j < len(xs))
            out["dist_crossover_log10"] = round(math.log10(d + 1), 3)
            out["near_crossover_flag"] = int(d <= 50000)
        return out


class BedMask:
    """Interval lookup over one or more BED files (0-based half-open) -> overlap flag / fraction."""

    def __init__(self, paths: Iterable[str]):
        self.iv: Dict[str, List[Tuple[int, int]]] = {}
        for p in paths:
            op = gzip.open if str(p).endswith(".gz") else open
            with op(p, "rt") as fh:
                for ln in fh:
                    if not ln.strip() or ln.startswith(("#", "track", "browser")):
                        continue
                    f = ln.split("\t")
                    self.iv.setdefault(f[0], []).append((int(f[1]), int(f[2])))
        self.starts: Dict[str, List[int]] = {}
        for c, lst in self.iv.items():
            lst.sort()
            merged: List[Tuple[int, int]] = []
            for s, e in lst:
                if merged and s <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], e))
                else:
                    merged.append((s, e))
            self.iv[c] = merged
            self.starts[c] = [s for s, _ in merged]

    def overlap_bp(self, chrom: str, start1: int, end1: int) -> int:
        lst = self.iv.get(chrom)
        if not lst:
            return 0
        s0, e0 = start1 - 1, end1
        i = bisect.bisect_right(self.starts[chrom], s0) - 1
        tot = 0
        for j in range(max(0, i), len(lst)):
            s, e = lst[j]
            if s >= e0:
                break
            tot += max(0, min(e, e0) - max(s, s0))
        return tot

    def flag(self, chrom: str, start1: int, end1: int, frac: float = 0.5) -> int:
        L = max(1, end1 - start1 + 1)
        ov = self.overlap_bp(chrom, start1, end1)
        return int(ov > 0) if L <= 50 else int(ov / L >= frac)


class SeqContext:
    """Homopolymer run, 21-bp Shannon entropy and CpG flag from a reference FASTA (pysam, lazy)."""

    def __init__(self, fasta: str):
        import pysam
        self.fa = pysam.FastaFile(fasta)

    def at(self, chrom: str, pos1: int, ref: str, alt: str) -> Dict[str, object]:
        seq = self.fa.fetch(chrom, max(0, pos1 - 11), pos1 + 10).upper()
        centre = min(10, pos1 - 1)
        base = seq[centre] if centre < len(seq) else "N"
        run = 1
        i = centre - 1
        while i >= 0 and seq[i] == base:
            run += 1; i -= 1
        i = centre + 1
        while i < len(seq) and seq[i] == base:
            run += 1; i += 1
        cnt = Counter(seq)
        n = len(seq) or 1
        ent = -sum((c / n) * math.log2(c / n) for c in cnt.values() if c)
        cpg = 0
        if len(ref) == 1 and len(alt) == 1:
            nxt = seq[centre + 1] if centre + 1 < len(seq) else "N"
            prv = seq[centre - 1] if centre >= 1 else "N"
            cpg = int((base == "C" and nxt == "G") or (base == "G" and prv == "C"))
        return {"homopolymer_run": run, "seq_entropy_21bp": round(ent, 3), "cpg_context": cpg}


# ----------------------------------------------------------------------------------------------
# extraction
# ----------------------------------------------------------------------------------------------
def caller_features(rec: CandidateRecord, sex: str, ev: Optional[Dict[str, object]] = None) -> Dict[str, object]:
    pl = rec.class_payload
    a = int(pl.get("allele_index", 1)) if rec.variant_class in ("SNV", "INDEL") else 1
    c_ad, f_ad, m_ad = _ints(rec.child_ad), _ints(rec.father_ad), _ints(rec.mother_ad)
    c_pl, f_pl, m_pl = _ints(rec.child_pl), _ints(rec.father_pl), _ints(rec.mother_pl)
    f: Dict[str, object] = {}
    f["child_AR"] = allele_ratio(c_ad, a); f["min_AR"], f["max_AR"] = _mm([allele_ratio(f_ad, a), allele_ratio(m_ad, a)])
    f["child_GQ"] = rec.caller_gq; f["min_GQ"], f["max_GQ"] = _mm([rec.father_gq, rec.mother_gq])
    f["child_DP"] = rec.caller_dp; f["min_DP"], f["max_DP"] = _mm([rec.father_dp, rec.mother_dp])
    for i, (cv, fv, mv) in enumerate(zip(pl_triplet(c_pl, a), pl_triplet(f_pl, a), pl_triplet(m_pl, a))):
        f["child_PL%d" % i] = cv; f["min_PL%d" % i], f["max_PL%d" % i] = _mm([fv, mv])
    f["child_AB"] = allele_balance(c_ad, a)
    f["indel_flag"] = int(rec.variant_class == "INDEL")
    sx = normalise_sex(sex)
    nonpar_x = rec.chrom in ("chrX", "X") and not any(lo <= rec.start <= hi for lo, hi in PAR_GRCH38)
    f["haploid_flag"] = int(sx == "M" and (nonpar_x or rec.chrom in ("chrY", "Y")))
    # QD, computed here rather than taken from the caller (P24 correction 2026-09-14, JS).
    # Raw site QUAL is unusable as a feature: in a joint callset it rises with the number of carriers, AND our real
    # candidates come from the per-FAMILY joint VCF while the synthetic ones come from the 105-sample cohort BCF, so the
    # same variant scores differently on the two sides of the label. GATK's QD (QUAL normalised by the depth of the
    # informative samples) is the frequency-independent form, but QD is GATK-only and DeepVariant/GLnexus, sawfish and
    # TRGT do not emit it. So it is derived from PER-SAMPLE evidence, which is on the same scale in both sources:
    #   numerator   PL[0] = the child's Phred evidence against hom-ref (the per-sample analogue of QUAL); GQ when the
    #               caller emits no PL (sawfish, TRGT)
    #   denominator the child's READABLE depth from our own six-haplotype matrix, falling back to the caller's DP
    f["site_qual"] = rec.caller_qual        # kept in the evidence tables for provenance; NOT a registry feature
    _num = _f(f.get("child_PL0"))
    if _num is None:
        _num = _f(f.get("child_GQ"))
    # the readable six-haplotype depth lives in the EVIDENCE row, not in this dict: reading it from `f` (as this did
    # until 2026-09-15) always found nothing and silently fell through to the caller's DP, which is not what P24 says.
    _ev = ev or {}
    _den = (_f(_ev.get("c_dp_hapA")) or 0.0) + (_f(_ev.get("c_dp_hapO")) or 0.0)
    if not _den:
        _den = _f(f.get("child_DP")) or 0.0
    f["qd_child"] = round(_num / _den, 4) if (_num is not None and _den) else None
    f["site_filter_fail"] = int(rec.caller_filter not in (".", "PASS", ""))
    f["multiallelic"] = int(int(pl.get("n_alts", 1) or 1) > 1)
    f["glnexus_rnc_child"] = pl.get("rnc")
    if rec.variant_class == "INDEL":
        f["indel_len"] = abs(len(rec.alt) - len(rec.ref))
    if rec.variant_class == "SV":
        f["svtype"] = pl.get("svtype")
        f["svlen_log10"] = round(math.log10(abs(int(pl["svlen"]))), 3) if pl.get("svlen") not in (None, 0) else None
        f["bp_homology_len"] = pl.get("homlen")
        f["child_sv_support"] = c_ad[1] if c_ad and len(c_ad) > 1 else None
        f["max_parent_sv_support"] = _mm([f_ad[1] if f_ad and len(f_ad) > 1 else None, m_ad[1] if m_ad and len(m_ad) > 1 else None])[1]
    if rec.variant_class == "TR":
        c_al = pl.get("child_AL") or []; idx = int(pl.get("outlier_allele_idx", 0))
        f["child_AL_expanded"] = c_al[idx] if idx < len(c_al) else None
        f["child_AL_other"] = min((v for i, v in enumerate(c_al) if i != idx), default=None)
        c_sd = _ints(pl.get("child_SD")); f["child_SD_expanded"] = c_sd[idx] if c_sd and idx < len(c_sd) else None
        f["min_parent_SD"] = _mm([min(_ints(pl.get("father_SD")) or [None]), min(_ints(pl.get("mother_SD")) or [None])])[0]
        f["delta_motif_units_nearest_parent"] = pl.get("delta_units")
        f["motif_len"] = pl.get("motif_unit_bp")
        f["n_motifs"] = len((pl.get("motifs") or "").split(",")) if pl.get("motifs") else None
        f["locus_len_ref"] = rec.end - rec.start + 1
        allr = (pl.get("child_ALLR") or "").split(",")
        try:
            lo, hi = allr[idx].split("-"); f["child_ALLR_width"] = int(hi) - int(lo)
        except (IndexError, ValueError):
            f["child_ALLR_width"] = None
    return f


def extract_child(evidence_tsv: str, candidates_tsv: str, registry: Registry, sex: str, out_tsv: str,
                  rf_out_tsv: Optional[str] = None, mask: Optional[BedMask] = None, seqctx: Optional[SeqContext] = None,
                  annot: Optional[Dict[str, dict]] = None, geom: Optional["Geometry"] = None) -> Dict[str, object]:
    cands: Dict[str, CandidateRecord] = {}
    for r in read_candidates(candidates_tsv):
        cands[r.variant_id] = r
    with open(evidence_tsv, newline="") as fh:
        ev_rows = list(csv.DictReader(fh, delimiter="\t"))
    if not ev_rows:
        return {"rows": 0}
    vclass_set = {r["variant_class"] for r in ev_rows}
    feats = sorted({f.name for vc in vclass_set for f in registry.for_class(vc)})
    cols = ID_COLUMNS + feats
    produced: Counter = Counter()
    # A column can be "produced" and still carry nothing: annotate.write_annot writes the sentinel "." for
    # "checked, no hit" (strchive_locus) and callers write "." for fields they do not emit, and a column that is
    # one value in every row is exactly what the presence-leak guard drops at training time. Counting distinct
    # values (capped at 2 - we only need "more than one") makes that visible in the summary instead of leaving it
    # to be found by hand, as glnexus_rnc_child and parent_rnc_any were.
    distinct: Dict[str, set] = {}
    n = 0
    rf_cols = None
    rf_fh = None
    if rf_out_tsv:
        seen = set(); rf_cols = []
        for vc in sorted(vclass_set):
            for c in registry.rf_matrix_columns(vc):
                if c not in seen:
                    seen.add(c); rf_cols.append(c)
        registry.assert_rf_safe(rf_cols)
        rf_fh = open(rf_out_tsv, "w", newline="")
        rf_fh.write("\t".join(["family_id", "sample_id", "variant_id", "variant_class"] + rf_cols) + "\n")
    with open(out_tsv, "w", newline="") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in ev_rows:
            rec = cands.get(r["variant_id"])
            row: Dict[str, object] = {k: r.get(k) for k in ID_COLUMNS}
            # D + C + transmission blocks: already in the evidence row under registry names
            for k in feats:
                if k in r and r[k] not in (None, ""):
                    row[k] = r[k]
            if geom is not None:
                for k, v in geom.at(r["chrom"], int(float(r["start"]))).items():
                    if k in registry.features and v is not None:
                        row[k] = v
            if annot is not None and r["variant_id"] in annot:
                for k, v in annot[r["variant_id"]].items():
                    if k in registry.features and v not in (None, ""):
                        row[k] = v
            if rec is not None:
                for k, v in caller_features(rec, sex, ev=r).items():
                    if k in registry.features:
                        row[k] = v
                if mask is not None:
                    row["segdup_overlap"] = mask.flag(rec.chrom, rec.start, rec.end)
                if seqctx is not None and rec.variant_class in ("SNV", "INDEL"):
                    for k, v in seqctx.at(rec.chrom, rec.start, rec.ref, rec.alt).items():
                        if k in registry.features:
                            row[k] = v
            for k in feats:
                v = row.get(k)
                if v not in (None, ""):
                    produced[k] += 1
                    d = distinct.setdefault(k, set())
                    if len(d) < 2:
                        d.add(str(v))
            fh.write("\t".join("" if row.get(c) is None else str(row.get(c)) for c in cols) + "\n")
            if rf_fh is not None:
                rf_fh.write("\t".join([str(row.get(k, "")) for k in ("family_id", "sample_id", "variant_id", "variant_class")]
                                      + ["" if row.get(k) is None else str(row.get(k)) for k in rf_cols]) + "\n")
            n += 1
    if rf_fh is not None:
        rf_fh.close()
    never = [k for k in feats if produced[k] == 0]
    constant = [k for k in feats if produced[k] and len(distinct.get(k, ())) < 2]
    return {"rows": n, "features_applicable": len(feats), "features_produced": len(feats) - len(never),
            "never_produced": never, "constant_columns": constant,
            "features_informative": len(feats) - len(never) - len(constant),
            "rf_columns": len(rf_cols) if rf_cols else None, "registry_sha256": registry.sha256}
