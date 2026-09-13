"""M2 candidate generation (DESIGN P5): every putative de novo event of every class, UNFILTERED, from the family
joint callsets. Three generators share one contract - yield CandidateRecords for one child - and differ only in
the class-specific rule for "the child carries something neither parent's genotype has":

  SNV/INDEL  family joint small-variant VCF (DeepVariant + GLnexus): the child carries an alt allele index that is
             in neither parent's called genotype; any GQ / DP / AB; parents must be CALLED (a missing parent is not
             evidence of absence). One row per such alt allele.
  SV         family joint sawfish VCF: the child's GT carries an alt allele and both parents are called 0/0.
             Payload carries SVTYPE, SVLEN, END, HOMLEN, IMPRECISE, MATEID, the caller id, and per-sample AD/CN.
  TR         family joint TRGT VCF: a child allele length (AL) that exceeds the LONGEST parental allele, or
             undercuts the SHORTEST, by >= `min_units` motif units (unit = shortest motif). Payload carries TRID,
             MOTIFS, STRUC, per-sample AL/SD/ALLR/MC and which child allele is the outlier. Both directions are kept
             (the existing pipeline looks at expansions only; contraction candidates are flagged `direction`).

Existing per-family candidate lists (slivar hiconf, denovo_sv.tsv, denovo_tr.tsv) are merged in afterwards by
`merge_lists` so `source_tier` / `source_list` record where each row was already seen (P18 needs this).
"""
from __future__ import annotations

import csv
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

from .io.vcf import Record, VcfReader, iter_records, _parse_gt
from .records import CandidateRecord


def _ival(d: Dict[str, str], k: str) -> Optional[int]:
    v = d.get(k)
    if v in (None, ".", ""):
        return None
    try:
        return int(v)
    except ValueError:
        try:
            return int(float(v))
        except ValueError:
            return None


def _called(alleles) -> bool:
    return alleles is not None and len(alleles) > 0


# ----------------------------------------------------------------------------------------------
# SNV / INDEL
# ----------------------------------------------------------------------------------------------
def snv_indel_candidates(vcf_path: str, family: str, child: str, father: str, mother: str,
                         caller: str = "deepvariant_glnexus") -> Iterator[CandidateRecord]:
    rd = VcfReader(vcf_path)
    for s in (child, father, mother):
        if s not in rd.samples:
            raise KeyError("sample not in %s" % vcf_path)
    for rec in iter_records(rd, (child, father, mother)):
        c_al, _ = rec.gt(child)
        f_al, _ = rec.gt(father)
        m_al, _ = rec.gt(mother)
        if not (_called(c_al) and _called(f_al) and _called(m_al)):
            continue
        parental = set(f_al) | set(m_al)
        novel = sorted(a for a in set(c_al) if a > 0 and a not in parental)
        if not novel:
            continue
        cs, fs, ms = rec.samples[child], rec.samples[father], rec.samples[mother]
        for a in novel:
            alt = rec.alts[a - 1]
            vclass = "SNV" if len(rec.ref) == 1 and len(alt) == 1 else "INDEL"
            end = rec.pos + max(len(rec.ref), 1) - 1
            yield CandidateRecord(
                family_id=family, sample_id=child, variant_id="%s:%d:%s:%s" % (rec.chrom, rec.pos, rec.ref, alt),
                chrom=rec.chrom, start=rec.pos, end=end, ref=rec.ref, alt=alt, variant_class=vclass, caller=caller,
                caller_gt=cs.get("GT", "."), caller_gq=_ival(cs, "GQ"), caller_dp=_ival(cs, "DP"), caller_qual=rec.qual,
                caller_filter=rec.filter, child_ad=cs.get("AD", "."), father_gt=fs.get("GT", "."), mother_gt=ms.get("GT", "."),
                father_gq=_ival(fs, "GQ"), mother_gq=_ival(ms, "GQ"), father_dp=_ival(fs, "DP"), mother_dp=_ival(ms, "DP"),
                father_ad=fs.get("AD", "."), mother_ad=ms.get("AD", "."),
                child_pl=cs.get("PL", "."), father_pl=fs.get("PL", "."), mother_pl=ms.get("PL", "."),
                class_payload={"allele_index": a, "n_alts": len(rec.alts), "rnc": cs.get("RNC", ".")})


# ----------------------------------------------------------------------------------------------
# SV (sawfish)
# ----------------------------------------------------------------------------------------------
def sv_candidates(vcf_path: str, family: str, child: str, father: str, mother: str, caller: str = "sawfish"
                  ) -> Iterator[CandidateRecord]:
    rd = VcfReader(vcf_path)
    for s in (child, father, mother):
        if s not in rd.samples:
            raise KeyError("sample not in %s" % vcf_path)
    for rec in iter_records(rd, (child, father, mother)):
        c_al, _ = rec.gt(child)
        f_al, _ = rec.gt(father)
        m_al, _ = rec.gt(mother)
        if not (_called(c_al) and _called(f_al) and _called(m_al)):
            continue
        if any(a > 0 for a in f_al) or any(a > 0 for a in m_al) or not any(a > 0 for a in c_al):
            continue
        cs, fs, ms = rec.samples[child], rec.samples[father], rec.samples[mother]
        svtype = rec.info.get("SVTYPE", "NA")
        svlen = _ival(rec.info, "SVLEN")
        end = _ival(rec.info, "END") or (rec.pos + abs(svlen) if svlen and svtype == "DEL" else rec.pos)
        alt = rec.alts[0]
        yield CandidateRecord(
            family_id=family, sample_id=child, variant_id=rec.id if rec.id not in (".", "") else "%s:%d:%s" % (rec.chrom, rec.pos, svtype),
            chrom=rec.chrom, start=rec.pos, end=end, ref=rec.ref[:50], alt=alt if len(alt) <= 50 else "<%s:%dbp>" % (svtype, len(alt)),
            variant_class="SV", caller=caller, caller_gt=cs.get("GT", "."), caller_gq=_ival(cs, "GQ"), caller_dp=None,
            caller_qual=rec.qual, caller_filter=rec.filter, child_ad=cs.get("AD", "."), father_gt=fs.get("GT", "."),
            mother_gt=ms.get("GT", "."), father_gq=_ival(fs, "GQ"), mother_gq=_ival(ms, "GQ"), father_ad=fs.get("AD", "."),
            mother_ad=ms.get("AD", "."), child_pl=cs.get("PL", "."), father_pl=fs.get("PL", "."), mother_pl=ms.get("PL", "."),
            class_payload={"svtype": svtype, "svlen": svlen, "end": end, "homlen": _ival(rec.info, "HOMLEN"),
                           "imprecise": rec.info.get("IMPRECISE") == "1", "mateid": rec.info.get("MATEID"),
                           "svclaim": rec.info.get("SVCLAIM"), "caller_id": rec.id, "child_cn": cs.get("CN"),
                           "father_cn": fs.get("CN"), "mother_cn": ms.get("CN"), "insseq_len": len(rec.info.get("INSSEQ", "")) or None})


# ----------------------------------------------------------------------------------------------
# TR (TRGT)
# ----------------------------------------------------------------------------------------------
def _al(d: Dict[str, str]) -> Optional[List[int]]:
    v = d.get("AL")
    if v in (None, ".", ""):
        return None
    out = []
    for p in v.split(","):
        if p in (".", ""):
            return None
        out.append(int(float(p)))
    return out


def tr_candidates(vcf_path: str, family: str, child: str, father: str, mother: str, min_units: int = 1,
                  min_bp: int = 1, caller: str = "trgt") -> Iterator[CandidateRecord]:
    rd = VcfReader(vcf_path)
    for s in (child, father, mother):
        if s not in rd.samples:
            raise KeyError("sample not in %s" % vcf_path)
    for rec in iter_records(rd, (child, father, mother)):
        cs, fs, ms = rec.samples[child], rec.samples[father], rec.samples[mother]
        c, f, m = _al(cs), _al(fs), _al(ms)
        if not (c and f and m):
            continue                                            # uncalled in anyone: not evaluable (as the existing rule)
        motifs = rec.info.get("MOTIFS", "")
        unit = max(1, min((len(x) for x in motifs.split(",") if x), default=1))
        margin = max(min_bp, min_units * unit)
        pmax, pmin = max(f + m), min(f + m)
        for idx, al in enumerate(c):
            if al >= pmax + margin:
                direction, delta = "expansion", al - pmax
            elif al <= pmin - margin:
                direction, delta = "contraction", pmin - al
            else:
                continue
            end = _ival(rec.info, "END") or rec.pos
            yield CandidateRecord(
                family_id=family, sample_id=child, variant_id=rec.info.get("TRID", "%s:%d" % (rec.chrom, rec.pos)),
                chrom=rec.chrom, start=rec.pos, end=end, ref="<TR>", alt="<AL=%d>" % al, variant_class="TR", caller=caller,
                caller_gt=cs.get("GT", "."), caller_gq=None, caller_dp=None, caller_qual=rec.qual, caller_filter=rec.filter,
                child_ad=cs.get("SD", "."), father_gt=fs.get("GT", "."), mother_gt=ms.get("GT", "."), father_ad=fs.get("SD", "."),
                mother_ad=ms.get("SD", "."),
                class_payload={"trid": rec.info.get("TRID"), "motifs": motifs, "struc": rec.info.get("STRUC"), "motif_unit_bp": unit,
                               "child_AL": c, "father_AL": f, "mother_AL": m, "child_SD": cs.get("SD"), "father_SD": fs.get("SD"),
                               "mother_SD": ms.get("SD"), "child_ALLR": cs.get("ALLR"), "child_MC": cs.get("MC"),
                               "outlier_allele_idx": idx, "direction": direction, "delta_bp": delta,
                               "delta_units": round(delta / unit, 2)})


# ----------------------------------------------------------------------------------------------
# merge with existing per-family lists (P18): mark rows already present in a list
# ----------------------------------------------------------------------------------------------
def merge_lists(records: Iterable[CandidateRecord], lists: Sequence[Tuple[str, str, str]]) -> List[CandidateRecord]:
    """lists: (list_name, tier HIGH|LOW, tsv path with chrom,pos,ref,alt columns for small variants or chrom,pos,svtype
    for SV or trid for TR). A record found in a list gets that tier and the list name appended to source_list."""
    keys: Dict[str, set] = {}
    tiers: Dict[str, str] = {}
    for name, tier, path in lists:
        tiers[name] = tier
        ks = set()
        with open(path, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if "trid" in r:
                    ks.add(("TR", r["trid"]))
                elif "svtype" in r:
                    ks.add(("SV", "%s:%s:%s" % (r["chrom"], r["pos"], r["svtype"])))
                else:
                    ks.add(("SMALL", "%s:%s:%s:%s" % (r["chrom"], r["pos"], r["ref"], r["alt"])))
        keys[name] = ks
    out = []
    for rec in records:
        k = ("TR", rec.variant_id) if rec.variant_class == "TR" else \
            ("SV", "%s:%d:%s" % (rec.chrom, rec.start, rec.class_payload.get("svtype"))) if rec.variant_class == "SV" else \
            ("SMALL", rec.variant_id)
        hits = [name for name, ks in keys.items() if k in ks]
        if hits:
            rec.source_list = ";".join([rec.source_list] + hits)
            rec.source_tier = "HIGH" if any(tiers[h] == "HIGH" for h in hits) else "LOW"
        out.append(rec)
    return out
