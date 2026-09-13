"""Spike-in harness (DESIGN R3, §5; README `sim/spike.py`): plant de novo variants of every class into the REAL
haplotagged reads of a trio and check that the six-haplotype review recovers class, parent of origin and the
transmitted-haplotype evidence. Truth is known by construction; nothing here touches folds or thresholds.

Scenarios (per class SNV, INDEL, SV, TR):
  G    germline DNM on child haplotype h          child tagged reads of h: fraction 1.0, untagged 0.5
  CM   child postzygotic mosaic                    child tagged reads of h: fraction c, untagged c/2
  PM   parental mosaic, transmitted                as G, plus the ORIGIN parent's TRANSMITTED haplotype at fraction m
  IM   inherited, missed in the parent             as G, plus the origin parent's transmitted haplotype at 1.0 (untagged 0.5)
Parent of origin of h and the parent's transmitted haplotype come from the M1 label tables at the site, exactly as
the review reads them; sites without labels are not planted (spike-ins test the review, not the phasing).

Three steps, so the expensive one (apply) can be inspected before the review runs:
  plan      choose sites that pass a read-level QC (unanimous readable base on every haplotype, CIGAR-clean window
            for indels/SVs, consistent repeat lengths for TR) and write plan.tsv + candidates.tsv + truth.tsv
  apply     write per-sample SLICE BAMs (all reads in +-pad of every site, edited where the plan says) + .bai
  evaluate  join truth with the reviewed (and likelihood-scored) evidence table -> recovery per class x scenario

The edits themselves are `sim/edit.py` (pure functions, unit-tested); pysam is imported lazily here.
Identifiers never appear in code; read names are left alone (the review hashes them).
"""
from __future__ import annotations

import csv
import json
import os
import random
import statistics as st
from collections import defaultdict
from dataclasses import asdict, dataclass, fields
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from ..records import CandidateRecord, write_candidates
from . import edit as E

SCENARIOS = ("G", "CM", "PM", "IM")
EXPECTED_CLASS = {"G": "germline_DNM_phased", "CM": "child_postzygotic_mosaic",
                  "PM": "parental_mosaic_transmitted", "IM": "inherited_missed_in_parent"}
SUBTYPES = {"SNV": ("SNV",), "INDEL": ("DEL", "INS"), "SV": ("DEL", "INS"), "TR": ("EXP",)}
LENGTHS = {("INDEL", "DEL"): (1, 2, 4, 8), ("INDEL", "INS"): (1, 3, 6, 12),
           ("SV", "DEL"): (200, 500, 1500), ("SV", "INS"): (100, 250, 400), ("TR", "EXP"): (3, 5, 8)}   # TR in motif units
FRACTIONS = {"CM": (0.15, 0.3), "PM": (0.1, 0.25)}
BASES = "ACGT"


@dataclass
class PlanRow:
    variant_id: str
    chrom: str
    pos: int                      # 1-based anchor (VCF POS); SNV: the base itself
    variant_class: str            # SNV | INDEL | SV | TR
    subtype: str                  # SNV | DEL | INS | EXP
    length: int                   # bp (TR: bp = units * motif_unit_bp)
    ref: str                      # SNV/INDEL: VCF-style REF; SV/TR: "."
    alt: str                      # SNV/INDEL: VCF-style ALT; INS/SV-INS: inserted sequence; TR: motif; SV-DEL: "."
    scenario: str
    child_hap: int
    child_frac: float
    parent: str                   # F | M | .
    parent_hap: int               # transmitted haplotype of the origin parent (0 when unused)
    parent_frac: float
    expected_poo: str             # paternal | maternal
    expected_class: str
    child_ps: int
    parent_ps: int
    trid: str = "."
    motif_unit_bp: int = 0
    child_base_al: str = "."      # TR: per-haplotype repeat lengths before planting, "h1,h2"
    father_al: str = "."
    mother_al: str = "."
    seed: int = 0
    locus_end: int = 0            # TR: 1-based locus end (insertion goes to the locus midpoint)


PLAN_COLS = [f.name for f in fields(PlanRow)]


# ----------------------------------------------------------------------------------------------
# helpers on reads (pysam objects in, plain data out)
# ----------------------------------------------------------------------------------------------
def aln_of(read) -> E.Aln:
    q = read.query_qualities
    return E.Aln(read.reference_start, list(read.cigartuples or []), read.query_sequence or "", list(q) if q is not None else None)


def write_back(read, a: E.Aln) -> None:
    read.cigartuples = a.cigar
    read.query_sequence = a.seq
    read.query_qualities = a.qual
    for tag in ("MD",):
        if read.has_tag(tag):
            read.set_tag(tag, None)


def _hp(read) -> Optional[int]:
    v = read.get_tag("HP") if read.has_tag("HP") else None
    return v if v in (1, 2) else None


def _ps(read) -> Optional[int]:
    return read.get_tag("PS") if read.has_tag("PS") else None


def _primary(read, min_mapq: int) -> bool:
    return not (read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate) and read.mapping_quality >= min_mapq


def clean_window(a: E.Aln, s0: int, e0: int) -> bool:
    """Aligned through [s0, e0) by M ops only (no indel, no clip inside)."""
    ref = a.reference_start
    if ref > s0 or a.reference_end < e0:
        return False
    for op, L in a.cigar:
        if op in (E.M, E.EQ, E.X):
            nxt = ref + L
            if ref <= s0 and nxt >= e0:
                return True
            ref = nxt
        elif op in (E.D, E.N):
            if ref < e0 and ref + L > s0:
                return False
            ref += L
        elif op == E.I:
            if s0 < ref < e0:
                return False
        if ref >= e0:
            break
    return False


def consensus_base(alns: Sequence[E.Aln], pos0: int) -> Tuple[Optional[str], int, float]:
    """(majority base, n readable, majority fraction) at a 0-based position."""
    c: Dict[str, int] = defaultdict(int)
    for a in alns:
        qi = E.query_index_at(a, pos0)
        if qi is not None:
            c[a.seq[qi].upper()] += 1
    n = sum(c.values())
    if not n:
        return None, 0, 0.0
    b = max(c, key=c.get)
    return b, n, c[b] / n


def consensus_seq(alns: Sequence[E.Aln], s0: int, e0: int) -> Optional[str]:
    out = []
    for p in range(s0, e0):
        b, n, f = consensus_base(alns, p)
        if b is None or f < 0.9 or b == "N":
            return None
        out.append(b)
    return "".join(out)


# ----------------------------------------------------------------------------------------------
# planning
# ----------------------------------------------------------------------------------------------
@dataclass
class SiteQC:
    k: int = 5                    # readable reads per haplotype per sample
    unanimity: float = 0.97
    clean_frac: float = 0.9       # fraction of reads (all samples) with an indel-free window
    margin: int = 30
    min_mapq: int = 20
    pad: int = 2500               # slice half-width around a site (>= review windows + read overhang)
    spacing: int = 60000          # min distance between planted sites (no read spans two)
    max_tries_per_site: int = 400


def _reads_at(bam, chrom: str, s0: int, e0: int, qc: SiteQC):
    return [r for r in bam.fetch(chrom, max(0, s0), e0) if _primary(r, qc.min_mapq) and r.reference_start <= s0 and (r.reference_end or 0) >= e0]


def _labels_at(labels, chrom: str, pos1: int, child_reads, parent_reads: Dict[str, list]):
    """(child_hap1_is, child_ps, {F: (transmitted, ps), M: ...}) from the majority PS of the spanning reads."""
    def maj_ps(reads):
        c: Dict[int, int] = defaultdict(int)
        for r in reads:
            p = _ps(r)
            if p is not None and _hp(r) is not None:
                c[p] += 1
        return max(c, key=c.get) if c else None
    cps = maj_ps(child_reads)
    h1 = labels.child_hap1_is(chrom, cps, pos1)
    par = {}
    for role in ("F", "M"):
        pps = maj_ps(parent_reads[role])
        par[role] = (labels.parent_transmitted(role, chrom, pps, pos1), pps)
    return h1, cps, par


def _scenario_grid(rng: random.Random, n_per: int) -> List[Tuple[str, str, str, int, float, float]]:
    """(class, subtype, scenario, length, child_frac, parent_frac) x n_per, shuffled."""
    grid = []
    for cls, subs in SUBTYPES.items():
        for sub in subs:
            lengths = LENGTHS.get((cls, sub), (1,))
            for sc in SCENARIOS:
                for i in range(n_per):
                    L = lengths[i % len(lengths)]
                    cf = 1.0 if sc in ("G", "PM", "IM") else FRACTIONS["CM"][i % 2]
                    pf = {"G": 0.0, "CM": 0.0, "PM": FRACTIONS["PM"][i % 2], "IM": 1.0}[sc]
                    grid.append((cls, sub, sc, L, cf, pf))
    rng.shuffle(grid)
    return grid


def plan_sites(bams: Dict[str, object], tr_bams: Dict[str, object], labels, regions: List[Tuple[str, int, int]],
               n_per: int, seed: int, family: str, child: str, qc: SiteQC, trgt_vcf: Optional[str] = None,
               mask=None, log=None, trgt_vcf_index: Optional[str] = None) -> Tuple[List[PlanRow], List[CandidateRecord], Dict[str, int]]:
    rng = random.Random(seed)
    grid = _scenario_grid(rng, n_per)
    used: Dict[str, List[int]] = defaultdict(list)
    plan: List[PlanRow] = []
    cands: List[CandidateRecord] = []
    stats: Dict[str, int] = defaultdict(int)
    tr_loci = _tr_loci(trgt_vcf, regions, child, trgt_vcf_index) if trgt_vcf else []
    rng.shuffle(tr_loci)
    tr_iter = iter(tr_loci)
    stats["planned"] = len(grid)
    for n_item, (cls, sub, sc, L, cf, pf) in enumerate(grid):
        placed = False
        for _ in range(qc.max_tries_per_site):
            stats["tries"] += 1
            if cls == "TR":
                locus = next(tr_iter, None)
                if locus is None:
                    break
                if any(abs(locus["start"] - u) < qc.spacing for u in used[locus["chrom"]]):
                    stats["skip_spacing"] += 1
                    continue
                site = _try_tr_site(locus, tr_bams, labels, rng, qc, sc, L, cf, pf, family, child, n_item, seed)
            else:
                chrom, rs, re_ = regions[rng.randrange(len(regions))]
                pos1 = rng.randrange(rs + qc.pad, re_ - qc.pad)
                if any(abs(pos1 - u) < qc.spacing for u in used[chrom]):
                    stats["skip_spacing"] += 1
                    continue
                if mask is not None and mask.overlap_bp(chrom, pos1, pos1 + max(1, L)) > 0:
                    stats["skip_mask"] += 1
                    continue
                site = _try_site(chrom, pos1, cls, sub, L, bams, labels, rng, qc, sc, cf, pf, family, child, n_item, seed, stats)
            if site is None:
                continue
            row, rec = site
            used[row.chrom].append(row.pos)
            plan.append(row); cands.append(rec)
            placed = True
            break
        stats["placed" if placed else "unplaced"] += 1
        if log and n_item % 25 == 0:
            log("spike plan: %d/%d items, %d placed" % (n_item + 1, len(grid), len(plan)))
    return plan, cands, dict(stats)


def _try_site(chrom, pos1, cls, sub, L, bams, labels, rng, qc, sc, cf, pf, family, child, n_item, seed, stats):
    pos0 = pos1 - 1
    s0, e0 = pos0 - qc.margin, pos0 + max(1, L) + qc.margin + 1
    reads = {role: _reads_at(b, chrom, s0, e0, qc) for role, b in bams.items()}
    alns = {role: [aln_of(r) for r in rs] for role, rs in reads.items()}
    # every haplotype of every sample: >= k readable reads, unanimous base at the anchor
    for role in ("C", "F", "M"):
        by_hp: Dict[int, List[E.Aln]] = defaultdict(list)
        for r, a in zip(reads[role], alns[role]):
            h = _hp(r)
            if h:
                by_hp[h].append(a)
        for h in (1, 2):
            b, n, f = consensus_base(by_hp[h], pos0)
            if n < qc.k or f < qc.unanimity or b is None or b == "N":
                stats["skip_hap_qc"] += 1
                return None
    allr = [a for role in alns for a in alns[role]]
    if cls != "SNV":
        clean = sum(1 for a in allr if clean_window(a, s0, e0))
        if not allr or clean / len(allr) < qc.clean_frac:
            stats["skip_window"] += 1
            return None
    h1, cps, par = _labels_at(labels, chrom, pos1, reads["C"], {"F": reads["F"], "M": reads["M"]})
    if h1 is None or cps is None:
        stats["skip_unoriented"] += 1
        return None
    child_hap = rng.choice((1, 2))
    origin = ("F" if h1 == "P" else "M") if child_hap == 1 else ("M" if h1 == "P" else "F")
    t_hap, pps = par[origin]
    if sc in ("PM", "IM") and (t_hap is None or pps is None):
        stats["skip_untransmitted"] += 1
        return None
    # the variant itself, from the reads' consensus (reference-free)
    ref_seq = consensus_seq(allr, pos0, pos0 + (L + 1 if (cls, sub) in (("INDEL", "DEL"), ("SV", "DEL")) else 1))
    if ref_seq is None:
        stats["skip_consensus"] += 1
        return None
    if cls == "SNV":
        alt = rng.choice([b for b in BASES if b != ref_seq[0]])
        ref, vid, end = ref_seq[0], "spike:%s:%d:%s:%s" % (chrom, pos1, ref_seq[0], alt), pos1
        payload = {"spike": True}
        rec_ref, rec_alt = ref, alt
    elif cls == "INDEL":
        if sub == "DEL":
            ref, alt = ref_seq, ref_seq[0]
        else:
            ins = "".join(rng.choice(BASES) for _ in range(L))
            ref, alt = ref_seq[0], ref_seq[0] + ins
        vid, end = "spike:%s:%d:%s:%s" % (chrom, pos1, ref, alt), pos1 + (L if sub == "DEL" else 0)
        payload = {"spike": True}
        rec_ref, rec_alt = ref, alt
    else:  # SV
        if sub == "DEL":
            svlen, end, alt = -L, pos1 + L, "<DEL>"
        else:
            svlen, end, alt = L, pos1, "".join(rng.choice(BASES) for _ in range(L))
        vid = "spike:%s:%d:%s:%d" % (chrom, pos1, sub, L)
        payload = {"svtype": sub, "svlen": svlen, "end": end, "homlen": None, "imprecise": False, "mateid": None, "svclaim": None,
                   "caller_id": vid, "child_cn": None, "father_cn": None, "mother_cn": None, "insseq_len": L if sub == "INS" else None, "spike": True}
        rec_ref, rec_alt = ref_seq[0], (alt if sub == "DEL" else "<INS:%dbp>" % L)
    row = PlanRow(variant_id=vid, chrom=chrom, pos=pos1, variant_class=cls, subtype=sub, length=L, ref=rec_ref if cls != "SV" else ".",
                  alt=(alt if cls != "SV" or sub == "INS" else "."), scenario=sc, child_hap=child_hap, child_frac=cf,
                  parent=origin if sc in ("PM", "IM") else ".", parent_hap=t_hap or 0, parent_frac=pf,
                  expected_poo="paternal" if origin == "F" else "maternal", expected_class=EXPECTED_CLASS[sc],
                  child_ps=cps, parent_ps=pps or 0, seed=seed)
    rec = CandidateRecord(family_id=family, sample_id=child, variant_id=vid, chrom=chrom, start=pos1, end=end, ref=rec_ref, alt=rec_alt,
                          variant_class=cls, caller="spike", class_payload=payload)
    return row, rec


def _tr_loci(trgt_vcf: str, regions, child: str, index: Optional[str] = None) -> List[dict]:
    import pysam
    out = []
    vf = pysam.VariantFile(trgt_vcf, index_filename=index) if index else pysam.VariantFile(trgt_vcf)
    if child not in list(vf.header.samples):
        return out
    for chrom, s, e in regions:
        try:
            it = vf.fetch(chrom, s, e)
        except ValueError:
            continue
        for rec in it:
            motifs = rec.info.get("MOTIFS")
            motifs = ",".join(motifs) if isinstance(motifs, tuple) else (motifs or "")
            unit = max(1, min((len(x) for x in motifs.split(",") if x), default=1))
            out.append({"trid": rec.info.get("TRID"), "chrom": chrom, "start": rec.pos, "end": rec.stop, "motifs": motifs, "unit": unit})
    vf.close()
    return out


def _try_tr_site(locus, tr_bams, labels, rng, qc, sc, n_units, cf, pf, family, child, n_item, seed):
    chrom, start1, end1, trid, unit = locus["chrom"], locus["start"], locus["end"], locus["trid"], locus["unit"]
    from ..evidence.readers import tr_read_length
    lengths: Dict[str, Dict[int, List[int]]] = {r: defaultdict(list) for r in ("C", "F", "M")}
    reads: Dict[str, list] = {}
    for role, bam in tr_bams.items():
        rs = []
        for r in bam.fetch(chrom, max(0, start1 - 1), end1 + 1):
            if not _primary(r, 0) or not r.has_tag("TR") or r.get_tag("TR") != trid:
                continue
            L = tr_read_length(r)
            if L is None:
                continue
            rs.append(r)
            h = _hp(r)
            if h:
                lengths[role][h].append(L)
        reads[role] = rs
    for role in ("C", "F", "M"):
        for h in (1, 2):
            v = lengths[role][h]
            if len(v) < qc.k or (len(v) > 1 and st.pstdev(v) > 2.0):
                return None
    al = {role: (int(st.median(lengths[role][1])), int(st.median(lengths[role][2]))) for role in ("C", "F", "M")}
    h1, cps, par = _labels_at(labels, chrom, start1, reads["C"], {"F": reads["F"], "M": reads["M"]})
    if h1 is None or cps is None:
        return None
    child_hap = rng.choice((1, 2))
    origin = ("F" if h1 == "P" else "M") if child_hap == 1 else ("M" if h1 == "P" else "F")
    t_hap, pps = par[origin]
    if sc in ("PM", "IM") and (t_hap is None or pps is None):
        return None
    delta = n_units * unit
    base = al["C"][child_hap - 1]
    planted = base + delta
    other = al["C"][2 - child_hap]
    # separable from every competing allele by more than the tolerance (support_tr caps it at half the gap)
    if min(abs(planted - o) for o in (other, *al["F"], *al["M"])) < 2 * unit + 2:
        return None
    motif = locus["motifs"].split(",")[0] or "A"
    vid = "spike:%s:EXP:%d" % (trid, delta)
    payload = {"trid": trid, "motifs": locus["motifs"], "struc": None, "motif_unit_bp": unit,
               "child_AL": [other, planted], "father_AL": list(al["F"]), "mother_AL": list(al["M"]),
               "outlier_allele_idx": 1, "direction": "expansion", "delta_bp": planted - max(al["F"] + al["M"]),
               "delta_units": round((planted - max(al["F"] + al["M"])) / unit, 2), "spike": True}
    row = PlanRow(variant_id=vid, chrom=chrom, pos=start1, variant_class="TR", subtype="EXP", length=delta, ref=".", alt=motif,
                  scenario=sc, child_hap=child_hap, child_frac=cf, parent=origin if sc in ("PM", "IM") else ".",
                  parent_hap=t_hap or 0, parent_frac=pf, expected_poo="paternal" if origin == "F" else "maternal",
                  expected_class=EXPECTED_CLASS[sc], child_ps=cps, parent_ps=pps or 0, trid=trid, motif_unit_bp=unit,
                  child_base_al="%d,%d" % al["C"], father_al="%d,%d" % al["F"], mother_al="%d,%d" % al["M"], seed=seed, locus_end=end1)
    rec = CandidateRecord(family_id=family, sample_id=child, variant_id=vid, chrom=chrom, start=start1, end=end1, ref="<TR>",
                          alt="<AL=%d>" % planted, variant_class="TR", caller="spike", class_payload=payload)
    return row, rec


# ----------------------------------------------------------------------------------------------
# applying the plan to reads
# ----------------------------------------------------------------------------------------------
def _edit_one(a: E.Aln, row: PlanRow, rng: random.Random) -> Optional[E.Aln]:
    pos0 = row.pos - 1
    if row.variant_class == "SNV":
        return E.apply_snv(a, pos0, row.alt)
    if row.variant_class == "TR":
        motif = row.alt if row.alt not in (".", "") else "A"
        ins = (motif * (row.length // len(motif) + 1))[:row.length]
        mid = (row.pos + max(row.pos, row.locus_end)) // 2
        return E.apply_insertion(a, mid - 1, ins, margin=10)
    if row.subtype == "DEL":
        return E.apply_deletion(a, row.pos, row.length)             # deleted bases are 0-based [pos, pos+L): right after the anchor
    ins = row.alt[1:] if row.variant_class == "INDEL" else row.alt
    return E.apply_insertion(a, pos0, ins)


def _decide(read, role: str, row: PlanRow, rng: random.Random) -> bool:
    """Does this read get the edit?"""
    h = _hp(read)
    if role == "C":
        if h == row.child_hap:
            return rng.random() < row.child_frac
        if h is None:
            return rng.random() < row.child_frac * 0.5
        return False
    if role == row.parent and row.parent_frac > 0:
        if h == row.parent_hap:
            return rng.random() < row.parent_frac
        if h is None:
            return rng.random() < row.parent_frac * 0.5
    return False


def apply_plan(plan: List[PlanRow], bams: Dict[str, object], out_dir: str, tag: str, qc: SiteQC, seed: int,
               tr: bool = False, log=None) -> Dict[str, int]:
    """Write <out_dir>/<role>.<tag>.bam with every primary/secondary/supplementary read in +-pad of each planted site of
    the matching class group (tr=False: SNV/INDEL/SV from the genome BAMs; tr=True: TR from the TRGT spanning BAMs),
    edited where the plan says. Sites are >= spacing apart so no read is written twice."""
    import pysam
    rows = [r for r in plan if (r.variant_class == "TR") == tr]
    stats: Dict[str, int] = defaultdict(int)
    rng = random.Random(seed + (1 if tr else 0))
    os.makedirs(out_dir, exist_ok=True)
    for role, bam in bams.items():
        out_path = os.path.join(out_dir, "%s.%s.bam" % (role, tag))
        rows.sort(key=lambda r: (bam.get_tid(r.chrom), r.pos))
        with pysam.AlignmentFile(out_path, "wb", template=bam) as out:
            for row in rows:
                s0, e0 = row.pos - 1 - qc.pad, row.pos - 1 + row.length + qc.pad
                for read in bam.fetch(row.chrom, max(0, s0), e0):
                    if read.is_unmapped:
                        continue
                    stats["reads_%s" % role] += 1
                    if tr and (not read.has_tag("TR") or read.get_tag("TR") != row.trid):
                        out.write(read); continue
                    if _primary(read, 0) and _decide(read, role, row, rng):
                        a = _edit_one(aln_of(read), row, rng)
                        if a is not None:
                            write_back(read, a)
                            stats["edited_%s" % role] += 1
                        else:
                            stats["edit_refused_%s" % role] += 1
                    out.write(read)
        pysam.index(out_path)
        if log:
            log("spike apply: %s -> %s (%d reads, %d edited)" % (role, os.path.basename(out_path), stats["reads_%s" % role], stats["edited_%s" % role]))
    return dict(stats)


# ----------------------------------------------------------------------------------------------
# tables
# ----------------------------------------------------------------------------------------------
def write_plan(rows: List[PlanRow], path: str) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=PLAN_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(asdict(r))


def read_plan(path: str) -> List[PlanRow]:
    out = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            kw = {}
            for f in fields(PlanRow):
                v = r[f.name]
                kw[f.name] = (int(v) if f.type == "int" else float(v) if f.type == "float" else v)
            out.append(PlanRow(**kw))
    return out


# ----------------------------------------------------------------------------------------------
# evaluation
# ----------------------------------------------------------------------------------------------
def evaluate(plan: List[PlanRow], evidence_rows: Iterable[dict], k: int = 5) -> Tuple[List[dict], List[dict]]:
    """Per planted site: observed? class recovered? parent of origin recovered? Then a summary per class x scenario.
    'observable' = all six haplotypes observed at k (hap_obs_k<k>); recovery is reported both overall and among
    observable sites, because a spike-in on an unobservable haplotype tests the phasing coverage, not the review."""
    ev = {r["variant_id"]: r for r in evidence_rows}
    present = {r.get("variant_class") for r in ev.values()}
    if present:
        plan = [p for p in plan if p.variant_class in present]       # one evidence table per class group
    post_key = {"G": "lik_post_germline", "CM": "lik_post_child_mosaic", "PM": "lik_post_parental_mosaic", "IM": "lik_post_inherited"}
    per_site = []
    for p in plan:
        r = ev.get(p.variant_id)
        d = dict(variant_id=p.variant_id, variant_class=p.variant_class, subtype=p.subtype, length=p.length, scenario=p.scenario,
                 child_frac=p.child_frac, parent_frac=p.parent_frac, expected_class=p.expected_class, expected_poo=p.expected_poo,
                 reviewed=int(r is not None))
        if r is None:
            per_site.append(d); continue
        obs = r.get("hap_obs_k%d" % k, "")
        cls = r.get("phase_class", "")
        d.update(phase_class=cls, parent_of_origin=r.get("parent_of_origin", ""), rule_score=r.get("rule_score", ""),
                 phase_score=r.get("phase_score", ""), flags=r.get("flags", ""), hap_obs=obs,
                 observable=int(obs not in ("", None) and int(float(obs)) == 6),
                 class_ok=int(cls == p.expected_class),
                 class_ok_lenient=int(cls == p.expected_class or (p.scenario == "G" and cls == "germline_DNM_unphased")),
                 poo_ok=int(r.get("parent_of_origin") == p.expected_poo),
                 post_expected=r.get(post_key[p.scenario], ""),
                 c_alt_hapA=r.get("c_alt_hapA", ""), c_dp_hapA=r.get("c_dp_hapA", ""), t_alt_reads=r.get("t_alt_reads", ""), t_dp=r.get("t_dp", ""))
        per_site.append(d)
    summary = []
    groups: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    for d in per_site:
        groups[(d["variant_class"], d["scenario"])].append(d)
    for (cls, sc), rows in sorted(groups.items()):
        rev = [d for d in rows if d["reviewed"]]
        obs = [d for d in rev if d.get("observable")]
        def frac(xs, key):
            return round(sum(d[key] for d in xs) / len(xs), 3) if xs else None
        def mean(xs, key):
            v = [float(d[key]) for d in xs if d.get(key) not in ("", None)]
            return round(st.mean(v), 3) if v else None
        called: Dict[str, int] = defaultdict(int)
        for d in rev:
            called[d.get("phase_class") or "."] += 1
        summary.append(dict(variant_class=cls, scenario=sc, expected_class=EXPECTED_CLASS[sc], n_planted=len(rows), n_reviewed=len(rev),
                            n_observable=len(obs), class_ok_all=frac(rev, "class_ok"), class_ok_lenient_all=frac(rev, "class_ok_lenient"),
                            class_ok_observable=frac(obs, "class_ok"), poo_ok_all=frac(rev, "poo_ok"), poo_ok_observable=frac(obs, "poo_ok"),
                            mean_phase_score=mean(rev, "phase_score"), mean_phase_score_observable=mean(obs, "phase_score"),
                            mean_post_expected=mean(rev, "post_expected"), mean_rule_score=mean(rev, "rule_score"),
                            called_classes=";".join("%s:%d" % kv for kv in sorted(called.items(), key=lambda kv: -kv[1]))))
    return per_site, summary


def write_rows(rows: List[dict], path: str) -> None:
    if not rows:
        open(path, "w").close(); return
    cols = list(rows[0].keys())
    for r in rows[1:]:
        for c in r:
            if c not in cols:
                cols.append(c)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: ("" if r.get(c) is None else r.get(c)) for c in cols})


def parse_regions(spec: str) -> List[Tuple[str, int, int]]:
    out = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        chrom, rng = part.split(":")
        s, e = rng.replace(",", "").split("-")
        out.append((chrom, int(s), int(e)))
    return out
