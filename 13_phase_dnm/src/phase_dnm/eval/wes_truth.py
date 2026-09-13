"""P27 arm 2 — WES-confirmed external truth for small variants (R10). The SPARK iWES v3 DeepVariant/GLnexus pVCF
(per chromosome, 142k samples, FORMAT GT RNC DP AD GQ PL) is queried at the cohort's exonic candidate sites for the
101 cohort samples it contains, and each candidate of a child gets a label from the WES trio genotypes:

    1  WES de novo:   child het (GQ >= min_gq, DP >= min_dp, alt reads >= 2), both parents hom-ref (GQ >= min_gq,
                      DP >= min_dp, alt reads == 0)
    0  WES refutes:   child hom-ref at DP >= neg_dp (the allele is not in the child's exome reads), OR a parent carries
                      the allele (GT has it, or alt reads >= 3) - inherited, not de novo
   -1  not evaluable: not in the capture target, sample absent from WES, low depth / low GQ, allele representation mismatch

"Exonic" = inside the capture target BED(s) (+/- 50 bp versions), which is where the exome has reads; everything else is
-1 by construction. Sites are matched on the canonical (chrom, pos, deleted, inserted) key (concordance.norm_allele) after
`bcftools norm -m -any`. Output: wes_truth/<child>.snv_indel.wes.tsv (variant_id, label, child_gt, child_gq, child_dp,
child_ad, father_gt, mother_gt, father_ad, mother_ad, reason) - consumed by eval/external.evaluate_labelled.
Identifiers are values, never code.
"""
from __future__ import annotations

import csv
import glob
import os
import subprocess
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Set, Tuple

from ..concordance import norm_allele
from ..records import read_candidates


# ----------------------------------------------------------------------------------------------
# targets and sites
# ----------------------------------------------------------------------------------------------
class Targets:
    def __init__(self, beds: Iterable[str]):
        import bisect as _b
        self._b = _b
        iv: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
        for bed in beds:
            op = open
            if bed.endswith(".gz"):
                import gzip
                op = lambda p: gzip.open(p, "rt")
            with op(bed) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 3 or f[0].startswith(("#", "track", "browser")):
                        continue
                    try:
                        iv[f[0]].append((int(f[1]), int(f[2])))
                    except ValueError:
                        continue
        self.iv: Dict[str, List[Tuple[int, int]]] = {}
        for c, lst in iv.items():
            lst.sort()
            merged: List[Tuple[int, int]] = []
            for s, e in lst:
                if merged and s <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], e))
                else:
                    merged.append((s, e))
            self.iv[c] = merged
        self.starts = {c: [s for s, _ in v] for c, v in self.iv.items()}

    def covers(self, chrom: str, pos1: int) -> bool:
        v = self.iv.get(chrom)
        if not v:
            return False
        i = self._b.bisect_right(self.starts[chrom], pos1 - 1) - 1
        return i >= 0 and v[i][0] <= pos1 - 1 < v[i][1]


def exonic_sites(cand_paths: Iterable[str], targets: Targets) -> Dict[str, Set[int]]:
    """chrom -> set of 1-based positions of candidates inside the capture target (all children pooled)."""
    out: Dict[str, Set[int]] = defaultdict(set)
    for p in cand_paths:
        for rec in read_candidates(p):
            if targets.covers(rec.chrom, rec.start):
                out[rec.chrom].add(rec.start)
    return out


def write_regions(sites: Dict[str, Set[int]], chrom: str, path: str, pad: int = 1) -> int:
    n = 0
    with open(path, "w") as fh:
        for pos in sorted(sites.get(chrom, ())):
            fh.write("%s\t%d\t%d\n" % (chrom, max(0, pos - 1 - pad), pos + pad))
            n += 1
    return n


# ----------------------------------------------------------------------------------------------
# extraction (bcftools) -> per-site trio genotypes for the cohort samples
# ----------------------------------------------------------------------------------------------
def extract_chrom(pvcf: str, regions_bed: str, samples_file: str, out_tsv: str, bcftools: str) -> int:
    """bcftools view -R regions -S samples | norm -m -any | query -> TSV: chrom pos ref alt then per sample GT:GQ:DP:AD."""
    cmd = ("%s view -R %s -S %s --force-samples -Ou %s | %s norm -m -any -Ou | "
           "%s query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT:%%GQ:%%DP:%%AD]\\n' > %s" % (bcftools, regions_bed, samples_file, pvcf, bcftools, bcftools, out_tsv))
    subprocess.run(cmd, shell=True, check=True)
    with open(out_tsv) as fh:
        return sum(1 for _ in fh)


def sample_order(pvcf: str, samples_file: str, bcftools: str) -> List[str]:
    q = subprocess.run("%s view -S %s --force-samples -h %s | tail -1" % (bcftools, samples_file, pvcf), shell=True, check=True, capture_output=True, text=True).stdout
    return q.rstrip("\n").split("\t")[9:]


# ----------------------------------------------------------------------------------------------
# labelling (pure)
# ----------------------------------------------------------------------------------------------
def _parse(cell: str) -> Tuple[Optional[Tuple[int, ...]], Optional[int], Optional[int], Optional[List[int]]]:
    f = cell.split(":")
    if len(f) < 4:
        return None, None, None, None
    gt = None if f[0] in (".", "./.", ".|.") else tuple(int(x) for x in f[0].replace("|", "/").split("/") if x != ".")
    def i(x):
        try:
            return int(x)
        except ValueError:
            return None
    ad = None if f[3] in (".", "") else [i(x) for x in f[3].split(",")]
    return gt, i(f[1]), i(f[2]), ad


def label_trio(child: str, father: str, mother: str, alt_index: int, min_gq: int = 20, min_dp: int = 10, neg_dp: int = 20) -> Tuple[int, str]:
    """(label, reason) from the three GT:GQ:DP:AD cells of a split (single-alt) record; alt_index is 1 after norm -m -any."""
    cg, cq, cd, cad = _parse(child); fg, fq, fd, fad = _parse(father); mg, mq, md, mad = _parse(mother)
    if cg is None or cd is None:
        return -1, "child_uncalled"
    c_alt = cad[alt_index] if cad and len(cad) > alt_index and cad[alt_index] is not None else None
    f_alt = fad[alt_index] if fad and len(fad) > alt_index and fad[alt_index] is not None else None
    m_alt = mad[alt_index] if mad and len(mad) > alt_index and mad[alt_index] is not None else None
    child_has = alt_index in cg
    # a parent carries the allele -> inherited, not de novo (label 0)
    for g, a, who in ((fg, f_alt, "father"), (mg, m_alt, "mother")):
        if (g is not None and alt_index in g) or (a is not None and a >= 3):
            return 0, "parent_carries_%s" % who
    if child_has:
        if (cq or 0) < min_gq or cd < min_dp or (c_alt is not None and c_alt < 2):
            return -1, "child_low_quality"
        for g, q, d, a, who in ((fg, fq, fd, f_alt, "father"), (mg, mq, md, m_alt, "mother")):
            if g is None or d is None or (q or 0) < min_gq or d < min_dp:
                return -1, "%s_low_quality" % who
            if a is not None and a > 0:
                return -1, "%s_alt_read" % who
        return 1, "wes_de_novo"
    # child does not carry the allele in WES
    if all(x == 0 for x in cg) and cd >= neg_dp and (cq or 0) >= min_gq:
        return 0, "child_hom_ref"
    return -1, "child_inconclusive"


def label_children(extract_tsv: str, order: List[str], trios: Dict[str, Tuple[str, str]], cand_paths: Dict[str, str],
                   out_dir: str, chrom_filter: Optional[str] = None) -> Dict[str, Dict[str, int]]:
    """For every child with a candidates table, label its candidates on this chromosome's extraction."""
    idx = {s: i for i, s in enumerate(order)}
    wes: Dict[Tuple, List[str]] = {}
    with open(extract_tsv) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            try:
                key = (f[0],) + norm_allele(int(f[1]), f[2], f[3])
            except ValueError:
                continue
            wes[key] = f[4:]
    os.makedirs(out_dir, exist_ok=True)
    stats: Dict[str, Dict[str, int]] = {}
    for child, (fa, mo) in trios.items():
        if child not in cand_paths or child not in idx or fa not in idx or mo not in idx:
            continue
        cnt: Dict[str, int] = defaultdict(int)
        out = os.path.join(out_dir, "%s.snv_indel.wes.%s.tsv" % (child, chrom_filter or "all"))
        with open(out, "w", newline="") as fh:
            fh.write("variant_id\tlabel\treason\tchild_cell\tfather_cell\tmother_cell\n")
            for rec in read_candidates(cand_paths[child]):
                if chrom_filter and rec.chrom != chrom_filter:
                    continue
                key = (rec.chrom,) + norm_allele(rec.start, rec.ref, rec.alt)
                cells = wes.get(key)
                if cells is None:
                    lab, why, c, f_, m = -1, "not_in_wes", "", "", ""
                else:
                    c, f_, m = cells[idx[child]], cells[idx[fa]], cells[idx[mo]]
                    lab, why = label_trio(c, f_, m, 1)
                cnt[why] += 1
                fh.write("%s\t%d\t%s\t%s\t%s\t%s\n" % (rec.variant_id, lab, why, c, f_, m))
        stats[child] = dict(cnt)
    return stats
