"""Cohort / population annotation of the candidate sites (DESIGN P26): the columns the H2/H3 heuristic arms and the
cohort features need, computed once per class group over the union of every child's candidates.

  gnomad_af                 site property; `slivar gnotate` with the lab's gnomAD zip on a sites-only VCF of the candidates
                            (the raw joint VCF carries only AF/AQ/AC/AN; the WDL's tertiary VCF is already filtered)
  cohort_AC_loo             alternate-allele count among the 65 unaffected founders EXCLUDING the founders of the child's
  pon_founder_recurrence_loo  family and of the parents' family (identical for a real trio; symmetric for a synthetic one);
                            genotypes from the cohort BCF (small variants) / cohort SV VCF (SVs, matched by sawfish id),
                            both normalised with `bcftools norm -m -any` so a split multi-allelic record matches
  sib_shared                1 when the candidate allele is carried by a sibling in the family joint VCF (the two quads)

Outputs annot/<child>.<class>.annot.tsv: variant_id, gnomad_af, cohort_AC_loo, cohort_AN_loo, pon_founder_recurrence_loo,
sib_shared. The feature extractor merges them by variant_id. bcftools / slivar run as subprocesses (paths from the env);
everything else is plain Python. Identifiers are values, never code.
"""
from __future__ import annotations

import csv
import gzip
import os
import subprocess
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Set, Tuple

from .records import CandidateRecord, read_candidates
from .concordance import norm_allele


def founders(manifest_rows: Iterable[dict]) -> Dict[str, str]:
    """unaffected founders -> family_id (the 65 of the panel of normals; affected founders excluded as in the PON)."""
    out = {}
    for r in manifest_rows:
        if r.get("father_id") in ("0", "", ".") and r.get("mother_id") in ("0", "", ".") and str(r.get("affected", "1")) != "2":
            out[r["sample_id"]] = r["family_id"]
    return out


# ----------------------------------------------------------------------------------------------
# sites
# ----------------------------------------------------------------------------------------------
def union_sites(cand_paths: Iterable[str], class_group: str) -> Tuple[List[Tuple[str, int, str, str]], Dict[str, Set[str]]]:
    """Distinct small-variant sites (chrom, pos, ref, alt) or SV ids across candidate tables; also which ids each child has."""
    sites: Set[Tuple[str, int, str, str]] = set()
    for p in cand_paths:
        for rec in read_candidates(p):
            if class_group == "snv_indel":
                sites.add((rec.chrom, rec.start, rec.ref, rec.alt))
            elif class_group == "sv":
                sites.add((rec.chrom, rec.start, rec.class_payload.get("caller_id") or rec.variant_id, rec.class_payload.get("svtype") or ""))
            else:
                # TR: the locus, keyed by TRID. Until 2026-09-15 this branch was absent, so every TR annotation column
                # came out empty and both the P24 rarity gate and the P15 population rules were silently vacuous for TR.
                sites.add((rec.chrom, rec.start, rec.class_payload.get("trid") or rec.variant_id, str(rec.end)))
    return sorted(sites), {}


def write_sites_vcf(sites: Iterable[Tuple[str, int, str, str]], path: str, contigs: Iterable[str] = (), header_lines: Iterable[str] = ()) -> int:
    n = 0
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for l in header_lines:
            fh.write(l.rstrip("\n") + "\n")
        for c in contigs:
            fh.write("##contig=<ID=%s>\n" % c)
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref, alt in sites:
            if alt.startswith("<") or len(ref) > 1000 or len(alt) > 1000:
                continue
            fh.write("%s\t%d\t.\t%s\t%s\t.\t.\t.\n" % (chrom, pos, ref, alt))
            n += 1
    return n


def write_sites_bed(sites: Iterable[Tuple[str, int, str, str]], path: str, pad: int = 1) -> int:
    n = 0
    with open(path, "w") as fh:
        for chrom, pos, ref, alt in sites:
            fh.write("%s\t%d\t%d\n" % (chrom, max(0, pos - 1 - pad), pos + max(1, len(ref)) + pad))
            n += 1
    return n


# ----------------------------------------------------------------------------------------------
# gnomAD via slivar gnotate
# ----------------------------------------------------------------------------------------------
def gnotate_command(sites_vcf: str, out_vcf: str, slivar_cmd: List[str], gnomad_zip: str) -> List[str]:
    """slivar >= 0.3 removed the `gnotate` sub-command; `slivar expr --gnotate` with no expression annotates every record
    (the WDL's own tertiary step uses the same zip through slivar expr)."""
    return slivar_cmd + ["expr", "--vcf", sites_vcf, "--gnotate", gnomad_zip, "--out-vcf", out_vcf]


def contig_lines(vcf_or_bcf: str, bcftools: str) -> List[str]:
    """##contig header lines of a callset (so a sites-only VCF we write is a valid input for slivar / htslib)."""
    try:
        h = subprocess.run([bcftools, "view", "-h", vcf_or_bcf], check=True, capture_output=True, text=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return []
    return [l for l in h.splitlines() if l.startswith("##contig=")]


def gnotate(sites_vcf: str, out_vcf: str, slivar_cmd: List[str], gnomad_zip: str) -> Dict[Tuple[str, int, str, str], float]:
    """Annotate the sites VCF with gnomAD via slivar and return {normalised site: gnomad_af}."""
    subprocess.run(gnotate_command(sites_vcf, out_vcf, slivar_cmd, gnomad_zip), check=True)
    out: Dict[Tuple[str, int, str, str], float] = {}
    op = gzip.open if out_vcf.endswith(".gz") else open
    with op(out_vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            info = dict(kv.split("=", 1) if "=" in kv else (kv, "1") for kv in f[7].split(";") if kv)
            af = info.get("gnomad_af")
            if af is None:
                continue
            try:
                p, a, b = norm_allele(int(f[1]), f[3], f[4])
                out[(f[0], p, a, b)] = float(af)
            except ValueError:
                continue
    return out


# ----------------------------------------------------------------------------------------------
# founder genotypes at the sites -> leave-one-family-out counts
# ----------------------------------------------------------------------------------------------
def founder_genotypes(cohort_vcf: str, sites_bed: str, founder_ids: List[str], bcftools: str, work: str) -> Tuple[List[str], Dict[Tuple[str, int, str, str], List[str]]]:
    """{normalised site: [GT per founder]} from the cohort callset restricted to the sites (split multi-allelics)."""
    samples = os.path.join(work, "founders.txt")
    with open(samples, "w") as fh:
        fh.write("\n".join(founder_ids) + "\n")
    q = subprocess.run("%s view -R %s -S %s --force-samples -Ou %s | %s norm -m -any -Ou | %s query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT]\\n'"
                       % (bcftools, sites_bed, samples, cohort_vcf, bcftools, bcftools), shell=True, check=True, capture_output=True, text=True)
    order = subprocess.run("%s query -l -S %s --force-samples %s 2>/dev/null || cat %s" % (bcftools, samples, cohort_vcf, samples), shell=True,
                           capture_output=True, text=True).stdout.split()
    if len(order) != len(founder_ids):
        order = founder_ids
    out: Dict[Tuple[str, int, str, str], List[str]] = {}
    for line in q.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 5:
            continue
        try:
            p, a, b = norm_allele(int(f[1]), f[2], f[3])
        except ValueError:
            continue
        out[(f[0], p, a, b)] = f[4:]
    return order, out


def _alt_count(gt: str) -> int:
    return sum(1 for x in gt.replace("|", "/").split("/") if x not in ("0", ".", ""))


def _called(gt: str) -> int:
    return sum(1 for x in gt.replace("|", "/").split("/") if x != ".")


def lofo_counts(gts: List[str], order: List[str], founder_family: Dict[str, str], exclude_families: Set[str]) -> Tuple[int, int, int]:
    """(AC, AN, n_founders_with_alt) over founders whose family is not excluded."""
    ac = an = nf = 0
    for sid, gt in zip(order, gts):
        if founder_family.get(sid) in exclude_families:
            continue
        a = _alt_count(gt)
        ac += a; an += _called(gt); nf += 1 if a > 0 else 0
    return ac, an, nf


# ----------------------------------------------------------------------------------------------
# sib-shared (quads)
# ----------------------------------------------------------------------------------------------
def sib_shared_sites(joint_vcf: str, child: str, sibs: List[str], sites: Set[Tuple[str, int, str, str]]) -> Set[Tuple[str, int, str, str]]:
    """Normalised sites (of the given set) where any sibling carries the candidate allele in the family joint VCF."""
    from .io.vcf import VcfReader, iter_records
    if not sibs:
        return set()
    out: Set[Tuple[str, int, str, str]] = set()
    rd = VcfReader(joint_vcf)
    want = {(c, p) for c, p, _, _ in sites}
    for rec in iter_records(rd, tuple(sibs)):
        if (rec.chrom, rec.pos) not in want:
            continue
        for ai, alt in enumerate(rec.alts, start=1):
            key = (rec.chrom,) + norm_allele(rec.pos, rec.ref, alt)
            if key not in sites:
                continue
            for s in sibs:
                gt = rec.samples.get(s, {}).get("GT", ".")
                if str(ai) in gt.replace("|", "/").split("/"):
                    out.add(key); break
    return out


# ----------------------------------------------------------------------------------------------
# assembly per child
# ----------------------------------------------------------------------------------------------
def strchive_hit(chrom: str, start: int, end: int, index: Optional[Dict[str, list]], slop: int = 50):
    """The STRchive locus overlapping this repeat, or None.

    The join is by coordinate, not by name: our TRIDs are the lab TRGT catalogue's `chrom_start_end_motif`, while
    STRchive names loci by disease and gene (`HD_HTT`), so the two identifier spaces never meet. `slop` absorbs the
    small boundary differences between the two catalogue builds."""
    if not index:
        return None
    for a, b, rec in index.get(chrom, ()):
        if a - slop < end and start < b + slop:
            return rec
    return None


def write_annot(cand_path: str, out_path: str, class_group: str, gnomad: Dict, fgt: Dict, order: List[str], founder_family: Dict[str, str],
                exclude_families: Set[str], sib: Set, tr_table: Optional[Dict[str, List[List[int]]]] = None,
                strchive: Optional[Dict[str, list]] = None) -> Dict[str, int]:
    n = n_g = n_ac = 0
    with open(out_path, "w", newline="") as fh:
        fh.write("variant_id\tgnomad_af\tcohort_AC_loo\tcohort_AN_loo\tpon_founder_recurrence_loo\tsib_shared"
                 "\ttrgt_pop_p99_distance\tstrchive_locus\n")
        for rec in read_candidates(cand_path):
            p99 = None
            strc = ""
            if class_group == "snv_indel":
                key = (rec.chrom,) + norm_allele(rec.start, rec.ref, rec.alt)
            elif class_group == "sv":
                key = (rec.chrom, rec.start, rec.class_payload.get("caller_id") or rec.variant_id, rec.class_payload.get("svtype") or "")
            else:
                key = None
            af = gnomad.get(key) if key is not None else None
            ac = an = nf = None
            if key is not None and key in fgt:
                ac, an, nf = lofo_counts(fgt[key], order, founder_family, exclude_families)
                n_ac += 1
            if class_group == "tr" and tr_table is not None:
                pl = rec.class_payload or {}
                trid = pl.get("trid") or rec.variant_id.rsplit(":a", 1)[0]
                als = pl.get("child_AL") or []
                idx = pl.get("outlier_allele_idx")
                child_len = None
                try:
                    child_len = int(als[int(idx)]) if (als and idx is not None) else None
                except (ValueError, IndexError, TypeError):
                    child_len = None
                ac, an, nf, p99 = tr_lofo_stats(child_len, int(pl.get("motif_unit_bp") or 1), tr_table.get(trid),
                                                order, founder_family, exclude_families)
                if ac is not None:
                    n_ac += 1
                hit = strchive_hit(rec.chrom, rec.start, rec.end or rec.start, strchive)
                strc = hit["id"] if hit else ("" if strchive is None else ".")
            if af is not None:
                n_g += 1
            fh.write("%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\n"
                     % (rec.variant_id, "" if af is None else af, "" if ac is None else ac, "" if an is None else an,
                        "" if nf is None else nf, int(key in sib) if key is not None else 0,
                        "" if p99 is None else p99, strc))
            n += 1
    return {"rows": n, "gnomad_annotated": n_g, "founder_counts": n_ac}


def tr_bed(sites: Iterable[Tuple[str, int, str, str]], path: str, pad: int = 50) -> int:
    """BED of the candidate TR loci (chrom, start, trid, end), for restricting the cohort TRGT query."""
    n = 0
    with open(path, "w") as fh:
        for chrom, pos, _trid, end in sorted(sites):
            try:
                e = int(end)
            except (TypeError, ValueError):
                e = pos
            fh.write("%s\t%d\t%d\n" % (chrom, max(0, pos - pad - 1), max(e, pos) + pad))
            n += 1
    return n


def tr_founder_alleles(cohort_trgt_vcf: str, sites_bed: str, founder_ids: List[str], bcftools: str, work: str,
                       min_spanning: int = 5) -> Tuple[List[str], Dict[str, List[List[int]]]]:
    """Per locus, the QC-passing allele lengths of every unaffected founder, in founder order.

    The founder reference of `02_tiering/tr_outliers.py`, ported: an allele counts only when its spanning-read depth
    SD is at least `min_spanning`, because allele purity is not usable on this multi-motif catalogue. Returns the
    founder order and {TRID: [[allele lengths] per founder]}."""
    os.makedirs(work, exist_ok=True)
    sfile = os.path.join(work, "tr_founders.txt")
    with open(sfile, "w") as fh:
        fh.write("\n".join(founder_ids) + "\n")
    out = os.path.join(work, "tr_founder_alleles.tsv")
    if not os.path.exists(out) or os.path.getsize(out) == 0:
        cmd = ("%s view -R %s -S %s %s | %s query -f '%%INFO/TRID[\t%%AL|%%SD]\n' > %s"
               % (bcftools, sites_bed, sfile, cohort_trgt_vcf, bcftools, out))
        rc = subprocess.call(cmd, shell=True)
        if rc != 0:
            raise RuntimeError("tr founder query failed (rc=%d)" % rc)
    order = list(founder_ids)
    table: Dict[str, List[List[int]]] = {}
    with open(out) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            trid = f[0]
            per: List[List[int]] = []
            for cell in f[1:]:
                al_sd = cell.split("|")
                als = (al_sd[0] if al_sd else ".").split(",")
                sds = (al_sd[1] if len(al_sd) > 1 else ".").split(",")
                keep: List[int] = []
                for i, a in enumerate(als):
                    if a in (".", ""):
                        continue
                    try:
                        sd = float(sds[i]) if i < len(sds) and sds[i] not in (".", "") else 0.0
                    except ValueError:
                        sd = 0.0
                    if sd >= min_spanning:
                        try:
                            keep.append(int(float(a)))
                        except ValueError:
                            pass
                per.append(keep)
            table[trid] = per
    return order, table


def _pct(sorted_vals: List[int], q: float) -> Optional[float]:
    """Nearest-rank percentile, as in 02_tiering/tr_outliers.py."""
    if not sorted_vals:
        return None
    import math
    k = max(1, int(math.ceil(q / 100.0 * len(sorted_vals))))
    return float(sorted_vals[min(k, len(sorted_vals)) - 1])


def tr_lofo_stats(child_len: Optional[int], unit_bp: int, per_founder: Optional[List[List[int]]], order: List[str],
                  founder_family: Dict[str, str], exclude_families: Set[str], tol_units: float = 1.0
                  ) -> Tuple[Optional[int], Optional[int], Optional[int], Optional[float]]:
    """(cohort_AC_loo, cohort_AN_loo, pon_founder_recurrence_loo, trgt_pop_p99_distance) for one candidate allele.

    Leave-one-family-out exactly as for the other classes (P26): founders of the child's family, and of the surrogate
    parents' family for a synthetic trio, are removed from the reference before anything is counted.
      cohort_AC_loo   founder alleles whose length matches the child's outlier allele within `tol_units` motif units
      pon_founder_recurrence_loo   founders carrying any allele at least as long (a recurrent expansion at this locus)
      trgt_pop_p99_distance        (child allele - founder 99th percentile) in motif units; negative = inside the range
    """
    if child_len is None or per_founder is None:
        return None, None, None, None
    pool: List[int] = []
    n_rec = 0
    unit = max(1, int(unit_bp or 1))
    tol = max(unit * tol_units, 1.0)
    for i, sample in enumerate(order):
        if i >= len(per_founder):
            break
        if founder_family.get(sample) in exclude_families:
            continue
        als = per_founder[i]
        pool.extend(als)
        if any(a >= child_len for a in als):
            n_rec += 1
    if not pool:
        return None, 0, None, None
    ac = sum(1 for a in pool if abs(a - child_len) <= tol)
    pool.sort()
    p99 = _pct(pool, 99)
    dist = round((child_len - p99) / unit, 2) if p99 is not None else None
    return ac, len(pool), n_rec, dist


def founder_genotypes_sv(cohort_sv_vcf: str, founder_ids: List[str], bcftools: str, work: str) -> Tuple[List[str], Dict[Tuple, List[str]]]:
    """{(chrom, pos, sawfish id, svtype): [GT per founder]} for every record of the cohort SV VCF (the SV candidate key)."""
    samples = os.path.join(work, "founders.txt")
    with open(samples, "w") as fh:
        fh.write("\n".join(founder_ids) + "\n")
    q = subprocess.run("%s view -S %s --force-samples -Ou %s | %s query -f '%%CHROM\\t%%POS\\t%%ID\\t%%INFO/SVTYPE[\\t%%GT]\\n'" % (bcftools, samples, cohort_sv_vcf, bcftools),
                       shell=True, check=True, capture_output=True, text=True)
    out: Dict[Tuple, List[str]] = {}
    for line in q.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 5:
            continue
        out[(f[0], int(f[1]), f[2], f[3])] = f[4:]
    return founder_ids, out
