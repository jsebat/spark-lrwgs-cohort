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
def write_annot(cand_path: str, out_path: str, class_group: str, gnomad: Dict, fgt: Dict, order: List[str], founder_family: Dict[str, str],
                exclude_families: Set[str], sib: Set) -> Dict[str, int]:
    n = n_g = n_ac = 0
    with open(out_path, "w", newline="") as fh:
        fh.write("variant_id\tgnomad_af\tcohort_AC_loo\tcohort_AN_loo\tpon_founder_recurrence_loo\tsib_shared\n")
        for rec in read_candidates(cand_path):
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
            if af is not None:
                n_g += 1
            fh.write("%s\t%s\t%s\t%s\t%s\t%d\n" % (rec.variant_id, "" if af is None else af, "" if ac is None else ac, "" if an is None else an,
                                                  "" if nf is None else nf, int(key in sib) if key is not None else 0))
            n += 1
    return {"rows": n, "gnomad_annotated": n_g, "founder_counts": n_ac}


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
