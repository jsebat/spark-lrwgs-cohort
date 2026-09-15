"""Genotype the spiked slices so planted candidates carry the same caller-derived block as real ones (P29).

A planted variant never went through a variant caller, so `spike plan` writes its candidate row with the whole genotype
block empty -- caller_gt / caller_gq / caller_dp / caller_qual, child_ad, child_pl and both parents' equivalents. That
block feeds eighteen registry features (child_AR/AB/DP/GQ/PL0-2, the min/max parent forms, qd_child), and the population
annotation feeds four more, so until 2026-09-15 the planted-truth arm scored the classifier on a matrix in which 42 of
the 72 SNV/indel columns were empty. The presence-leak guard silently drops an all-empty column at training time, but at
scoring time the frozen model still expects it and reads a constant: the external arm has therefore been measuring the
classifier with over half its evidence removed, on exactly the arm used to validate it. This module closes that gap.

Small variants are genotyped for real: bcftools mpileup + call over the three spiked slice BAMs at the planted sites.
That is a caller run on the edited reads, so DP, AD, PL and GQ mean the same thing they mean in the real arm.

SVs and TRs cannot be genotyped that way -- bcftools has no model for a 35 kb deletion or a repeat expansion, and the
real arm gets those fields from sawfish and TRGT, which cannot be run on a slice. For those two classes the allelic
depths come from the planter's own ledger: `apply_plan` records, per site and per role, how many reads it edited (the
alt reads), how many it left alone (the ref reads) and how many it removed. Those counts are exact, where a caller's
would carry error, so this is an OPTIMISTIC bound on what the classifier could see. That direction is deliberate: the
classifier is the arm on trial in P28, and an optimistic bound cannot manufacture the conclusion that it does worse.
Genotype likelihoods are then the standard binomial model over those counts (`_pl_from_counts`), not a caller's own.

Population annotation is by construction, not by lookup: a planted variant is private to the child by definition, so
gnomad_af is absent, cohort_AC_loo and pon_founder_recurrence_loo are 0 and sib_shared is 0. Writing those as empty
would have told the classifier "unknown" where the truth is "private", which is the single most informative state a
de novo candidate can be in.
"""
from __future__ import annotations

import csv
import math
import os
import subprocess
import tempfile
from typing import Dict, List, Optional, Tuple

ROLE_COLS = {"child": ("caller_gt", "caller_gq", "caller_dp", "child_ad", "child_pl"),
             "father": ("father_gt", "father_gq", "father_dp", "father_ad", "father_pl"),
             "mother": ("mother_gt", "mother_gq", "mother_dp", "mother_ad", "mother_pl")}
ERR = 0.02          # per-read error rate of the binomial genotype model (HiFi; the same order sawfish assumes)


# ----------------------------------------------------------------------------------------------
# a genotype likelihood from allelic counts
# ----------------------------------------------------------------------------------------------
def _pl_from_counts(n_ref: int, n_alt: int, err: float = ERR) -> Tuple[str, str, str]:
    """(GT, GQ, PL) under the standard biallelic binomial model: P(read is alt | 0/0, 0/1, 1/1) = err, 0.5, 1-err.

    This is the genotype model, not a caller: it says what the allelic counts imply, with no mapping or realignment
    term. It is used only for the planted SV and TR rows, whose counts come from the planter's ledger."""
    n = n_ref + n_alt
    if n <= 0:
        return ".", "", "."
    p = (err, 0.5, 1.0 - err)
    ll = [n_alt * math.log10(max(x, 1e-12)) + n_ref * math.log10(max(1.0 - x, 1e-12)) for x in p]
    best = max(range(3), key=lambda i: ll[i])
    pl = [min(999, int(round(-10.0 * (ll[i] - ll[best])))) for i in range(3)]
    second = min(pl[i] for i in range(3) if i != best)
    gt = ("0/0", "0/1", "1/1")[best]
    return gt, str(min(99, second)), ",".join(str(x) for x in pl)


# ----------------------------------------------------------------------------------------------
# the planter's ledger
# ----------------------------------------------------------------------------------------------
def read_ledger(path: str) -> Dict[Tuple[str, str], Dict[str, int]]:
    """{(variant_id, role): {alt, ref, dropped}} from apply.ledger.tsv."""
    out: Dict[Tuple[str, str], Dict[str, int]] = {}
    if not os.path.exists(path):
        return out
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[(r["variant_id"], r["role"])] = {k: int(r.get(k) or 0) for k in ("alt", "ref", "dropped")}
    return out


# ----------------------------------------------------------------------------------------------
# bcftools on the spiked slices
# ----------------------------------------------------------------------------------------------
def _bcftools_call(sites: List[Tuple[str, int]], bams: Dict[str, str], reference: str, bcftools: str,
                   work: str, min_mapq: int = 5, min_bq: int = 1, max_depth: int = 500,
                   log=None) -> Dict[Tuple[str, int], Dict[str, Dict[str, str]]]:
    """{(chrom, pos): {role: {GT, GQ, DP, AD, PL, QUAL}}} by genotyping the three slice BAMs at the planted sites.

    Roles are taken from the ORDER the BAMs are passed, never from the read groups: the slices carry the real sample
    names, which must not leave the cluster in any table this module writes."""
    if not sites:
        return {}
    order = [r for r in ("child", "father", "mother") if r in bams]
    reg = os.path.join(work, "sites.tsv")
    with open(reg, "w", newline="\n") as fh:
        for chrom, pos in sorted(set(sites)):
            fh.write("%s\t%d\n" % (chrom, pos))
    vcf = os.path.join(work, "spiked.vcf")
    mp = [bcftools, "mpileup", "-f", reference, "-R", reg, "-a", "AD,DP", "-q", str(min_mapq), "-Q", str(min_bq),
          "-d", str(max_depth), "-Ou"] + [bams[r] for r in order]
    cl = [bcftools, "call", "-m", "-A", "-f", "GQ", "-Ov", "-o", vcf]   # -f GQ: bcftools omits GQ unless asked
    try:
        p1 = subprocess.Popen(mp, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p2 = subprocess.Popen(cl, stdin=p1.stdout, stderr=subprocess.PIPE)
        p1.stdout.close()
        _, e2 = p2.communicate()
        p1.wait()
        if p2.returncode != 0:
            raise RuntimeError((e2 or b"").decode()[:400])
    except (OSError, RuntimeError) as e:
        if log:
            log("spike genotype: bcftools failed (%s); the small-variant genotype block stays empty" % e)
        return {}
    out: Dict[Tuple[str, int], Dict[str, Dict[str, str]]] = {}
    with open(vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 + len(order):
                continue
            keys = f[8].split(":")
            per = {}
            for i, role in enumerate(order):
                vals = f[9 + i].split(":")
                d = dict(zip(keys, vals))
                per[role] = {"GT": d.get("GT", "."), "GQ": d.get("GQ", ""), "DP": d.get("DP", ""),
                             "AD": d.get("AD", "."), "PL": d.get("PL", "."), "QUAL": f[5]}
            out[(f[0], int(f[1]))] = per
    if log:
        log("spike genotype: bcftools genotyped %d of %d planted small-variant sites" % (len(out), len(set(sites))))
    return out


# ----------------------------------------------------------------------------------------------
# rewrite the candidates table
# ----------------------------------------------------------------------------------------------
def fill_candidates(candidates_path: str, out_path: str, bams: Dict[str, str], reference: Optional[str],
                    ledger_path: str, bcftools: str = "bcftools", log=None) -> Dict[str, int]:
    """Rewrite a spike candidates table with the genotype block filled. Returns per-class counts of what was filled."""
    with open(candidates_path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = list(rd.fieldnames or [])
        rows = list(rd)
    if not rows:
        return {"rows": 0}
    ledger = read_ledger(ledger_path)
    small = [(r["chrom"], int(r["start"])) for r in rows if r.get("variant_class") in ("SNV", "INDEL")]
    gt: Dict[Tuple[str, int], Dict[str, Dict[str, str]]] = {}
    if small and reference and os.path.exists(reference):
        work = tempfile.mkdtemp(prefix="spikegt.", dir=os.path.dirname(os.path.abspath(out_path)))
        try:
            gt = _bcftools_call(small, bams, reference, bcftools, work, log=log)
        finally:
            for fn in os.listdir(work):
                os.remove(os.path.join(work, fn))
            os.rmdir(work)
    elif small and log:
        log("spike genotype: no reference FASTA, so %d planted small variants keep an empty genotype block" % len(small))

    n = {"small_genotyped": 0, "sv_from_ledger": 0, "tr_from_ledger": 0, "rows": len(rows)}
    for r in rows:
        vc = r.get("variant_class")
        if vc in ("SNV", "INDEL"):
            per = gt.get((r["chrom"], int(r["start"])))
            if not per:
                continue
            for role, (c_gt, c_gq, c_dp, c_ad, c_pl) in ROLE_COLS.items():
                d = per.get(role)
                if not d:
                    continue
                r[c_gt], r[c_gq], r[c_dp], r[c_ad], r[c_pl] = d["GT"], d["GQ"], d["DP"], d["AD"], d["PL"]
            r["caller_qual"] = per.get("child", {}).get("QUAL", "")
            r["caller"] = "spike+bcftools"
            n["small_genotyped"] += 1
        else:
            filled = False
            for role, (c_gt, c_gq, c_dp, c_ad, c_pl) in ROLE_COLS.items():
                led = ledger.get((r["variant_id"], role))
                if not led:
                    continue
                n_alt, n_ref = led["alt"], led["ref"]
                g, q, p = _pl_from_counts(n_ref, n_alt)
                r[c_gt], r[c_gq], r[c_dp], r[c_ad], r[c_pl] = g, q, str(n_ref + n_alt), "%d,%d" % (n_ref, n_alt), p
                filled = True
            if filled:
                r["caller"] = "spike+ledger"
                n["sv_from_ledger" if vc == "SV" else "tr_from_ledger"] += 1
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    if log:
        log("spike genotype: %(rows)d rows -> %(small_genotyped)d small variants genotyped by bcftools, "
            "%(sv_from_ledger)d SV and %(tr_from_ledger)d TR filled from the planter ledger" % n)
    return n


def write_private_annot(candidates_path: str, out_path: str, log=None) -> int:
    """The annotation table a planted variant would have if it were looked up: it is private to the child, by
    construction. Absent gnomAD AF, zero leave-one-family-out cohort AC and founder recurrence, not sibling-shared.
    Written in the same schema `annotate write_annot` uses, so the features step joins it without a special case."""
    with open(candidates_path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    with open(out_path, "w", newline="\n") as fh:
        fh.write("variant_id\tgnomad_af\tcohort_AC_loo\tcohort_AN_loo\tpon_founder_recurrence_loo\tsib_shared"
                 "\ttrgt_pop_p99_distance\tstrchive_locus\n")
        for r in rows:
            # gnomad_af -1: the negative sentinel the rarity gate already reads as "absent from gnomAD"
            # (nested_cv._private_ids keeps a row when af < af_max OR af < 0). Leaving it empty would leave the
            # column unproduced, which the presence gate correctly refuses.
            # cohort_AN_loo: the planted site was not queried, so the denominator is left empty rather than invented.
            fh.write("%s\t-1\t0\t\t0\t0\t\t\n" % r["variant_id"])
    if log:
        log("spike genotype: wrote the by-construction private annotation for %d planted candidates" % len(rows))
    return len(rows)
