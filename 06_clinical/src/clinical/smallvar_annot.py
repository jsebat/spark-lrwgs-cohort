"""Consequence annotation for SNV/indel candidates, from a VEP+LOFTEE-annotated cohort VCF.

The de novo module scores ~34,000 SNV/indel candidates per child, most of them nowhere near a gene of interest. Arm B
needs consequences only for candidates inside a panel gene, so the candidates are first restricted to panel-gene spans
(GeneModel.genes_spanning) and only those are looked up, by tabix fetch, in the annotated VCF. The annotated VCF is the
03_tiering output (VEP 115 + LOFTEE + dbNSFP gnomAD 4.1 AF, one file per chromosome, `--pick_allele` so one CSQ per
allele); any VCF with a CSQ INFO field of that shape works.

Missense tiers (James's forward-selection thresholds; 06_clinical/legacy/miss_tier.py) need dbNSFP rankscores, which
live in a parquet store. They are attached when pyarrow can read it and left blank otherwise -- the impact is still
reported as missense, only the tier is missing, and the report says so.
"""
from __future__ import annotations

import csv
import glob
import os
import sys
from typing import Dict, Iterable, List, Optional, Set, Tuple

from .genes import GeneModel

TAB = chr(9)
Key = Tuple[str, int, str, str]

MISSENSE_THRESHOLDS = (("ClinPred_rankscore", 0.4298), ("AlphaMissense_rankscore", 0.9603),
                       ("popEVE_converted_rankscore", 0.9209), ("MPC_rankscore", 0.8947))
MISSENSE_TIER = {4: "miss_t1", 3: "miss_t2", 2: "miss_t3", 1: "miss_t4"}


def candidate_keys(final_rows: Iterable[Dict[str, str]], gm: GeneModel, panel_genes: Set[str]) -> Dict[Key, str]:
    """(chrom, pos, ref, alt) -> first panel gene whose span contains the site; only candidates inside a panel gene."""
    out: Dict[Key, str] = {}
    for r in final_rows:
        try:
            pos = int(float(r["start"]))
        except (KeyError, ValueError):
            continue
        genes = gm.genes_spanning(r["chrom"], pos) & panel_genes
        if genes:
            out[(r["chrom"], pos, r["ref"], r["alt"])] = sorted(genes)[0]
    return out


def _csq_fields(header) -> List[str]:
    d = header.info["CSQ"].description
    return d.split("Format: ")[-1].strip().split("|")


def annotate_from_vep(keys: Dict[Key, str], vep_glob: str) -> Tuple[Dict[Key, Dict[str, str]], List[Key]]:
    """Look each key up in the per-chromosome VEP VCFs. Returns (annotations, keys not found)."""
    import pysam
    by_chrom: Dict[str, List[Key]] = {}
    for k in keys:
        by_chrom.setdefault(k[0], []).append(k)
    files = {os.path.basename(p).split(".")[0]: p for p in glob.glob(vep_glob)}
    out: Dict[Key, Dict[str, str]] = {}
    missing: List[Key] = []
    for chrom, ks in sorted(by_chrom.items()):
        path = files.get(chrom)
        if not path:
            missing.extend(ks)
            continue
        vf = pysam.VariantFile(path)
        fields = _csq_fields(vf.header)
        for k in ks:
            _, pos, ref, alt = k
            hit = None
            try:
                it = vf.fetch(chrom, pos - 1, pos)
            except ValueError:
                it = []
            for v in it:
                if v.pos != pos or v.ref != ref or alt not in (v.alts or ()):
                    continue
                for c in v.info.get("CSQ", ()):
                    d = dict(zip(fields, c.split("|")))
                    # --pick_allele: one CSQ per allele; match by uploaded allele when several are present
                    ua = d.get("UPLOADED_ALLELE") or d.get("Allele")
                    if ua and ua not in (alt, alt[1:] if len(alt) > 1 else alt, "-"):
                        # VEP trims the shared anchor base for indels; accept the CSQ when the allele number matches
                        if d.get("ALLELE_NUM") and int(d["ALLELE_NUM"]) != (v.alts.index(alt) + 1):
                            continue
                    hit = d
                    break
                if hit:
                    break
            if hit:
                hit["gnomad_af"] = hit.get("gnomAD4.1_joint_AF", "")
                out[k] = hit
            else:
                missing.append(k)
        vf.close()
    return out, missing


def attach_missense_tiers(annot: Dict[Key, Dict[str, str]], parquet_dir: str) -> int:
    """n_flag / tier from dbNSFP rankscores; parquet '#chr' has no 'chr' prefix (legacy/miss_tier.py join fix)."""
    miss = [k for k, a in annot.items() if "missense" in (a.get("Consequence") or "")]
    if not miss or not parquet_dir or not os.path.isdir(parquet_dir):
        return 0
    try:
        import pyarrow.parquet as pq
        import pyarrow.compute as pc
    except ImportError:
        sys.stderr.write("missense tiers skipped: pyarrow not available\n")
        return 0
    n = 0
    by_chrom: Dict[str, List[Key]] = {}
    for k in miss:
        by_chrom.setdefault(k[0], []).append(k)
    cols = ["#chr", "pos(1-based)", "ref", "alt"] + [c for c, _ in MISSENSE_THRESHOLDS]
    for chrom, ks in by_chrom.items():
        pf = os.path.join(parquet_dir, "%s.parquet" % chrom)
        if not os.path.exists(pf):
            continue
        # the parquet stores pos(1-based) and #chr as STRINGS ('25244326', '2'); an integer filter raises
        # ArrowTypeError, and an exception here must not cost the VEP annotation already in hand
        want = {str(k[1]) for k in ks}
        try:
            t = pq.read_table(pf, columns=cols, filters=[("pos(1-based)", "in", sorted(want))])
            rows = t.to_pylist()
        except Exception as e:  # noqa: BLE001 - reported, annotation continues without tiers for this chromosome
            sys.stderr.write("missense tiers skipped for %s: %s\n" % (chrom, str(e)[:160]))
            continue
        idx = {}
        for r in rows:
            try:
                idx[(int(r["pos(1-based)"]), r["ref"], r["alt"])] = r
            except (TypeError, ValueError):
                pass
        for k in ks:
            r = idx.get((k[1], k[2], k[3]))
            if not r:
                continue
            nf = 0
            for c, thr in MISSENSE_THRESHOLDS:
                try:
                    if r[c] is not None and float(r[c]) >= thr:
                        nf += 1
                except (TypeError, ValueError):
                    pass
            annot[k]["missense_nflag"] = str(nf)
            annot[k]["missense_tier"] = MISSENSE_TIER.get(nf, "miss_none")
            n += 1
    return n


def write_annot(annot: Dict[Key, Dict[str, str]], missing: List[Key], path: str):
    cols = ["chrom", "pos", "ref", "alt", "SYMBOL", "Gene", "Consequence", "LoF", "LoF_filter", "gnomad_af", "missense_nflag", "missense_tier"]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter=TAB, lineterminator="\n")
        w.writerow(cols)
        for k, a in sorted(annot.items()):
            w.writerow([k[0], k[1], k[2], k[3]] + [a.get(c, "") for c in cols[4:]])
    with open(path + ".missing", "w") as fh:
        for k in missing:
            fh.write(TAB.join(str(x) for x in k) + "\n")


def read_annot(path: str) -> Dict[Key, Dict[str, str]]:
    out: Dict[Key, Dict[str, str]] = {}
    if not os.path.exists(path):
        return out
    for r in csv.DictReader(open(path, newline=""), delimiter=TAB):
        out[(r["chrom"], int(r["pos"]), r["ref"], r["alt"])] = r
    return out
