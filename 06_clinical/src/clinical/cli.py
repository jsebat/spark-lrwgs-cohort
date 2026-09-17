"""clinical: the targeted clinical arm of the de novo workflow (DESIGN P31-P33 in 05_denovo).

  clinical annotate-smallvar   consequences for SNV/indel candidates inside panel genes (VEP+LOFTEE VCFs, dbNSFP)
  clinical run                 build arm A and arm B tables and the per-child report
  clinical acceptance          the module must recover the cohort's two known pathogenic deletions (MECP2, DNMT3A)

Identifiers are arguments and data; none is written into code.
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import sys

from . import arms as A
from . import smallvar_annot as S
from .genes import GeneModel, Panel, load_shet

TAB = chr(9)


def _common(ap: argparse.ArgumentParser):
    ap.add_argument("--final-dir", required=True, help="05_denovo FINAL_DIR (*.{snv_indel,sv,tr}.dnm.tsv)")
    ap.add_argument("--evidence-dir", required=True, help="05_denovo EVIDENCE_DIR (raw six-haplotype SV evidence)")
    ap.add_argument("--gff3", required=True, help="Ensembl GFF3 gene model (exon/CDS/UTR)")
    ap.add_argument("--gene-sets", required=True, help="directory with sfari_hc.txt sfari_all.txt ddg2p_confident.txt ddg2p_all.txt")
    ap.add_argument("--genebayes", default="", help="GeneBayes Supplementary_Table_1.tsv (s_het by Ensembl id); reported, never gates")
    ap.add_argument("--out-dir", required=True)


def cmd_annotate(a) -> int:
    gm = GeneModel.from_gff3(a.gff3)
    panel = Panel.load(a.gene_sets)
    rows = A.load_final(a.final_dir, "snv_indel")
    keys = S.candidate_keys(rows, gm, panel.union)
    sys.stderr.write("SNV/indel candidates: %d; inside a panel gene span: %d\n" % (len(rows), len(keys)))
    annot, missing = S.annotate_from_vep(keys, a.vep_glob)
    n_miss = S.attach_missense_tiers(annot, a.dbnsfp_parquet) if a.dbnsfp_parquet else 0
    os.makedirs(a.out_dir, exist_ok=True)
    out = os.path.join(a.out_dir, "smallvar_annotation.tsv")
    S.write_annot(annot, missing, out)
    sys.stderr.write("annotated %d, not found in the VEP VCFs %d, missense tiers attached %d -> %s\n" % (len(annot), len(missing), n_miss, out))
    return 0


def cmd_run(a) -> int:
    gm = GeneModel.from_gff3(a.gff3)
    panel = Panel.load(a.gene_sets)
    shet = load_shet(a.genebayes, a.gff3) if a.genebayes else {}
    sys.stderr.write("gene model: %d symbols; panel union: %d genes; s_het: %d genes\n" % (len(gm.span), len(panel.union), len(shet)))
    os.makedirs(a.out_dir, exist_ok=True)

    sv = A.load_final(a.final_dir, "snv_indel" if False else "sv")
    ev = A.load_sv_evidence(a.evidence_dir)
    sys.stderr.write("SV candidates: %d; evidence rows: %d\n" % (len(sv), len(ev)))
    arm_a = A.arm_a(sv, ev, gm, panel, shet, min_svlen=a.min_svlen)
    A.write_rows(arm_a, os.path.join(a.out_dir, "arm_A.rare_large_sv.tsv"))

    b = A.arm_b_sv(sv, ev, gm, panel, shet)
    annot = S.read_annot(os.path.join(a.out_dir, "smallvar_annotation.tsv"))
    if annot:
        b += A.arm_b_smallvar(A.load_final(a.final_dir, "snv_indel"), annot, panel, shet, gm)
    else:
        sys.stderr.write("NOTE: no smallvar_annotation.tsv in out-dir; run `clinical annotate-smallvar` first for the SNV/indel arm\n")
    tr = A.load_final(a.final_dir, "tr")
    if tr:
        b += A.arm_b_tr(tr, gm, panel, shet)
    b.sort(key=lambda r: (r["sample_id"], {"PASS": 0, "REVIEW": 1, "FAIL": 2}[r["quality_verdict"]], r["panel_rank"], r["chrom"], int(r["start"])))
    A.write_rows(b, os.path.join(a.out_dir, "arm_B.clinically_led.tsv"))

    # per-child report: quality PASS/REVIEW rows of both arms, clinically led rows first
    rep = os.path.join(a.out_dir, "clinical_report.md")
    with open(rep, "w") as fh:
        fh.write("# Targeted clinical arm — per-child report\n\n")
        fh.write("discovery_mode = targeted_clinical for every row below. Nothing here enters a recall, FDR or rate "
                 "estimate (DESIGN P31). Quality and clinical relevance are recorded independently (P32).\n\n")
        fh.write("| arm | rows | PASS | REVIEW | FAIL |\n|---|---|---|---|---|\n")
        for name, rows in (("A rare+large SV", arm_a), ("B clinically led", b)):
            c = collections.Counter(r["quality_verdict"] for r in rows)
            fh.write("| %s | %d | %d | %d | %d |\n" % (name, len(rows), c["PASS"], c["REVIEW"], c["FAIL"]))
        fh.write("\n")
        by_child = collections.defaultdict(list)
        for r in b + arm_a:
            if r["quality_verdict"] in ("PASS", "REVIEW") and not str(r["clinical_relevance"]).startswith("NONE"):
                by_child[(r["family_id"], r["sample_id"])].append(r)
        for (fam, sid), rows in sorted(by_child.items()):
            fh.write("## %s / %s\n\n" % (fam, sid))
            fh.write("| arm | variant | gene | impact | gene sets | s_het | quality | clinical relevance | caller |\n|---|---|---|---|---|---|---|---|---|\n")
            for r in sorted(rows, key=lambda r: (str(r["clinical_relevance"])[:4], r["panel_rank"])):
                var = "%s:%s-%s %s %s" % (r["chrom"], r["start"], r["end"], r.get("svtype") or r["variant_class"], r.get("svlen") or "")
                fh.write("| %s | %s | %s | %s | %s | %s | %s (%s) | %s | %s/%s |\n" % (
                    r["arm"][:1], var, r["gene"], r["impact"], r["gene_sets"], r["s_het"] or "NA", r["quality_verdict"],
                    r["quality_criteria"], r["clinical_relevance"], r["dnm_call"], r["decision_reason"]))
            fh.write("\n")
    summary = dict(arm_A_rows=len(arm_a), arm_B_rows=len(b),
                   arm_A_verdicts=dict(collections.Counter(r["quality_verdict"] for r in arm_a)),
                   arm_B_verdicts=dict(collections.Counter(r["quality_verdict"] for r in b)),
                   children_with_reportable=len(by_child))
    json.dump(summary, open(os.path.join(a.out_dir, "summary.json"), "w"), indent=1)
    print(json.dumps(summary, indent=1))
    return 0


def cmd_acceptance(a) -> int:
    """The module must find the cohort's two known pathogenic de novo deletions in arm B, and grade them PASS or
    REVIEW. These are an ACCEPTANCE test, not a validation (P31): the depth-path constants were chosen with them in
    view. Expectations are given as data (a TSV of gene, sample, chrom, start, minimum verdict), never in code."""
    want = list(csv.DictReader(open(a.expect, newline=""), delimiter=TAB))
    rows = list(csv.DictReader(open(os.path.join(a.out_dir, "arm_B.clinically_led.tsv"), newline=""), delimiter=TAB))
    order = {"FAIL": 0, "REVIEW": 1, "PASS": 2}
    ok = True
    for w in want:
        hit = [r for r in rows if r["gene"] == w["gene"] and r["sample_id"] == w["sample_id"] and r["chrom"] == w["chrom"] and r["start"] == w["start"]]
        if not hit:
            print("FAIL  %-8s %s %s:%s  not in arm B" % (w["gene"], w["sample_id"], w["chrom"], w["start"]))
            ok = False
            continue
        r = hit[0]
        good = order[r["quality_verdict"]] >= order[w.get("min_verdict", "REVIEW")]
        ok &= good
        print("%s  %-8s %s %s:%s-%s %s %s  quality=%s (%s)  clinical=%s  caller=%s/%s" % (
            "ok  " if good else "FAIL", w["gene"], w["sample_id"], r["chrom"], r["start"], r["end"], r["svtype"], r["svlen"],
            r["quality_verdict"], r["quality_criteria"], r["clinical_relevance"], r["dnm_call"], r["decision_reason"]))
        if a.verbose:
            print("      reasons: %s" % r["quality_reasons"])
    print("ACCEPTANCE:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="clinical", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("annotate-smallvar"); _common(p)
    p.add_argument("--vep-glob", required=True, help="per-chromosome VEP+LOFTEE VCFs, e.g. '/path/vep/chr*.vep.vcf.gz'")
    p.add_argument("--dbnsfp-parquet", default="", help="dbNSFP parquet directory (<chrom>.parquet) for missense tiers")
    p.set_defaults(fn=cmd_annotate)
    p = sub.add_parser("run"); _common(p)
    p.add_argument("--min-svlen", type=int, default=5000, help="arm A size floor (P32: 5 kb)")
    p.set_defaults(fn=cmd_run)
    p = sub.add_parser("acceptance")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--expect", required=True, help="TSV: gene sample_id chrom start min_verdict")
    p.add_argument("--verbose", action="store_true")
    p.set_defaults(fn=cmd_acceptance)
    a = ap.parse_args(argv)
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
