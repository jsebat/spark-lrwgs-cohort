"""P26 annotation: founder set, leave-one-family-out counts, site keys, and the per-child table (bcftools/slivar stubbed)."""
import csv

from phase_dnm import annotate as AN
from phase_dnm.records import CandidateRecord, write_candidates


def test_founders_are_unaffected_parents_only():
    rows = [dict(sample_id="f1", family_id="A", father_id="0", mother_id="0", affected="1"),
            dict(sample_id="m1", family_id="A", father_id="0", mother_id="0", affected="2"),        # affected founder: excluded (as in the PON)
            dict(sample_id="k1", family_id="A", father_id="f1", mother_id="m1", affected="2"),
            dict(sample_id="f2", family_id="B", father_id="0", mother_id="0", affected="1")]
    assert AN.founders(rows) == {"f1": "A", "f2": "B"}


def test_lofo_counts_exclude_both_families():
    order = ["f1", "m1", "f2", "m2", "f3"]
    fam = {"f1": "A", "m1": "A", "f2": "B", "m2": "B", "f3": "C"}
    gts = ["0/1", "0/0", "1/1", "./.", "0/1"]
    assert AN.lofo_counts(gts, order, fam, set()) == (4, 8, 3)
    assert AN.lofo_counts(gts, order, fam, {"A"}) == (3, 4, 2)                # a real trio: its own founders out (the uncalled one adds 0 to AN)
    assert AN.lofo_counts(gts, order, fam, {"A", "B"}) == (1, 2, 1)           # a synthetic trio: both families out
    assert AN.lofo_counts(["0|1", "1|1"], ["f1", "m1"], fam, set()) == (3, 4, 2)


def test_union_sites_and_bed(tmp_path):
    recs = [CandidateRecord("fam", "kid", "chr1:100:A:G", "chr1", 100, 100, "A", "G", "SNV", "dv"),
            CandidateRecord("fam", "kid", "chr1:200:GAT:G", "chr1", 200, 202, "GAT", "G", "INDEL", "dv"),
            CandidateRecord("fam", "kid2", "chr1:100:A:G", "chr1", 100, 100, "A", "G", "SNV", "dv")]
    p1 = str(tmp_path / "a.tsv"); p2 = str(tmp_path / "b.tsv")
    write_candidates(recs[:2], p1); write_candidates(recs[2:], p2)
    sites, _ = AN.union_sites([p1, p2], "snv_indel")
    assert sites == [("chr1", 100, "A", "G"), ("chr1", 200, "GAT", "G")]
    bed = str(tmp_path / "s.bed")
    assert AN.write_sites_bed(sites, bed) == 2
    lines = open(bed).read().splitlines()
    assert lines[1].split("\t") == ["chr1", "198", "204"]
    vcf = str(tmp_path / "s.vcf")
    assert AN.write_sites_vcf(sites, vcf) == 2 and open(vcf).read().splitlines()[-1].startswith("chr1\t200\t.\tGAT\tG")


def test_write_annot_joins_by_normalised_key(tmp_path):
    recs = [CandidateRecord("fam", "kid", "v1", "chr1", 100, 100, "A", "G", "SNV", "dv"),
            CandidateRecord("fam", "kid", "v2", "chr1", 200, 202, "GAT", "G", "INDEL", "dv"),
            CandidateRecord("fam", "kid", "v3", "chr1", 300, 300, "C", "T", "SNV", "dv")]
    cp = str(tmp_path / "c.tsv"); write_candidates(recs, cp)
    gnomad = {("chr1", 100, "A", "G"): 0.0004, ("chr1", 201, "AT", ""): 0.2}       # the indel in canonical form
    order = ["f1", "m1", "f2"]; fam = {"f1": "A", "m1": "A", "f2": "B"}
    fgt = {("chr1", 100, "A", "G"): ["0/1", "0/0", "0/1"], ("chr1", 201, "AT", ""): ["1/1", "0/1", "0/0"]}
    out = str(tmp_path / "annot.tsv")
    st = AN.write_annot(cp, out, "snv_indel", gnomad, fgt, order, fam, {"A"}, sib={("chr1", 300, "C", "T")})
    rows = {r["variant_id"]: r for r in csv.DictReader(open(out), delimiter="\t")}
    assert st == {"rows": 3, "gnomad_annotated": 2, "founder_counts": 2}
    assert rows["v1"]["gnomad_af"] == "0.0004" and rows["v1"]["cohort_AC_loo"] == "1" and rows["v1"]["cohort_AN_loo"] == "2" and rows["v1"]["pon_founder_recurrence_loo"] == "1"
    assert rows["v2"]["gnomad_af"] == "0.2" and rows["v2"]["cohort_AC_loo"] == "0"             # family A excluded -> only f2 (0/0)
    assert rows["v3"]["gnomad_af"] == "" and rows["v3"]["cohort_AC_loo"] == "" and rows["v3"]["sib_shared"] == "1"
    assert rows["v1"]["sib_shared"] == "0"


def test_gnotate_uses_slivar_expr_and_contig_headers(tmp_path):
    cmd = AN.gnotate_command("s.vcf", "o.vcf", ["singularity", "exec", "-B", "/x:/x", "slivar.sif", "slivar"], "gnomad.zip")
    assert cmd[:6] == ["singularity", "exec", "-B", "/x:/x", "slivar.sif", "slivar"] and cmd[6] == "expr"
    assert "--gnotate" in cmd and "--out-vcf" in cmd and cmd[cmd.index("--vcf") + 1] == "s.vcf"
    vcf = str(tmp_path / "s.vcf")
    AN.write_sites_vcf([("chr1", 5, "A", "G")], vcf, header_lines=["##contig=<ID=chr1,length=248956422>"])
    lines = open(vcf).read().splitlines()
    assert lines[0] == "##fileformat=VCFv4.2" and lines[1].startswith("##contig=<ID=chr1") and lines[2].startswith("#CHROM")
