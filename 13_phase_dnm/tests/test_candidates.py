import gzip
import os

from phase_dnm import candidates as C
from phase_dnm.records import CandidateRecord, read_candidates, write_candidates

HDR = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=1000000>",
       '##FORMAT=<ID=GT,Number=1,Type=String,Description="">', '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="">',
       '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="">', '##FORMAT=<ID=AL,Number=.,Type=Integer,Description="">',
       '##FORMAT=<ID=SD,Number=.,Type=Integer,Description="">', '##FORMAT=<ID=ALLR,Number=.,Type=String,Description="">',
       '##FORMAT=<ID=MC,Number=.,Type=String,Description="">', '##FORMAT=<ID=CN,Number=1,Type=Float,Description="">',
       '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="">', '##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="">',
       '##INFO=<ID=END,Number=1,Type=Integer,Description="">', '##INFO=<ID=HOMLEN,Number=A,Type=Integer,Description="">',
       '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="">', '##INFO=<ID=TRID,Number=1,Type=String,Description="">',
       '##INFO=<ID=MOTIFS,Number=.,Type=String,Description="">', '##INFO=<ID=STRUC,Number=1,Type=String,Description="">',
       "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tkid\tdad\tmom"]


def _write(path, rows):
    with gzip.open(path, "wt") as fh:
        fh.write("\n".join(HDR + rows) + "\n")


# ---------------------------------------------------------------------------- SNV / INDEL on the minitrio joint VCF
def test_snv_indel_candidates_recover_planted_dnms(minitrio):
    c, f, m = minitrio["ids"]
    recs = list(C.snv_indel_candidates(minitrio["paths"]["joint"], "minifam", c, f, m))
    assert recs, "no candidates"
    got = {(r.chrom, r.start) for r in recs}
    planted = {(ch, pos) for ch, pos, ref, alt, hap in minitrio["truth"].dnms}
    # every planted DNM is a candidate unless a simulated genotype error hit it (rare); allow at most one miss
    assert len(planted - got) <= 1, "planted DNMs missing from the unfiltered candidate set: %s" % (planted - got)
    # the set is UNFILTERED: genotype errors make many more candidates than the planted 6
    assert len(recs) >= len(planted)
    assert {r.variant_class for r in recs} <= {"SNV", "INDEL"}
    assert all(r.source_tier == "UNFILTERED" and r.father_gt != "." and r.mother_gt != "." for r in recs)


# ---------------------------------------------------------------------------- SV (sawfish-shaped)
def test_sv_candidates(tmp_path):
    p = tmp_path / "sv.vcf.gz"
    rows = [
        # de novo DEL: kid 0/1, parents 0/0
        "chr1\t1000\tsawfish:0:1:0:0\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;SVLEN=-5000;END=6000;HOMLEN=3\tGT:GQ:AD:CN\t0/1:40:10,8:1.5\t0/0:50:20,0:2.0\t0/0:45:18,0:2.0",
        # inherited INS: kid 0/1, dad 0/1 -> not a candidate
        "chr1\t20000\tsawfish:0:2:0:0\tN\tACGTACGTAC\t60\tPASS\tSVTYPE=INS;SVLEN=10\tGT:GQ:AD\t0/1:40:10,8\t0/1:50:12,9\t0/0:45:18,0",
        # parent missing -> not a candidate (absence of evidence)
        "chr1\t30000\tsawfish:0:3:0:0\tN\t<DUP>\t60\tPASS\tSVTYPE=DUP;SVLEN=2000;END=32000\tGT:GQ:AD\t0/1:40:10,8\t./.:.:.\t0/0:45:18,0",
        # hom-alt child, parents 0/0: still a candidate (impossible biologically -> artefact class, kept unfiltered)
        "chr1\t40000\tsawfish:0:4:0:0\tN\t<INV>\t60\tPASS\tSVTYPE=INV;SVLEN=800;END=40800;IMPRECISE\tGT:GQ:AD\t1/1:20:0,6\t0/0:50:20,0\t0/0:45:18,0",
    ]
    _write(p, rows)
    recs = list(C.sv_candidates(str(p), "fam", "kid", "dad", "mom"))
    ids = [r.variant_id for r in recs]
    assert ids == ["sawfish:0:1:0:0", "sawfish:0:4:0:0"], ids
    d = recs[0]
    assert d.variant_class == "SV" and d.class_payload["svtype"] == "DEL" and d.class_payload["svlen"] == -5000
    assert d.start == 1000 and d.end == 6000 and d.class_payload["homlen"] == 3 and d.class_payload["child_cn"] == "1.5"
    assert recs[1].class_payload["imprecise"] is True and recs[1].caller_gt == "1/1"


# ---------------------------------------------------------------------------- TR (TRGT-shaped)
def test_tr_candidates(tmp_path):
    p = tmp_path / "tr.vcf.gz"
    rows = [
        # expansion: kid 30/120 vs parents max 36, motif CAG (unit 3): 120 >= 36+3 -> candidate, allele idx 1
        "chr1\t5000\t.\tCAGCAG\t<TR>\t.\tPASS\tTRID=locus_A;END=5100;MOTIFS=CAG;STRUC=(CAG)n\tGT:AL:SD:ALLR:MC\t0/1:30,120:10,8:29-31,110-125:10,40\t0/0:30,36:12,11:.:.\t0/0:33,30:9,10:.:."
        .replace("0/1:30,120", "1/2:30,120"),
        # within parental range: no candidate
        "chr1\t6000\t.\tA\t<TR>\t.\tPASS\tTRID=locus_B;END=6050;MOTIFS=A;STRUC=(A)n\tGT:AL:SD\t0/1:20,25:10,8\t0/1:20,26:12,11\t0/0:22,24:9,10",
        # contraction: kid 10 vs parents min 40, motif ACGT (unit 4): 10 <= 40-4 -> candidate flagged contraction
        "chr1\t7000\t.\tACGT\t<TR>\t.\tPASS\tTRID=locus_C;END=7080;MOTIFS=ACGT;STRUC=(ACGT)n\tGT:AL:SD\t0/1:10,44:10,8\t0/0:40,44:12,11\t0/0:41,48:9,10",
        # uncalled parent: not evaluable
        "chr1\t8000\t.\tA\t<TR>\t.\tPASS\tTRID=locus_D;END=8050;MOTIFS=A;STRUC=(A)n\tGT:AL:SD\t0/1:20,90:10,8\t./.:.:.\t0/0:22,24:9,10",
    ]
    _write(p, rows)
    recs = list(C.tr_candidates(str(p), "fam", "kid", "dad", "mom", min_units=1))
    assert [r.variant_id for r in recs] == ["locus_A:a1", "locus_C:a0"], [r.variant_id for r in recs]   # TRID:a<allele index>: unique per allele
    a, c = recs
    assert a.class_payload["direction"] == "expansion" and a.class_payload["outlier_allele_idx"] == 1
    assert a.class_payload["delta_bp"] == 84 and a.class_payload["motif_unit_bp"] == 3 and a.alt == "<AL=120>"
    assert c.class_payload["direction"] == "contraction" and c.class_payload["delta_bp"] == 30
    # a stricter margin removes the contraction (30 bp < 8 units*4=32)
    assert [r.variant_id for r in C.tr_candidates(str(p), "fam", "kid", "dad", "mom", min_units=8)] == ["locus_A:a1"]


# ---------------------------------------------------------------------------- shared record round-trip and list merge
def test_record_roundtrip_and_merge(tmp_path):
    recs = [CandidateRecord("fam", "kid", "chr1:100:A:G", "chr1", 100, 100, "A", "G", "SNV", "deepvariant_glnexus", caller_gt="0/1",
                            caller_gq=30, class_payload={"allele_index": 1}),
            CandidateRecord("fam", "kid", "locus_A", "chr1", 5000, 5100, "<TR>", "<AL=120>", "TR", "trgt", class_payload={"trid": "locus_A"}),
            CandidateRecord("fam", "kid", "sawfish:0:1:0:0", "chr1", 1000, 6000, "N", "<DEL>", "SV", "sawfish", class_payload={"svtype": "DEL"})]
    out = tmp_path / "cand.tsv"
    assert write_candidates(recs, str(out)) == 3
    back = list(read_candidates(str(out)))
    assert [r.key for r in back] == [r.key for r in recs] and back[1].class_payload["trid"] == "locus_A" and back[0].caller_gq == 30
    hi = tmp_path / "hiconf.tsv"; hi.write_text("chrom\tpos\tref\talt\nchr1\t100\tA\tG\n")
    trl = tmp_path / "tr.tsv"; trl.write_text("trid\tgain_bp\nlocus_A\t84\n")
    merged = C.merge_lists(back, [("slivar_hiconf", "HIGH", str(hi)), ("denovo_tr", "LOW", str(trl))])
    assert merged[0].source_tier == "HIGH" and "slivar_hiconf" in merged[0].source_list
    assert merged[1].source_tier == "LOW" and "denovo_tr" in merged[1].source_list
    assert merged[2].source_tier == "UNFILTERED" and merged[2].source_list == "joint_vcf"
