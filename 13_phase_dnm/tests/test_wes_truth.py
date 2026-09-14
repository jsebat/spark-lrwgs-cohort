"""WES external truth: target lookup, site collection and the trio labelling rules (pure parts)."""
from phase_dnm.eval import wes_truth as W
from phase_dnm.records import CandidateRecord, write_candidates


def test_targets_and_exonic_sites(tmp_path):
    bed = tmp_path / "t.bed"
    bed.write_text("chr1\t100\t200\nchr1\t150\t300\nchr2\t10\t20\n")
    t = W.Targets([str(bed)])
    assert t.iv["chr1"] == [(100, 300)]                       # merged
    assert t.covers("chr1", 101) and t.covers("chr1", 300) and not t.covers("chr1", 301) and not t.covers("chr1", 100)
    assert t.covers("chr2", 11) and not t.covers("chr3", 11)
    recs = [CandidateRecord("f", "k", "a", "chr1", 150, 150, "A", "G", "SNV", "dv"), CandidateRecord("f", "k", "b", "chr1", 500, 500, "A", "G", "SNV", "dv"),
            CandidateRecord("f", "k", "c", "chr2", 15, 15, "A", "G", "SNV", "dv")]
    p = str(tmp_path / "c.tsv"); write_candidates(recs, p)
    sites = W.exonic_sites([p], t)
    assert sites == {"chr1": {150}, "chr2": {15}}
    out = str(tmp_path / "r.bed")
    assert W.write_regions(sites, "chr1", out) == 1 and open(out).read() == "chr1\t148\t151\n"


def test_label_trio_rules():
    # WES de novo: child het with alt reads, parents clean hom-ref
    assert W.label_trio("0/1:45:30:15,15", "0/0:50:28:28,0", "0/0:48:31:31,0", 1) == (1, "wes_de_novo")
    assert W.label_trio("1/1:27:13:0,13", "0/0:50:46:46,0", "0/0:50:50:50,0", 1) == (-1, "child_hom_alt")     # not a single DNM (audit 2026-09-14)
    # a parent carries -> inherited (label 0), even if the child looks de novo
    assert W.label_trio("0/1:45:30:15,15", "0/1:50:28:14,14", "0/0:48:31:31,0", 1) == (0, "parent_carries_father")
    assert W.label_trio("0/1:45:30:15,15", "0/0:50:28:25,3", "0/0:48:31:31,0", 1) == (0, "parent_carries_father")   # 3 alt reads in a "hom-ref" parent
    # one parental alt read: not evaluable (the slivar rule's AD[1]==0 territory)
    assert W.label_trio("0/1:45:30:15,15", "0/0:50:28:27,1", "0/0:48:31:31,0", 1) == (-1, "father_alt_read")
    # child hom-ref with enough depth refutes the candidate
    assert W.label_trio("0/0:50:40:40,0", "0/0:50:28:28,0", "0/0:48:31:31,0", 1) == (0, "child_hom_ref")
    assert W.label_trio("0/0:50:12:12,0", "0/0:50:28:28,0", "0/0:48:31:31,0", 1) == (-1, "child_inconclusive")   # DP 12 < 20
    # low-quality child / parent -> not evaluable
    assert W.label_trio("0/1:12:30:15,15", "0/0:50:28:28,0", "0/0:48:31:31,0", 1) == (-1, "child_low_quality")
    assert W.label_trio("0/1:45:30:15,15", "0/0:9:28:28,0", "0/0:48:31:31,0", 1) == (-1, "father_low_quality")
    assert W.label_trio("./.:.:.:.", "0/0:50:28:28,0", "0/0:48:31:31,0", 1) == (-1, "child_uncalled")


def test_label_children_joins_by_canonical_key(tmp_path):
    ex = tmp_path / "chr1.tsv"
    # split multi-allelic record with padded alleles: GAT>GT is the same deletion as AT>T at 101
    ex.write_text("chr1\t100\tGAT\tGT\t0/1:45:30:15,15\t0/0:50:28:28,0\t0/0:48:31:31,0\n"
                  "chr1\t300\tC\tT\t0/0:50:40:40,0\t0/0:50:28:28,0\t0/0:48:31:31,0\n")
    recs = [CandidateRecord("f", "kid", "v1", "chr1", 101, 102, "AT", "T", "INDEL", "dv"), CandidateRecord("f", "kid", "v2", "chr1", 300, 300, "C", "T", "SNV", "dv"),
            CandidateRecord("f", "kid", "v3", "chr1", 900, 900, "A", "G", "SNV", "dv")]
    cp = str(tmp_path / "kid.tsv"); write_candidates(recs, cp)
    st = W.label_children(str(ex), ["kid", "dad", "mom"], {"kid": ("dad", "mom")}, {"kid": cp}, str(tmp_path / "out"), chrom_filter="chr1")
    assert st["kid"] == {"wes_de_novo": 1, "child_hom_ref": 1, "not_in_wes": 1}
    rows = open(tmp_path / "out" / "kid.snv_indel.wes.chr1.tsv").read().splitlines()
    assert rows[1].startswith("v1\t1\twes_de_novo") and rows[2].startswith("v2\t0\tchild_hom_ref") and rows[3].startswith("v3\t-1\tnot_in_wes")
