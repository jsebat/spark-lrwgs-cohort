import csv
import os

import pytest

from phase_dnm.evidence.review import evidence_columns
from phase_dnm.features import extract as X
from phase_dnm.features.registry import Registry, RfSafetyError
from phase_dnm.records import CandidateRecord, write_candidates


def test_registry_loads_and_gates():
    reg = Registry()
    assert reg.version and len(reg.features) > 100
    unsafe = set(reg.unsafe_columns())
    for name in ("parent_of_origin", "t_alt_reads", "u_alt_reads", "dist_crossover_log10", "sib_shared", "graphtyper_outcome"):
        assert name in unsafe, name
    snv_rf = reg.rf_matrix_columns("SNV")
    assert "child_AR" in snv_rf and "c_alt_hap_frac" in snv_rf and "hap_obs_k5" in snv_rf
    assert not (set(snv_rf) & unsafe)
    reg.assert_rf_safe(snv_rf)
    with pytest.raises(RfSafetyError):
        reg.assert_rf_safe(snv_rf + ["t_alt_reads"])
    with pytest.raises(RfSafetyError):
        reg.assert_rf_safe(["not_a_registered_feature"])
    m = reg.manifest()
    assert len(m["features_sha256"]) == 64 and m["n_rf_safe"] > 50


def test_synthdnm_universal_derivations():
    rec = CandidateRecord("fam", "kid", "chr1:100:A:G", "chr1", 100, 100, "A", "G", "SNV", "deepvariant_glnexus",
                          caller_gt="0/1", caller_gq=45, caller_dp=30, child_ad="14,16", father_gt="0/0", mother_gt="0/0",
                          father_gq=50, mother_gq=40, father_dp=28, mother_dp=25, father_ad="28,0", mother_ad="25,0",
                          child_pl="60,0,70", father_pl="0,80,900", mother_pl="0,70,800", class_payload={"allele_index": 1, "n_alts": 1})
    f = X.caller_features(rec, "M")
    assert f["child_AR"] == 14 / 17 and f["min_AR"] == 25.0 and f["max_AR"] == 28.0
    assert (f["child_GQ"], f["min_GQ"], f["max_GQ"]) == (45, 40, 50)
    assert (f["child_PL0"], f["child_PL1"], f["child_PL2"]) == (60, 0, 70)
    assert (f["min_PL1"], f["max_PL1"]) == (70, 80)
    assert abs(f["child_AB"] - 16 / 30) < 1e-9 and f["indel_flag"] == 0 and f["haploid_flag"] == 0
    # male chrX outside the PAR is haploid; PAR is not
    x = CandidateRecord("fam", "kid", "chrX:50000000:A:G", "chrX", 50_000_000, 50_000_000, "A", "G", "SNV", "dv")
    assert X.caller_features(x, "1")["haploid_flag"] == 1
    par = CandidateRecord("fam", "kid", "chrX:100000:A:G", "chrX", 100_000, 100_000, "A", "G", "SNV", "dv")
    assert X.caller_features(par, "1")["haploid_flag"] == 0 and X.caller_features(x, "2")["haploid_flag"] == 0
    # multi-allelic second allele: PL indices 0, 3, 5
    r2 = CandidateRecord("fam", "kid", "chr1:200:A:T", "chr1", 200, 200, "A", "T", "SNV", "dv", child_ad="10,0,9",
                         child_pl="50,300,400,0,350,60", class_payload={"allele_index": 2, "n_alts": 2})
    f2 = X.caller_features(r2, "F")
    assert (f2["child_PL0"], f2["child_PL1"], f2["child_PL2"]) == (50, 0, 60) and f2["multiallelic"] == 1
    assert f2["child_AR"] == 10 / 10


def test_class_specific_caller_features():
    sv = CandidateRecord("fam", "kid", "sawfish:0:1:0:0", "chr1", 1000, 6000, "N", "<DEL>", "SV", "sawfish", child_ad="10,8",
                         father_ad="20,0", mother_ad="18,1", class_payload={"svtype": "DEL", "svlen": -5000, "homlen": 3})
    f = X.caller_features(sv, "F")
    assert f["svtype"] == "DEL" and abs(f["svlen_log10"] - 3.699) < 1e-3 and f["child_sv_support"] == 8 and f["max_parent_sv_support"] == 1
    tr = CandidateRecord("fam", "kid", "locus_A", "chr1", 5000, 5100, "<TR>", "<AL=120>", "TR", "trgt",
                         class_payload={"child_AL": [30, 120], "outlier_allele_idx": 1, "child_SD": "10,8", "father_SD": "12,11",
                                        "mother_SD": "9,10", "delta_units": 28.0, "motif_unit_bp": 3, "motifs": "CAG",
                                        "child_ALLR": "29-31,110-125"})
    f = X.caller_features(tr, "F")
    assert f["child_AL_expanded"] == 120 and f["child_AL_other"] == 30 and f["child_SD_expanded"] == 8
    assert f["min_parent_SD"] == 9 and f["child_ALLR_width"] == 15 and f["motif_len"] == 3 and f["locus_len_ref"] == 101


def test_bed_mask():
    import tempfile
    d = tempfile.mkdtemp()
    p = os.path.join(d, "mask.bed")
    open(p, "w").write("chr1\t1000\t2000\nchr1\t1500\t3000\nchr2\t10\t20\n")
    m = X.BedMask([p])
    assert m.iv["chr1"] == [(1000, 3000)]                    # merged
    assert m.flag("chr1", 1500, 1500) == 1 and m.flag("chr1", 5000, 5000) == 0
    assert m.flag("chr1", 2500, 4000) == 0                    # 500/1501 < 0.5 of a >50 bp interval
    assert m.flag("chr1", 2000, 3500) == 1                    # 1000/1501 >= 0.5


def test_extract_child_end_to_end(tmp_path):
    reg = Registry()
    recs = [CandidateRecord("fam", "kid", "chr1:100:A:G", "chr1", 100, 100, "A", "G", "SNV", "dv", caller_gt="0/1", caller_gq=45,
                            caller_dp=30, child_ad="14,16", father_gq=50, mother_gq=40, father_ad="28,0", mother_ad="25,0",
                            child_pl="60,0,70", father_pl="0,80,900", mother_pl="0,70,800", class_payload={"allele_index": 1, "n_alts": 1})]
    cpath = tmp_path / "kid.snv_indel.candidates.tsv"; write_candidates(recs, str(cpath))
    cols = evidence_columns((3, 5))
    ev = {c: "" for c in cols}
    ev.update(family_id="fam", sample_id="kid", variant_id="chr1:100:A:G", chrom="chr1", start=100, end=100, ref="A", alt="G",
              variant_class="SNV", caller="dv", phase_class="germline_DNM_phased", rule_score=6, flags="", parent_of_origin="paternal",
              poo_reason="OK", poo_confidence=1.0, c_alt_hapA=10, c_alt_hapO=0, c_alt_hap_frac=1.0, c_alt_confined=1, hap_obs_k5=6,
              hap_obs_k3=6, p_max_alt_any_hap=0, t_alt_reads=0, t_dp=11, u_alt_reads=0, u_dp=10)
    epath = tmp_path / "kid.snv_indel.evidence.tsv"
    with open(epath, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerow(ev)
    mask = tmp_path / "mask.bed"; mask.write_text("chr1\t50\t150\n")
    out = tmp_path / "kid.snv_indel.features.tsv"; rf = tmp_path / "kid.snv_indel.features.rf.tsv"
    summ = X.extract_child(str(epath), str(cpath), reg, "M", str(out), str(rf), mask=X.BedMask([str(mask)]))
    assert summ["rows"] == 1 and summ["features_produced"] > 25
    row = next(csv.DictReader(open(out), delimiter="\t"))
    assert row["child_AR"].startswith("0.8235") and row["hap_obs_k5"] == "6" and row["segdup_overlap"] == "1"
    assert row["t_alt_reads"] == "0" and row["parent_of_origin"] == "paternal"      # full table keeps them
    rf_row = next(csv.DictReader(open(rf), delimiter="\t"))
    assert "t_alt_reads" not in rf_row and "parent_of_origin" not in rf_row and "poo_confidence" not in rf_row
    assert rf_row["c_alt_hap_frac"] == "1.0" and rf_row["child_GQ"] == "45"
    # everything never produced is reported, not silently absent
    assert "gnomad_af" in summ["never_produced"] and "mappability_k100" in summ["never_produced"]
