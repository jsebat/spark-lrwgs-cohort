"""Module 3: the P15 decision in both modes, class-specific columns, table assembly, VCF INFO, and the P18 concordance."""
import csv
import json
import os

from phase_dnm import concordance as C
from phase_dnm import integrate as I
from phase_dnm.io import vcfinfo as V

P = I.FinalParams(tau={"SNV": 0.8, "SV": 0.7}, tau_rescue={"SNV": 0.4}, phase_only_min_score=0.9)


def row(cls="germline_DNM_phased", hap=6, ps=0.98, vclass="SNV", **kw):
    r = dict(phase_class=cls, hap_obs_k5=str(hap), phase_score=str(ps), variant_class=vclass)
    r.update(kw)
    return r


def test_decide_rf_mode():
    assert I.decide(row(), P, 0.9) == dict(dnm_call="YES", call_mode="rf+phase", decision_reason="RF", mosaic_flag=0)
    # rescue: below tau, above tau_rescue, phased germline with six haplotypes observed
    assert I.decide(row(), P, 0.5)["decision_reason"] == "RESCUED"
    assert I.decide(row(hap=5), P, 0.5)["dnm_call"] == "NO" and I.decide(row(hap=5), P, 0.5)["decision_reason"] == "BELOW_TAU"
    assert I.decide(row(cls="germline_DNM_unphased"), P, 0.5)["dnm_call"] == "NO"
    # demotion beats any probability
    d = I.decide(row(cls="inherited_missed_in_parent"), P, 0.99)
    assert d["dnm_call"] == "NO" and d["decision_reason"] == "DEMOTED:inherited_missed_in_parent"
    assert I.decide(row(cls="phase_conflict_artifact"), P, 0.99)["dnm_call"] == "NO"
    # mosaics are never YES, flagged
    d = I.decide(row(cls="child_postzygotic_mosaic"), P, 0.99)
    assert d["dnm_call"] == "NO" and d["mosaic_flag"] == 1 and d["decision_reason"].startswith("MOSAIC:")
    # a class without tau falls back to the provisional rule even when rf_prob exists
    assert I.decide(row(vclass="TR"), P, 0.99)["call_mode"] == "phase_only"


def test_decide_provisional_mode():
    assert I.decide(row(), P, None) == dict(dnm_call="YES", call_mode="phase_only", decision_reason="PHASE_ONLY", mosaic_flag=0)
    assert I.decide(row(ps=0.5), P, None)["decision_reason"] == "LOW_POSTERIOR"
    assert I.decide(row(hap=4), P, None)["decision_reason"] == "HAP_UNOBSERVED"
    assert I.decide(row(cls="inconclusive"), P, None)["decision_reason"] == "NOT_PHASED_GERMLINE"
    assert I.decide(row(cls="germline_DNM_unphased", ps=0.99), P, None)["dnm_call"] == "NO"


def test_class_columns_tr_and_sv():
    pl = json.dumps({"trid": "chr1_100_200_CAG", "motifs": "CAG,CCG", "child_AL": [30, 60], "father_AL": [30, 33], "mother_AL": [28, 30], "father_gt": "0|1", "mother_gt": "1|0"})
    r = dict(variant_class="TR", class_payload=pl, caller_gt="0|1", child_hap1_is="M", F_transmitted_hap="2", M_transmitted_hap="1")
    c = I.class_columns(r)
    assert c["trid"] == "chr1_100_200_CAG" and c["motif"] == "CAG"
    assert c["child_AL_mat"] == 30 and c["child_AL_pat"] == 60          # hap1 is maternal: hap1 allele (idx 0 = 30) is maternal
    assert c["father_AL_T"] == 33 and c["father_AL_U"] == 30            # father transmitted hap2 -> allele idx 1
    assert c["mother_AL_T"] == 30 and c["mother_AL_U"] == 28            # mother gt 1|0, transmitted hap1 -> allele idx 1 = 30
    r2 = dict(variant_class="TR", class_payload=pl, caller_gt="0/1", child_hap1_is="M", F_transmitted_hap="2", M_transmitted_hap="1")
    assert I.class_columns(r2)["child_AL_pat"] is None                    # unphased genotype: not resolvable
    sv = dict(variant_class="SV", class_payload=json.dumps({"svtype": "DEL", "svlen": -500, "imprecise": False}))
    assert I.class_columns(sv) == dict(svtype="DEL", svlen=-500, bp_precision="PRECISE")


def _write(path, rows, cols):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})


def test_integrate_table_and_vcf(tmp_path):
    ev_cols = ["family_id", "sample_id", "variant_id", "chrom", "start", "end", "ref", "alt", "variant_class", "caller", "caller_gt", "caller_qual",
               "source_tier", "source_list", "mask_overlap", "phase_class", "rule_score", "flags", "parent_of_origin", "poo_reason", "poo_confidence",
               "child_hap1_is", "F_transmitted_hap", "M_transmitted_hap", "hap_obs_k3", "hap_obs_k5", "c_alt_hap_frac", "c_alt_hapO", "t_alt_reads",
               "u_alt_reads", "t_dp", "phase_score", "lik_best_alternative", "lik_log10lr_germline", "thresholds_version", "class_payload"]
    base = dict(family_id="fam", sample_id="kid", chrom="chr1", ref="A", alt="G", variant_class="SNV", caller="dv", caller_gt="0/1", caller_qual="50",
                source_tier="UNFILTERED", source_list="joint_vcf", mask_overlap="0", rule_score="6", flags="", poo_reason="OK", poo_confidence="1.0",
                child_hap1_is="P", F_transmitted_hap="1", M_transmitted_hap="2", hap_obs_k3="6", c_alt_hapO="0", u_alt_reads="0", t_dp="10",
                lik_best_alternative="parental_mosaic", lik_log10lr_germline="0.8", thresholds_version="0.2.0", class_payload="{}")
    ev = [dict(base, variant_id="v1", start="100", end="100", phase_class="germline_DNM_phased", hap_obs_k5="6", c_alt_hap_frac="1.0", t_alt_reads="0", parent_of_origin="paternal", phase_score="0.98"),
          dict(base, variant_id="v2", start="200", end="200", phase_class="inherited_missed_in_parent", hap_obs_k5="6", c_alt_hap_frac="1.0", t_alt_reads="6", parent_of_origin="maternal", phase_score="0.0"),
          dict(base, variant_id="v3", start="300", end="300", phase_class="germline_DNM_unphased", hap_obs_k5="4", c_alt_hap_frac="0.9", t_alt_reads="", parent_of_origin="undetermined", phase_score="0.95"),
          dict(base, variant_id="v4", start="400", end="400", phase_class="child_postzygotic_mosaic", hap_obs_k5="6", c_alt_hap_frac="0.4", t_alt_reads="0", parent_of_origin="paternal", phase_score="0.0")]
    evp = str(tmp_path / "ev.tsv"); _write(evp, ev, ev_cols)
    fe_cols = ["family_id", "sample_id", "variant_id", "chrom", "start", "end", "ref", "alt", "variant_class", "caller", "child_GQ", "AR_child", "c_alt_hap_frac"]
    fe = [dict(family_id="fam", sample_id="kid", variant_id=v, child_GQ=str(30 + i), AR_child="0.02", c_alt_hap_frac="1.0") for i, v in enumerate(["v1", "v2", "v3", "v4"])]
    fep = str(tmp_path / "fe.tsv"); _write(fep, fe, fe_cols)
    out = str(tmp_path / "final.tsv")
    summ = I.integrate_table(evp, fep, out, I.FinalParams(), "snv_indel")
    assert summ["counts"] == {"rows": 4, "YES": 1, "NO": 3, "mosaic": 1}
    assert summ["reasons"] == {"PHASE_ONLY": 1, "DEMOTED:inherited_missed_in_parent": 1, "NOT_PHASED_GERMLINE": 1, "MOSAIC:child_postzygotic_mosaic": 1}
    rows = V.rows_from_final(out)
    cols = list(rows[0].keys())
    assert cols[:8] == I.FINAL_CORE[:8] and "child_GQ" in cols and "AR_child" in cols
    assert cols.index("rf_prob") < cols.index("phase_class") < cols.index("child_GQ")
    by = {r["variant_id"]: r for r in rows}
    assert by["v1"]["dnm_call"] == "YES" and by["v1"]["call_mode"] == "phase_only" and by["v1"]["child_GQ"] == "30"
    assert by["v2"]["dnm_call"] == "NO" and by["v2"]["decision_reason"].startswith("DEMOTED") and by["v2"]["transmitted_parent_alt_reads"] == "6"
    assert by["v3"]["decision_reason"] == "NOT_PHASED_GERMLINE" and by["v4"]["mosaic_flag"] == "1"
    assert by["v1"]["child_alt_hap_frac"] == "1.0" and by["v1"]["child_alt_other_hap"] == "0" and by["v1"]["rf_prob"] == ""
    # with M4 probabilities the same table flips to rf+phase
    rf = {"v1": 0.2, "v3": 0.95}
    summ2 = I.integrate_table(evp, fep, str(tmp_path / "f2.tsv"), I.FinalParams(tau={"SNV": 0.8}, tau_rescue={"SNV": 0.1}), "snv_indel", rf_probs=rf)
    by2 = {r["variant_id"]: r for r in V.rows_from_final(str(tmp_path / "f2.tsv"))}
    assert by2["v1"]["decision_reason"] == "RESCUED" and by2["v3"]["decision_reason"] == "RF" and by2["v3"]["dnm_call"] == "YES"
    assert by2["v2"]["decision_reason"].startswith("DEMOTED")
    # VCF
    info = V.info_string(by["v1"])
    assert "PDNM_CALL=YES" in info and "PDNM_POO=paternal" in info and "PDNM_SCORE=0.98" in info and "PDNM_MOSAIC" not in info
    assert "PDNM_MOSAIC" in V.info_string(by["v4"])
    vcf = str(tmp_path / "f.vcf")
    assert V.write_sites_vcf(rows, vcf) == 4
    lines = open(vcf).read().splitlines()
    assert lines[0] == "##fileformat=VCFv4.3" and any(l.startswith("##INFO=<ID=PDNM_CALL") for l in lines)
    body = [l for l in lines if not l.startswith("#")]
    assert body[0].split("\t")[:5] == ["chr1", "100", "v1", "A", "G"] and "PDNM_CALL=YES" in body[0]


def test_concordance_small_and_sv(tmp_path):
    b = tmp_path / "baselines"; b.mkdir()
    with open(b / "denovo_tiered.tsv", "w") as fh:
        fh.write("family\tproband\tchrom\tpos\tref\talt\tgene\ttier\n")
        fh.write("fam\tkid\tchr1\t100\tA\tG\tX\tt1\n")          # concordant
        fh.write("fam\tkid\tchr1\t200\tC\tT\tY\tt2\n")          # original only (module: demoted)
        fh.write("fam\tkid\tchr1\t900\tG\tA\tZ\tt1\n")          # never a candidate in the module
    with open(b / "denovo_sv_all.bed", "w") as fh:
        fh.write("chr1\t1000\t2000\tfam\tkid\tDEL\t1000\n")
        fh.write("chr2\t5000\t5002\tfam\tkid\tINS\t300\n")
    final_small = [dict(sample_id="kid", variant_id="a", variant_class="SNV", chrom="chr1", start="100", end="100", ref="A", alt="G", dnm_call="YES", phase_class="germline_DNM_phased", decision_reason="PHASE_ONLY", parent_of_origin="paternal", source_tier="UNFILTERED", mask_overlap="0"),
                   dict(sample_id="kid", variant_id="b", variant_class="SNV", chrom="chr1", start="200", end="200", ref="C", alt="T", dnm_call="NO", phase_class="inherited_missed_in_parent", decision_reason="DEMOTED:inherited_missed_in_parent", parent_of_origin="maternal", source_tier="UNFILTERED", mask_overlap="0"),
                   dict(sample_id="kid", variant_id="c", variant_class="SNV", chrom="chr1", start="300", end="300", ref="T", alt="C", dnm_call="YES", phase_class="germline_DNM_phased", decision_reason="PHASE_ONLY", parent_of_origin="maternal", source_tier="UNFILTERED", mask_overlap="1"),
                   dict(sample_id="kid", variant_id="d", variant_class="SNV", chrom="chr1", start="400", end="400", ref="T", alt="C", dnm_call="NO", phase_class="inconclusive", decision_reason="NOT_PHASED_GERMLINE", parent_of_origin="undetermined", source_tier="UNFILTERED", mask_overlap="0")]
    per_row, per_prob, s = C.concordance(final_small, "snv_indel", str(b))
    st = {d["variant_id"]: d["status"] for d in per_row}
    assert st == {"a": "concordant_YES", "b": "original_only", "c": "module_only", "d": "neither"}
    assert s["concordant_YES"] == 1 and s["original_only"] == 1 and s["module_only"] == 1 and s["original_unseen_as_candidate"] == 1
    assert s["original_only_by_phase_class"] == {"inherited_missed_in_parent": 1} and s["module_only_in_mask"] == 1
    assert s["module_YES_paternal_fraction"] == 0.5 and per_prob[0]["module_YES"] == 2 and per_prob[0]["original_unseen"] == 1
    assert s["per_proband_original_median"] == 3 and s["per_proband_module_YES_median"] == 2
    final_sv = [dict(sample_id="kid", variant_id="s1", variant_class="SV", chrom="chr1", start="1010", end="1990", svtype="DEL", dnm_call="YES", phase_class="germline_DNM_phased", decision_reason="PHASE_ONLY", parent_of_origin="paternal", source_tier="UNFILTERED", mask_overlap="0"),
                dict(sample_id="kid", variant_id="s2", variant_class="SV", chrom="chr2", start="5100", end="5100", svtype="INS", dnm_call="NO", phase_class="inconclusive", decision_reason="NOT_PHASED_GERMLINE", parent_of_origin="undetermined", source_tier="UNFILTERED", mask_overlap="0"),
                dict(sample_id="kid", variant_id="s3", variant_class="SV", chrom="chr3", start="1", end="500", svtype="DEL", dnm_call="YES", phase_class="germline_DNM_phased", decision_reason="PHASE_ONLY", parent_of_origin="maternal", source_tier="UNFILTERED", mask_overlap="0")]
    per_row, per_prob, s = C.concordance(final_sv, "sv", str(b))
    st = {d["variant_id"]: d["status"] for d in per_row}
    assert st == {"s1": "concordant_YES", "s2": "original_only", "s3": "module_only"} and s["original_unseen_as_candidate"] == 0


def test_small_variant_matching_is_representation_independent():
    assert C.norm_allele(100, "GAT", "GT") == (101, "A", "")          # padded deletion -> first changed base, deleted seq
    assert C.norm_allele(101, "AT", "T") == (101, "A", "")            # anchor-after form of the same deletion
    assert C.norm_allele(99, "GA", "GAT") == (101, "", "T")           # padded insertion
    assert C.norm_allele(100, "A", "AT") == (101, "", "T")
    assert C.norm_allele(100, "A", "G") == (100, "A", "G")            # SNV untouched
    small = {("kid", "chr1", 101, "A", ""): {"tier": "t1"}}
    assert C.small_match(small, "kid", "chr1", 100, "GAT", "GT") == ("kid", "chr1", 101, "A", "")      # exact after normalisation
    assert C.small_match(small, "kid", "chr1", 105, "CA", "C") == ("kid", "chr1", 101, "A", "")        # same 1-bp deletion, shifted (repeat)
    assert C.small_match(small, "kid", "chr1", 105, "CAA", "C") is None                                # different length change
    assert C.small_match(small, "kid", "chr1", 101, "A", "T") is None                                   # SNVs: exact only
    # a baseline insertion truncated at 30 characters matches the full-length module allele at the same position
    ins = "GGCCGAGGCGGGTGGATCATGAGGTCAGGAGATCGAGACCAACCTGG"
    small2 = {("kid", "chr2", 501, "", ins[:29]): {"tier": ""}}
    assert C.small_match(small2, "kid", "chr2", 500, "C", "C" + ins) == ("kid", "chr2", 501, "", ins[:29])
    assert C.small_match(small2, "kid", "chr2", 500, "C", "C" + "ACGT" * 12) is None                 # different insertion, same length class
