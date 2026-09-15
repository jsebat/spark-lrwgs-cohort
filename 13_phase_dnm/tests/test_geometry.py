"""Phase-segment geometry and crossover-distance features from the M1 tables; registry drops of the B-block duplicates."""
import math

from phase_dnm.features import extract as X
from phase_dnm.features.registry import Registry


def test_geometry_features(tmp_path):
    ori = tmp_path / "o.tsv"
    ori.write_text("chrom\tphase_block_id\tstart\tend\torientation\tswitch_pos\n"
                   "chr1\t1000\t1000\t101000\tHAP1_PAT\t0\n"
                   "chr1\t200000\t200000\t210000\tHAP1_MAT\t0\n")
    cp = tmp_path / "c.tsv"
    cp.write_text("parent\tchrom\tleft_pos\tright_pos\tstatus\n"
                  "F\tchr1\t150000\t150100\tCROSSOVER\n"
                  "M\tchr1\t900000\t900200\tSWITCH_ERROR\n")
    g = X.Geometry(str(ori), str(cp))
    a = g.at("chr1", 51000)                                  # middle of a 100 kb segment
    assert a["c_block_len_log10"] == round(math.log10(100001), 3) and a["c_dist_block_edge_log10"] == round(math.log10(50001), 3)
    assert a["dist_crossover_log10"] == round(math.log10(99050 + 1), 3) and a["near_crossover_flag"] == 0
    b = g.at("chr1", 160000)                                 # outside any segment, 10 kb from the crossover
    assert b["c_block_len_log10"] is None and b["c_dist_block_edge_log10"] == -1.0 and b["near_crossover_flag"] == 1
    c = g.at("chr1", 200000)                                 # at a segment start
    assert c["c_dist_block_edge_log10"] == 0.0
    assert g.at("chr2", 5) == {}                             # no tables for this contig
    assert X.Geometry(None, None).at("chr1", 1) == {}


def test_registry_drops_superseded_b_block():
    reg = Registry()
    for name in ("child_alt_mapq_mean", "child_alt_readpos_frac_median", "child_alt_rq_mean"):
        assert name not in reg.features
    for name in ("c_alt_mapq0_frac", "c_alt_supp_frac", "c_alt_readlen_median", "c_alt_rq_mean", "c_block_len_log10", "c_dist_block_edge_log10"):
        assert name in reg.features and reg.features[name].rf_safe
    assert not reg.features["dist_crossover_log10"].rf_safe and not reg.features["near_crossover_flag"].rf_safe
