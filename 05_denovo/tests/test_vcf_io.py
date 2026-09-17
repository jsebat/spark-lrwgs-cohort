from phase_dnm.io.vcf import VcfReader, iter_trio, read_manifest, trio_of, normalise_sex, _parse_gt


def test_parse_gt():
    assert _parse_gt("0|1") == ((0, 1), True)
    assert _parse_gt("1/0") == ((1, 0), False)
    assert _parse_gt("1") == ((1,), False)          # hemizygous
    assert _parse_gt("./.") == (None, False)
    assert _parse_gt(".|1") == (None, False)
    assert _parse_gt("0/2") == ((0, 2), False)      # multi-allelic index


def test_reader_header_and_sites(minitrio):
    r = VcfReader(minitrio["paths"]["child"])
    assert r.contigs == ["chr21", "chr22"]
    assert r.samples == [minitrio["ids"][0]]
    sites = list(r.sites())
    assert len(sites) > 2000
    phased = [s for s in sites if s.phased]
    assert phased and all(s.ps is not None for s in phased)
    assert all(s.gq is not None and s.dp is not None for s in sites)
    # positions non-decreasing within contig
    last = ("", -1)
    for s in sites:
        assert (s.chrom, s.pos) > last or s.chrom != last[0]
        last = (s.chrom, s.pos)


def test_trio_merge_matches_every_child_site(minitrio):
    p = minitrio["paths"]
    trio = list(iter_trio(VcfReader(p["child"]), VcfReader(p["father"]), VcfReader(p["mother"])))
    n_child = sum(1 for _ in VcfReader(p["child"]).sites())
    assert len(trio) == n_child
    # the simulator writes every site for every sample, so parents are never missing
    assert all(t.father is not None and t.mother is not None for t in trio)
    assert all(t.father.pos == t.pos and t.mother.pos == t.pos for t in trio)


def test_trio_merge_from_joint_vcf(minitrio):
    p = minitrio["paths"]
    c, f, m = minitrio["ids"]
    trio = list(iter_trio(VcfReader(p["joint"]), VcfReader(p["joint"]), VcfReader(p["joint"]), c, f, m))
    assert len(trio) == sum(1 for _ in VcfReader(p["child"]).sites())


def test_manifest_and_duo_rule(minitrio, tmp_path):
    man = read_manifest(minitrio["paths"]["manifest"])
    c, f, m = minitrio["ids"]
    assert trio_of(man, c) == (f, m)
    assert normalise_sex(man[c]["sex"]) == "F"
    # a duo (father '0') must be refused, DESIGN P19
    man[c]["father_id"] = "0"
    import pytest
    with pytest.raises(ValueError):
        trio_of(man, c)
