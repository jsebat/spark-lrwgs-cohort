"""Folds are family-grouped and balanced; swaps never cross a fold; the whole construction is deterministic by seed."""
import pytest

from phase_dnm.train import folds as F
from phase_dnm.train import swap as S


def manifest():
    rows = []
    for i in range(33):
        fam = "fam%02d" % i
        rows += [dict(sample_id="f%02d" % i, family_id=fam, father_id="0", mother_id="0", role="founder"),
                 dict(sample_id="m%02d" % i, family_id=fam, father_id="0", mother_id="0", role="founder"),
                 dict(sample_id="c%02d" % i, family_id=fam, father_id="f%02d" % i, mother_id="m%02d" % i, role="offspring")]
    rows.append(dict(sample_id="c00b", family_id="fam00", father_id="f00", mother_id="m00", role="offspring"))       # a quad
    rows += [dict(sample_id="mduo", family_id="duo", father_id="0", mother_id="0", role="founder"),
             dict(sample_id="cduo", family_id="duo", father_id="0", mother_id="mduo", role="offspring")]             # a duo: excluded
    return rows


def test_families_from_manifest_excludes_duos_and_keeps_quads():
    fams = F.families_from_manifest(manifest(), sources={"fam05": "blood"})
    assert len(fams) == 33 and all(f.family_id != "duo" for f in fams)
    q = [f for f in fams if f.family_id == "fam00"][0]
    assert q.children == ("c00", "c00b")
    assert [f for f in fams if f.source == "blood"][0].family_id == "fam05"


def test_outer_folds_balanced_deterministic_and_blood_not_alone():
    fams = F.families_from_manifest(manifest(), sources={"fam05": "blood"})
    a = F.outer_folds(fams, n_folds=5, seed=3)
    b = F.outer_folds(fams, n_folds=5, seed=3)
    assert a == b
    sizes = [sum(1 for v in a.values() if v == k) for k in range(5)]
    assert sum(sizes) == 33 and max(sizes) - min(sizes) <= 1
    assert a["fam05"] == 0 and sizes[0] >= 3
    assert F.outer_folds(fams, n_folds=5, seed=4) != a                     # a different seed is a different split
    with pytest.raises(ValueError):
        F.outer_folds(fams, n_folds=20, seed=0)                            # too many folds for 33 families


def test_swap_is_fold_closed_and_covers_every_family():
    fams = F.families_from_manifest(manifest())
    assign = F.outer_folds(fams, n_folds=5, seed=1)
    parents = {"fam%02d" % i: ("f%02d" % i, "m%02d" % i) for i in range(33)}
    table = S.pairings(fams, assign, parents, seed=1)
    S.check_closed(table, assign)
    donors = {r["child_family"] for r in table}
    assert donors == {f.family_id for f in fams}
    # the quad's two children go to the same surrogate parents
    q = [r for r in table if r["child_family"] == "fam00"]
    assert len({(r["surrogate_father"], r["surrogate_mother"]) for r in q}) == 1 and len(q) == 2
    # no child is paired with its own parents; ids carry no sample identifiers
    assert all(r["surrogate_father"] != "f" + r["child"][1:3] for r in table if r["child"].startswith("c") and len(r["child"]) == 3)
    assert all(r["synthetic_id"].startswith("syn_") and r["child"] not in r["synthetic_id"] for r in table)
    assert table == S.pairings(fams, assign, parents, seed=1)               # deterministic
    with pytest.raises(ValueError):
        bad = dict(table[0]); bad["parent_family"] = [f for f in assign if assign[f] != assign[bad["child_family"]]][0]
        S.check_closed([bad], assign)
