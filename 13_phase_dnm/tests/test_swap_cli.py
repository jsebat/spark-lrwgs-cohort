"""The swap CLI writes fold tables, synthetic-trio tables and per-trio manifests that the downstream steps consume unchanged."""
import csv
import os
import subprocess
import sys

from phase_dnm.io.vcf import read_manifest, trio_of

COLS = ["sample_id", "family_id", "father_id", "mother_id", "sex", "affected", "role", "lr_haplotagged_bam", "lr_family_sv_vcf",
        "lr_family_smallvar_vcf", "lr_family_trgt_vcf"]


def write_manifest(path, n_fam=12):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for i in range(n_fam):
            fam = "FX%02d" % i
            for sid, fa, mo, sex, role in (("p%02d" % i, "0", "0", "1", "founder"), ("q%02d" % i, "0", "0", "2", "founder"),
                                           ("k%02d" % i, "p%02d" % i, "q%02d" % i, "1", "offspring")):
                w.writerow(dict(sample_id=sid, family_id=fam, father_id=fa, mother_id=mo, sex=sex, affected="2" if role == "offspring" else "1",
                                role=role, lr_haplotagged_bam=".", lr_family_sv_vcf=".", lr_family_smallvar_vcf=".", lr_family_trgt_vcf="."))


def test_swap_cli_outputs(tmp_path):
    man = str(tmp_path / "manifest.tsv"); write_manifest(man)
    out = str(tmp_path / "train")
    env = dict(os.environ, PYTHONPATH=os.path.join(os.path.dirname(__file__), "..", "src"))
    r = subprocess.run([sys.executable, "-m", "phase_dnm.cli", "swap", "--manifest", man, "--out-dir", out, "--n-folds", "3", "--seeds", "0,1",
                        "--blood-family-prefix", "FX05"], capture_output=True, text=True, env=env)
    assert r.returncode == 0, r.stderr
    folds = list(csv.DictReader(open(os.path.join(out, "folds.seed0.tsv")), delimiter="\t"))
    assert len(folds) == 12 and {f["outer_fold"] for f in folds} == {"0", "1", "2"}
    assert [f for f in folds if f["family_id"] == "FX05"][0]["outer_fold"] == "0"
    syn = list(csv.DictReader(open(os.path.join(out, "synthetic_trios.seed1.tsv")), delimiter="\t"))
    assert len(syn) == 12 and all(r_["child_family"] != r_["parent_family"] for r_ in syn)
    fold_of = {f["family_id"]: f["outer_fold"] for f in csv.DictReader(open(os.path.join(out, "folds.seed1.tsv")), delimiter="\t")}
    assert all(fold_of[r_["child_family"]] == fold_of[r_["parent_family"]] for r_ in syn)
    # the synthetic manifest resolves the child's trio to the SURROGATE parents through the ordinary helpers
    m = read_manifest(os.path.join(out, "manifests", syn[0]["synthetic_id"] + ".manifest.tsv"))
    assert trio_of(m, syn[0]["child"]) == (syn[0]["surrogate_father"], syn[0]["surrogate_mother"])
    assert m[syn[0]["child"]]["family_id"] == syn[0]["synthetic_id"] and len(m) == 3


def test_swap_cli_blood_family_by_sample_prefix(tmp_path):
    """The cohort's blood family has REACH sample ids but an F0xxx family id: the sample-id prefix must find it."""
    man = str(tmp_path / "manifest.tsv"); write_manifest(man)
    rows = list(csv.DictReader(open(man), delimiter="\t"))
    for r_ in rows:
        if r_["family_id"] == "FX07":
            r_["sample_id"] = "REACH" + r_["sample_id"]
            for k in ("father_id", "mother_id"):
                if r_[k] != "0":
                    r_[k] = "REACH" + r_[k]
    with open(man, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    out = str(tmp_path / "train")
    env = dict(os.environ, PYTHONPATH=os.path.join(os.path.dirname(__file__), "..", "src"))
    env.pop("BLOOD_FAMILY_PREFIX", None); env.pop("BLOOD_SAMPLE_PREFIX", None)
    r = subprocess.run([sys.executable, "-m", "phase_dnm.cli", "swap", "--manifest", man, "--out-dir", out, "--n-folds", "3", "--seeds", "0"],
                       capture_output=True, text=True, env=env)
    assert r.returncode == 0, r.stderr
    assert "blood families: 1 (FX07)" in r.stderr + r.stdout
    folds = list(csv.DictReader(open(os.path.join(out, "folds.seed0.tsv")), delimiter="\t"))
    assert [f for f in folds if f["family_id"] == "FX07"][0]["outer_fold"] == "0"
