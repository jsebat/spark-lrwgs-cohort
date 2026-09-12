import csv
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))          # make_minitrio_vcf importable
from make_minitrio_vcf import Sim, simulate, CHILD, FATHER, MOTHER, FAMILY  # noqa: E402


@pytest.fixture(scope="session")
def minitrio(tmp_path_factory):
    d = tmp_path_factory.mktemp("minitrio")
    truth = simulate(Sim(seed=7), str(d))
    paths = {
        "dir": str(d),
        "child": os.path.join(d, "%s.%s.joint.GRCh38.small_variants.phased.vcf.gz" % (CHILD, FAMILY)),
        "father": os.path.join(d, "%s.%s.joint.GRCh38.small_variants.phased.vcf.gz" % (FATHER, FAMILY)),
        "mother": os.path.join(d, "%s.%s.joint.GRCh38.small_variants.phased.vcf.gz" % (MOTHER, FAMILY)),
        "joint": os.path.join(d, "%s.joint.GRCh38.small_variants.vcf.gz" % FAMILY),
        "manifest": os.path.join(d, "manifest.tsv"),
    }
    return {"paths": paths, "truth": truth, "ids": (CHILD, FATHER, MOTHER)}


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))
