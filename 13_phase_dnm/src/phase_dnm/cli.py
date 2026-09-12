"""phase-dnm command line. One sub-command per sub-module; identifiers are always arguments."""
from __future__ import annotations

import argparse
import os
import sys

import yaml

from . import __version__
from .io.vcf import VcfReader, iter_trio, read_manifest, trio_of
from .phasing import orient as O

PLANNED = {
    "transmission": "M1b, week 2", "haplotag": "M1c, week 3", "phase-qc": "M1d, week 3",
    "candidates": "M2, week 4", "review": "M2, weeks 4-5", "features": "weeks 6-8",
    "spike": "weeks 6-8", "integrate": "M3, weeks 9-10", "train": "M4, weeks 11-13", "classify": "M4",
}


def load_thresholds(path: str | None) -> dict:
    if path is None:
        path = os.path.join(os.path.dirname(__file__), "..", "..", "config", "thresholds.yaml")
    with open(path) as fh:
        return yaml.safe_load(fh)


def cmd_orient(a: argparse.Namespace) -> int:
    thr = load_thresholds(a.thresholds)
    p = O.OrientParams(**{k: v for k, v in thr["orient"].items() if k in O.OrientParams.__dataclass_fields__})
    if a.manifest:
        man = read_manifest(a.manifest)
        father, mother = trio_of(man, a.child)
        sex = man[a.child]["sex"]
        if a.father and a.father != father or a.mother and a.mother != mother:
            sys.exit("parents given on the command line disagree with the manifest for %s" % a.child)
    else:
        if not (a.father and a.mother and a.sex):
            sys.exit("without --manifest, --father, --mother and --sex are required")
        father, mother, sex = a.father, a.mother, a.sex
    if a.joint_vcf:
        readers = [VcfReader(a.joint_vcf) for _ in range(3)]
        names = (a.child, father, mother)
    else:
        if not (a.child_vcf and a.father_vcf and a.mother_vcf):
            sys.exit("give --joint-vcf or all of --child-vcf --father-vcf --mother-vcf")
        readers = [VcfReader(a.child_vcf), VcfReader(a.father_vcf), VcfReader(a.mother_vcf)]
        names = tuple(None if len(r.samples) == 1 else s for r, s in zip(readers, (a.child, father, mother)))
    trio = iter_trio(readers[0], readers[1], readers[2], *names)
    blocks, stats = O.orient_child(trio, sex, p)
    os.makedirs(a.out_dir, exist_ok=True)
    stem = os.path.join(a.out_dir, a.child)
    O.write_orientation(blocks, stem + ".orientation.tsv")
    O.write_dissent(blocks, stem + ".orientation.dissent.tsv")
    summ = O.summarise(stats, blocks, p)
    summ["child"], summ["father"], summ["mother"], summ["sex"] = a.child, father, mother, sex
    summ["thresholds_version"] = thr.get("version")
    O.write_summary(summ, stem + ".orientation.summary.json")
    t = summ["total"]
    sys.stderr.write("orient %s: %d blocks (%d split at located switches) -> %d segments, %d oriented "
                     "(%.1f%% of block bp ambiguous); informative %d, Mendelian-inconsistent %d, dissenting votes %d\n"
                     % (a.child, summ["n_blocks"], summ["n_blocks_split"], summ["n_segments"],
                        summ["n_segments_oriented"], 100 * (summ["frac_bp_ambiguous"] or 0),
                        t.get("n_informative", 0), t.get("n_mendel_inconsistent", 0), t.get("n_dissent", 0)))
    return 0


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="phase-dnm", description=__doc__)
    ap.add_argument("--version", action="version", version=__version__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    o = sub.add_parser("orient", help="M1a: orient the child's phase blocks to parent of origin")
    o.add_argument("--child", required=True, help="child sample id")
    o.add_argument("--father"), o.add_argument("--mother"), o.add_argument("--sex", help="child sex 1/2/M/F")
    o.add_argument("--manifest", help="cohort manifest TSV; supplies parents and sex")
    o.add_argument("--child-vcf"), o.add_argument("--father-vcf"), o.add_argument("--mother-vcf")
    o.add_argument("--joint-vcf", help="alternatively one multi-sample phased VCF")
    o.add_argument("--out-dir", required=True)
    o.add_argument("--thresholds", help="config/thresholds.yaml (default: the module's)")
    o.set_defaults(func=cmd_orient)

    for name, when in PLANNED.items():
        s = sub.add_parser(name, help="planned (%s)" % when)
        s.set_defaults(func=lambda a, n=name, w=when: sys.exit("%s is not implemented yet (%s)" % (n, w)))
    return ap


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
