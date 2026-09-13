"""phase-dnm command line. One sub-command per sub-module; identifiers are always arguments."""
from __future__ import annotations

import argparse
import json
import os
import sys

import yaml

from . import __version__
from .io.vcf import VcfReader, iter_trio, read_manifest, trio_of
from .phasing import orient as O
from .phasing import transmission as T

PLANNED = {
    "haplotag": "M1c export to BAM (optional, IGV); labels are the orientation/transmission tables",
    "candidates": "M2, week 4", "review": "M2, weeks 4-5", "features": "weeks 6-8",
    "spike": "weeks 6-8", "integrate": "M3, weeks 9-10", "train": "M4, weeks 11-13", "classify": "M4",
}


def load_thresholds(path: str | None) -> dict:
    if path is None:
        path = os.path.join(os.path.dirname(__file__), "..", "..", "config", "thresholds.yaml")
    with open(path) as fh:
        return yaml.safe_load(fh)


def _resolve_trio(a: argparse.Namespace):
    """(father, mother, sex, trio_iterator) from the manifest or explicit flags; shared by orient and transmission."""
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
    return father, mother, sex, iter_trio(readers[0], readers[1], readers[2], *names)


def cmd_orient(a: argparse.Namespace) -> int:
    thr = load_thresholds(a.thresholds)
    p = O.OrientParams(**{k: v for k, v in thr["orient"].items() if k in O.OrientParams.__dataclass_fields__})
    father, mother, sex, trio = _resolve_trio(a)
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


def cmd_transmission(a: argparse.Namespace) -> int:
    thr = load_thresholds(a.thresholds)
    p = T.TransmissionParams(**{k: v for k, v in thr["transmission"].items() if k in T.TransmissionParams.__dataclass_fields__})
    father, mother, sex, trio = _resolve_trio(a)
    stem = os.path.join(a.out_dir, a.child)
    orient_path = a.orientation or (stem + ".orientation.tsv")
    if not os.path.exists(orient_path):
        sys.exit("orientation table not found (%s); run `phase-dnm orient` first" % orient_path)
    orientation = T.Orientation.from_tsv(orient_path)
    segments, changes, stats = T.build_transmission(trio, orientation, sex, p)
    os.makedirs(a.out_dir, exist_ok=True)
    T.write_segments(segments, stem + ".transmission.tsv")
    T.write_changes(changes, stem + ".changepoints.tsv")
    summ = T.summarise(stats, p)
    summ.update(child=a.child, father=father, mother=mother, sex=sex, thresholds_version=thr.get("version"))
    T.write_summary(summ, stem + ".transmission.summary.json")
    for parent in ("F", "M"):
        pp = summ["per_parent"].get(parent, {})
        sys.stderr.write("transmission %s parent %s: %d blocks, %d resolved segments (%.1f%% of parent hets), "
                         "%d change points (crossover-or-switch candidates), informative %d, Mendelian-inconsistent %d\n"
                         % (a.child, parent, pp.get("n_blocks", 0), pp.get("n_segments_resolved", 0),
                            100 * (pp.get("frac_het_resolved") or 0), pp.get("n_change_points", 0),
                            pp.get("n_informative", 0), pp.get("n_mendel_inconsistent", 0)))
    return 0


def cmd_hapdepth(a: argparse.Namespace) -> int:
    from .phasing.hapdepth import hapdepth  # pysam, lazy
    region = None
    if a.region:
        chrom, span = a.region.split(":")
        s, e = span.replace(",", "").split("-")
        region = (chrom, int(s), int(e))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    summ = hapdepth(a.bam, a.bai, a.out, a.summary, bin_size=a.bin_size, min_mapq=a.min_mapq, region=region)
    tot = sum(v["n_reads"] for v in summ.values())
    sys.stderr.write("hapdepth: %d contigs/regions, %d primary reads, written %s\n" % (len(summ), tot, a.out))
    return 0


def cmd_xo_reads(a: argparse.Namespace) -> int:
    from .phasing.xo_reads import classify  # pysam, lazy
    thr = load_thresholds(a.thresholds)
    xo = thr.get("xo_reads", {})
    parent_bam = {}
    parent_vcf = {}
    if a.father_bam:
        parent_bam["F"] = (a.father_bam, a.father_bai)
        parent_vcf["F"] = (a.father_vcf, a.father_vcf_index, a.father)
    if a.mother_bam:
        parent_bam["M"] = (a.mother_bam, a.mother_bai)
        parent_vcf["M"] = (a.mother_vcf, a.mother_vcf_index, a.mother)
    if not parent_bam:
        sys.exit("give at least one parent's --*-bam and --*-vcf")
    counts = classify(a.changepoints, a.out, parent_bam, parent_vcf,
                      min_spanning=xo.get("min_spanning_reads", 3), min_mapq=xo.get("min_mapq", 20),
                      crossover_max_disc=xo.get("crossover_max_disc", 0.2), switch_min_disc=xo.get("switch_min_disc", 0.8),
                      child_orientation_tsv=a.child_orientation)
    sys.stderr.write("xo-reads: %s -> %s\n" % (a.changepoints, " ".join("%s=%d" % kv for kv in sorted(counts.items()))))
    return 0


def cmd_phase_qc(a: argparse.Namespace) -> int:
    from .phasing import qc as Q
    thr = load_thresholds(a.thresholds)
    gates = thr.get("qc_gates", {})
    man = read_manifest(a.manifest)
    families = [a.family] if a.family else sorted({r["family_id"] for r in man.values()})
    rows = []
    for fam in families:
        fam_dir = os.path.join(a.phase_dir, fam)
        if not os.path.isdir(fam_dir):
            continue
        for sid, r in man.items():
            if r["family_id"] != fam or r["role"] != "offspring":
                continue
            try:
                father, mother = trio_of(man, sid)
            except ValueError:
                continue                                        # duo (P19)
            q = Q.child_qc(fam_dir, sid, father, mother, a.hapdepth_dir, gates, thr.get("version"))
            with open(os.path.join(fam_dir, sid + ".phase_qc.json"), "w") as fh:
                json.dump(q, fh, indent=1, sort_keys=True)
            rows.append(Q.cohort_row(fam, q))
    if not rows:
        sys.exit("no children with M1 outputs under %s" % a.phase_dir)
    out = a.out or os.path.join(a.phase_dir, "cohort_phase_qc.tsv")
    Q.write_cohort_table(rows, out)
    summ = Q.cohort_summary(rows)
    with open(out.replace(".tsv", ".summary.json"), "w") as fh:
        json.dump(summ, fh, indent=1, sort_keys=True)
    sys.stderr.write("phase-qc: %d children, %d PASS; flags seen: %s; wrote %s\n"
                     % (summ["n_children"], summ["n_pass"], ", ".join(summ["flags"]) or "none", out))
    return 0


def _trio_args(sp: argparse.ArgumentParser):
    sp.add_argument("--child", required=True, help="child sample id")
    sp.add_argument("--father"), sp.add_argument("--mother"), sp.add_argument("--sex", help="child sex 1/2/M/F")
    sp.add_argument("--manifest", help="cohort manifest TSV; supplies parents and sex")
    sp.add_argument("--child-vcf"), sp.add_argument("--father-vcf"), sp.add_argument("--mother-vcf")
    sp.add_argument("--joint-vcf", help="alternatively one multi-sample phased VCF")
    sp.add_argument("--out-dir", required=True)
    sp.add_argument("--thresholds", help="config/thresholds.yaml (default: the module's)")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="phase-dnm", description=__doc__)
    ap.add_argument("--version", action="version", version=__version__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    o = sub.add_parser("orient", help="M1a: orient the child's phase blocks to parent of origin")
    _trio_args(o)
    o.set_defaults(func=cmd_orient)

    t = sub.add_parser("transmission", help="M1b: transmitted/untransmitted map per parent + change points")
    _trio_args(t)
    t.add_argument("--orientation", help="<child>.orientation.tsv (default: <out-dir>/<child>.orientation.tsv)")
    t.set_defaults(func=cmd_transmission)

    h = sub.add_parser("hapdepth", help="M1: per-haplotype depth in fixed bins from one haplotagged BAM (pysam)")
    h.add_argument("--bam", required=True), h.add_argument("--bai", help="index if not beside the BAM")
    h.add_argument("--out", required=True, help="<sample>.hapdepth.tsv.gz")
    h.add_argument("--summary", help="<sample>.hapdepth.summary.json")
    h.add_argument("--bin-size", type=int, default=1000), h.add_argument("--min-mapq", type=int, default=20)
    h.add_argument("--region", help="chrom:start-end (smoke tests); default all primary contigs")
    h.set_defaults(func=cmd_hapdepth)

    x = sub.add_parser("xo-reads", help="M1b2: classify change points as CROSSOVER / SWITCH_ERROR from parent reads (pysam)")
    x.add_argument("--changepoints", required=True, help="<child>.changepoints.tsv from `transmission`")
    x.add_argument("--out", required=True, help="<child>.changepoints.resolved.tsv")
    x.add_argument("--father"), x.add_argument("--mother")
    x.add_argument("--father-bam"), x.add_argument("--father-bai"), x.add_argument("--father-vcf"), x.add_argument("--father-vcf-index")
    x.add_argument("--mother-bam"), x.add_argument("--mother-bai"), x.add_argument("--mother-vcf"), x.add_argument("--mother-vcf-index")
    x.add_argument("--child-orientation", help="<child>.orientation.tsv: marks intervals containing a located child switch")
    x.add_argument("--thresholds")
    x.set_defaults(func=cmd_xo_reads)

    qc = sub.add_parser("phase-qc", help="M1d: per-child phase_qc.json and the cohort QC table with gates")
    qc.add_argument("--manifest", required=True)
    qc.add_argument("--phase-dir", required=True, help="$PHASE_DIR holding <FAMILY>/ sub-directories")
    qc.add_argument("--hapdepth-dir", help="$PHASE_DIR/hapdepth (optional)")
    qc.add_argument("--family", help="one family only (default: every family present under --phase-dir)")
    qc.add_argument("--out", help="cohort table path (default: <phase-dir>/cohort_phase_qc.tsv)")
    qc.add_argument("--thresholds")
    qc.set_defaults(func=cmd_phase_qc)

    for name, when in PLANNED.items():
        s = sub.add_parser(name, help="planned (%s)" % when)
        s.set_defaults(func=lambda a, n=name, w=when: sys.exit("%s is not implemented yet (%s)" % (n, w)))
    return ap


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
