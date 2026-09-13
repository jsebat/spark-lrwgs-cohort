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
 "train": "M4, weeks 11-13", "classify": "M4",
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


def cmd_candidates(a: argparse.Namespace) -> int:
    from . import candidates as C
    from .records import write_candidates
    man = read_manifest(a.manifest)
    father, mother = trio_of(man, a.child)
    family = man[a.child]["family_id"]
    thr = load_thresholds(a.thresholds)
    gen = {"snv_indel": lambda: C.snv_indel_candidates(a.vcf, family, a.child, father, mother),
           "sv": lambda: C.sv_candidates(a.vcf, family, a.child, father, mother),
           "tr": lambda: C.tr_candidates(a.vcf, family, a.child, father, mother,
                                          min_units=thr.get("candidates", {}).get("tr_min_units", 1),
                                          min_bp=thr.get("candidates", {}).get("tr_min_bp", 1))}[a.variant_class]
    recs = list(gen())
    if a.max_rows:
        # positives for M4 are thinned BEFORE the BAM pass (P23): reservoir per variant_class, seeded
        import random
        rng = random.Random(a.seed)
        kept = {}
        for r in recs:
            b = kept.setdefault(r.variant_class, [])
            b.append(r)
        recs = []
        for vc, b in sorted(kept.items()):
            recs += b if len(b) <= a.max_rows else rng.sample(b, a.max_rows)
        recs.sort(key=lambda r: (r.chrom, r.start))
    lists = []
    for spec in a.list or []:
        name, tier, path = spec.split(":", 2)
        lists.append((name, tier.upper(), path))
    if lists:
        recs = C.merge_lists(recs, lists)
    os.makedirs(a.out_dir, exist_ok=True)
    out = os.path.join(a.out_dir, "%s.%s.candidates.tsv" % (a.child, a.variant_class))
    n = write_candidates(recs, out)
    by_tier = {}
    for r in recs:
        by_tier[r.source_tier] = by_tier.get(r.source_tier, 0) + 1
    sys.stderr.write("candidates %s %s: %d rows (%s) -> %s\n" % (a.child, a.variant_class, n,
                     ", ".join("%s=%d" % kv for kv in sorted(by_tier.items())), out))
    return 0


def _m2_params(thr: dict):
    """HapParams / ClassParams from thresholds.yaml (shared by review and reclassify so both layers agree)."""
    from .evidence import hapmatrix as H
    hm, pc = thr.get("hapmatrix", {}), thr.get("phase_class", {})
    hp = H.HapParams(k=tuple(hm.get("k", [3, 5])), error_reads=hm.get("error_reads", 1), error_frac=hm.get("error_frac", 0.05),
                     min_mapq=hm.get("min_mapq", 20), amb_flag_frac=hm.get("amb_flag_frac", 0.3))
    cp = H.ClassParams(**{k: v for k, v in pc.items() if k in H.ClassParams.__dataclass_fields__}, working_k=max(hp.k))
    return hp, cp


def cmd_reclassify(a: argparse.Namespace) -> int:
    """Re-run the rule layer (features + transmission + P8 class) from an existing evidence table, no BAM access.
    Input is the immutable review output; the output replaces the working evidence table downstream."""
    from .evidence.review import reclassify_table
    thr = load_thresholds(a.thresholds)
    hp, cp = _m2_params(thr)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    counts = reclassify_table(a.evidence, a.out, hp, cp, thresholds_version=thr.get("version"))
    sys.stderr.write("reclassify %s -> %s: %s\n" % (os.path.basename(a.evidence), a.out,
                     " ".join("%s=%d" % kv for kv in sorted(counts.items()))))
    return 0


def cmd_spike(a: argparse.Namespace) -> int:
    """Spike-in harness: plan | apply | evaluate (sim/spike.py)."""
    from .sim import spike as SP
    log = lambda m: sys.stderr.write(m + "\n")
    if a.step == "evaluate":
        import csv
        plan = SP.read_plan(a.plan)
        with open(a.evidence, newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        per_site, summary = SP.evaluate(plan, rows, k=a.k)
        SP.write_rows(per_site, a.out.replace(".tsv", ".per_site.tsv"))
        SP.write_rows(summary, a.out)
        for s_ in summary:
            log("spike %-5s %-3s planted=%3d reviewed=%3d observable=%3d class_ok(obs)=%s poo_ok(obs)=%s post_expected=%s called=%s" % (
                s_["variant_class"], s_["scenario"], s_["n_planted"], s_["n_reviewed"], s_["n_observable"],
                s_["class_ok_observable"], s_["poo_ok_observable"], s_["mean_post_expected"], s_["called_classes"]))
        return 0
    import pysam
    from .evidence import hapmatrix as H
    thr = load_thresholds(a.thresholds)
    hm = thr.get("hapmatrix", {})
    qc = SP.SiteQC(k=max(hm.get("k", [3, 5])), min_mapq=hm.get("min_mapq", 20), pad=a.pad, spacing=a.spacing)
    op = lambda p_, i: pysam.AlignmentFile(p_, "rb", index_filename=i) if i else pysam.AlignmentFile(p_, "rb")
    bams = {"C": op(a.child_bam, a.child_bai), "F": op(a.father_bam, a.father_bai), "M": op(a.mother_bam, a.mother_bai)}
    tr_bams = {"C": op(a.child_tr_bam, a.child_tr_bai), "F": op(a.father_tr_bam, a.father_tr_bai), "M": op(a.mother_tr_bam, a.mother_tr_bai)} if a.child_tr_bam else {}
    os.makedirs(a.out_dir, exist_ok=True)
    if a.step == "plan":
        from .features.extract import BedMask
        labels = H.LabelTables.load(a.orientation, a.transmission, a.changepoints)
        mask = BedMask(a.mask) if a.mask else None
        plan, cands, stats = SP.plan_sites(bams, tr_bams, labels, SP.parse_regions(a.regions), a.n_per_scenario, a.seed,
                                           a.family, a.child, qc, trgt_vcf=a.trgt_vcf, mask=mask, log=log, trgt_vcf_index=a.trgt_vcf_index)
        SP.write_plan(plan, os.path.join(a.out_dir, "plan.tsv"))
        from .records import write_candidates
        for cls_group, classes in (("snv_indel", ("SNV", "INDEL")), ("sv", ("SV",)), ("tr", ("TR",))):
            write_candidates([c for c in cands if c.variant_class in classes], os.path.join(a.out_dir, "%s.candidates.tsv" % cls_group))
        with open(os.path.join(a.out_dir, "plan.stats.json"), "w") as fh:
            json.dump(stats, fh, indent=1, sort_keys=True)
        log("spike plan: %d sites planted of %d planned; %s" % (len(plan), stats.get("planned", 0), " ".join("%s=%d" % kv for kv in sorted(stats.items()))))
    elif a.step == "apply":
        plan = SP.read_plan(a.plan)
        st1 = SP.apply_plan(plan, bams, a.out_dir, "spiked", qc, a.seed, tr=False, log=log)
        st2 = SP.apply_plan(plan, tr_bams, a.out_dir, "spiked.tr", qc, a.seed, tr=True, log=log) if tr_bams else {}
        with open(os.path.join(a.out_dir, "apply.stats.json"), "w") as fh:
            json.dump({"genome": st1, "tr": st2}, fh, indent=1, sort_keys=True)
    for b in list(bams.values()) + list(tr_bams.values()):
        b.close()
    return 0


def _final_params(thr: dict):
    from .integrate import FinalParams
    f = thr.get("final", {})
    return FinalParams(tau=dict(f.get("tau", {}) or {}), tau_rescue=dict(f.get("tau_rescue", {}) or {}),
                       phase_only_min_score=f.get("phase_only_min_score", 0.9), require_hap_obs=f.get("require_hap_obs", 6))


def cmd_integrate(a: argparse.Namespace) -> int:
    """M3: final unfiltered table for one child and class group (+ sites VCF, + Parquet when pyarrow exists)."""
    from . import integrate as I
    from .io import vcfinfo as V
    thr = load_thresholds(a.thresholds)
    p = _final_params(thr)
    rf = I.load_rf_probs(a.rf_probs)
    summ = I.integrate_table(a.evidence, a.features, a.out, p, a.class_group, rf_probs=rf)
    summ["thresholds_version"] = thr.get("version")
    summ["call_mode"] = "rf+phase" if rf else "phase_only"
    pq = I.write_parquet(a.out)
    summ["parquet"] = pq
    if a.vcf_out:
        rows = V.rows_from_final(a.out)
        summ["vcf_records"] = V.write_sites_vcf(rows, a.vcf_out)
    with open(a.out.replace(".tsv", ".summary.json"), "w") as fh:
        json.dump(summ, fh, indent=1, sort_keys=True)
    c = summ["counts"]
    sys.stderr.write("integrate %s [%s]: %d rows, YES=%d NO=%d mosaic=%d; %s; %d columns (%d features)%s\n" % (
        os.path.basename(a.evidence), summ["call_mode"], c["rows"], c["YES"], c["NO"], c["mosaic"],
        " ".join("%s=%d" % kv for kv in sorted(summ["reasons"].items())), summ["columns"], summ["features"],
        "" if pq else "; no pyarrow: TSV only"))
    return 0


def cmd_concordance(a: argparse.Namespace) -> int:
    """M3/P18: concordance of the final tables (one class group, all families) with the original pipeline's de novo set."""
    import csv, glob
    from . import concordance as C
    rows = []
    for f in sorted(glob.glob(os.path.join(a.final_dir, "*.%s.dnm.tsv" % a.class_group))):
        with open(f, newline="") as fh:
            rows.extend(csv.DictReader(fh, delimiter="\t"))
    per_row, per_proband, summ = C.concordance(rows, a.class_group, a.baselines, sv_bp_tol=a.sv_bp_tol, sv_recip=a.sv_recip)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    def dump(path, rs):
        if not rs:
            open(path, "w").close(); return
        cols = []
        for r in rs:
            for k in r:
                if k not in cols:
                    cols.append(k)
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
            w.writeheader()
            for r in rs:
                w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in cols})
    dump(a.out, [d for d in per_row if d["status"] != "neither"] if a.drop_neither else per_row)
    dump(a.out.replace(".tsv", ".per_proband.tsv"), per_proband)
    with open(a.out.replace(".tsv", ".summary.json"), "w") as fh:
        json.dump(summ, fh, indent=1, sort_keys=True)
    sys.stderr.write("concordance %s: rows=%d probands=%d concordant_YES=%d original_only=%d module_only=%d original_unseen=%d; "
                     "per-proband original median %s -> module YES median %s; module YES paternal fraction %s\n" % (
                         a.class_group, summ["n_rows"], summ["n_probands"], summ["concordant_YES"], summ["original_only"], summ["module_only"],
                         summ["original_unseen_as_candidate"], summ["per_proband_original_median"], summ["per_proband_module_YES_median"],
                         summ["module_YES_paternal_fraction"]))
    return 0


def cmd_swap(a: argparse.Namespace) -> int:
    """M4: swap-closed family folds + within-fold pedigree swaps for each seed; per-synthetic-trio manifests (P12, P23)."""
    import csv
    from .train import folds as FO
    from .train import swap as SW
    with open(a.manifest, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    cols = list(rows[0].keys())
    by_id = {r["sample_id"]: r for r in rows}
    sources = {r["family_id"]: "blood" for r in rows if a.blood_family_prefix and r["family_id"].startswith(a.blood_family_prefix)}
    fams = FO.families_from_manifest(rows, sources)
    parents = {}
    for f in fams:
        c0 = by_id[f.children[0]]
        parents[f.family_id] = (c0["father_id"], c0["mother_id"])
    os.makedirs(os.path.join(a.out_dir, "manifests"), exist_ok=True)
    n_syn = 0
    for seed in [int(x) for x in a.seeds.split(",")]:
        assign = FO.outer_folds(fams, n_folds=a.n_folds, seed=seed)
        table = SW.pairings(fams, assign, parents, seed=seed)
        SW.check_closed(table, assign)
        with open(os.path.join(a.out_dir, "folds.seed%d.tsv" % seed), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["family_id", "outer_fold", "seed"], delimiter="\t", lineterminator="\n")
            w.writeheader()
            for r in FO.fold_table(assign, seed):
                w.writerow(r)
        with open(os.path.join(a.out_dir, "synthetic_trios.seed%d.tsv" % seed), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(table[0].keys()), delimiter="\t", lineterminator="\n")
            w.writeheader()
            for r in table:
                w.writerow(r)
        for r in table:
            child = dict(by_id[r["child"]]); child.update(family_id=r["synthetic_id"], father_id=r["surrogate_father"], mother_id=r["surrogate_mother"])
            fa = dict(by_id[r["surrogate_father"]]); fa.update(family_id=r["synthetic_id"])
            mo = dict(by_id[r["surrogate_mother"]]); mo.update(family_id=r["synthetic_id"])
            with open(os.path.join(a.out_dir, "manifests", "%s.manifest.tsv" % r["synthetic_id"]), "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
                w.writeheader()
                for row in (child, fa, mo):
                    w.writerow(row)
            n_syn += 1
        sizes = [sum(1 for v in assign.values() if v == k) for k in range(a.n_folds)]
        sys.stderr.write("swap seed %d: %d families in %d folds (sizes %s), %d synthetic trios\n" % (seed, len(fams), a.n_folds, sizes, len(table)))
    sys.stderr.write("swap: %d synthetic-trio manifests -> %s\n" % (n_syn, os.path.join(a.out_dir, "manifests")))
    return 0


def cmd_review(a: argparse.Namespace) -> int:
    from .evidence import hapmatrix as H
    from .evidence.readers import TrioBams
    from .evidence.review import review_child
    thr = load_thresholds(a.thresholds)
    hp, cp = _m2_params(thr)
    labels = H.LabelTables.load(a.orientation, a.transmission, a.changepoints)
    bams = TrioBams({"C": (a.child_bam, a.child_bai), "F": (a.father_bam, a.father_bai), "M": (a.mother_bam, a.mother_bai)},
                    {"C": (a.child_tr_bam, a.child_tr_bai), "F": (a.father_tr_bam, a.father_tr_bai), "M": (a.mother_tr_bam, a.mother_tr_bai)}
                    if a.child_tr_bam else None)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    counts = review_child(a.candidates, a.out, bams, labels, hp, cp, salt=a.salt or a.out, supporting_json=a.sv_supporting_reads,
                          reads_jsonl_gz=a.reads_out, max_per_class=a.max_per_class, thresholds_version=thr.get("version"),
                          class_filter=a.only_class)
    bams.close()
    sys.stderr.write("review %s -> %s: %s\n" % (os.path.basename(a.candidates), a.out,
                     " ".join("%s=%d" % kv for kv in sorted(counts.items()))))
    return 0


def cmd_features(a: argparse.Namespace) -> int:
    from .features import extract as X
    from .features.registry import Registry
    reg = Registry(a.registry)
    man = read_manifest(a.manifest)
    sex = man[a.child]["sex"]
    mask = X.BedMask(a.mask) if a.mask else None
    seqctx = X.SeqContext(a.reference) if a.reference else None
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    summ = X.extract_child(a.evidence, a.candidates, reg, sex, a.out, a.rf_out, mask=mask, seqctx=seqctx)
    summ["registry"] = reg.manifest()
    with open(a.out.replace(".tsv", ".summary.json"), "w") as fh:
        json.dump(summ, fh, indent=1, sort_keys=True)
    sys.stderr.write("features %s: %d rows, %d/%d applicable features produced; never produced: %s\n"
                     % (os.path.basename(a.evidence), summ["rows"], summ.get("features_produced", 0), summ.get("features_applicable", 0),
                        ", ".join(summ.get("never_produced", [])[:12]) + (" ..." if len(summ.get("never_produced", [])) > 12 else "")))
    return 0


def cmd_likelihood(a: argparse.Namespace) -> int:
    """Post-step over an evidence table: adds lik_post_* / phase_score columns (P9)."""
    import csv
    from .classify import likelihood as L
    thr = load_thresholds(a.thresholds).get("likelihood", {})
    prm = L.LikParams(eps=thr.get("eps", 0.01), delta=thr.get("delta", 0.03), grid=thr.get("grid", 40), prior=thr.get("prior"))
    if a.eps is not None:
        prm.eps = a.eps
    n = 0
    agree = 0
    with open(a.evidence, newline="") as fh, open(a.out, "w", newline="") as out:
        rd = csv.DictReader(fh, delimiter="\t")
        extra = ["lik_post_" + h for h in L.HYPS] + ["lik_best_alternative", "lik_log10lr_germline", "phase_score"]
        w = csv.DictWriter(out, fieldnames=rd.fieldnames + [c for c in extra if c not in rd.fieldnames], delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rd:
            r.update({k: ("" if v is None else v) for k, v in L.score_evidence_row(r, prm).items()})
            n += 1
            if r.get("phase_class") in ("germline_DNM_phased",) and r.get("phase_score") not in ("", None) and float(r["phase_score"]) > 0.5:
                agree += 1
            w.writerow(r)
    sys.stderr.write("likelihood: %d rows scored -> %s (eps=%.3f delta=%.3f)\n" % (n, a.out, prm.eps, prm.delta))
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

    cd = sub.add_parser("candidates", help="M2: unfiltered candidate records for one child and one class (P5)")
    cd.add_argument("--child", required=True), cd.add_argument("--manifest", required=True)
    cd.add_argument("--class", dest="variant_class", required=True, choices=["snv_indel", "sv", "tr"])
    cd.add_argument("--vcf", required=True, help="family joint VCF of that class (GLnexus / sawfish / TRGT)")
    cd.add_argument("--list", action="append", help="existing list to cross-reference: NAME:TIER:path.tsv (repeatable)")
    cd.add_argument("--out-dir", required=True), cd.add_argument("--thresholds")
    cd.add_argument("--max-rows", type=int, help="M4 synthetic trios: keep at most this many rows per variant class (seeded reservoir)")
    cd.add_argument("--seed", type=int, default=0)
    cd.set_defaults(func=cmd_candidates)

    rv = sub.add_parser("review", help="M2: six-haplotype evidence + phase class for every candidate of one child (pysam)")
    rv.add_argument("--candidates", required=True, help="<child>.<class>.candidates.tsv")
    rv.add_argument("--out", required=True, help="<child>.<class>.evidence.tsv")
    rv.add_argument("--reads-out", help="<child>.<class>.reads.jsonl.gz (per-read observations, hashed ids)")
    rv.add_argument("--orientation"), rv.add_argument("--transmission"), rv.add_argument("--changepoints")
    for role in ("child", "father", "mother"):
        rv.add_argument("--%s-bam" % role, required=True), rv.add_argument("--%s-bai" % role)
        rv.add_argument("--%s-tr-bam" % role), rv.add_argument("--%s-tr-bai" % role)
    rv.add_argument("--sv-supporting-reads", help="sawfish supporting_reads.json.gz of the family")
    rv.add_argument("--max-per-class", type=int, help="cap rows per class (smoke tests)")
    rv.add_argument("--only-class", action="append", choices=["SNV", "INDEL", "SV", "TR"])
    rv.add_argument("--salt", help="salt for read-id hashing (default: the output path)")
    rv.add_argument("--thresholds")
    rv.set_defaults(func=cmd_review)

    fe = sub.add_parser("features", help="registry-driven feature vectors from an evidence table (+ rf_safe classifier matrix)")
    fe.add_argument("--evidence", required=True), fe.add_argument("--candidates", required=True)
    fe.add_argument("--child", required=True), fe.add_argument("--manifest", required=True)
    fe.add_argument("--out", required=True, help="<child>.<class>.features.tsv")
    fe.add_argument("--rf-out", help="<child>.<class>.features.rf.tsv (rf_safe columns only; refused if unsafe)")
    fe.add_argument("--mask", action="append", help="BED file(s) of the lab region mask (flag, never a filter)")
    fe.add_argument("--reference", help="reference FASTA for sequence-context features (pysam)")
    fe.add_argument("--registry", help="config/features.yaml (default: the module's)")
    fe.set_defaults(func=cmd_features)

    sk = sub.add_parser("spike", help="spike-in harness: plan sites, edit haplotagged reads into slice BAMs, evaluate recovery")
    sk.add_argument("step", choices=["plan", "apply", "evaluate"])
    sk.add_argument("--out-dir", help="plan/apply: directory for plan.tsv, candidates, slice BAMs")
    sk.add_argument("--family"), sk.add_argument("--child")
    sk.add_argument("--regions", help="plan: chrom:start-end[,...] to sample sites from")
    sk.add_argument("--n-per-scenario", type=int, default=12, help="plan: sites per class x subtype x scenario")
    sk.add_argument("--seed", type=int, default=13)
    sk.add_argument("--pad", type=int, default=2500), sk.add_argument("--spacing", type=int, default=60000)
    sk.add_argument("--mask", action="append", help="plan: BED(s) to avoid (spike-ins are planted OUTSIDE the mask)")
    sk.add_argument("--trgt-vcf", help="plan: TRGT VCF (tabix-indexed) for TR loci and motifs"), sk.add_argument("--trgt-vcf-index")
    sk.add_argument("--orientation"), sk.add_argument("--transmission"), sk.add_argument("--changepoints")
    for r in ("child", "father", "mother"):
        sk.add_argument("--%s-bam" % r), sk.add_argument("--%s-bai" % r), sk.add_argument("--%s-tr-bam" % r), sk.add_argument("--%s-tr-bai" % r)
    sk.add_argument("--plan", help="apply/evaluate: plan.tsv")
    sk.add_argument("--evidence", help="evaluate: the reviewed (+likelihood) evidence table of the spiked candidates")
    sk.add_argument("--out", help="evaluate: summary TSV (per-site table written next to it)")
    sk.add_argument("--k", type=int, default=5), sk.add_argument("--thresholds")
    sk.set_defaults(func=cmd_spike)
    ig = sub.add_parser("integrate", help="M3: final unfiltered table per child and class group (P15 decision, class columns, feature vector)")
    ig.add_argument("--evidence", required=True, help="<child>.<class>.evidence.lik.tsv"), ig.add_argument("--features", help="<child>.<class>.features.tsv")
    ig.add_argument("--class-group", required=True, choices=["snv_indel", "sv", "tr"])
    ig.add_argument("--out", required=True, help="final/<FAMILY>.<child>.<class>.dnm.tsv"), ig.add_argument("--vcf-out")
    ig.add_argument("--rf-probs", help="M4 output: TSV with variant_id, rf_prob (absent -> provisional phase_only mode)")
    ig.add_argument("--thresholds")
    ig.set_defaults(func=cmd_integrate)
    cc = sub.add_parser("concordance", help="M3/P18: final calls vs the original pipeline's de novo set, one class group")
    cc.add_argument("--final-dir", required=True), cc.add_argument("--baselines", required=True)
    cc.add_argument("--class-group", required=True, choices=["snv_indel", "sv", "tr"]), cc.add_argument("--out", required=True)
    cc.add_argument("--sv-bp-tol", type=int, default=500), cc.add_argument("--sv-recip", type=float, default=0.5)
    cc.add_argument("--drop-neither", action="store_true", help="omit rows that are neither called nor in the baseline from the per-row table")
    cc.set_defaults(func=cmd_concordance)
    sw = sub.add_parser("swap", help="M4: swap-closed family folds and within-fold pedigree swaps (synthetic trios) per seed")
    sw.add_argument("--manifest", required=True), sw.add_argument("--out-dir", required=True)
    sw.add_argument("--n-folds", type=int, default=5), sw.add_argument("--seeds", default="0,1,2,3,4")
    sw.add_argument("--blood-family-prefix", default="REACH", help="family-id prefix of the blood-derived family (kept with company in fold 0)")
    sw.set_defaults(func=cmd_swap)
    rc = sub.add_parser("reclassify", help="re-run the P8 rule layer from an existing evidence table (no BAMs)")
    rc.add_argument("--evidence", required=True, help="the review output (immutable)"), rc.add_argument("--out", required=True)
    rc.add_argument("--thresholds")
    rc.set_defaults(func=cmd_reclassify)
    lk = sub.add_parser("likelihood", help="P9: posterior over germline/inherited/mosaic/artefact from an evidence table")
    lk.add_argument("--evidence", required=True), lk.add_argument("--out", required=True)
    lk.add_argument("--eps", type=float, help="override the per-read false-alt rate")
    lk.add_argument("--thresholds")
    lk.set_defaults(func=cmd_likelihood)

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
