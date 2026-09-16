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
 "classify": "M4",
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
    if a.step == "genotype":
        # P29: fill the caller-derived block of the planted candidates. Without this the planted-truth arm scores the
        # classifier on a matrix whose genotype and annotation columns are empty -- 42 of 72 for SNV/indel.
        from .sim import genotype as GT
        bams = {"child": a.child_bam, "father": a.father_bam, "mother": a.mother_bam}
        bams = {k: v for k, v in bams.items() if v}
        led = os.path.join(a.out_dir, "apply.ledger.tsv")
        n_any = 0
        for cls_group in ("snv_indel", "sv", "tr"):
            cand = os.path.join(a.out_dir, "%s.candidates.tsv" % cls_group)
            if not os.path.exists(cand) or os.path.getsize(cand) == 0:
                continue
            st = GT.fill_candidates(cand, cand + ".tmp", bams, a.reference, led, bcftools=a.bcftools, log=log)
            if st.get("rows"):
                os.replace(cand + ".tmp", cand)
                GT.write_private_annot(cand, os.path.join(a.out_dir, "%s.annot.tsv" % cls_group), log=log)
                n_any += st["rows"]
            elif os.path.exists(cand + ".tmp"):
                os.remove(cand + ".tmp")
        log("spike genotype: %d planted candidate rows now carry the caller and annotation blocks" % n_any)
        return 0
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
        led = os.path.join(a.out_dir, "apply.ledger.tsv")
        if os.path.exists(led):
            os.remove(led)                                          # the ledger is appended to, once per class group
        st1 = SP.apply_plan(plan, bams, a.out_dir, "spiked", qc, a.seed, tr=False, log=log, ledger_path=led)
        st2 = SP.apply_plan(plan, tr_bams, a.out_dir, "spiked.tr", qc, a.seed, tr=True, log=log, ledger_path=led) if tr_bams else {}
        with open(os.path.join(a.out_dir, "apply.stats.json"), "w") as fh:
            json.dump({"genome": st1, "tr": st2}, fh, indent=1, sort_keys=True)
    for b in list(bams.values()) + list(tr_bams.values()):
        b.close()
    return 0


def cmd_audit(a: argparse.Namespace) -> int:
    """P30: compare every declaration against the artefact that is supposed to satisfy it. Non-zero on any FAIL."""
    from . import audit as AU
    from .features.registry import Registry
    reg = Registry(a.registry)
    groups = tuple(a.class_group) if a.class_group else ("snv_indel", "sv", "tr")
    return AU.run(a.evidence_dir, reg, harness_dir=a.harness_dir, spike=a.spike, class_groups=groups,
                  log=lambda m: sys.stderr.write(m + chr(10)))


def _final_params(thr: dict):
    from .integrate import FinalParams
    f = thr.get("final", {})
    from .integrate import DEFAULT_RULES
    rules = dict(DEFAULT_RULES); rules.update(f.get("rules") or {})
    return FinalParams(tau=dict(f.get("tau", {}) or {}), tau_rescue=dict(f.get("tau_rescue", {}) or {}), tau_tier2=dict(f.get("tau_tier2", {}) or {}),
                       phase_only_min_score=f.get("phase_only_min_score", 0.9), require_hap_obs=f.get("require_hap_obs", 6),
                       tr_rescue_min_units=float(f.get("tr_rescue_min_units", 3)), rules=rules, apply_rules=bool(f.get("apply_rules", True)),
                       sv_depth_rule_tier1=bool(f.get("sv_depth_rule_tier1", True)),
                       sv_depth_min_rule_score=float(f.get("sv_depth_min_rule_score", 6)))


def cmd_integrate(a: argparse.Namespace) -> int:
    """M3: final unfiltered table for one child and class group (+ sites VCF, + Parquet when pyarrow exists)."""
    from . import integrate as I
    from .io import vcfinfo as V
    thr = load_thresholds(a.thresholds)
    p = _final_params(thr)
    score_col = a.score_column
    if a.tau_json and os.path.exists(a.tau_json):
        with open(a.tau_json) as fh:
            tj = json.load(fh)
        if score_col == "auto":
            score_col = tj.get("score_column", "rf_prob")
        fq = thr.get("final", {})
        for vc in ({"snv_indel": ("SNV", "INDEL"), "sv": ("SV",), "tr": ("TR",)}[a.class_group]):
            if score_col == "rf_q":
                # thresholds.yaml final.tau_q / tau_q_tier2 win over the tau json's uniform value; a key per variant class (SNV, INDEL)
                # overrides the class-group key (snv_indel) - indels get their own operating point (JS 2026-09-14)
                def _pick(d, default):
                    d = d or {}
                    v = d.get(vc, d.get(a.class_group, default))
                    return None if v is None else float(v)
                p.tau[vc] = _pick(fq.get("tau_q"), tj.get("tau_q", 0.999))
                t2 = _pick(fq.get("tau_q_tier2"), fq.get("tau_q_rescue", tj.get("tau_q_rescue", 0.99)))
                p.tau_tier2[vc] = t2; p.tau_rescue[vc] = t2
            else:
                if tj.get("tau") is not None:
                    p.tau[vc] = float(tj["tau"])
                if tj.get("tau_rescue") is not None:
                    p.tau_rescue[vc] = float(tj["tau_rescue"])
    if score_col == "auto":
        score_col = "rf_prob"
    rf = I.load_rf_probs(a.rf_probs, column=score_col)
    summ = I.integrate_table(a.evidence, a.features, a.out, p, a.class_group, rf_probs=rf)
    summ["thresholds_version"] = thr.get("version")
    summ["call_mode"] = "rf+phase" if rf else "phase_only"
    summ["score_column"] = score_col if rf else None
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
    # blood-derived family (fold stratification, DESIGN P12): found by SAMPLE-id prefix (the REACH quad's family id follows the
    # cohort's F0xxx pattern, so a family-id prefix alone matched nothing — fixed 2026-09-13) or by family-id prefix; both from
    # the env file when set (BLOOD_SAMPLE_PREFIX / BLOOD_FAMILY_PREFIX).
    sources = {}
    for r in rows:
        if (a.blood_sample_prefix and r["sample_id"].startswith(a.blood_sample_prefix)) or (
            a.blood_family_prefix and r["family_id"].startswith(a.blood_family_prefix)
        ):
            sources[r["family_id"]] = "blood"
    sys.stderr.write("SWAP blood families: %d (%s)\n" % (len(sources), ", ".join(sorted(sources)) or "none"))
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


def cmd_train(a: argparse.Namespace) -> int:
    """M4: the one harness - nested CV per class group over seeds, ablations, baselines, heuristic arms, rf_probs, tau, frozen model."""
    from .eval import harness as HZ
    from .features.registry import Registry
    man = read_manifest(a.manifest)
    fam_of = {sid: r["family_id"] for sid, r in man.items()}
    reg = Registry(a.registry)
    log = lambda m: sys.stderr.write(m + "\n")
    seeds = [int(x) for x in a.seeds.split(",")]
    classes = {"snv_indel": ("SNV", "INDEL"), "sv": ("SV",), "tr": ("TR",)}[a.class_group]
    allowed = {c for vc in classes for c in reg.rf_matrix_columns(vc)}
    rep = HZ.run_class(a.class_group, a.evidence_dir, a.train_dir, a.folds_dir, seeds, a.out_dir, fam_of, a.max_real_per_child,
                       reg.manifest(), log=log, baselines=not a.no_baselines, freeze=not a.no_freeze, allowed_cols=allowed)
    for arm, sm in sorted(rep["summary"].items()):
        log("SUMMARY %-22s roc_auc %s [%s-%s] pr_auc %s%s" % (arm, None if sm["roc_auc_mean"] is None else round(sm["roc_auc_mean"], 4),
            None if sm["roc_auc_min"] is None else round(sm["roc_auc_min"], 4), None if sm["roc_auc_max"] is None else round(sm["roc_auc_max"], 4),
            None if sm["pr_auc_mean"] is None else round(sm["pr_auc_mean"], 4),
            "" if sm.get("op_tpr_mean") is None else "  op tpr %.3f fpr %.4f" % (sm["op_tpr_mean"], sm["op_fpr_mean"])))
    log("tau: %s" % json.dumps(rep["tau"]))
    return 0


def cmd_annotate(a: argparse.Namespace) -> int:
    """P26: gnomAD AF (slivar gnotate), leave-one-family-out founder counts, sib-shared, for one class group.
    Real trios: all candidate tables matching --cand-glob; synthetic trio: one table with --exclude-families."""
    import csv, glob, shlex
    from . import annotate as AN
    from .io.vcf import read_manifest
    man = read_manifest(a.manifest)
    rows = list(man.values())
    founder_family = AN.founders(rows)
    paths = sorted(glob.glob(a.cand_glob)) if a.cand_glob else [a.candidates]
    if not paths:
        sys.stderr.write("annotate: no candidate tables\n"); return 2
    os.makedirs(a.out_dir, exist_ok=True); os.makedirs(a.work, exist_ok=True)
    sites, _ = AN.union_sites(paths, a.class_group)
    log = lambda m: sys.stderr.write(m + "\n")
    log("annotate %s: %d candidate tables, %d distinct sites" % (a.class_group, len(paths), len(sites)))
    gnomad = {}
    if a.class_group == "snv_indel" and a.gnomad_zip and a.slivar_cmd:
        sv_path = os.path.join(a.work, "sites.%s.vcf" % a.class_group)
        n = AN.write_sites_vcf(sites, sv_path, header_lines=AN.contig_lines(a.cohort_vcf, a.bcftools) if a.cohort_vcf else ())
        gnomad = AN.gnotate(sv_path, os.path.join(a.work, "sites.%s.gnotate.vcf" % a.class_group), shlex.split(a.slivar_cmd), a.gnomad_zip)
        log("gnotate: %d sites written, %d with gnomad_af" % (n, len(gnomad)))
    fgt, order = {}, []
    if a.cohort_vcf and a.class_group in ("snv_indel", "sv"):
        bed = os.path.join(a.work, "sites.%s.bed" % a.class_group)
        AN.write_sites_bed(sites, bed)
        if a.class_group == "snv_indel":
            order, fgt = AN.founder_genotypes(a.cohort_vcf, bed, sorted(founder_family), a.bcftools, a.work)
        else:
            order, fgt = AN.founder_genotypes_sv(a.cohort_vcf, sorted(founder_family), a.bcftools, a.work) if hasattr(AN, "founder_genotypes_sv") else ([], {})
        log("founder genotypes at %d sites (%d founders)" % (len(fgt), len(order)))
    tr_table = None
    strchive = None
    if a.class_group == "tr" and a.cohort_vcf:
        # The TR branch did not exist until 2026-09-15, so every TR annotation column came out empty and both the P24
        # rarity gate and the P15 population rules were vacuous for this class.
        bed = os.path.join(a.work, "sites.tr.bed")
        n_bed = AN.tr_bed(sites, bed)
        order, tr_table = AN.tr_founder_alleles(a.cohort_vcf, bed, sorted(founder_family), a.bcftools, a.work,
                                                min_spanning=a.tr_min_spanning)
        log("TR founder reference: %d candidate loci, %d loci with founder alleles, %d founders" % (n_bed, len(tr_table), len(order)))
        if a.strchive:
            import json as _json
            try:
                cat = _json.load(open(a.strchive))
                recs = cat if isinstance(cat, list) else cat.get("loci", [])
                strchive = {}
                for x in recs:
                    if not isinstance(x, dict):
                        continue
                    ch, s0, e0 = x.get("chrom"), x.get("start_hg38"), x.get("stop_hg38")
                    if not ch or s0 is None or e0 is None:
                        continue
                    strchive.setdefault(ch, []).append((int(s0), int(e0), {
                        "id": str(x.get("id") or ""), "gene": str(x.get("gene") or ""),
                        "disease": str(x.get("disease") or ""), "inheritance": ",".join(x.get("inheritance") or []),
                        "pathogenic_min": x.get("pathogenic_min"), "benign_max": x.get("benign_max"),
                        "motif_len": x.get("motif_len")}))
                for ch in strchive:
                    strchive[ch].sort()
                log("STRchive catalogue: %d known pathogenic loci over %d contigs"
                    % (sum(len(v) for v in strchive.values()), len(strchive)))
            except (OSError, ValueError) as e:
                log("STRchive catalogue could not be read (%s); strchive_locus stays empty" % e)
        else:
            log("no --strchive catalogue given: strchive_locus stays EMPTY (known pathogenic repeat loci are not flagged)")
    excl_extra = set(x for x in (a.exclude_families or "").split(",") if x)
    total = {"rows": 0, "gnomad_annotated": 0, "founder_counts": 0}
    for cp in paths:
        child = os.path.basename(cp).split(".")[0]
        fam = man[child]["family_id"] if child in man else None
        excl = set(excl_extra)
        if fam:
            excl.add(fam)
        sib: set = set()
        if a.class_group == "snv_indel" and a.joint_vcf_pattern and fam and not excl_extra:
            sibs = [r["sample_id"] for r in rows if r["family_id"] == fam and r.get("role") == "offspring" and r["sample_id"] != child
                    and r.get("father_id") in man and r.get("mother_id") in man]
            jv = glob.glob(a.joint_vcf_pattern.replace("{FAMILY}", fam))
            if sibs and jv:
                keys = {(r_.chrom,) + AN.norm_allele(r_.start, r_.ref, r_.alt) for r_ in AN.read_candidates(cp)}
                sib = AN.sib_shared_sites(jv[0], child, sibs, keys)
        out = os.path.join(a.out_dir, "%s.%s.annot.tsv" % (child, a.class_group))
        st = AN.write_annot(cp, out, a.class_group, gnomad, fgt, order, founder_family, excl, sib,
                            tr_table=tr_table, strchive=strchive)
        for k in total:
            total[k] += st[k]
    log("annotate %s done: %s -> %s" % (a.class_group, " ".join("%s=%d" % kv for kv in total.items()), a.out_dir))
    return 0


def cmd_external(a: argparse.Namespace) -> int:
    """P27: external-truth arm (spike-ins) scored by the held-out fold models of one seed."""
    import csv
    from .eval import external as EX
    man = read_manifest(a.manifest)
    fam_of = {sid: r["family_id"] for sid, r in man.items()}
    fold_of = {r["family_id"]: int(r["outer_fold"]) for r in csv.DictReader(open(os.path.join(a.folds_dir, "folds.seed%d.tsv" % a.seed)), delimiter="\t")}
    taus = {}
    tj = os.path.join(a.harness_dir, "tau.%s.json" % a.class_group)
    if os.path.exists(tj):
        t = json.load(open(tj))
        if t.get("tau") is not None:
            taus[a.class_group] = float(t["tau"])
        for k in ("score_column", "tau_q", "tau_q_rescue"):
            if k in t:
                taus[k] = t[k]
    log = lambda m: sys.stderr.write(m + "\n")
    if a.truth == "wes":
        rep = EX.evaluate_labelled(a.evidence_dir, a.harness_dir, a.class_group, a.seed, fam_of, fold_of, taus, a.labels_dir, log=log)
    else:
        rep = EX.evaluate_spikes(a.evidence_dir, a.harness_dir, a.class_group, a.seed, fam_of, fold_of, taus, max_real_per_child=a.max_real_per_child, log=log)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as fh:
        json.dump(rep, fh, indent=1, sort_keys=True, default=str)
    for arm, r in sorted(rep.get("arms", {}).items()):
        log("EXTERNAL %-22s roc_auc %s pr_auc %s%s" % (arm, None if r["roc_auc"] is None else round(r["roc_auc"], 4), None if r["pr_auc"] is None else round(r["pr_auc"], 4),
            "" if "tpr" not in r else "  op tpr %.3f fpr %.4f" % (r["tpr"], r["fpr"])))
    for k in ("recall_at_tau_by_class", "mosaic_sensitivity_at_tau_CM", "mosaic_sensitivity_at_tau_PM", "tau_q_for_recall", "recall_at_tau", "neg_pass_at_tau", "n_pos", "n_rows"):
        if k in rep:
            log("EXTERNAL %s: %s" % (k, json.dumps(rep[k])))
    return 0


def cmd_rescore(a: argparse.Namespace) -> int:
    """M4 post-step: fold-quantile score rf_q for every real candidate from the saved fold models (train/rescore.py)."""
    from .train import rescore as RS
    man = read_manifest(a.manifest)
    fam_of = {sid: r["family_id"] for sid, r in man.items()}
    log = lambda m: sys.stderr.write(m + "\n")
    rep = RS.rescore_class(a.class_group, a.harness_dir, a.evidence_dir, a.folds_dir, [int(x) for x in a.seeds.split(",")], fam_of,
                           max_ref_per_child=a.max_real_per_child, log=log)
    log("rescore %s: %s" % (a.class_group, json.dumps(rep)))
    return 0


def cmd_score(a: argparse.Namespace) -> int:
    """P21 transfer path: score a cohort with a FROZEN model (models/<class>.xgb.json); rf_q against the scored cohort itself."""
    from .train import score as SC
    man = read_manifest(a.manifest)
    ids = set(man)
    kids = sorted({(r["family_id"], sid) for sid, r in man.items()
                   if r.get("role") == "offspring" and r.get("father_id") in ids and r.get("mother_id") in ids})
    log = lambda m: sys.stderr.write(m + "\n")
    rep = SC.score_class(a.class_group, a.model, a.evidence_dir, a.out_dir, kids, manifest_path=a.model_manifest, log=log)
    log("score %s: %s" % (a.class_group, json.dumps(rep)))
    return 0


def cmd_attribution(a: argparse.Namespace) -> int:
    """P14: TreeSHAP share by feature family on held-out rows (real, synthetic, classifier-called real rows)."""
    from .eval import attribution as AT
    man = read_manifest(a.manifest)
    fam_of = {sid: r["family_id"] for sid, r in man.items()}
    log = lambda m: sys.stderr.write(m + "\n")
    rep = AT.run(a.class_group, a.harness_dir, a.evidence_dir, a.train_dir, a.folds_dir, a.seed, fam_of, log=log)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as fh:
        json.dump(rep, fh, indent=1, sort_keys=True)
    for kind, d in rep["share_by_family"].items():
        log("ATTRIBUTION %-12s %s (n=%d)" % (kind, json.dumps(d), rep["n_rows"].get(kind, 0)))
    log("ATTRIBUTION top features: %s" % json.dumps(rep["top_features"]))
    return 0


def cmd_wes_truth(a: argparse.Namespace) -> int:
    """P27 arm 2: WES-confirmed labels for exonic small-variant candidates from the iWES pVCF (per chromosome)."""
    import glob
    from .eval import wes_truth as W
    man = read_manifest(a.manifest)
    cand_paths = {os.path.basename(p).split(".")[0]: p for p in glob.glob(a.cand_glob)}
    log = lambda m: sys.stderr.write(m + "\n")
    os.makedirs(a.out_dir, exist_ok=True)
    if a.step in ("sites", "extract"):
        targets = W.Targets(a.target_bed)
        sites = W.exonic_sites(cand_paths.values(), targets)
        n = W.write_regions(sites, a.chrom, os.path.join(a.out_dir, "sites.%s.bed" % a.chrom))
        log("wes-truth %s: %d exonic candidate positions" % (a.chrom, n))
        if a.step == "sites" or n == 0:
            return 0
        samples = os.path.join(a.out_dir, "cohort_samples.txt")
        with open(samples, "w") as fh:
            fh.write("\n".join(sorted(man)) + "\n")
        pvcf = a.pvcf_pattern.replace("{CHROM}", a.chrom)
        n_rec = W.extract_chrom(pvcf, os.path.join(a.out_dir, "sites.%s.bed" % a.chrom), samples, os.path.join(a.out_dir, "wes.%s.tsv" % a.chrom), a.bcftools)
        order = W.sample_order(pvcf, samples, a.bcftools)
        with open(os.path.join(a.out_dir, "wes.%s.samples.txt" % a.chrom), "w") as fh:
            fh.write("\n".join(order) + "\n")
        log("wes-truth %s: %d records extracted for %d cohort samples present" % (a.chrom, n_rec, len(order)))
        return 0
    # label
    order = open(os.path.join(a.out_dir, "wes.%s.samples.txt" % a.chrom)).read().split()
    trios = {}
    for sid, r in man.items():
        if r.get("role") == "offspring" and r.get("father_id") in man and r.get("mother_id") in man:
            trios[sid] = (r["father_id"], r["mother_id"])
    st = W.label_children(os.path.join(a.out_dir, "wes.%s.tsv" % a.chrom), order, trios, cand_paths, os.path.join(a.out_dir, "labels"), chrom_filter=a.chrom)
    tot = {}
    for c in st.values():
        for k, v in c.items():
            tot[k] = tot.get(k, 0) + v
    log("wes-truth label %s: %d children; %s" % (a.chrom, len(st), " ".join("%s=%d" % kv for kv in sorted(tot.items()))))
    return 0


def cmd_wes_merge(a: argparse.Namespace) -> int:
    """Concatenate per-chromosome label files into one <child>.snv_indel.wes.tsv per child (labels != -1 summarised)."""
    import glob, csv
    from collections import Counter
    labels_dir = os.path.join(a.out_dir, "labels")
    by_child = {}
    for f in sorted(glob.glob(os.path.join(labels_dir, "*.snv_indel.wes.chr*.tsv"))):
        by_child.setdefault(os.path.basename(f).split(".")[0], []).append(f)
    tot = Counter()
    for child, files in by_child.items():
        out = os.path.join(labels_dir, "%s.snv_indel.wes.tsv" % child)
        with open(out, "w", newline="") as fh:
            w = None
            for f in files:
                with open(f, newline="") as ih:
                    rd = csv.DictReader(ih, delimiter="\t")
                    if w is None:
                        w = csv.DictWriter(fh, fieldnames=rd.fieldnames, delimiter="\t", lineterminator="\n"); w.writeheader()
                    for r in rd:
                        w.writerow(r); tot[r["label"]] += 1
    sys.stderr.write("wes-merge: %d children; labels %s\n" % (len(by_child), dict(tot)))
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
                          class_filter=a.only_class, smallvar_vcf=a.smallvar_vcf, child_sample=a.child_sample)
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
    annot = None
    if a.annot and os.path.exists(a.annot):
        import csv
        with open(a.annot, newline="") as fh:
            annot = {r["variant_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    geom = X.Geometry(a.orientation, a.changepoints) if (a.orientation or a.changepoints) else None
    summ = X.extract_child(a.evidence, a.candidates, reg, sex, a.out, a.rf_out, mask=mask, seqctx=seqctx, annot=annot, geom=geom)
    summ["registry"] = reg.manifest()
    with open(a.out.replace(".tsv", ".summary.json"), "w") as fh:
        json.dump(summ, fh, indent=1, sort_keys=True)
    never = summ.get("never_produced", [])
    if never and a.require_all_features:
        sys.stderr.write("ERROR: %d registry features applicable to this class were never produced: %s\n"
                         "       A declared feature that is always empty is silently dropped by the presence-leak guard and the\n"
                         "       classifier never sees it. Implement it or set its registry status to drop. (--no-require-all-features\n"
                         "       to override.)\n" % (len(never), ", ".join(never)))
        return 4
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
    rv.add_argument("--smallvar-vcf", help="the family's PHASED small-variant VCF; SV heterozygous-SNV persistence inside the interval vs the flanks (P17)")
    rv.add_argument("--child-sample", help="sample id of the child in --smallvar-vcf (default: the candidate's sample_id)")
    rv.add_argument("--max-per-class", type=int, help="cap rows per class (smoke tests)")
    rv.add_argument("--only-class", action="append", choices=["SNV", "INDEL", "SV", "TR"])
    rv.add_argument("--salt", help="salt for read-id hashing (default: the output path)")
    rv.add_argument("--thresholds")
    rv.set_defaults(func=cmd_review)

    fe = sub.add_parser("features", help="registry-driven feature vectors from an evidence table (+ rf_safe classifier matrix)")
    fe.add_argument("--evidence", required=True), fe.add_argument("--candidates", required=True)
    fe.add_argument("--child", required=True), fe.add_argument("--manifest", required=True)
    fe.add_argument("--out", required=True, help="<child>.<class>.features.tsv")
    fe.add_argument("--annot", help="annot/<child>.<class>.annot.tsv from `annotate` (gnomad_af, cohort_AC_loo, pon_founder_recurrence_loo, sib_shared)")
    fe.add_argument("--orientation", help="M1 <child>.orientation.tsv: segment geometry features (child-only, rf_safe)")
    fe.add_argument("--changepoints", help="M1 <child>.changepoints.resolved.tsv: crossover distance (rf_safe: false)")
    fe.add_argument("--rf-out", help="<child>.<class>.features.rf.tsv (rf_safe columns only; refused if unsafe)")
    fe.add_argument("--mask", action="append", help="BED file(s) of the lab region mask (flag, never a filter)")
    fe.add_argument("--reference", help="reference FASTA for sequence-context features (pysam)")
    fe.add_argument("--registry", help="config/features.yaml (default: the module's)")
    fe.add_argument("--require-all-features", dest="require_all_features", action="store_true", default=True,
                    help="fail when a registry feature applicable to the class is never produced (a declared-but-empty column is "
                         "dropped by the presence-leak guard and never reaches the classifier; this gate would have caught the SV "
                         "depth block, declared in 2026-09 and unimplemented until 2026-09-14)")
    fe.add_argument("--no-require-all-features", dest="require_all_features", action="store_false")
    fe.set_defaults(func=cmd_features)

    sk = sub.add_parser("spike", help="spike-in harness: plan sites, edit haplotagged reads into slice BAMs, evaluate recovery")
    sk.add_argument("step", choices=["plan", "apply", "genotype", "evaluate"])
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
    sk.add_argument("--reference", help="genotype: reference FASTA for bcftools mpileup over the spiked slices")
    sk.add_argument("--bcftools", default="bcftools", help="genotype: bcftools executable")
    sk.set_defaults(func=cmd_spike)
    au = sub.add_parser("audit", help="P30: assert that what the config declares is what the data contains")
    au.add_argument("--evidence-dir", required=True)
    au.add_argument("--harness-dir", help="a harness output dir; enables the train/score column check")
    au.add_argument("--registry", help="config/features.yaml (default: the module's)")
    au.add_argument("--class-group", action="append", choices=["snv_indel", "sv", "tr"],
                    help="restrict to these class groups (repeatable; default all three)")
    au.add_argument("--no-spike", dest="spike", action="store_false", default=True)
    au.set_defaults(func=cmd_audit)
    ig = sub.add_parser("integrate", help="M3: final unfiltered table per child and class group (P15 decision, class columns, feature vector)")
    ig.add_argument("--evidence", required=True, help="<child>.<class>.evidence.lik.tsv"), ig.add_argument("--features", help="<child>.<class>.features.tsv")
    ig.add_argument("--class-group", required=True, choices=["snv_indel", "sv", "tr"])
    ig.add_argument("--out", required=True, help="final/<FAMILY>.<child>.<class>.dnm.tsv"), ig.add_argument("--vcf-out")
    ig.add_argument("--rf-probs", help="M4 output: TSV with variant_id, rf_prob (absent -> provisional phase_only mode)")
    ig.add_argument("--tau-json", help="M4 output tau.<class>.json: tau / tau_rescue for this class group (overrides thresholds.yaml final:)")
    ig.add_argument("--score-column", default="auto", choices=["auto", "rf_prob", "rf_q"], help="which rf_probs column drives the decision (auto: the tau json's score_column)")
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
    sw.add_argument("--blood-sample-prefix", default=os.environ.get("BLOOD_SAMPLE_PREFIX", "REACH"),
                    help="sample-id prefix identifying the blood-derived family (kept with company in fold 0); env BLOOD_SAMPLE_PREFIX")
    sw.add_argument("--blood-family-prefix", default=os.environ.get("BLOOD_FAMILY_PREFIX", ""),
                    help="family-id prefix identifying the blood-derived family (alternative to the sample prefix); env BLOOD_FAMILY_PREFIX")
    sw.set_defaults(func=cmd_swap)
    tr = sub.add_parser("train", help="M4: nested CV harness per class group (RF arms, ablations, baselines, heuristic sweeps), rf_probs, tau, frozen model")
    tr.add_argument("--class-group", required=True, choices=["snv_indel", "sv", "tr"])
    tr.add_argument("--manifest", required=True), tr.add_argument("--evidence-dir", required=True), tr.add_argument("--train-dir", required=True)
    tr.add_argument("--folds-dir", required=True), tr.add_argument("--out-dir", required=True)
    tr.add_argument("--seeds", default="0"), tr.add_argument("--max-real-per-child", type=int, default=20000)
    tr.add_argument("--registry"), tr.add_argument("--no-baselines", action="store_true"), tr.add_argument("--no-freeze", action="store_true")
    tr.set_defaults(func=cmd_train)
    an = sub.add_parser("annotate", help="P26: gnomAD AF, leave-one-family-out founder counts, sib-shared for one class group")
    an.add_argument("--class-group", required=True, choices=["snv_indel", "sv", "tr"]), an.add_argument("--manifest", required=True)
    an.add_argument("--cand-glob", help="real trios: glob of <child>.<class>.candidates.tsv tables"), an.add_argument("--candidates", help="one table (synthetic trio)")
    an.add_argument("--exclude-families", help="synthetic trio: comma-separated families to exclude from the founder panel (child's and parents')")
    an.add_argument("--out-dir", required=True), an.add_argument("--work", required=True)
    an.add_argument("--cohort-vcf", help="cohort BCF (snv_indel) or cohort SV VCF (sv)"), an.add_argument("--bcftools", default="bcftools")
    an.add_argument("--slivar-cmd", help="command prefix, e.g. 'singularity exec -B /expanse:/expanse slivar.sif slivar'"), an.add_argument("--gnomad-zip")
    an.add_argument("--joint-vcf-pattern", help="family joint small-variant VCF glob with {FAMILY} (sib-shared in quads)")
    an.add_argument("--tr-min-spanning", type=int, default=5, help="TR: minimum spanning-read depth (SD) for a founder allele to enter the reference (02_tiering/tr_outliers.py uses 5)")
    an.add_argument("--strchive", help="TR: STRchive catalogue JSON; sets strchive_locus for known pathogenic repeat loci")
    an.set_defaults(func=cmd_annotate)
    ex = sub.add_parser("external", help="P27: external-truth arm - spike-ins scored by held-out fold models, same arms as the harness")
    ex.add_argument("--class-group", required=True, choices=["snv_indel", "sv", "tr"]), ex.add_argument("--manifest", required=True)
    ex.add_argument("--evidence-dir", required=True), ex.add_argument("--harness-dir", required=True), ex.add_argument("--folds-dir", required=True)
    ex.add_argument("--seed", type=int, default=0), ex.add_argument("--max-real-per-child", type=int, default=20000), ex.add_argument("--out", required=True)
    ex.add_argument("--truth", default="spike", choices=["spike", "wes"]), ex.add_argument("--labels-dir", help="wes: directory of <child>.snv_indel.wes.tsv label files")
    ex.set_defaults(func=cmd_external)
    rs = sub.add_parser("rescore", help="M4 post-step: fold-quantile score rf_q from saved fold models; writes rf_probs/ (rf_prob + rf_q) and tau_q")
    rs.add_argument("--class-group", required=True, choices=["snv_indel", "sv", "tr"]), rs.add_argument("--manifest", required=True)
    rs.add_argument("--evidence-dir", required=True), rs.add_argument("--harness-dir", required=True), rs.add_argument("--folds-dir", required=True)
    rs.add_argument("--seeds", default="0"), rs.add_argument("--max-real-per-child", type=int, default=20000)
    rs.set_defaults(func=cmd_rescore)
    sc = sub.add_parser("score", help="P21: score a cohort with a FROZEN model (models/<class>.xgb.json) -> rf_probs/ (rf_prob, rf_q) + tau json; then `integrate`")
    sc.add_argument("--class-group", required=True, choices=["snv_indel", "sv", "tr"]), sc.add_argument("--manifest", required=True)
    sc.add_argument("--evidence-dir", required=True), sc.add_argument("--out-dir", required=True, help="written: rf_probs/<child>.<class>.rf_probs.tsv, tau.<class>.json")
    sc.add_argument("--model", required=True, help="frozen model, e.g. models/snv_indel.xgb.json"), sc.add_argument("--model-manifest", help="default: <model>.training_manifest.json beside it")
    sc.set_defaults(func=cmd_score)
    at = sub.add_parser("attribution", help="P14: TreeSHAP share by feature family (A caller / B context / C reads / D phase) on held-out rows")
    at.add_argument("--class-group", required=True, choices=["snv_indel", "sv", "tr"]), at.add_argument("--manifest", required=True)
    at.add_argument("--evidence-dir", required=True), at.add_argument("--train-dir", required=True), at.add_argument("--harness-dir", required=True)
    at.add_argument("--folds-dir", required=True), at.add_argument("--seed", type=int, default=0), at.add_argument("--out", required=True)
    at.set_defaults(func=cmd_attribution)
    wt = sub.add_parser("wes-truth", help="P27 arm 2: exonic candidate sites -> iWES pVCF extraction -> trio labels, one chromosome")
    wt.add_argument("step", choices=["sites", "extract", "label"]), wt.add_argument("--chrom", required=True)
    wt.add_argument("--manifest", required=True), wt.add_argument("--cand-glob", required=True, help="<child>.snv_indel.candidates.tsv tables")
    wt.add_argument("--target-bed", action="append", default=[], help="capture target BED(s); repeatable")
    wt.add_argument("--pvcf-pattern", help="per-chromosome pVCF path with {CHROM}"), wt.add_argument("--bcftools", default="bcftools")
    wt.add_argument("--out-dir", required=True)
    wt.set_defaults(func=cmd_wes_truth)
    wm = sub.add_parser("wes-merge", help="merge per-chromosome WES label files per child")
    wm.add_argument("--out-dir", required=True)
    wm.set_defaults(func=cmd_wes_merge)
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
