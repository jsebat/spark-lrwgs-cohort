#!/usr/bin/env python3
"""Generate per-family HiFi-human-WGS-WDL inputs.json from a cohort manifest.

Manifest = TSV, one row per individual, columns:
  cohort family_id sample_id sex affected father_id mother_id hifi_reads
    sex        MALE | FEMALE | (blank -> omitted; pipeline infers)
    affected   true | false
    father_id  sample_id of father, or blank
    mother_id  sample_id of mother, or blank
    hifi_reads one or more paths, comma-separated

Usage:
  make_inputs.py --manifest samples.tsv --outdir config/inputs [--family FAM01]
                 [--config-dir DIR] [--check-files/--no-check-files]
"""
import argparse, csv, json, os, sys, collections

REQUIRED = ["cohort","family_id","sample_id","sex","affected","father_id","mother_id","hifi_reads"]

def parse_bool(v, where):
    s = (v or "").strip().lower()
    if s in ("true","yes","1","y"):  return True
    if s in ("false","no","0","n",""): return False
    sys.exit(f"ERROR {where}: affected must be true/false, got {v!r}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--config-dir",
        default=os.environ.get("UPSTREAM_CONFIG_DIR", "config"))
    ap.add_argument("--family", action="append",
        help="only this family_id (repeatable); default all")
    ap.add_argument("--check-files", dest="check", action="store_true", default=True)
    ap.add_argument("--no-check-files", dest="check", action="store_false")
    a = ap.parse_args()

    with open(a.manifest) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        sys.exit("ERROR: manifest is empty")
    missing = [c for c in REQUIRED if c not in rows[0]]
    if missing:
        sys.exit(f"ERROR: manifest missing columns: {missing}")

    fams = collections.OrderedDict()
    for r in rows:
        fams.setdefault(r["family_id"].strip(), []).append(r)
    if a.family:
        want = set(a.family)
        unknown = want - set(fams)
        if unknown: sys.exit(f"ERROR: family_id not in manifest: {sorted(unknown)}")
        fams = {k: v for k, v in fams.items() if k in want}

    os.makedirs(a.outdir, exist_ok=True)
    # absolute: driver_family.sb changes directory before `miniwdl run`, so a relative path that passed this check
    # would not resolve there. The files are config/*.template.tsv with ${RESOURCES_ROOT} filled in, saved as *.hpc.tsv.
    ref  = os.path.abspath(os.path.join(a.config_dir, "GRCh38.ref_map.v3p1p0.hpc.tsv"))
    tert = os.path.abspath(os.path.join(a.config_dir, "GRCh38.tertiary_map.v3p1p0.hpc.tsv"))
    for p in (ref, tert):
        if not os.path.exists(p):
            sys.exit(f"ERROR: map file not found: {p}; fill in the matching config/*.template.tsv and save it under this name")

    problems, written = [], []
    for fid, members in fams.items():
        ids = {m["sample_id"].strip() for m in members}
        if len(ids) != len(members):
            problems.append(f"{fid}: duplicate sample_id")
        samples = []
        for m in members:
            sid = m["sample_id"].strip()
            where = f"{fid}/{sid}"
            reads = [p.strip() for p in m["hifi_reads"].split(",") if p.strip()]
            if not reads:
                problems.append(f"{where}: no hifi_reads"); continue
            if a.check:
                for p in reads:
                    if not os.path.exists(p): problems.append(f"{where}: missing {p}")
            s = {"sample_id": sid, "hifi_reads": reads,
                 "affected": parse_bool(m["affected"], where)}
            sex = (m["sex"] or "").strip().upper()
            if sex:
                if sex not in ("MALE","FEMALE"):
                    problems.append(f"{where}: sex must be MALE/FEMALE, got {sex!r}")
                else: s["sex"] = sex
            fa, mo = (m["father_id"] or "").strip(), (m["mother_id"] or "").strip()
            # a parent named but not sequenced in this family cannot be used
            if fa:
                if fa in ids: s["father_id"] = fa
                else: problems.append(f"{where}: father_id {fa} not a member of {fid}")
            if mo:
                if mo in ids: s["mother_id"] = mo
                else: problems.append(f"{where}: mother_id {mo} not a member of {fid}")
            samples.append(s)

        doc = {"humanwgs_family.family": {"family_id": fid, "samples": samples},
               "humanwgs_family.ref_map_file": ref,
               "humanwgs_family.tertiary_map_file": tert,
               "humanwgs_family.backend": "HPC",
               "humanwgs_family.preemptible": False}
        out = os.path.join(a.outdir, f"{fid}.inputs.json")
        with open(out, "w") as fh:
            json.dump(doc, fh, indent=2); fh.write("\n")
        written.append((fid, len(samples), sum(1 for s in samples if s["affected"]),
                        sum(1 for s in samples if "father_id" in s or "mother_id" in s)))

    w = max((len(f) for f, *_ in written), default=8)
    print(f"{'family':<{w}}  samples  affected  with_parents")
    for fid, n, naff, nkid in written:
        print(f"{fid:<{w}}  {n:>7}  {naff:>8}  {nkid:>12}")
    print(f"\nwrote {len(written)} inputs.json -> {a.outdir}")
    if problems:
        print(f"\n!! {len(problems)} PROBLEM(S):", file=sys.stderr)
        for p in problems[:40]: print("   " + p, file=sys.stderr)
        if len(problems) > 40: print(f"   ... and {len(problems)-40} more", file=sys.stderr)
        sys.exit(1)
    print("all checks passed")

if __name__ == "__main__":
    main()
