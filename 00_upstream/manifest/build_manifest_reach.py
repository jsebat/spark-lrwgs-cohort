#!/usr/bin/env python3
"""Emit a cohort manifest (samples.tsv) for REACH families from the
REACH long-read metadata manifest.

Usage:
  build_manifest_reach.py --reads-root DIR [--metadata FILE] [--only-family F] [--out FILE]

--reads-root is the directory holding the delivered SMRT Link tree, e.g.
  <staging root>            (staging)
  <permanent root>/<FAMILY>/raw   (permanent)
Reads are located by <movie>.hifi_reads.<barcode>.bam anywhere beneath it, so
the same command works against either tier.
"""
import argparse, csv, glob, os, sys, collections

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reads-root", required=True)
    ap.add_argument("--metadata",
        default=os.environ.get("REACH_METADATA_TSV"), help="REACH long-read manifest TSV (env REACH_METADATA_TSV)")
    ap.add_argument("--only-family", action="append")
    ap.add_argument("--out", default="-")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.metadata), delimiter="\t"))

    # index every hifi_reads bam under reads-root by its barcode
    bybc = collections.defaultdict(list)
    for p in glob.glob(os.path.join(a.reads_root, "**", "*.hifi_reads.*.bam"), recursive=True):
        base = os.path.basename(p)
        if ".unassigned." in base: continue
        bc = base.rsplit(".hifi_reads.", 1)[1][:-4]      # strip trailing .bam
        bybc[bc].append(p)

    # relationship -> parent lookup, per family
    fams = collections.defaultdict(list)
    for r in rows: fams[r["family"]].append(r)

    out = sys.stdout if a.out == "-" else open(a.out, "w")
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(["cohort","family_id","sample_id","sex","affected",
                "father_id","mother_id","hifi_reads"])
    n = 0
    missing = []
    for fid, members in sorted(fams.items()):
        if a.only_family and fid not in a.only_family: continue
        byrole = {m["relationship"]: m["reach_id"] for m in members}
        fa = byrole.get("Dad", ""); mo = byrole.get("Mom", "")
        for m in members:
            bams = sorted(bybc.get(m["lr_barcode"], []))
            if not bams:
                missing.append(f'{m["reach_id"]} (barcode {m["lr_barcode"]})'); continue
            offspring = m["relationship"] in ("Proband", "Sibling")
            w.writerow(["REACH", fid, m["reach_id"],
                        (m.get("sex") or "").upper(),
                        (m.get("affected") or ("true" if m["relationship"] == "Proband" else "false")).lower(),   # relationship is the fallback, not "false" for everyone
                        fa if offspring else "",
                        mo if offspring else "",
                        ",".join(bams)])
            n += 1
    if out is not sys.stdout: out.close()
    print(f"# {n} samples written", file=sys.stderr)
    if missing:
        print(f"# NO READS FOUND for: {', '.join(missing)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
