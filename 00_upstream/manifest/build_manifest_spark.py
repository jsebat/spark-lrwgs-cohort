#!/usr/bin/env python3
"""Emit a cohort manifest (samples.tsv rows) for the SPARK PacBio pilot.

Reads SPARK's own metadata and discovers the per-movie uBAMs under --reads-root,
so the same command works whether reads live on lustre staging or the filer.

Usage:
  build_manifest_spark.py --reads-root DIR [--metadata DIR] [--min-members N]
                          [--only-family SFID] [--out samples.tsv]
"""
import argparse, csv, glob, os, sys, collections

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reads-root", required=True,
        help="dir containing phase1/<SPID>/*.hifi_reads.bam")
    ap.add_argument("--metadata",
        default=os.environ.get("SPARK_METADATA_DIR"), help="directory with SPARK individual metadata (env SPARK_METADATA_DIR)")
    ap.add_argument("--only-family", action="append")
    ap.add_argument("--min-members", type=int, default=0)
    ap.add_argument("--out", default="-")
    a = ap.parse_args()

    ind = list(csv.DictReader(
        open(os.path.join(a.metadata,"SPARK.PacBioPilot.2026_06.individual_metadata.tsv")),
        delimiter="\t"))

    reads = {}
    for spid_dir in glob.glob(os.path.join(a.reads_root, "*", "*")):
        spid = os.path.basename(spid_dir)
        bams = sorted(glob.glob(os.path.join(spid_dir, "*.hifi_reads.bam")))
        if bams: reads.setdefault(spid, []).extend(bams)

    fams = collections.defaultdict(list)
    for r in ind: fams[r["sfid"]].append(r)

    SEX = {"male":"MALE","female":"FEMALE"}
    out = sys.stdout if a.out == "-" else open(a.out,"w")
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(["cohort","family_id","sample_id","sex","affected",
                "father_id","mother_id","hifi_reads"])
    nfam = nsamp = 0
    skipped = []
    for fid, members in sorted(fams.items()):
        if a.only_family and fid not in a.only_family: continue
        seq = [m for m in members if m["spid"] in reads]
        if len(seq) < a.min_members or not seq:
            skipped.append((fid, len(members), len(seq))); continue
        ids = {m["spid"] for m in seq}
        for m in seq:
            fa = m["father_spid"] if m["father_spid"] in ids else ""
            mo = m["mother_spid"] if m["mother_spid"] in ids else ""
            w.writerow(["SPARK", fid, m["spid"],
                        SEX.get(m["sex"].strip().lower(),""),
                        "true" if m["asd"].strip()=="True" else "false",
                        fa, mo, ",".join(reads[m["spid"]])])
            nsamp += 1
        nfam += 1
    if out is not sys.stdout: out.close()
    print(f"# {nfam} families, {nsamp} samples", file=sys.stderr)
    if skipped:
        print(f"# skipped {len(skipped)} families below --min-members", file=sys.stderr)

if __name__ == "__main__":
    main()
