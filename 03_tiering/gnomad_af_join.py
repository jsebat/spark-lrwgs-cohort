#!/usr/bin/env python3
"""Replace the gnomAD AF column of lof_sites.tsv with the gnomAD v4.1 GENOMES frequency from a slivar-gnotated VCF.

Why: lr_vep.sb took gnomAD AF from the dbNSFP VEP plugin, which covers coding SNVs only. Every indel LoF -- the
frameshifts that are the bulk of HC LoF -- had an empty AF, and the rare filter downstream read "empty" as "absent
from gnomAD = rare", so common frameshifts entered lof_tiered.rare.tsv, the burden counts, the TDT and the inherited
tables (review 2026-09-16, T2). The lab's slivar zip is built from gnomad.genomes.v4.1.sites (SNVs and indels; on
this cohort's HC-LoF set 656 of 1,077 indels and 365 of 500 SNVs have a frequency, the rest are genuinely absent).

Columns 1-10 are unchanged in position (lof_tier.sb reads them by index); the dbNSFP value moves to a trailing
`gnomad_af_dbnsfp` column for provenance. slivar writes -1 for a site absent from gnomAD; that is written as "." here,
which the consumers already treat as absent (= rare).

usage: gnomad_af_join.py lof_sites.tsv gnotated.vcf out.tsv
"""
import csv
import gzip
import sys


def main():
    sites, vcf, out = sys.argv[1:4]
    af = {}
    op = gzip.open if vcf.endswith(".gz") else open
    with op(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            v = info.get("gnomad_af", "-1")
            af[(f[0], f[1], f[3], f[4])] = "." if v in ("-1", "-1.0", ".", "") else v
    n = hit = absent = 0
    with open(sites, newline="") as fh, open(out, "w", newline="") as oh:
        rd = csv.reader(fh, delimiter="\t")
        hdr = next(rd)
        assert hdr[:4] == ["chrom", "pos", "ref", "alt"] and hdr[7] == "gnomad_af", hdr
        w = csv.writer(oh, delimiter="\t", lineterminator="\n")
        w.writerow(hdr + ["gnomad_af_dbnsfp"])
        for r in rd:
            if len(r) < 10:
                continue
            n += 1
            key = (r[0], r[1], r[2], r[3])
            old = r[7]
            if key in af:
                r[7] = af[key]
                hit += af[key] != "."
                absent += af[key] == "."
            else:
                r[7] = "."
                absent += 1
            w.writerow(r + [old])
    sys.stderr.write("gnomad_af_join: %d sites; gnomAD genomes v4.1 frequency for %d, absent %d\n" % (n, hit, absent))


if __name__ == "__main__":
    main()
