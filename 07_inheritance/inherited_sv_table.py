#!/usr/bin/env python
"""Rare inherited structural variants disrupting coding sequence, transmitted to affected children.

The companion to the small-variant table. The structural work in this project was de novo only; nothing
enumerated the inherited side, which is what this does, using the same gene sets and the same coding
annotation the de novo structural analysis used.

Inclusion
  1. Deletion, duplication, inversion or insertion from the cohort structural call set. Breakends are counted
     but not interpreted, because a single breakend record does not define a gene-disrupting event without
     resolving its mate.
  2. Overlapping coding sequence of a gene in one of the panels: SFARI autism genes, DDG2P developmental
     disorder genes, or genes under strong selection against heterozygous loss (shet >= 0.1).
  3. Rare in this cohort: carrier allele count at or below COHORT_AC_MAX across 105 genomes. There is no
     population frequency filter available for these calls, so cohort frequency does all the work, which
     makes it a weaker filter than the small-variant equivalent and that must be said plainly.
  4. Transmitted: carried by an affected child and by at least one observed carrier parent.

Genotype quality is carried through per carrier so the table can be read quality-first, as the small-variant
table is.
"""
import csv, gzip, io, pathlib, subprocess, sys
from collections import defaultdict

F = "/expanse/projects/sebat1/longread_cohort_2026/freeze1/freeze1.cohort.sv.vcf.gz"
T = pathlib.Path("/expanse/lustre/projects/ddp195/jsebat/longread-autism/tiering")
D = T / "denovo_sv"
OUT = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else str(T / "inherited_sv_constrained.tsv"))
BCF = str(pathlib.Path.home() / "micromamba/envs/dnmt3a-py/bin/bcftools")
COHORT_AC_MAX = 4
SHET_MIN = 0.10
KEEP_TYPES = {"DEL", "DUP", "INV", "INS"}

# ---------------------------------------------------------------- pedigree and gene sets
ped = {}
for line in open(T / "_all.ped"):
    p = line.split()
    if len(p) >= 6:
        ped[p[1]] = dict(fam=p[0], pat=p[2], mat=p[3], aff=p[5])
cases = {s for s, v in ped.items() if v["pat"] != "0" and v["mat"] != "0" and v["aff"] == "2"}

def load(fn):
    return {l.strip() for l in open(D / fn) if l.strip()}

sfari, ddg2p = load("sfari_all.set"), load("ddg2p_all.set")
sfari_hc, ddg2p_conf = load("sfari_hc.set"), load("ddg2p_conf.set")
shet = {}
for line in open(D / "shet_by_symbol.tsv"):
    p = line.split()
    if len(p) >= 2:
        try:
            shet[p[0]] = float(p[1])
        except ValueError:
            pass
constrained = {g for g, v in shet.items() if v >= SHET_MIN}
panel = sfari | ddg2p | constrained
print(f"cases {len(cases)}; SFARI {len(sfari)}, DDG2P {len(ddg2p)}, shet>={SHET_MIN} {len(constrained)}; "
      f"panel union {len(panel)}", file=sys.stderr)

# ---------------------------------------------------------------- panel CDS intervals
cds = defaultdict(list)
n_cds = 0
for line in open(D / "gencode_cds.bed"):
    p = line.rstrip("\n").split("\t")
    if len(p) < 4:
        continue
    for g in [x for x in p[3].split(",") if x in panel]:       # merged CDS carry "GENEA,GENEB" (T12)
        cds[p[0]].append((int(p[1]), int(p[2]), g))
        n_cds += 1
for c in cds:
    cds[c].sort()
print(f"panel CDS intervals {n_cds:,} on {len(cds)} contigs", file=sys.stderr)

bed = OUT.parent / "_panel_cds.bed"
with open(bed, "w") as fh:
    for c, ivs in cds.items():
        for a, b, g in ivs:
            fh.write(f"{c}\t{a}\t{b}\n")
subprocess.run(f"sort -k1,1 -k2,2n {bed} -o {bed}", shell=True, check=True)

def genes_hit(chrom, start, end):
    out = set()
    for a, b, g in cds.get(chrom, ()):
        if a < end and b > start:
            out.add(g)
        elif a >= end:
            break
    return out

# ---------------------------------------------------------------- stream the call set over panel regions
q = (f"{BCF} query -R {bed} "
     f"-f '%CHROM\\t%POS\\t%INFO/END\\t%INFO/SVTYPE\\t%INFO/SVLEN\\t%FILTER[\\t%SAMPLE|%GT|%GQ]\\n' {F}")
samples = subprocess.run(f"{BCF} query -l {F}", shell=True, capture_output=True, text=True).stdout.split()
print(f"samples {len(samples)}", file=sys.stderr)
CARRIER = {"0/1", "1/0", "1/1", "0|1", "1|0", "1|1"}
rows, n_rec, n_bnd, n_type, n_rare, n_gene = [], 0, 0, 0, 0, 0
proc = subprocess.Popen(q, shell=True, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
for line in proc.stdout:
    p = line.rstrip("\n").split("\t")
    if len(p) < 7:
        continue
    n_rec += 1
    chrom, pos, end, svtype, svlen, filt = p[0], p[1], p[2], p[3], p[4], p[5]
    if svtype == "BND":
        n_bnd += 1
        continue
    if svtype not in KEEP_TYPES:
        continue
    n_type += 1
    gt = {}
    for cell in p[6:]:
        b = cell.split("|")
        if len(b) >= 2:
            gt[b[0]] = dict(gt=b[1], gq=(int(b[2]) if len(b) > 2 and b[2].lstrip("-").isdigit() else None))
    carriers = [s for s, d in gt.items() if d["gt"] in CARRIER]
    ac = sum(2 if gt[s]["gt"] in ("1/1", "1|1") else 1 for s in carriers)
    if ac == 0 or ac > COHORT_AC_MAX:
        continue
    n_rare += 1
    s0 = int(pos)
    s1 = int(end) if end not in (".", "") else s0 + 1
    if s1 < s0:
        s0, s1 = s1, s0
    hits = genes_hit(chrom, s0, s1)
    if not hits:
        continue
    n_gene += 1
    for child in carriers:
        if child not in cases:
            continue
        pat, mat = ped[child]["pat"], ped[child]["mat"]
        fp = gt.get(pat, {}).get("gt") in CARRIER
        fm = gt.get(mat, {}).get("gt") in CARRIER
        if not (fp or fm):
            continue
        src = "father" if fp and not fm else ("mother" if fm and not fp else "both")
        for g in sorted(hits):
            rows.append(dict(
                gene=g, svtype=svtype, variant=f"{chrom}:{s0}-{s1}",
                length=(abs(int(svlen)) if svlen not in (".", "") and svlen.lstrip("-").isdigit() else s1 - s0),
                filter=filt, cohort_ac=ac, shet=f"{shet[g]:.3f}" if g in shet else "",
                sfari="hc" if g in sfari_hc else ("yes" if g in sfari else ""),
                ddg2p="conf" if g in ddg2p_conf else ("yes" if g in ddg2p else ""),
                constrained="yes" if g in constrained else "no",
                proband=child, family=ped[child]["fam"], transmitted_from=src,
                child_gt=gt[child]["gt"], child_gq=gt[child]["gq"]))
proc.wait()
print(f"\nrecords over panel CDS {n_rec:,}; breakends skipped {n_bnd:,}; typed {n_type:,} -> "
      f"rare in cohort {n_rare:,} -> hitting a panel gene's CDS {n_gene:,}", file=sys.stderr)

rows.sort(key=lambda r: (-float(r["shet"] or 0), r["gene"]))
if rows:
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
print(f"\n{len(rows)} carrier-gene rows, {len({r['variant'] for r in rows})} distinct SVs, "
      f"{len({r['gene'] for r in rows})} genes, {len({r['proband'] for r in rows})} probands")
H = "%-12s %-5s %-26s %9s %5s %7s %6s %6s %-12s %-8s %4s"
print("\n" + H % ("gene", "type", "variant", "length", "AC", "shet", "SFARI", "DDG2P", "proband", "from", "GQ"))
print("-" * 120)
for r in rows[:60]:
    print(H % (r["gene"], r["svtype"], r["variant"][:26], f"{r['length']:,}", r["cohort_ac"], r["shet"],
               r["sfari"], r["ddg2p"], r["proband"], r["transmitted_from"], r["child_gq"]))
print(f"\nwrote {OUT}")
