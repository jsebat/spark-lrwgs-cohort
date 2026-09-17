#!/usr/bin/env python
"""Rare inherited tandem-repeat expansions in panel genes, transmitted to affected children.

Third of the three inherited tables, alongside small variants and structural variants. A repeat expansion is
not a presence/absence call, so "rare" and "disruptive" have to be defined against the cohort itself.

Definitions, applied per locus
  outlier      the child's longer allele exceeds the cohort 99th percentile for that locus AND exceeds the
               cohort median by at least EXCESS_BP base pairs and by at least EXCESS_FOLD times. Both an
               absolute and a relative floor are needed, because a locus with a median of 20 bp and one with a
               median of 2,000 bp do not behave alike.
  rare         at most COHORT_MAX_CARRIERS individuals in the cohort carry an allele that long. This is the
               analogue of cohort allele count in the other two tables.
  transmitted  a parent carries an allele at least PARENT_FRAC of the child's, so the long allele was passed
               on rather than arising de novo or by genotyping error.
  in a gene    the locus lies within the span of a panel gene, taken as first to last coding base, so
               intronic and untranslated repeats are retained. That matters: the known pathogenic expansions
               are mostly not in coding sequence.

Spanning-read depth is carried through, because a repeat genotype supported by few spanning reads is the
dominant failure mode for long alleles.
"""
import csv, pathlib, statistics as st, subprocess, sys
from collections import defaultdict

V = "/expanse/projects/sebat1/longread_cohort_2026/freeze1/freeze1.cohort.trgt.vcf.gz"
T = pathlib.Path("/expanse/lustre/projects/ddp195/jsebat/longread-autism/tiering")
D = T / "denovo_sv"
OUT = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else str(T / "inherited_tr_constrained.tsv"))
BCF = str(pathlib.Path.home() / "micromamba/envs/dnmt3a-py/bin/bcftools")
EXCESS_BP, EXCESS_FOLD = 50, 1.5
COHORT_MAX_CARRIERS = 4
PARENT_FRAC = 0.80
SHET_MIN = 0.10
MIN_SD = 3

ped = {}
for line in open(T / "_all.ped"):
    p = line.split()
    if len(p) >= 6:
        ped[p[1]] = dict(fam=p[0], pat=p[2], mat=p[3], aff=p[5])
cases = {s for s, v in ped.items() if v["pat"] != "0" and v["mat"] != "0" and v["aff"] == "2"}

load = lambda fn: {l.strip() for l in open(D / fn) if l.strip()}
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

span = {}
for line in open(D / "gencode_cds.bed"):
    p = line.rstrip("\n").split("\t")
    if len(p) < 4 or p[3] not in panel:
        continue
    g, c, a, b = p[3], p[0], int(p[1]), int(p[2])
    if g in span:
        s = span[g]
        if s[0] == c:
            span[g] = (c, min(s[1], a), max(s[2], b))
    else:
        span[g] = (c, a, b)
print(f"panel genes {len(panel)}, with a span {len(span)}; cases {len(cases)}", file=sys.stderr)

bed = OUT.parent / "_panel_span.bed"
by_chrom = defaultdict(list)
with open(bed, "w") as fh:
    for g, (c, a, b) in span.items():
        fh.write(f"{c}\t{a}\t{b}\n")
        by_chrom[c].append((a, b, g))
for c in by_chrom:
    by_chrom[c].sort()
subprocess.run(f"sort -k1,1 -k2,2n {bed} -o {bed}", shell=True, check=True)

def genes_at(chrom, pos):
    return [g for a, b, g in by_chrom.get(chrom, ()) if a <= pos <= b]

samples = subprocess.run(f"{BCF} query -l {V}", shell=True, capture_output=True, text=True).stdout.split()
q = f"{BCF} query -R {bed} -f '%CHROM\\t%POS\\t%INFO/TRID\\t%INFO/MOTIFS[\\t%AL|%SD]\\n' {V}"
proc = subprocess.Popen(q, shell=True, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)

rows, n_loci, n_out, n_rare = [], 0, 0, 0
for line in proc.stdout:
    p = line.rstrip("\n").split("\t")
    if len(p) < 5:
        continue
    n_loci += 1
    chrom, pos, trid, motifs = p[0], int(p[1]), p[2], p[3]
    longest, sd = {}, {}
    for s, cell in zip(samples, p[4:]):
        a, _, d = cell.partition("|")
        als = [int(x) for x in a.split(",") if x.lstrip("-").isdigit()]
        sds = [int(x) for x in d.split(",") if x.lstrip("-").isdigit()]
        if als:
            longest[s] = max(als)
            sd[s] = max(sds) if sds else 0
    if len(longest) < 60:
        continue
    vals = sorted(longest.values())
    med = st.median(vals)
    p99 = vals[min(len(vals) - 1, int(0.99 * len(vals)))]
    for child in cases:
        L = longest.get(child)
        if L is None or sd.get(child, 0) < MIN_SD:
            continue
        if not (L > p99 and L - med >= EXCESS_BP and (med == 0 or L >= EXCESS_FOLD * med)):
            continue
        n_out += 1
        carriers = [s for s, v in longest.items() if v >= L * 0.95]
        if len(carriers) > COHORT_MAX_CARRIERS:
            continue
        n_rare += 1
        pat, mat = ped[child]["pat"], ped[child]["mat"]
        pl, ml = longest.get(pat), longest.get(mat)
        fp = pl is not None and pl >= PARENT_FRAC * L
        fm = ml is not None and ml >= PARENT_FRAC * L
        if not (fp or fm):
            continue
        src = "father" if fp and not fm else ("mother" if fm and not fp else "both")
        for g in genes_at(chrom, pos):
            rows.append(dict(
                gene=g, locus=f"{chrom}:{pos}", trid=trid, motif=motifs.split(",")[0][:24],
                child_bp=L, cohort_median_bp=int(med), cohort_p99_bp=int(p99),
                fold=round(L / med, 2) if med else "", carriers=len(carriers),
                shet=f"{shet[g]:.3f}" if g in shet else "",
                sfari="hc" if g in sfari_hc else ("yes" if g in sfari else ""),
                ddg2p="conf" if g in ddg2p_conf else ("yes" if g in ddg2p else ""),
                proband=child, family=ped[child]["fam"], transmitted_from=src,
                parent_bp=(pl if fp else ml), child_spanning_reads=sd.get(child, 0)))
proc.wait()
print(f"\nloci in panel gene spans {n_loci:,}; child outlier calls {n_out:,} -> rare in cohort {n_rare:,}",
      file=sys.stderr)

rows.sort(key=lambda r: (-float(r["shet"] or 0), -r["child_bp"]))
if rows:
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
print(f"\n{len(rows)} carrier-gene rows, {len({r['locus'] for r in rows})} loci, "
      f"{len({r['gene'] for r in rows})} genes, {len({r['proband'] for r in rows})} probands")
H = "%-12s %-20s %-14s %8s %8s %6s %5s %6s %6s %-12s %-8s %4s"
print("\n" + H % ("gene", "locus", "motif", "child bp", "median", "fold", "carr", "shet", "SFARI",
                  "proband", "from", "SD"))
print("-" * 130)
for r in rows[:50]:
    print(H % (r["gene"], r["locus"], r["motif"][:14], f"{r['child_bp']:,}", f"{r['cohort_median_bp']:,}",
               r["fold"], r["carriers"], r["shet"], r["sfari"], r["proband"], r["transmitted_from"],
               r["child_spanning_reads"]))
print(f"\nwrote {OUT}")
