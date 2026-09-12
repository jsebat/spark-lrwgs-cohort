#!/usr/bin/env python
"""Rare inherited loss-of-function variants transmitted to affected children, in genes that matter.

The cohort report states that no inherited variant met reportable thresholds and that rare loss-of-function
burden was null. Neither statement lets a reader see what was actually there. This enumerates it.

Inclusion, applied in this order
  1. Loss of function, from the tiering call set: frameshift, stop gained, splice acceptor or donor, start
     lost, with LOFTEE annotation retained.
  2. Rare in the population: gnomAD allele frequency below 1e-4, or absent from gnomAD entirely.
  3. Rare in this cohort: cohort allele count at or below COHORT_AC_MAX, counted across all 105 genomes.
     This is the filter that removes the recurrent artifacts, which are the dominant failure mode here.
  4. Transmitted to a case: carried by an affected child AND by at least one of that child's parents, with
     the parent genotype actually observed rather than assumed.
  5. In a gene of interest: SFARI autism gene, or loss-of-function intolerant by gnomAD constraint.

Reported per variant and per carrier, with the transmitting parent named, so the table can be read as a list
of candidates rather than as a summary.
"""
import csv, gzip, io, math, pathlib, sys
import numpy as np

T = pathlib.Path("/expanse/lustre/projects/ddp195/jsebat/longread-autism/tiering")
G = pathlib.Path.home() / "genelists"
OUT = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else str(T / "inherited_lof_constrained.tsv"))
AF_MAX = 1e-4
COHORT_AC_MAX = 4
LOEUF_MAX = 0.35

# ---------------------------------------------------------------- pedigree
ped = {}
for line in open(T / "_all.ped"):
    p = line.split()
    if len(p) < 6:
        continue
    ped[p[1]] = dict(fam=p[0], iid=p[1], pat=p[2], mat=p[3], sex=p[4], aff=p[5])
kids = [s for s, v in ped.items() if v["pat"] != "0" and v["mat"] != "0"]
cases = [s for s in kids if ped[s]["aff"] == "2"]
print(f"pedigree {len(ped)} individuals, {len(kids)} offspring, {len(cases)} affected offspring", file=sys.stderr)

# ---------------------------------------------------------------- gene sets
sfari = {}
with open(G / "sfari_genes.csv", newline="", encoding="utf-8", errors="replace") as fh:
    for r in csv.DictReader(fh):
        sym = (r.get("gene-symbol") or "").strip()
        if sym:
            sfari[sym] = (r.get("gene-score") or "").strip() or "S"
loeuf = {}
with open(G / "gnomad_constraint.tsv") as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    for r in rd:
        if str(r.get("mane_select", "")).lower() not in ("true", "1"):
            continue
        try:
            v = float(r["lof.oe_ci.upper"])
        except (KeyError, ValueError):
            continue
        g = r["gene"]
        if g not in loeuf or v < loeuf[g]:
            loeuf[g] = v
intolerant = {g for g, v in loeuf.items() if v < LOEUF_MAX}
print(f"SFARI genes {len(sfari):,}; LoF-intolerant (LOEUF<{LOEUF_MAX}) {len(intolerant):,}", file=sys.stderr)

# ---------------------------------------------------------------- genotypes, keyed by site
gt = {}
with open(T / "_gt.tsv") as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < 3:
            continue
        key = (p[0], p[1])
        d = {}
        for cell in p[2:]:
            bits = cell.split("|")
            if len(bits) >= 2:
                d[bits[0]] = dict(gt=bits[1],
                                  dp=int(bits[2]) if len(bits) > 2 and bits[2].isdigit() else None,
                                  gq=int(bits[3]) if len(bits) > 3 and bits[3].isdigit() else None)
        gt[key] = d
print(f"genotype rows {len(gt):,}", file=sys.stderr)

CARRIER = {"0/1", "1/0", "1/1", "0|1", "1|0", "1|1"}
LOF_CSQ = ("frameshift", "stop_gained", "splice_acceptor", "splice_donor", "start_lost")

rows = []
n_lof = n_rare = n_cohort = n_gene = 0
with open(T / "lof_tiered.rare.tsv") as fh:
    for v in csv.DictReader(fh, delimiter="\t"):
        csq = v["consequence"]
        if not any(c in csq for c in LOF_CSQ):
            continue
        n_lof += 1
        af_raw = v.get("gnomad_af", ".")
        try:
            af = float(af_raw)
        except (TypeError, ValueError):
            af = 0.0                      # '.' means absent from gnomAD
        if af > AF_MAX:
            continue
        n_rare += 1
        key = (v["chrom"], v["pos"])
        g = gt.get(key, {})
        carriers = [s for s, d in g.items() if d["gt"] in CARRIER]
        ac = sum(2 if g[s]["gt"] in ("1/1", "1|1") else 1 for s in carriers)
        if ac > COHORT_AC_MAX or ac == 0:
            continue
        n_cohort += 1
        sym = v["symbol"]
        in_sfari = sym in sfari
        in_int = sym in intolerant
        if not (in_sfari or in_int):
            continue
        n_gene += 1
        for child in carriers:
            if child not in cases:
                continue
            pat, mat = ped[child]["pat"], ped[child]["mat"]
            from_pat = g.get(pat, {}).get("gt") in CARRIER
            from_mat = g.get(mat, {}).get("gt") in CARRIER
            if not (from_pat or from_mat):
                continue                  # not transmitted, or parent missing
            src = "father" if from_pat and not from_mat else ("mother" if from_mat and not from_pat else "both")
            rows.append(dict(
                gene=sym, variant=f"{v['chrom']}:{v['pos']}:{v['ref']}>{v['alt']}",
                consequence=csq.split("&")[0], loftee=v.get("lof", "."), tier=v.get("lof_tier", ""),
                gnomad_af=("absent" if af_raw in (".", "", None) else f"{af:.2e}"),
                cohort_ac=ac, s_het=v.get("s_het", ""),
                loeuf=f"{loeuf[sym]:.2f}" if sym in loeuf else "",
                sfari_score=sfari.get(sym, ""), lof_intolerant="yes" if in_int else "no",
                proband=child, family=ped[child]["fam"], transmitted_from=src,
                child_gt=g[child]["gt"], child_dp=g[child]["dp"], child_gq=g[child]["gq"]))

print(f"\nfiltering: {n_lof:,} LoF -> {n_rare:,} rare in gnomAD -> {n_cohort:,} rare in cohort "
      f"(AC<={COHORT_AC_MAX}) -> {n_gene:,} in a gene of interest", file=sys.stderr)

rows.sort(key=lambda r: (r["tier"], float(r["loeuf"] or 9), r["gene"]))
cols = list(rows[0].keys()) if rows else []
with open(OUT, "w", newline="") as fh:
    if rows:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
print(f"\n{len(rows)} carrier rows, {len({r['variant'] for r in rows})} distinct variants, "
      f"{len({r['gene'] for r in rows})} genes, {len({r['proband'] for r in rows})} probands")

H = "%-12s %-30s %-20s %-6s %-7s %-10s %4s %6s %6s %-12s %-8s"
print("\n" + H % ("gene", "variant", "consequence", "LOFTEE", "tier", "gnomAD AF", "AC", "LOEUF", "SFARI",
                  "proband", "from"))
print("-" * 136)
for r in rows:
    print(H % (r["gene"], r["variant"][:30], r["consequence"][:20], r["loftee"], r["tier"],
               r["gnomad_af"], r["cohort_ac"], r["loeuf"], r["sfari_score"], r["proband"],
               r["transmitted_from"]))
print(f"\nwrote {OUT}")
