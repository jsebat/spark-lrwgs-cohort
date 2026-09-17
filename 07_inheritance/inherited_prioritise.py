#!/usr/bin/env python
"""Prioritise the inherited small-variant, structural and repeat findings into three tiers.

On inheritance from an unaffected parent. For a highly penetrant Mendelian syndrome this argues against a
diagnosis. For autism it does not: inherited risk variants carried by unaffected parents are the norm, because
penetrance is incomplete, expressivity is variable, and females carrying the same variant are affected less
often than males. So transmission by an unaffected parent is NOT treated as evidence against a contributory
role here. It is treated as evidence against a fully penetrant syndromic diagnosis, which is a different claim.

What does raise confidence
  co-segregation       the transmitting parent is themselves affected
  X-linked in a male   a maternally transmitted X variant in an affected son
  mechanism match      the variant class matches the established disease mechanism for that gene

Tiers
  1  plausibly contributory to neurodevelopmental disorder or autism
  2  incidental clinical findings, split into dominant actionable and recessive carrier status
  3  loss of function or expansion in a constrained gene with no established disease association

Quality gates are applied first and reported, because a tier assignment resting on an unreliable genotype is
worse than no assignment. Gene-disease curation is explicit per variant rather than derived from a list match.
"""
import csv, io, pathlib, sys

T = pathlib.Path("/expanse/lustre/projects/ddp195/jsebat/longread-autism/tiering")
OUT = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else str(T / "inherited_prioritised.tsv"))

ped = {}
for line in open(T / "_all.ped"):
    p = line.split()
    if len(p) >= 6:
        ped[p[1]] = dict(fam=p[0], pat=p[2], mat=p[3],
                         sex=("M" if p[4] == "1" else "F"), aff=(p[5] == "2"))

def parent_id(child, src):
    return ped[child]["pat"] if src == "father" else ped[child]["mat"]

def parent_affected(child, src):
    if src == "both":
        return ped[ped[child]["pat"]]["aff"] or ped[ped[child]["mat"]]["aff"]
    pid = parent_id(child, src)
    return ped.get(pid, {}).get("aff", False)

# ---- curated gene-disease layer. Each entry: tier hint, mechanism, and the reason, written out.
NDD = {   # established neurodevelopmental gene where loss or expansion is a plausible mechanism
    "DLGAP2": ("promoter VNTR expansion; high-confidence autism gene", "expansion"),
    "ZC4H2":  ("X-linked; Wieacker-Wolff syndrome", "expansion"),
    "EHMT1":  ("Kleefstra syndrome; haploinsufficiency", "expansion"),
    "KDM5A":  ("high-confidence autism gene and confirmed developmental disorder gene", "duplication"),
    "PCDH9":  ("high-confidence autism gene", "deletion"),
    "KDM2B":  ("SFARI category 1 autism gene; haploinsufficiency", "lof"),
    "SHANK1": ("high-confidence autism gene", "expansion"),
    "CACNA1C": ("high-confidence autism gene; Timothy syndrome is gain of function", "expansion"),
    "FOXP2":  ("speech and language disorder; haploinsufficiency", "expansion"),
    "BCL11A": ("high-confidence autism gene; intellectual disability", "expansion"),
    "BRSK2":  ("high-confidence autism gene", "expansion"),
    "TCF12":  ("craniosynostosis and developmental disorder", "expansion"),
    "GRID2":  ("cerebellar ataxia, recessive", "expansion"),
    "MACF1":  ("lissencephaly; confirmed developmental disorder gene", "insertion"),
    "RAD21":  ("Cornelia de Lange spectrum; confirmed", "inversion"),
    "STAG1":  ("intellectual disability; confirmed", "expansion"),
    "NCKAP1": ("high-confidence autism gene", "expansion"),
    "DIP2C":  ("high-confidence autism gene", "expansion"),
}
INCIDENTAL_DOM = {
    "FGF14": ("GAA expansion above ~250 units causes spinocerebellar ataxia type 27B, adult onset dominant",
              "expansion"),
    "DPYD":  ("loss of function predicts severe fluoropyrimidine toxicity; pharmacogenomic, actionable",
              "lof"),
    "COL5A1": ("classical Ehlers-Danlos syndrome by haploinsufficiency", "lof"),
}
INCIDENTAL_CARRIER = {
    "USH2A": "Usher syndrome type 2, recessive",
    "CEP290": "Joubert syndrome and Leber congenital amaurosis, recessive",
    "PCCB": "propionic acidemia, recessive",
    "AHI1": "Joubert syndrome, recessive",
    "DNAH10": "primary ciliary dyskinesia spectrum, recessive",
    "ASPM": "primary microcephaly, recessive",
}

rows = []

def add(cls, r, gene, proband, src, qual_ok, qual_note, extra):
    co = parent_affected(proband, src)
    xl = r.get("chrom", "").startswith("chrX") or str(r.get("locus", "")).startswith("chrX") \
         or str(r.get("variant", "")).startswith("chrX")
    male = ped[proband]["sex"] == "M"
    if gene in INCIDENTAL_DOM:
        tier, why = 2, "incidental, dominant: " + INCIDENTAL_DOM[gene][0]
    elif gene in INCIDENTAL_CARRIER:
        tier, why = 2, "incidental, carrier only: " + INCIDENTAL_CARRIER[gene]
    elif gene in NDD:
        base = NDD[gene][0]
        if co:
            tier, why = 1, f"{base}; co-segregates with an affected transmitting parent"
        elif xl and male and src == "mother":
            tier, why = 1, f"{base}; X-linked, maternally transmitted to an affected son"
        else:
            tier, why = 1, f"{base}; inherited from an unaffected parent, which is expected for autism risk alleles"
    else:
        tier, why = 3, "loss of function or expansion in a constrained gene, no established disease link"
    rows.append(dict(tier=tier, klass=cls, gene=gene, proband=proband, proband_sex=ped[proband]["sex"],
                     transmitted_from=src, parent_affected="yes" if co else "no",
                     quality="pass" if qual_ok else "FAIL", quality_note=qual_note,
                     rationale=why, **extra))

# ---- small variants
for r in csv.DictReader(io.open(T / "inherited_lof_constrained.tsv", encoding="utf-8"), delimiter="\t"):
    ok = int(r["child_dp"]) >= 10 and int(r["child_gq"]) >= 20
    add("small variant", dict(chrom=r["variant"].split(":")[0]), r["gene"], r["proband"],
        r["transmitted_from"], ok, f"DP {r['child_dp']}, GQ {r['child_gq']}",
        dict(detail=f"{r['consequence']} {r['variant']}", metric=f"LOEUF {r['loeuf']}", cohort=r["cohort_ac"]))
# ---- structural
for r in csv.DictReader(io.open(T / "inherited_sv_constrained.tsv", encoding="utf-8"), delimiter="\t"):
    L = int(r["length"])
    if L > 5_000_000:
        continue
    ok = not (r["child_gq"] and int(r["child_gq"]) < 20)
    add("structural", dict(variant=r["variant"]), r["gene"], r["proband"], r["transmitted_from"], ok,
        f"GQ {r['child_gq']}", dict(detail=f"{r['svtype']} {L:,} bp {r['variant']}",
                                    metric=f"shet {r['shet']}", cohort=r["cohort_ac"]))
# ---- repeats
for r in csv.DictReader(io.open(T / "inherited_tr_constrained.tsv", encoding="utf-8"), delimiter="\t"):
    sd = int(r["child_spanning_reads"])
    simple = len(r["motif"]) <= 2
    ok = sd >= 8 and not simple
    note = f"{sd} spanning reads" + (", simple motif" if simple else "")
    add("repeat", dict(locus=r["locus"]), r["gene"], r["proband"], r["transmitted_from"], ok, note,
        dict(detail=f"{r['motif']} {r['child_bp']} bp vs median {r['cohort_median_bp']} ({r['fold']}x)",
             metric=f"shet {r['shet']}", cohort=r["carriers"]))

rows.sort(key=lambda r: (r["tier"], r["quality"] == "FAIL", r["gene"]))
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)

for t in (1, 2, 3):
    sel = [r for r in rows if r["tier"] == t]
    npass = sum(1 for r in sel if r["quality"] == "pass")
    print(f"\n{'='*110}\nTIER {t}: {len(sel)} rows, {npass} passing quality")
    for r in sel:
        if t == 3 and r["quality"] != "pass":
            continue
        print(f"  [{r['quality']:4}] {r['gene']:<10} {r['klass']:<14} {r['proband']:<12} "
              f"{r['proband_sex']} from {r['transmitted_from']:<7} "
              f"{'AFFECTED parent' if r['parent_affected']=='yes' else '':<16} {r['detail'][:46]}")
        if t <= 2:
            print(f"           {r['rationale'][:104]}")
print(f"\nwrote {OUT}")
