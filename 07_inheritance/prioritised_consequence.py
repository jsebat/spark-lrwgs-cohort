#!/usr/bin/env python
"""Assign an honest predicted functional consequence to every prioritised variant.

The prioritised table carried variant class but not predicted consequence, and the three classes are not
equivalent. Only the small variants are loss-of-function predictions. The others were selected for overlapping
coding sequence or for being length outliers, which is a much weaker statement, and the table should say so.

  small variants   already loss of function by construction: frameshift, stop gained, splice acceptor or donor,
                   start lost, and every one LOFTEE high confidence. Nothing further to compute.
  structural       selected for overlapping coding sequence, which is not the same as disrupting the gene.
                   Computed here per variant: what fraction of the gene's coding bases the event covers, and
                   whether it is contained inside the gene or spans past it. Deletions removing coding bases
                   are loss of function; duplications are dosage changes whose effect depends on whether the
                   whole gene is duplicated or only part of it; inversions depend on where the breakpoints
                   fall; insertions are frameshifts only if their length is not a multiple of three.
  repeat           NOT a loss-of-function prediction at all. A length outlier has no predicted consequence
                   unless the repeat is a known pathogenic locus. Computed here: whether the repeat overlaps
                   coding sequence or is intronic or untranslated.
"""
import csv, io, pathlib, sys
from collections import defaultdict

T = pathlib.Path("/expanse/lustre/projects/ddp195/jsebat/longread-autism/tiering")
D = T / "denovo_sv"
OUT = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else str(T / "prioritised_with_consequence.tsv"))

cds = defaultdict(list)
gene_cds_bp = defaultdict(int)
gene_span = {}
for line in open(D / "gencode_cds.bed"):
    p = line.rstrip("\n").split("\t")
    if len(p) < 4:
        continue
    c, a, b = p[0], int(p[1]), int(p[2])
    for g in p[3].split(","):                       # merged CDS carry "GENEA,GENEB" (T12)
        cds[(g, c)].append((a, b))
        gene_cds_bp[g] += b - a
        if g in gene_span and gene_span[g][0] == c:
            gene_span[g] = (c, min(gene_span[g][1], a), max(gene_span[g][2], b))
        else:
            gene_span.setdefault(g, (c, a, b))

def cds_overlap(gene, chrom, s, e):
    return sum(min(b, e) - max(a, s) for a, b in cds.get((gene, chrom), ()) if min(b, e) > max(a, s))

rows = list(csv.DictReader(io.open(T / "inherited_prioritised.tsv", encoding="utf-8"), delimiter="\t"))
out = []
for r in rows:
    g, d = r["gene"], r["detail"]
    csq, note = "", ""
    if r["klass"] == "small variant":
        csq = d.split()[0].replace("_variant", "").replace("_", " ")
        note = "LOFTEE high confidence"
    elif r["klass"] == "structural":
        parts = d.split()
        svtype, iv = parts[0], parts[-1]
        chrom, rng = iv.split(":")
        s, e = (int(x) for x in rng.split("-"))
        ov = cds_overlap(g, chrom, s, e)
        tot = gene_cds_bp.get(g, 0)
        frac = ov / tot if tot else 0
        sp = gene_span.get(g)
        spans_gene = bool(sp and s <= sp[1] and e >= sp[2])
        # the insertion LENGTH is the "N bp" token of detail (inherited_prioritise writes f"{svtype} {L:,} bp {variant}");
        # e - s is the 0/1 bp POS-END interval and called every coding insertion's frame from that (T8)
        length = int(parts[1].replace(",", "")) if len(parts) > 2 and parts[2] == "bp" and parts[1].replace(",", "").isdigit() else e - s
        if svtype == "DEL":
            csq = "whole-gene deletion" if spans_gene else "partial coding deletion"
            note = f"removes {ov:,} of {tot:,} coding bp ({frac:.0%})"
        elif svtype == "DUP":
            csq = "whole-gene duplication" if spans_gene else "partial duplication, uncertain"
            note = (f"covers {ov:,} of {tot:,} coding bp ({frac:.0%}); "
                    + ("dosage gain" if spans_gene else "may disrupt the transcript"))
        elif svtype == "INV":
            csq = "inversion spanning the gene" if spans_gene else "inversion with a breakpoint in the gene"
            note = f"{ov:,} coding bp within the inverted segment"
        elif svtype == "INS":
            csq = ("coding insertion, frameshift" if ov > 0 and length % 3 else
                   "coding insertion, in frame" if ov > 0 else "insertion, non-coding")
            note = f"{length} bp, {'not ' if length % 3 else ''}a multiple of three"
        else:
            csq, note = svtype.lower(), ""
    else:
        parts = d.split()
        locus = r.get("detail", "")
        chrom = None
        for src in (r.get("detail", ""), ):
            pass
        # locus is carried in the original repeat table; recover from the prioritised detail is not possible,
        # so fall back to the repeat table
        csq, note = "repeat length outlier", "no loss-of-function prediction; consequence unknown"
    out.append(dict(r, predicted_consequence=csq, consequence_note=note))

# repeats: recover coordinates from the repeat table to say coding vs non-coding
tr = {}
for r in csv.DictReader(io.open(T / "inherited_tr_constrained.tsv", encoding="utf-8"), delimiter="\t"):
    tr[(r["gene"], r["proband"], r["child_bp"])] = r
for r in out:
    if r["klass"] != "repeat":
        continue
    bp = r["detail"].split()[1]
    t = tr.get((r["gene"], r["proband"], bp))
    if not t:
        continue
    chrom, pos = t["locus"].split(":")
    pos = int(pos)
    ov = cds_overlap(r["gene"], chrom, pos, pos + int(bp))
    r["predicted_consequence"] = "coding repeat expansion" if ov > 0 else "non-coding repeat expansion"
    r["consequence_note"] = ("overlaps coding sequence" if ov > 0 else
                             "intronic or untranslated; no loss-of-function prediction")

with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(out)

sel = [r for r in out if r["tier"] in ("1", "2")]
print(f"{len(sel)} prioritised rows (tiers 1 and 2)\n")
from collections import Counter
for k in ("small variant", "structural", "repeat"):
    sub = [r for r in sel if r["klass"] == k]
    print(f"{k} ({len(sub)}):")
    for c, n in Counter(r["predicted_consequence"] for r in sub).most_common():
        print(f"    {n:3d}  {c}")
print("\n=== every tier 1 and 2 row with its predicted consequence ===")
H = "%-10s %-12s %-30s %-34s %-12s"
print(H % ("gene", "class", "predicted consequence", "note", "proband"))
print("-" * 104)
for r in sorted(sel, key=lambda r: (r["tier"], r["klass"], r["gene"])):
    print(H % (r["gene"], r["klass"][:12], r["predicted_consequence"][:30], r["consequence_note"][:34],
               r["proband"]))
print(f"\nwrote {OUT}")
