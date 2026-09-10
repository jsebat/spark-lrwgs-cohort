#!/usr/bin/env python
"""Fetch canonical (Ensembl canonical / MANE) transcript models for the flag genes -> resources/flag_genes/.

Writes, for each gene in config.flag_genes:
  resources/flag_genes/<GENE>.<transcript>.json   exons, CDS coordinates, CDS sequence, strand (Ensembl REST, GRCh38)
and one combined BED:
  resources/flag_genes/flag_genes.cds_exons.hg38.bed   chrom start end gene transcript exon_rank   (0-based half-open, CDS part of each exon)
plus SOURCE.md with the Ensembl release and date. Public data only; safe to commit.
"""
import argparse
import datetime as dt
import json
import pathlib
import urllib.request

import yaml


def get(rest, path):
    req = urllib.request.Request(rest + path, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--out", default="resources/flag_genes")
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config))
    rest = cfg["ensembl"]["rest_url"]
    out = pathlib.Path(a.out); out.mkdir(parents=True, exist_ok=True)
    rel = get(rest, "/info/data")["releases"]
    bed, src = [], [f"# flag gene models\n\nEnsembl REST {rest}, release {rel}, GRCh38, fetched {dt.date.today()}.\n",
                    "| gene | transcript | canonical flag | RefSeq | strand | CDS len | n exons |", "|---|---|---|---|---|---|---|"]
    for gene in cfg["flag_genes"]:
        g = get(rest, f"/lookup/symbol/homo_sapiens/{gene}?expand=1")
        canon = [t for t in g["Transcript"] if t.get("is_canonical") == 1] or [t for t in g["Transcript"] if t["biotype"] == "protein_coding"]
        t = canon[0]
        tid = t["id"]
        tr = get(rest, f"/lookup/id/{tid}?expand=1")
        cds = get(rest, f"/sequence/id/{tid}?type=cds")["seq"]
        xr = [x["display_id"] for x in get(rest, f"/xrefs/id/{tid}?external_db=RefSeq_mRNA")]
        strand = tr["strand"]
        cds_lo, cds_hi = tr["Translation"]["start"], tr["Translation"]["end"]
        exons = sorted(tr["Exon"], key=lambda e: e["start"], reverse=(strand == -1))
        model = dict(gene=gene, transcript=f"{tid}.{tr.get('version')}", refseq=xr, chrom="chr" + tr["seq_region_name"], strand=strand,
                     gene_start=tr["start"], gene_end=tr["end"], cds_genomic_lo=cds_lo, cds_genomic_hi=cds_hi, cds_seq=cds,
                     exons=[dict(rank=i + 1, id=e["id"], start=e["start"], end=e["end"]) for i, e in enumerate(exons)],
                     ensembl_release=rel)
        (out / f"{gene}.{tid}.json").write_text(json.dumps(model, indent=1))
        for i, e in enumerate(exons):
            lo, hi = max(e["start"], cds_lo), min(e["end"], cds_hi)
            if lo <= hi:
                bed.append((model["chrom"], lo - 1, hi, gene, tid, i + 1))
        src.append(f"| {gene} | {model['transcript']} | {t.get('is_canonical')} | {', '.join(xr) or '-'} | {strand} | {len(cds)} | {len(exons)} |")
        print(f"{gene}: {model['transcript']} {model['chrom']}:{tr['start']}-{tr['end']} strand {strand}, CDS {len(cds)} bp, {len(exons)} exons")
    bed.sort()
    (out / "flag_genes.cds_exons.hg38.bed").write_text("".join("\t".join(map(str, r)) + "\n" for r in bed))
    (out / "SOURCE.md").write_text("\n".join(src) + "\n")


if __name__ == "__main__":
    main()
