#!/usr/bin/env python
"""Step 01 - annotate the proband DNMT3A deletion against the canonical transcript.

Data source: Ensembl REST (GRCh38). Every response is cached verbatim under
config.ensembl.cache_dir for provenance; the Ensembl release is recorded in the output.

Reports, per coordinate convention (both if config says "unknown"):
  * exons overlapped (transcript numbering, 5'->3') and exonic bp removed from each
  * total exonic bp removed -> in-frame vs frameshift
  * CDS / codon range removed; predicted protein under two splicing scenarios:
      (a) partial exons fuse (deletion removes the intron and both flanking splice sites)
      (b) all touched exons skipped
    with first premature stop, protein length, and NMD prediction (>50 nt rule)
  * protein domains (Pfam / SMART / PROSITE / CDD) overlapping the lost codons
  * distance to codon R882 and the domain containing it
"""
import argparse
import datetime as dt
import hashlib
import json
import pathlib
import urllib.request

import yaml

STOPS = {"TAA", "TAG", "TGA"}
R882_CODON = 882
DOMAIN_DBS = ("Pfam", "Smart", "PROSITE_profiles", "CDD")


def fetch(rest, path, cache_dir):
    """GET JSON from Ensembl REST with an on-disk cache of the verbatim response."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    slug = path.strip("/").replace("/", "_").replace("?", "_").replace("=", "-")[:80]
    fn = cache_dir / f"{hashlib.sha1(path.encode()).hexdigest()[:12]}_{slug}.json"
    if fn.exists():
        return json.loads(fn.read_text())
    req = urllib.request.Request(rest + path, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        raw = r.read()
    fn.write_bytes(raw)
    return json.loads(raw)


class Transcript:
    """Minus-strand transcript with genomic <-> transcript <-> CDS coordinate maps (all 1-based)."""

    def __init__(self, tr, cds):
        assert tr["strand"] == -1, "script assumes a minus-strand transcript (DNMT3A)"
        self.tr = tr
        self.cds = cds
        self.exons = sorted(tr["Exon"], key=lambda e: -e["start"])  # 5'->3' == descending genomic
        self.cds_start_g = tr["Translation"]["end"]  # A of ATG has the highest coordinate
        self.cds_end_g = tr["Translation"]["start"]
        self.tcds_start = self.g2t(self.cds_start_g)
        self.tcds_end = self.g2t(self.cds_end_g)
        assert self.tcds_end - self.tcds_start + 1 == len(cds), "CDS length mismatch vs exon structure"

    def g2t(self, pos):
        t = 0
        for e in self.exons:
            if e["start"] <= pos <= e["end"]:
                return t + (e["end"] - pos) + 1
            t += e["end"] - e["start"] + 1
        return None

    def t2g(self, t):
        acc = 0
        for e in self.exons:
            n = e["end"] - e["start"] + 1
            if t <= acc + n:
                return e["end"] - (t - acc - 1)
            acc += n
        return None

    def g2c(self, pos):
        t = self.g2t(pos)
        return None if t is None else t - self.tcds_start + 1

    @staticmethod
    def codon_of_c(c):
        return (c - 1) // 3 + 1

    def exon_table(self):
        rows = []
        for i, e in enumerate(self.exons, 1):
            c1 = self.g2t(e["end"]) - self.tcds_start + 1
            c2 = self.g2t(e["start"]) - self.tcds_start + 1
            coding = not (c2 < 1 or c1 > len(self.cds))
            rows.append(dict(
                tx_exon=i, ensembl_exon=e["id"], start=e["start"], end=e["end"],
                length=e["end"] - e["start"] + 1, cds_from=c1, cds_to=c2, coding=coding,
                codon_from=self.codon_of_c(max(c1, 1)) if coding else None,
                codon_to=self.codon_of_c(min(c2, len(self.cds))) if coding else None,
                phase_start=(c1 - 1) % 3 if coding else None,
            ))
        return rows

    def last_junction_c(self):
        """CDS coordinate of the first base of the last exon (the last exon-exon junction)."""
        return self.g2t(self.exons[-1]["end"]) - self.tcds_start + 1


def translate_scan(new_cds):
    codons = [new_cds[i:i + 3] for i in range(0, len(new_cds) - 2, 3)]
    for i, c in enumerate(codons):
        if c in STOPS:
            return i  # protein length in aa
    return None


def scenario(T, removed_c, label):
    """Predict the protein after removing a set of CDS positions."""
    removed_c = set(removed_c)
    new = "".join(T.cds[i - 1] for i in range(1, len(T.cds) + 1) if i not in removed_c)
    plen = translate_scan(new)
    wt_len = len(T.cds) // 3 - 1
    out = dict(scenario=label, cds_bp_removed=len(removed_c),
               frame="in-frame" if len(removed_c) % 3 == 0 else "frameshift",
               new_cds_len=len(new), protein_length=plen, wt_protein_length=wt_len)
    if plen is None:
        out["note"] = "no stop codon found in new reading frame (read-through into 3' UTR not modelled)"
        return out
    stop_c_new = plen * 3 + 1
    junction_new = T.last_junction_c() - sum(1 for c in removed_c if c < T.last_junction_c())
    dist = junction_new - stop_c_new
    premature = stop_c_new != len(new) - 2  # stop is not the native terminal codon
    out.update(stop_codon_new_cds_pos=stop_c_new, last_junction_new_cds_pos=junction_new,
               premature_stop=premature, stop_upstream_of_last_junction_nt=dist,
               nmd_predicted=bool(premature and dist > 50))
    first_codon = T.codon_of_c(min(removed_c))
    if premature and len(removed_c) % 3 != 0:
        out.update(first_altered_codon=first_codon, aberrant_residues=plen - (first_codon - 1))
    elif len(removed_c) % 3 == 0:
        out.update(residues_deleted=len(removed_c) // 3,
                   deleted_codons=f"{first_codon}-{T.codon_of_c(max(removed_c))}")
    return out


def domains_hit(pfeat, a, b):
    hits = [dict(db=f["type"], id=f["id"], description=f.get("description"), start=f["start"], end=f["end"])
            for f in pfeat if f.get("type") in DOMAIN_DBS and not (f["end"] < a or f["start"] > b)]
    return sorted(hits, key=lambda x: (x["db"], x["start"]))


def analyse(T, pfeat, chrom, a, b, convention):
    """a, b: 1-based inclusive genomic interval of deleted bases."""
    ex_rows = T.exon_table()
    removed_c = sorted(c for c in (T.g2c(p) for p in range(a, b + 1)) if c is not None and 1 <= c <= len(T.cds))
    per_exon = []
    for r in ex_rows:
        ov = max(0, min(b, r["end"]) - max(a, r["start"]) + 1)
        if not ov:
            continue
        if ov == r["length"]:
            side = "whole exon"
        elif a <= r["start"]:
            side = "3' (donor) end"      # low-coordinate side == transcript 3' on minus strand
        elif b >= r["end"]:
            side = "5' (acceptor) end"
        else:
            side = "internal"
        per_exon.append(dict(tx_exon=r["tx_exon"], ensembl_exon=r["ensembl_exon"], exon_start=r["start"],
                             exon_end=r["end"], exon_length=r["length"], bp_deleted=ov, side_lost=side))
    res = dict(convention=convention, deleted_interval_1based=f"{chrom}:{a}-{b}", deleted_bp=b - a + 1,
               exons_overlapped=per_exon, exonic_bp_deleted=len(removed_c),
               intronic_bp_deleted=(b - a + 1) - len(removed_c),
               frame="in-frame" if len(removed_c) % 3 == 0 else "frameshift",
               cds_positions_deleted=f"c.{min(removed_c)}_{max(removed_c)}" if removed_c else None,
               codons_affected=f"p.{T.codon_of_c(min(removed_c))}-{T.codon_of_c(max(removed_c))}" if removed_c else None)
    if not per_exon:
        res["splice_sites_lost"] = "none (intronic-only deletion)"
        return res
    touched = sorted({e["tx_exon"] for e in per_exon})
    if len(touched) == 2 and touched[1] == touched[0] + 1:
        res["splice_sites_lost"] = (f"donor of exon {touched[0]}, acceptor of exon {touched[1]}, "
                                    f"and the entire intervening intron")
    else:
        res["splice_sites_lost"] = "see exons_overlapped"
    res["scenario_fused_partial_exons"] = scenario(T, removed_c, "partial exons fuse; remaining splice sites used")
    skip_c = set()
    for r in ex_rows:
        if r["tx_exon"] in touched and r["coding"]:
            skip_c |= set(range(max(r["cds_from"], 1), min(r["cds_to"], len(T.cds)) + 1))
    res["scenario_exon_skipping"] = scenario(T, skip_c, "exon(s) %s skipped" % ",".join(map(str, touched)))
    c1, c2 = T.codon_of_c(min(removed_c)), T.codon_of_c(max(removed_c))
    res["domains_overlapping_deleted_codons"] = domains_hit(pfeat, c1, c2)
    res["codons_from_deletion_to_R882"] = R882_CODON - c1
    return res


def fmt_domains(ds):
    return "; ".join(f"{d['db']} {d['id']} {d['description']} (aa {d['start']}-{d['end']})" for d in ds) or "none annotated"


def render_md(out):
    m = out["meta"]
    v = m["variant_as_configured"]
    L = ["# Step 01 - proband DNMT3A deletion annotation", ""]
    L.append(f"Generated {m['generated']}. Ensembl REST release {m['ensembl_release']}, assembly {m['assembly']}. "
             f"Transcript {m['transcript']} ({m['transcript_name']}; RefSeq {', '.join(m['refseq_xrefs']) or 'n/a'}), "
             f"protein {m['protein']} ({m['protein_length_aa']} aa, {m['n_exons']} exons, minus strand, gene span {m['gene_span']}).")
    L.append("")
    L.append(f"**Variant as configured:** {v['chrom']}:{v['start']}-{v['end']}, {v['zygosity']}, "
             f"coordinate convention `{v['coordinate_convention']}`, inheritance `{v.get('inheritance', 'unknown')}`.")
    L += ["", "## Exon overlap and reading frame", ""]
    for r in out["analyses"]:
        L.append(f"### Convention `{r['convention']}` -> deleted bases {r['deleted_interval_1based']} ({r['deleted_bp']} bp)")
        L.append("")
        L.append("| tx exon | Ensembl exon | exon (hg38, 1-based) | exon len | bp deleted | part lost |")
        L.append("|---|---|---|---|---|---|")
        for e in r["exons_overlapped"]:
            L.append(f"| {e['tx_exon']} | {e['ensembl_exon']} | {e['exon_start']}-{e['exon_end']} | {e['exon_length']} | {e['bp_deleted']} | {e['side_lost']} |")
        L.append("")
        L.append(f"- Exonic bp deleted: **{r['exonic_bp_deleted']}** ({r['intronic_bp_deleted']} bp intronic) -> **{r['frame'].upper()}** at the DNA level")
        L.append(f"- CDS positions removed: {r['cds_positions_deleted']}; codons {r['codons_affected']}")
        L.append(f"- Splice sites lost: {r.get('splice_sites_lost')}")
        for key in ("scenario_fused_partial_exons", "scenario_exon_skipping"):
            s = r.get(key)
            if not s:
                continue
            L.append(f"- **Scenario: {s['scenario']}** -> {s['frame']}, {s['cds_bp_removed']} CDS bp removed, "
                     f"predicted protein {s['protein_length']} aa (WT {s['wt_protein_length']}).")
            if s.get("premature_stop"):
                nmd = "predicted (>50 nt rule)" if s["nmd_predicted"] else "NOT predicted (escapes the 50-nt rule)"
                L.append(f"    - premature stop {s['stop_upstream_of_last_junction_nt']} nt upstream of the last exon-exon junction -> NMD {nmd}; "
                         f"aberrant residues from codon {s.get('first_altered_codon')}: {s.get('aberrant_residues')}")
            elif s["frame"] == "in-frame":
                L.append(f"    - in-frame loss of {s['residues_deleted']} residues (codons {s['deleted_codons']}); no premature stop")
        L.append("- Protein domains overlapping the deleted codons: " + fmt_domains(r.get("domains_overlapping_deleted_codons", [])))
        L.append(f"- Codons from first deleted codon to R882: {r.get('codons_from_deletion_to_R882')}")
        L.append("")
    R = m["R882"]
    L += ["## R882 reference point", "",
          f"Codon 882 = `{R['codon']}` at {R['genomic']} (tx exon {R['tx_exon']}); Ensembl domain hits at this residue: {fmt_domains(R['domains'])}. "
          "(The Ensembl Pfam PF00145 hit stops at aa 766; R882 sits in the C-terminal part of the catalytic MTase region, downstream of every hit above.)", ""]
    L += ["## Protein domain annotations (Ensembl protein features)", "", "| db | id | description | aa start | aa end |", "|---|---|---|---|---|"]
    for d in m["domains_all"]:
        L.append(f"| {d['db']} | {d['id']} | {d['description']} | {d['start']} | {d['end']} |")
    L += ["", "## Full exon table (transcript order, 5'->3')", "",
          "| tx exon | Ensembl exon | hg38 start | hg38 end | len | CDS from | CDS to | codons | phase |",
          "|---|---|---|---|---|---|---|---|---|"]
    for e in out["exons"]:
        cod = f"{e['codon_from']}-{e['codon_to']}" if e["coding"] else "UTR"
        ph = e["phase_start"] if e["coding"] else ""
        L.append(f"| {e['tx_exon']} | {e['ensembl_exon']} | {e['start']} | {e['end']} | {e['length']} | {e['cds_from']} | {e['cds_to']} | {cod} | {ph} |")
    L += ["", "## Checkpoint status (PLAN step 01)", ""]
    any_exonic = any(r["exonic_bp_deleted"] > 0 for r in out["analyses"])
    L.append("Exon-disrupting: **%s**. Frame and NMD calls are given per coordinate convention above; "
             "interpretation against the LOF hypothesis is recorded in the step log after review with Jonathan." % ("YES" if any_exonic else "NO (intronic-only)"))
    L.append("")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", required=True)
    ap.add_argument("--out-md", required=True)
    ap.add_argument("--out-json", required=True)
    args = ap.parse_args()
    cfg = yaml.safe_load(open(args.config))
    v = cfg["proband_variant"]
    rest = cfg["ensembl"]["rest_url"]
    cache = pathlib.Path(cfg["ensembl"].get("cache_dir", "results/01_variant/ensembl_cache"))
    tid = v["transcript"]

    tr = fetch(rest, f"/lookup/id/{tid}?expand=1", cache)
    cds = fetch(rest, f"/sequence/id/{tid}?type=cds", cache)["seq"]
    pid = tr["Translation"]["id"]
    pfeat = fetch(rest, f"/overlap/translation/{pid}?feature=protein_feature", cache)
    rel = fetch(rest, "/info/data", cache)
    xr = [x["display_id"] for x in fetch(rest, f"/xrefs/id/{tid}?external_db=RefSeq_mRNA", cache)]
    T = Transcript(tr, cds)

    chrom, s, e = v["chrom"], int(v["start"]), int(v["end"])
    conv = str(v.get("coordinate_convention", "unknown")).lower()
    convs = {"1-based-inclusive": (s, e), "bed": (s + 1, e)}
    todo = convs if conv == "unknown" else {conv: convs[conv]}
    analyses = [analyse(T, pfeat, chrom, a, b, name) for name, (a, b) in todo.items()]

    c882 = T.tcds_start + (R882_CODON - 1) * 3
    r882_g = sorted([T.t2g(c882), T.t2g(c882 + 2)])
    r882_exon = next(r["tx_exon"] for r in T.exon_table() if r["start"] <= r882_g[0] <= r["end"])
    meta = dict(
        generated=dt.datetime.now().isoformat(timespec="seconds"), ensembl_release=rel.get("releases"),
        assembly=tr.get("assembly_name"), transcript=f"{tid}.{tr.get('version')}", transcript_name=tr.get("display_name"),
        refseq_xrefs=xr, protein=pid, protein_length_aa=len(cds) // 3 - 1, n_exons=len(T.exons), strand=tr["strand"],
        gene_span=f"{chrom}:{tr['start']}-{tr['end']}", cds_len=len(cds),
        R882=dict(codon=cds[(R882_CODON - 1) * 3:R882_CODON * 3], genomic=f"{chrom}:{r882_g[0]}-{r882_g[1]}",
                  tx_exon=r882_exon, domains=domains_hit(pfeat, R882_CODON, R882_CODON)),
        domains_all=domains_hit(pfeat, 1, len(cds) // 3),
        variant_as_configured=v,
    )
    out = dict(meta=meta, exons=T.exon_table(), analyses=analyses)
    pathlib.Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    pathlib.Path(args.out_json).write_text(json.dumps(out, indent=2))
    pathlib.Path(args.out_md).write_text(render_md(out), encoding="utf-8")
    print(f"wrote {args.out_md} and {args.out_json}")


if __name__ == "__main__":
    main()
