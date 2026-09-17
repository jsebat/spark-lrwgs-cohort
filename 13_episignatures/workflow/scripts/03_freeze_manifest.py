#!/usr/bin/env python
"""Step 03 - build the sample manifest from SPARK metadata, pre-register variant flags, freeze it.

Runs on Expanse only (reads controlled-access metadata and cohort variant calls). Never prints sample IDs;
the step log reports counts only.

Manifest columns (PLAN.md + extras):
  sample_id family_id role age sex mean_depth flag_gene flag_note
  asd pilot_phase extraction_method n_movies bam_source bam_paths is_proband proband_family qc_note

Flags (config.flag_genes; canonical transcript models in resources/flag_genes/):
  * SNV/indel from the glnexus cohort BCF: variants inside CDS exons +/- SPLICE_PAD bp with a non-silent consequence
    (missense / nonsense / frameshift / in-frame indel / splice-region / start-loss / stop-loss), PASS or unfiltered,
    cohort allele count <= max_cohort_ac (rarity is cohort-internal: no gnomAD available offline - stated in the log).
  * SVs from the sawfish joint VCF: PASS records overlapping any CDS exon of a flag gene.
  * The configured proband deletion identifies the proband; its parents get proband_family = True.

Outputs (results/03_manifest/, gitignored except STEP_LOG.md):
  samples.tsv (also copied to config/samples.tsv), samples.frozen.tsv, samples.frozen.sha256, STEP_LOG.md
"""
import argparse
import datetime as dt
import glob
import hashlib
import json
import os
import pathlib
import shutil
import subprocess
import sys

import pandas as pd
import yaml

SPLICE_PAD = 8
BASES = "TCAG"
AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AAS[i * 16 + j * 4 + k] for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}
COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def rc(s):
    return s.translate(COMP)[::-1]


class GeneModel:
    def __init__(self, js):
        self.g = js
        self.chrom, self.strand = js["chrom"], js["strand"]
        self.cds = js["cds_seq"]
        self.exons = js["exons"]  # transcript order
        self.lo, self.hi = js["cds_genomic_lo"], js["cds_genomic_hi"]
        # genomic position -> CDS index (0-based), CDS portion only
        self.g2c = {}
        c = 0
        for e in self.exons:
            rng = range(e["end"], e["start"] - 1, -1) if self.strand == -1 else range(e["start"], e["end"] + 1)
            for p in rng:
                if self.lo <= p <= self.hi:
                    self.g2c[p] = c
                    c += 1
        assert c == len(self.cds), f"{js['gene']}: CDS map {c} != CDS len {len(self.cds)}"

    def cds_exon_intervals(self):
        for e in self.exons:
            lo, hi = max(e["start"], self.lo), min(e["end"], self.hi)
            if lo <= hi:
                yield lo, hi

    def consequence(self, pos, ref, alt):
        """Classify a VCF record (1-based pos) against this model. Returns (class, detail) or None if outside CDS+-pad."""
        ref_span = range(pos, pos + len(ref))
        in_cds = [p for p in ref_span if p in self.g2c]
        near = any(lo - SPLICE_PAD <= p <= hi + SPLICE_PAD for p in ref_span for lo, hi in self.cds_exon_intervals())
        if not near:
            return None
        if not in_cds:
            # intronic within pad of a CDS exon boundary: splice region if within 1-2 bp core, else splice-region
            core = any(p in (lo - 1, lo - 2, hi + 1, hi + 2) for p in ref_span for lo, hi in self.cds_exon_intervals())
            return ("splice_donor_acceptor" if core else "splice_region", f"{pos}{ref}>{alt}")
        if len(ref) == len(alt) == 1:
            ci = self.g2c[pos]
            codon_i = ci // 3
            codon = self.cds[codon_i * 3: codon_i * 3 + 3]
            base = alt if self.strand == 1 else rc(alt)
            off = ci % 3
            newc = codon[:off] + base + codon[off + 1:]
            a0, a1 = CODON.get(codon, "?"), CODON.get(newc, "?")
            if a0 == a1:
                cls = "synonymous"
            elif a1 == "*":
                cls = "nonsense"
            elif a0 == "*":
                cls = "stop_loss"
            elif codon_i == 0:
                cls = "start_loss"
            else:
                cls = "missense"
            return (cls, f"p.{a0}{codon_i + 1}{a1}")
        dlen = len(alt) - len(ref)
        if len(in_cds) < len(ref_span):  # indel spanning an exon boundary
            return ("splice_indel", f"{pos}{ref}>{alt}")
        if dlen % 3 == 0:
            return ("inframe_indel", f"c.{min(self.g2c[p] for p in in_cds) + 1}{'del' if dlen < 0 else 'ins'}{abs(dlen)}")
        return ("frameshift", f"c.{min(self.g2c[p] for p in in_cds) + 1}fs")


NONSILENT = {"missense", "nonsense", "stop_loss", "start_loss", "frameshift", "inframe_indel", "splice_donor_acceptor", "splice_indel"}


def run(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"command failed: {cmd}\n{r.stderr[-2000:]}")
    return r.stdout


def load_metadata(cfg):
    st = pathlib.Path(cfg["paths"]["staging_dir"])
    ind = pd.read_csv(st / cfg["paths"]["metadata"]["individuals"], sep="\t", dtype=str)
    smp = pd.read_csv(st / cfg["paths"]["metadata"]["samples"], sep="\t", dtype=str)
    smp = smp[smp["sequenced"].astype(str).str.lower() == "true"].copy()
    smp["total_coverage"] = pd.to_numeric(smp["total_coverage"], errors="coerce")
    # one row per spid: keep the highest-coverage sequenced sample
    smp = smp.sort_values("total_coverage", ascending=False).drop_duplicates("spid")
    m = ind.merge(smp, on="spid", how="inner", suffixes=("", "_s"))
    spids = set(m["spid"])
    is_child = m["mother_spid"].isin(spids) | m["father_spid"].isin(spids)
    is_parent = m["spid"].isin(set(m["mother_spid"]) | set(m["father_spid"]))
    m["role"] = "other"
    m.loc[is_parent & (m["sex"].str.lower() == "female"), "role"] = "mother"
    m.loc[is_parent & (m["sex"].str.lower() == "male"), "role"] = "father"
    m.loc[is_child, "role"] = "child"  # a child who is also a parent (multi-gen) stays 'child' for the primary null
    return m


def resolve_bams(m, cfg):
    merged = pathlib.Path(cfg["paths"]["merged_bam_dir"])
    bam_dir = pathlib.Path(cfg["paths"]["bam_dir"])
    src, paths, nmov = [], [], []
    for _, r in m.iterrows():
        mb = merged / f"{r['spid']}.GRCh38.hifi_reads.bam"
        movies = sorted(glob.glob(str(bam_dir / "*" / r["spid"] / "*.bam")))
        nmov.append(len(movies))
        if mb.exists() and (mb.with_suffix(".bam.bai").exists() or pathlib.Path(str(mb) + ".bai").exists()):
            src.append("merged"); paths.append(str(mb))
        elif movies:
            src.append("movies"); paths.append(";".join(movies))
        else:
            src.append("none"); paths.append("")
    m["bam_source"], m["bam_paths"], m["n_movies"] = src, paths, nmov
    return m


def flag_small_variants(m, cfg, models, bed, max_ac):
    bcf = cfg["paths"]["variant_calls"]["small_variants_bcf"]
    pad_bed = bed.with_name("flag_genes.cds_exons.pad8.bed")
    b = pd.read_csv(bed, sep="\t", header=None)
    b[1] = (b[1] - SPLICE_PAD).clip(lower=0); b[2] = b[2] + SPLICE_PAD
    b.to_csv(pad_bed, sep="\t", header=False, index=False)
    q = run(f"bcftools view -R {pad_bed} -i 'FILTER=\"PASS\" || FILTER=\".\"' -Ou {bcf} | "
            "bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%AC\\t%AN[\\t%SAMPLE=%GT]\\n'")
    hits, n_sites, n_nonsilent = {}, 0, 0
    for line in q.splitlines():
        f = line.split("\t")
        chrom, pos, ref, alts = f[0], int(f[1]), f[2], f[3].split(",")
        n_sites += 1
        # glnexus BCFs carry no INFO/AC: count alt alleles from the genotypes themselves
        gts = [sg.split("=")[1].replace("|", "/").split("/") for sg in f[6:]]
        acs = [sum(a.count(str(ai + 1)) for a in gts) for ai in range(len(alts))]
        for ai, alt in enumerate(alts):
            if alt in ("*", "<NON_REF>"):
                continue
            for gm in models:
                if gm.chrom != chrom:
                    continue
                res = gm.consequence(pos, ref, alt)
                if not res or res[0] not in NONSILENT:
                    continue
                n_nonsilent += 1
                this_ac = acs[ai]
                if this_ac > max_ac:
                    continue
                for sg in f[6:]:
                    s, gt = sg.split("=")
                    if str(ai + 1) in gt.replace("|", "/").split("/"):
                        hits.setdefault(s, []).append(f"{gm.g['gene']}:{res[0]}:{res[1]}:AC={this_ac}")
    return hits, n_sites, n_nonsilent


def flag_svs(m, cfg, models, bed, pv):
    vcf = cfg["paths"]["variant_calls"]["sv_vcf"]
    q = run(f"bcftools view -R {bed} -f PASS -Ou {vcf} | "
            "bcftools query -f '%CHROM\\t%POS\\t%INFO/END\\t%INFO/SVTYPE\\t%INFO/SVLEN[\\t%SAMPLE=%GT]\\n'")
    hits, proband, n = {}, set(), 0
    for line in q.splitlines():
        f = line.split("\t")
        chrom, pos, end, svtype = f[0], int(f[1]), f[2], f[3]
        end = int(end) if end not in (".", "") else pos
        genes = sorted({gm.g["gene"] for gm in models if gm.chrom == chrom and any(lo <= end and hi >= pos for lo, hi in gm.cds_exon_intervals())})
        if not genes:
            continue
        n += 1
        is_pv = chrom == pv["chrom"] and pos == int(pv["start"]) and end == int(pv["end"])
        for sg in f[5:]:
            s, gt = sg.split("=")
            if "1" in gt.replace("|", "/").split("/"):
                hits.setdefault(s, []).append(f"{'/'.join(genes)}:SV_{svtype}:{chrom}:{pos}-{end}")
                if is_pv:
                    proband.add(s)
    return hits, proband, n


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", required=True)
    ap.add_argument("--max-cohort-ac", type=int, default=6, help="max cohort allele count for a small variant to count as rare")
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config))
    out = pathlib.Path(cfg["paths"]["output_root"], "03_manifest"); out.mkdir(parents=True, exist_ok=True)
    gdir = pathlib.Path(cfg["paths"]["resources_dir"], "flag_genes")
    models = [GeneModel(json.load(open(p))) for p in sorted(gdir.glob("*.ENST*.json"))]
    bed = gdir / "flag_genes.cds_exons.hg38.bed"

    m = load_metadata(cfg)
    m = resolve_bams(m, cfg)
    small, n_sites, n_nonsilent = flag_small_variants(m, cfg, models, bed, a.max_cohort_ac)
    svs, proband_ids, n_sv = flag_svs(m, cfg, models, bed, cfg["proband_variant"])

    m["is_proband"] = m["spid"].isin(proband_ids)
    prob_fams = set(m.loc[m["is_proband"], "sfid"])
    m["proband_family"] = m["sfid"].isin(prob_fams)
    m["flag_gene"], m["flag_note"] = "", ""
    for i, r in m.iterrows():
        notes = small.get(r["spid"], []) + svs.get(r["spid"], [])
        if notes:
            m.at[i, "flag_gene"] = ";".join(sorted({n.split(":")[0] for n in notes}))
            m.at[i, "flag_note"] = ";".join(notes)
        if r["is_proband"]:
            m.at[i, "flag_note"] = (m.at[i, "flag_note"] + ";" if m.at[i, "flag_note"] else "") + "DNMT3A proband (configured deletion)"
        elif r["sfid"] in prob_fams:
            m.at[i, "flag_note"] = (m.at[i, "flag_note"] + ";" if m.at[i, "flag_note"] else "") + "parent/sib of DNMT3A proband"
    m["qc_note"] = ""
    m.loc[m["bam_source"] == "none", "qc_note"] = "no BAM found"

    cols = ["spid", "sfid", "role", "age_y", "sex", "total_coverage", "flag_gene", "flag_note", "asd", "pilot_phase",
            "extraction_method", "n_movies", "bam_source", "bam_paths", "is_proband", "proband_family", "qc_note"]
    man = m[cols].rename(columns={"spid": "sample_id", "sfid": "family_id", "age_y": "age", "total_coverage": "mean_depth"})
    man = man.sort_values(["family_id", "role", "sample_id"])
    frozen = out / "samples.frozen.tsv"
    candidate = out / "samples.candidate.tsv"
    man.to_csv(candidate, sep="\t", index=False)
    if frozen.exists():
        old = hashlib.sha256(frozen.read_bytes()).hexdigest()
        new = hashlib.sha256(candidate.read_bytes()).hexdigest()
        if old != new:
            # the freeze is the freeze: samples.tsv and config/samples.tsv (what the Snakefile and step 05 read) must not
            # silently move away from samples.frozen.tsv, which is what happened before (X26)
            print(f"WARNING: frozen manifest exists and differs (sha256 {old[:12]} vs {new[:12]}); samples.tsv and "
                  f"{cfg['paths']['samples_tsv']} are LEFT AS FROZEN. The new manifest is in {candidate}; delete the frozen file deliberately to re-freeze.")
            return
    shutil.copy(candidate, out / "samples.tsv")
    shutil.copy(candidate, cfg["paths"]["samples_tsv"])
    if not frozen.exists():
        shutil.copy(candidate, frozen)
    sha = hashlib.sha256(frozen.read_bytes()).hexdigest()
    (out / "samples.frozen.sha256").write_text(f"{sha}  samples.frozen.tsv\n")

    ch = man[man["role"] == "child"]; par = man[man["role"].isin(["mother", "father"])]
    flagged = man[man["flag_gene"] != ""]
    log = [f"# Step 03 - manifest freeze ({dt.datetime.now().isoformat(timespec='seconds')})", "",
           f"Frozen manifest sha256: `{sha}`", "",
           "| metric | n |", "|---|---|",
           f"| individuals in manifest (sequenced) | {len(man)} |",
           f"| families | {man['family_id'].nunique()} |",
           f"| children (>=1 sequenced parent) | {len(ch)} |",
           f"| children with ASD | {int((ch['asd'].astype(str).str.lower() == 'true').sum())} |",
           f"| parents | {len(par)} |",
           f"| role = other | {int((man['role'] == 'other').sum())} |",
           f"| BAM source merged / movies / none | {int((man['bam_source']=='merged').sum())} / {int((man['bam_source']=='movies').sum())} / {int((man['bam_source']=='none').sum())} |",
           f"| proband rows (configured deletion carriers) | {int(man['is_proband'].sum())} |",
           f"| proband-family rows (incl. proband) | {int(man['proband_family'].sum())} |",
           f"| flagged samples (any flag gene) | {len(flagged)} |",
           f"| flagged children / parents | {int((flagged['role']=='child').sum())} / {int(flagged['role'].isin(['mother','father']).sum())} |",
           f"| small-variant sites screened in flag-gene CDS+-{SPLICE_PAD}bp | {n_sites} |",
           f"| non-silent alt alleles (any AC) | {n_nonsilent} |",
           f"| PASS SVs overlapping flag-gene CDS | {n_sv} |", "",
           "Flag counts by gene (samples): " + ", ".join(f"{g}={int(man['flag_gene'].str.contains(g).sum())}" for g in cfg["flag_genes"]), "",
           f"Rarity criterion for small variants: cohort allele count <= {a.max_cohort_ac} (no population frequency available offline). "
           "Consequence classes counted: " + ", ".join(sorted(NONSILENT)) + ". Synonymous and deep-intronic variants ignored.",
           "mean_depth = SPARK metadata total_coverage; step 04 re-measures depth from the BAMs.",
           "Sample identifiers are confined to results/03_manifest/ and config/samples.tsv (both gitignored)."]
    (out / "STEP_LOG.md").write_text("\n".join(log) + "\n")
    print("\n".join(log))


if __name__ == "__main__":
    main()
