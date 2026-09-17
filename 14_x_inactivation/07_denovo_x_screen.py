#!/usr/bin/env python3
"""Does any female's X-inactivation skew have a genetic cause on the X itself?

Two mechanisms are worth screening for. Cell selection: a variant that impairs proliferation or survival makes
cells that silenced the healthy X lose out, and the tissue drifts toward silencing the mutant X. Primary
skewing: a variant at the X-inactivation centre biases the choice itself, which is why XIST is searched by
coordinate here even though it carries no constraint score and appears on no disease gene list.

Three screens, because each answers a different question and each has a different blind spot.

  denovo   X-linked de novo loss-of-function and missense in candidate genes, trio daughters only.
           Read the raw de novo count in the output before reading the filtered one. Coding sequence is about
           one per cent of the X, so the expected number of de novo coding variants is well under one per
           cohort of this size: an empty result here is a statement about power, not about biology.

  burden   Rare X-linked damaging variants in candidate genes counted in every female and tested against
           skew. This is the screen with real power, and the only one available for mothers, whose parents are
           not sequenced and who therefore have no de novo status at all. Counts are compared against rare
           synonymous variants on the same chromosome in the same individual, which absorbs the dominant
           confounder: someone whose ancestry is less well represented in the reference panel carries more
           rare variants of every class and would otherwise look enriched for candidates.

  xic      Any rare variant in the X-inactivation centre, no consequence filter.

Annotation is read from a sites-only VEP VCF and joined to the call set on (chrom, pos, ref, alt). Verify that
join on a handful of sites before trusting an empty result: if the two files disagree on allele representation
the screen returns zero for a reason that has nothing to do with the data.
"""
import argparse
import csv
import gzip
import os
import subprocess
import sys
from collections import Counter, defaultdict

LOF_CSQ = {"frameshift_variant", "stop_gained", "splice_acceptor_variant", "splice_donor_variant",
           "start_lost", "stop_lost", "transcript_ablation"}
MIS_CSQ = {"missense_variant", "inframe_deletion", "inframe_insertion"}
XIC_GRCH38 = (73_790_000, 73_880_000)      # XIST / TSIX
XIC_GRCH37 = (72_990_000, 73_080_000)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ped", required=True)
    p.add_argument("--vep", required=True, help="sites-only VEP-annotated chrX VCF (bgzipped)")
    p.add_argument("--calls", required=True, help="cohort call set with genotypes (BCF/VCF, indexed)")
    p.add_argument("--skew", required=True, help="output of 04_xci_skew.py")
    p.add_argument("--gene-sets", required=True,
                   help="directory of newline-delimited gene symbol lists (*.set) defining candidate genes")
    p.add_argument("--constraint", default=None,
                   help="optional TSV with columns 'gene' and a constraint metric column")
    p.add_argument("--constraint-col", default="lof.oe_ci.upper",
                   help="column name for the constraint metric (default LOEUF)")
    p.add_argument("--constraint-max", type=float, default=0.35,
                   help="gene is constrained below this value (default 0.35)")
    p.add_argument("--out", required=True, help="output directory")
    p.add_argument("--chrom", default="chrX")
    p.add_argument("--build", choices=["GRCh38", "GRCh37"], default="GRCh38")
    p.add_argument("--max-af", type=float, default=0.001, help="gnomAD AF ceiling (default 0.001)")
    p.add_argument("--max-carriers", type=int, default=3,
                   help="variant must be seen in at most this many cohort females (default 3)")
    p.add_argument("--min-gq", type=int, default=20)
    p.add_argument("--min-dp", type=int, default=10)
    p.add_argument("--marked-skew", type=float, default=0.80,
                   help="threshold defining markedly skewed for the group comparison (default 0.80)")
    p.add_argument("--bcftools", default="bcftools")
    return p.parse_args()


def read_ped(path):
    ped = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) >= 6:
                ped[f[1]] = dict(fam=f[0], pat=f[2], mat=f[3], sex=f[4], aff=(f[5] == "2"))
    return ped


def load_gene_sets(d):
    sets = {}
    for path in sorted(os.listdir(d)):
        if not path.endswith(".set"):
            continue
        name = os.path.splitext(path)[0]
        with open(os.path.join(d, path)) as fh:
            sets[name] = {ln.strip() for ln in fh if ln.strip()}
    return sets


def load_constraint(path, col, ceiling):
    if not path:
        return set()
    out = set()
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            gi, ci = header.index("gene"), header.index(col)
        except ValueError:
            sys.stderr.write("constraint file lacks 'gene' or '%s'; ignoring it\n" % col)
            return out
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > max(gi, ci):
                try:
                    if float(f[ci]) < ceiling:
                        out.add(f[gi])
                except ValueError:
                    continue
    return out


def scan_vep(path, candidate, max_af, xic):
    """classify each site as candidate-damaging, synonymous control, or inside the XIC"""
    fields = None
    klass, annot = {}, {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ"):
                fields = line.split("Format: ")[1].split('">')[0].split("|")
                continue
            if line.startswith("#") or "CSQ=" not in line:
                continue
            f = line.rstrip("\n").split("\t")
            key = (f[0], f[1], f[3], f[4])
            pos = int(f[1])
            info = f[7].split("CSQ=")[1].split(";")[0]
            best, syn = None, False
            for tr in info.split(","):
                d = dict(zip(fields, tr.split("|")))
                sym = d.get("SYMBOL", "")
                cons = set(d.get("Consequence", "").split("&"))
                af_raw = next((d[k] for k in d if k.startswith("gnomAD") and k.endswith("AF")), "")
                try:
                    af = float(af_raw) if af_raw not in ("", ".") else 0.0
                except ValueError:
                    af = 0.0
                if "synonymous_variant" in cons and af <= max_af:
                    syn = True
                if not sym or af > max_af:
                    continue
                kind = "LOF" if cons & LOF_CSQ else ("missense" if cons & MIS_CSQ else None)
                if not kind or sym not in candidate:
                    continue
                rec = dict(gene=sym, kind=kind, consequence=d.get("Consequence", ""), af=af,
                           canonical=(d.get("CANONICAL") == "YES"))
                if (best is None or (rec["canonical"] and not best["canonical"])
                        or (rec["kind"] == "LOF" and best["kind"] != "LOF")):
                    best = rec
            if best:
                klass[key] = "candidate"
                annot[key] = best
            elif syn:
                klass[key] = "synonymous"
            if xic[0] <= pos <= xic[1]:
                klass.setdefault(key, "xic")
                annot.setdefault(key, dict(gene="(XIC)", kind="XIC",
                                           consequence="", af=0.0, canonical=False))
    return klass, annot


def main():
    a = parse_args()
    os.makedirs(a.out, exist_ok=True)
    xic = XIC_GRCH38 if a.build == "GRCh38" else XIC_GRCH37
    ped = read_ped(a.ped)

    with open(a.skew) as fh:
        skew_rows = list(csv.DictReader(fh, delimiter="\t"))
    skew = {r["sample"]: float(r["skew"]) for r in skew_rows}
    role = {r["sample"]: r["role"] for r in skew_rows}
    females = sorted(skew)
    if not females:
        sys.exit("no females in %s" % a.skew)

    sets = load_gene_sets(a.gene_sets)
    candidate = set().union(*sets.values()) if sets else set()
    candidate |= load_constraint(a.constraint, a.constraint_col, a.constraint_max)
    sys.stderr.write("candidate genes: %d from %d set(s)%s\n"
                     % (len(candidate), len(sets), " plus constraint" if a.constraint else ""))

    klass, annot = scan_vep(a.vep, candidate, a.max_af, xic)
    counts = Counter(klass.values())
    sys.stderr.write("annotated sites: %s\n" % dict(counts))
    if not counts.get("candidate"):
        sys.stderr.write("WARNING: no candidate sites found; check gene-set symbols match the VEP SYMBOL "
                         "field\n")

    regions = os.path.join(a.out, "_sites.txt")
    with open(regions, "w") as fh:
        for chrom, pos, _, _ in sorted(klass, key=lambda k: int(k[1])):
            fh.write("%s\t%s\n" % (chrom, pos))

    cmd = [a.bcftools, "query", "-R", regions, "-s", ",".join(females),
           "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT:%GQ:%DP:%AD]\n", a.calls]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.exit("bcftools query failed: %s" % proc.stderr.strip()[-400:])
    idx = {s: i for i, s in enumerate(females)}

    def cell(text):
        p = text.split(":")
        gt = p[0].replace("|", "/")

        def num(x):
            try:
                return int(x)
            except ValueError:
                return -1
        gq = num(p[1]) if len(p) > 1 else -1
        dp = num(p[2]) if len(p) > 2 else -1
        ad = [num(x) for x in p[3].split(",")] if len(p) > 3 else []
        return gt, gq, dp, ad

    per_sample = defaultdict(lambda: Counter())
    hits, denovo, xic_hits = [], [], []
    joined = 0
    kids = {s for s in females if ped.get(s, {}).get("pat", "0") in ped
            and ped.get(s, {}).get("mat", "0") in ped}

    # parents are needed for the de novo test but not scored themselves
    parent_cols = sorted({ped[s]["pat"] for s in kids} | {ped[s]["mat"] for s in kids})
    pcmd = [a.bcftools, "query", "-R", regions, "-s", ",".join(parent_cols),
            "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT:%GQ:%DP:%AD]\n", a.calls] if parent_cols else None
    parent_at = {}
    if pcmd:
        pr = subprocess.run(pcmd, capture_output=True, text=True)
        if pr.returncode == 0:
            pidx = {s: i for i, s in enumerate(parent_cols)}
            for line in pr.stdout.splitlines():
                f = line.rstrip("\n").split("\t")
                parent_at[(f[0], f[1], f[2], f[3])] = (f[4:], pidx)

    for line in proc.stdout.splitlines():
        f = line.rstrip("\n").split("\t")
        key = (f[0], f[1], f[2], f[3])
        k = klass.get(key)
        if not k or "," in f[3]:
            continue
        joined += 1
        cells = f[4:]
        carriers = []
        for s in females:
            gt, gq, dp, ad = cell(cells[idx[s]])
            if gt not in ("0/1", "1/0", "1/1") or gq < a.min_gq or dp < a.min_dp:
                continue
            carriers.append((s, gt, ad))
        if not carriers or len(carriers) > a.max_carriers:
            continue
        for s, gt, ad in carriers:
            per_sample[s][k] += 1
            ann = annot.get(key, {})
            row = dict(sample=s, skew=skew[s], role=role[s], chrom=f[0], pos=int(f[1]),
                       ref=f[2], alt=f[3], gene=ann.get("gene", ""), kind=ann.get("kind", ""),
                       consequence=ann.get("consequence", ""), gnomad_af=ann.get("af", ""),
                       n_female_carriers=len(carriers), gt=gt)
            if k == "candidate":
                hits.append(row)
            elif k == "xic":
                xic_hits.append(row)
            # de novo: both parents called homozygous reference with support
            if k == "candidate" and s in kids and key in parent_at:
                pcells, pidx = parent_at[key]
                try:
                    dg, dq, dd, da = cell(pcells[pidx[ped[s]["pat"]]])
                    mg, mq, md, ma = cell(pcells[pidx[ped[s]["mat"]]])
                except (IndexError, KeyError):
                    continue
                if (dg == "0/0" and dq >= a.min_gq and dd >= a.min_dp
                        and (len(da) < 2 or da[1] <= 1)
                        and mg == "0/0" and mq >= a.min_gq and md >= a.min_dp
                        and (len(ma) < 2 or ma[1] <= 1)
                        and len(ad) >= 2 and ad[1] >= 5
                        and ad[1] / max(sum(ad), 1) >= 0.25):
                    denovo.append(row)

    sys.stderr.write("sites joined to the call set: %d of %d annotated\n" % (joined, len(klass)))
    if joined == 0:
        sys.stderr.write("WARNING: nothing joined. The VEP file and the call set probably disagree on "
                         "allele representation; check REF/ALT on a few sites before trusting this.\n")

    def write(name, rows):
        path = os.path.join(a.out, name)
        if not rows:
            open(path, "w").write("")
            return
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    write("candidate_variants.tsv", sorted(hits, key=lambda r: -r["skew"]))
    write("denovo_candidates.tsv", sorted(denovo, key=lambda r: -r["skew"]))
    write("xic_variants.tsv", sorted(xic_hits, key=lambda r: -r["skew"]))

    counts_path = os.path.join(a.out, "burden_counts.tsv")
    with open(counts_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "role", "skew", "n_candidate", "n_synonymous"])
        for s in females:
            w.writerow([s, role[s], "%.4f" % skew[s],
                        per_sample[s]["candidate"], per_sample[s]["synonymous"]])

    print("=== de novo screen (trio daughters: %d) ===" % len(kids))
    print("  X-linked de novo LOF/missense in candidate genes: %d" % len(denovo))
    for r in denovo:
        print("    %s skew %.3f  %s %s %s" % (r["sample"], r["skew"], r["gene"], r["kind"],
                                              r["consequence"][:30]))
    print("  NOTE: coding sequence is ~1%% of the X, so the expected yield is well under one variant.")
    print("  An empty result here bounds power, it does not exclude a genetic cause.")

    print("\n=== burden vs skew (all %d females) ===" % len(females))
    try:
        import numpy as np
        from scipy import stats
        sk = np.array([skew[s] for s in females])
        cand = np.array([per_sample[s]["candidate"] for s in females], float)
        syn = np.array([per_sample[s]["synonymous"] for s in females], float)
        ratio = cand / np.maximum(syn, 1)
        for label, v in (("candidate count", cand), ("candidate/synonymous", ratio),
                         ("synonymous (control)", syn)):
            rho, p = stats.spearmanr(sk, v)
            print("  %-24s vs skew: Spearman rho = %+.3f, p = %.3f" % (label, rho, p))
        hi = sk >= a.marked_skew
        if hi.sum() and (~hi).sum():
            print("  markedly skewed (>=%.2f, n=%d) vs rest (n=%d):"
                  % (a.marked_skew, hi.sum(), (~hi).sum()))
            for label, v in (("candidate count", cand), ("candidate/synonymous", ratio)):
                u, p = stats.mannwhitneyu(v[hi], v[~hi])
                print("    %-22s %.2f vs %.2f, Mann-Whitney p = %.3f"
                      % (label, v[hi].mean(), v[~hi].mean(), p))
    except ImportError:
        print("  (install numpy and scipy for the association test; counts are in %s)" % counts_path)

    print("\n=== X-inactivation centre ===")
    print("  rare variants in any female: %d" % len(xic_hits))
    print("\nwrote %s/" % a.out)


if __name__ == "__main__":
    main()
