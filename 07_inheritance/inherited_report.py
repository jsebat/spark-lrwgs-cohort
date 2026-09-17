"""Rare INHERITED clinically reportable variants, long-read cohort.

Main result 2 (inherited arm). Each inheritance model is matched to the gene's
DDG2P allelic requirement -- a dominant call must be in a monoallelic gene, a
recessive call in a biallelic gene. Reporting a homozygous variant in a
dominant-only gene, or vice versa, is not clinically defensible.

Models:
  DOMINANT_INHERITED  proband het, gene monoallelic_autosomal, transmitted from
                      a parent. Parents are unaffected, so this is reduced
                      penetrance and is reported as such.
  HOMOZYGOUS          proband hom-alt, gene biallelic_autosomal, both parents het
  COMPOUND_HET        two rare variants, gene biallelic_autosomal, one from each parent
  X_LINKED            male proband hemizygous, gene monoallelic_X_hemizygous,
                      mother carrier
Filters: DDG2P confidence definitive|strong; rare; coding/impactful; genotype QC
(GQ>=20, DP>=10); lab segdup/repeat mask; unaffected sibling must not share the
same qualifying genotype.
"""
import csv
import glob
import gzip
import os
import subprocess
import sys
import bisect
import collections

R = "/expanse/lustre/projects/ddp195/jsebat/longread-autism"
T = "/expanse/lustre/projects/ddp195/jsebat/longread-autism/tiering"
GS = "/expanse/projects/sebat1/s3/data/sebat/nf_rare_spark_wes/resources/gene_sets"
MASKF = os.path.join(T, "mask.segdup_repeat.bed.gz")
SHET = os.path.join(T, "denovo_sv", "shet_by_symbol.tsv")
SIF = subprocess.Popen("ls /expanse/projects/sebat1/jsebat/singularity_cache/miniwdl/*pb_wdl_base*.sif",
                       shell=True, stdout=subprocess.PIPE).communicate()[0].decode().split()[0]

AF_DOM = 1e-4
AF_REC = 5e-3
NHOM_REC = 2
GQ_MIN, DP_MIN = 20, 10
LOF = ("stop_gained", "frameshift", "splice_acceptor", "splice_donor", "start_lost")
OK_CONF = ("definitive", "strong")


def sx(cmd):
    p = subprocess.Popen(["singularity", "exec", "-B", "/expanse:/expanse", SIF, "bash", "-c", cmd],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    o, e = p.communicate()
    if p.returncode != 0:
        sys.stderr.write("[sx] rc=%s%s%s%s" % (p.returncode, chr(10), e.decode()[:400], chr(10)))
    return o.decode()


# ---- DDG2P
g2p = {}
with open(os.path.join(GS, "DDG2P.csv"), newline="", encoding="utf-8", errors="replace") as fh:
    for r in csv.DictReader(fh):
        g = (r.get("gene symbol") or "").strip()
        conf = (r.get("confidence") or "").strip().lower()
        ar = (r.get("allelic requirement") or "").strip().lower()
        if not g or conf not in OK_CONF:
            continue
        e = g2p.setdefault(g, {"ar": set(), "dis": set(), "conf": set(), "mech": set()})
        e["ar"].add(ar)
        if r.get("disease name"):
            e["dis"].add(r["disease name"].strip())
        e["conf"].add(conf)
        if r.get("molecular mechanism"):
            e["mech"].add(r["molecular mechanism"].strip())
print("DDG2P genes (definitive|strong): %d" % len(g2p))

shet = {}
if os.path.exists(SHET):
    for ln in open(SHET):
        f = ln.rstrip(chr(10)).split(chr(9))
        if len(f) >= 2:
            try:
                shet[f[0]] = float(f[1])
            except ValueError:
                pass

# ---- mask
mask = {}
with gzip.open(MASKF, "rt") as fh:
    for ln in fh:
        f = ln.rstrip(chr(10)).split(chr(9))
        if len(f) >= 3:
            mask.setdefault(f[0], []).append((int(f[1]), int(f[2])))
starts = {}
for c in mask:
    mask[c].sort()
    starts[c] = [a for a, b in mask[c]]


def is_masked(c, p):
    iv = mask.get(c)
    if not iv:
        return False
    i = bisect.bisect_right(starts[c], p - 1) - 1
    if 0 <= i < len(iv):
        a, b = iv[i]
        return a < p <= b
    return False


def impact(csq):
    if any(k in csq for k in LOF):
        return "LoF"
    if "missense" in csq:
        return "missense"
    if "inframe" in csq:
        return "inframe"
    return None


rows = []
fams = 0
for ped in sorted(glob.glob(os.path.join(R, "run_*", "analysis", "*.ped"))):
    if os.environ.get("EXCLUDE_FAMILY", "__none__") in ped:
        continue
    fam = os.path.basename(ped)[:-4]
    members, aff, sex, pat, mat = [], {}, {}, {}, {}
    for ln in open(ped):
        f = ln.split()
        if len(f) < 6:
            continue
        members.append(f[1])
        pat[f[1]], mat[f[1]] = f[2], f[3]
        sex[f[1]] = f[4]
        aff[f[1]] = (f[5] == "2")
    probands = [s for s in members if aff.get(s) and pat[s] != "0"]
    if not probands:
        continue
    od = sorted(glob.glob(os.path.join(R, "run_" + fam, "*_humanwgs_family", "out")))
    if not od:
        continue
    vcfs = glob.glob(os.path.join(od[-1], "tertiary_small_variant_filtered_vcf", "*.vcf.gz"))
    if not vcfs:
        continue
    V = vcfs[0]
    fams += 1
    fmt = ("%CHROM" + chr(9) + "%POS" + chr(9) + "%REF" + chr(9) + "%ALT" + chr(9) +
           "%INFO/gnomad_af" + chr(9) + "%INFO/gnomad_nhomalt" + chr(9) + "%INFO/BCSQ" +
           "[" + chr(9) + "%SAMPLE|%GT|%GQ|%DP]" + chr(10))
    q = ("bcftools query -i 'INFO/gnomad_af<%g' -f '%s' '%s'" % (AF_REC, fmt, V))
    for ln in sx(q).splitlines():
        f = ln.split(chr(9))
        if len(f) < 8:
            continue
        chrom, pos, ref, alt = f[0], int(f[1]), f[2], f[3]
        try:
            af = float(f[4])
        except ValueError:
            af = -1.0
        try:
            nhom = int(f[5])
        except ValueError:
            nhom = -1
        bcsq = f[6]
        if is_masked(chrom, pos):
            continue
        gene = ""
        csq = ""
        for blk in bcsq.split(","):
            p = blk.split("|")
            if len(p) >= 2 and p[1]:
                imp = impact(p[0])
                if imp:
                    gene, csq = p[1], p[0]
                    break
        if not gene or gene not in g2p:
            continue
        imp = impact(csq)
        if imp not in ("LoF", "missense", "inframe"):
            continue
        gt = {}
        for cell in f[7:]:
            pp = cell.split("|")
            if len(pp) >= 4:
                try:
                    gq, dp = int(pp[2]), int(pp[3])
                except ValueError:
                    gq, dp = 0, 0
                gt[pp[0]] = (pp[1], gq, dp)

        def carrier(s):
            v = gt.get(s)
            if not v:
                return None
            g, gq, dp = v
            if g in ("./.", ".|.", "."):
                return None
            if gq < GQ_MIN or dp < DP_MIN:
                return None
            n = g.replace("|", "/").split("/")
            return n

        ars = g2p[gene]["ar"]
        for pb in probands:
            n = carrier(pb)
            if not n or "1" not in n:
                continue
            mo, fa = mat[pb], pat[pb]
            nm, nf = carrier(mo), carrier(fa)
            homalt = all(x == "1" for x in n) and len(n) >= 2
            het = ("1" in n and "0" in n)
            model = None
            inh = ""
            if chrom == "chrX" and sex.get(pb) == "1" and any("x_hemizygous" in a for a in ars):
                if nm and "1" in nm and (not nf or "1" not in nf):
                    model, inh = ("X_LINKED", "maternal") if af < AF_DOM else (None, "")
            elif homalt and "biallelic_autosomal" in ars and af < AF_REC and nhom <= NHOM_REC:
                if nm and "1" in nm and nf and "1" in nf:
                    model, inh = "HOMOZYGOUS", "both parents"
            elif het and "monoallelic_autosomal" in ars and af < AF_DOM:
                if nm and "1" in nm and (not nf or "1" not in nf):
                    model, inh = "DOMINANT_INHERITED", "maternal"
                elif nf and "1" in nf and (not nm or "1" not in nm):
                    model, inh = "DOMINANT_INHERITED", "paternal"
            if not model:
                continue
            # unaffected sibling must not share it
            shared = ""
            for s in members:
                if s in (pb, mo, fa) or aff.get(s):
                    continue
                ns = carrier(s)
                if ns and "1" in ns:
                    if model == "HOMOZYGOUS" and not all(x == "1" for x in ns):
                        continue
                    shared = "shared_with_unaffected_sib(%s)" % s
            rows.append([fam, pb, chrom, str(pos), ref[:20], alt[:20], gene, csq, imp,
                         model, inh, ("%.3g" % af) if af >= 0 else "absent",
                         ("%.4g" % shet[gene]) if gene in shet else "NA",
                         ";".join(sorted(g2p[gene]["conf"])),
                         ";".join(sorted(ars)),
                         ";".join(sorted(g2p[gene]["dis"]))[:60], shared])

print("families scanned: %d" % fams)
out = os.path.join(T, "inherited_reportable.tsv")
cols = ["family", "proband", "chrom", "pos", "ref", "alt", "gene", "consequence", "impact",
        "model", "inherited_from", "gnomad_af", "s_het", "ddg2p_confidence",
        "allelic_requirement", "disease", "flag"]
with open(out, "w") as fh:
    fh.write(chr(9).join(cols) + chr(10))
    for r in sorted(rows, key=lambda x: (x[8] != "LoF", x[9], x[0])):
        fh.write(chr(9).join(r) + chr(10))
print("")
print("wrote %s (%d rows)" % (out, len(rows)))
print("=" * 70)
c = collections.Counter((r[9], r[8]) for r in rows)
for k, v in sorted(c.items()):
    print("  %-20s %-9s %d" % (k[0], k[1], v))
print("  clean (not shared with an unaffected sib): %d" % sum(1 for r in rows if not r[16]))
print("=" * 70)
for r in sorted(rows, key=lambda x: (x[8] != "LoF", x[9])):
    if r[16]:
        continue
    print("  %-11s %-11s %-9s %-7s %-9s %-11s %-9s af=%-9s %s" %
          (r[0], r[1], r[6], r[8], r[9][:9], r[10], r[13][:9], r[11], r[15][:34]))
