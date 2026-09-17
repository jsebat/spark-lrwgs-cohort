"""Tier the phase-aware module's de novo call set with the cohort's own tiering resources, and say what is NEW.

Nothing is re-annotated: gene, consequence, LOFTEE and gnomAD AF all come from tiering/vep/chr*.vep.vcf.gz, the same
annotation behind the cohort report's denovo_tiered.tsv, pulled by tabix at the module's call positions. LoF tiers use
GeneBayes s_het on the cohort's thresholds (t1 >= 0.18, t2 0.03-0.18, t3 otherwise). SVs and TRs are tiered by CDS
overlap against gencode_cds.bed and the same SFARI / DDG2P / constrained gene sets the cohort SV tiering used.

A call is NEW when the original (heuristic) pipeline's tables do not contain it, matched on (proband, chrom, pos) with a
window for SVs and TRs, whose breakpoints differ between callers and runs.

Usage: tier_module.py FINAL_DIR TIERING_DIR OUT_DIR [--tabix PATH] [--sv-window N]
"""
import sys, os, glob, csv, gzip, subprocess, collections, argparse

ap = argparse.ArgumentParser()
ap.add_argument("final_dir"); ap.add_argument("tiering_dir"); ap.add_argument("out_dir")
ap.add_argument("--tabix", default="tabix")
ap.add_argument("--sv-window", type=int, default=1000, help="bp slop when matching an SV/TR to the original call set")
A = ap.parse_args()
T, OUT = A.tiering_dir, A.out_dir
DSV = os.path.join(T, "denovo_sv")
os.makedirs(OUT, exist_ok=True)
SHET_T1, SHET_T2 = 0.18, 0.03
log = lambda m: sys.stderr.write(m + "\n")

# ------------------------------------------------------------------ resources
def load_set(p):
    return set(x.strip() for x in open(p)) if os.path.exists(p) else set()

sfari_hc = load_set(os.path.join(DSV, "sfari_hc.set"))
sfari_all = load_set(os.path.join(DSV, "sfari_all.set"))
ddg2p_conf = load_set(os.path.join(DSV, "ddg2p_conf.set"))
ddg2p_all = load_set(os.path.join(DSV, "ddg2p_all.set"))
NDD = sfari_hc | sfari_all | ddg2p_conf | ddg2p_all
shet = {}
p = os.path.join(DSV, "shet_by_symbol.tsv")
for ln in (open(p) if os.path.exists(p) else []):
    f = ln.rstrip("\n").split("\t")
    if len(f) >= 2:
        try:
            shet[f[0]] = float(f[1])
        except ValueError:
            pass
CONSTRAINED = set(g for g, v in shet.items() if v >= SHET_T1)
log("gene sets: SFARI hc %d, SFARI all %d, DDG2P conf %d, DDG2P all %d; s_het %d genes, %d constrained (>=%.2f)"
    % (len(sfari_hc), len(sfari_all), len(ddg2p_conf), len(ddg2p_all), len(shet), len(CONSTRAINED), SHET_T1))

def gene_flags(g):
    out = []
    if g in sfari_hc:
        out.append("SFARI_hc")
    elif g in sfari_all:
        out.append("SFARI")
    if g in ddg2p_conf:
        out.append("DDG2P_conf")
    elif g in ddg2p_all:
        out.append("DDG2P")
    if g in CONSTRAINED:
        out.append("CONSTRAINED")
    return ";".join(out)

# ------------------------------------------------------------------ the module's call set
calls = []
for cls in ("snv_indel", "sv", "tr"):
    for fp in sorted(glob.glob(os.path.join(A.final_dir, "*.%s.dnm.tsv" % cls))):
        with open(fp, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("dnm_call") not in ("YES", "CANDIDATE"):
                    continue
                r["_cls"] = cls
                calls.append(r)
log("module calls (tier 1 + tier 2): %d  %s" % (len(calls), dict(collections.Counter(c["_cls"] for c in calls))))

# ------------------------------------------------------------------ VEP by tabix at the call positions
small = [c for c in calls if c["_cls"] == "snv_indel"]
vep = {}
csq_fmt = []
by_chrom = collections.defaultdict(list)
for c in small:
    by_chrom[c["chrom"]].append(int(c["start"]))
for chrom, poss in sorted(by_chrom.items()):
    vp = os.path.join(T, "vep", "%s.vep.vcf.gz" % chrom)
    if not os.path.exists(vp):
        log("  no VEP for %s" % chrom)
        continue
    if not csq_fmt:
        with gzip.open(vp, "rt") as fh:
            for line in fh:
                if line.startswith("##INFO=<ID=CSQ"):
                    csq_fmt = line.split("Format: ")[-1].rstrip().rstrip('">').split("|")
                    break
                if not line.startswith("#"):
                    break
    reg = os.path.join(OUT, "_reg.%s.bed" % chrom)
    with open(reg, "w", newline="\n") as fh:
        for p0 in sorted(set(poss)):
            fh.write("%s\t%d\t%d\n" % (chrom, max(0, p0 - 1), p0 + 1))
    try:
        txt = subprocess.run([A.tabix, "-R", reg, vp], capture_output=True, text=True, check=True).stdout
    except (subprocess.CalledProcessError, OSError) as e:
        log("  tabix failed on %s: %s" % (chrom, e))
        continue
    for line in txt.splitlines():
        f = line.split("\t", 8)
        if len(f) < 8:
            continue
        csq = ""
        for kv in f[7].split(";"):
            if kv.startswith("CSQ="):
                csq = kv[4:]
        if not csq:
            continue
        best, best_rank = None, -1
        for tr in csq.split(","):
            d = dict(zip(csq_fmt, tr.split("|")))
            # the cohort tiering reports the MANE Select transcript where there is one, else the canonical;
            # a high-confidence LoF call is never lost to transcript choice
            rank = 4 if d.get("LoF") == "HC" else 3 if d.get("MANE_SELECT") else 2 if d.get("CANONICAL") == "YES" else 1
            if rank > best_rank:
                best, best_rank = d, rank
        if best is not None:
            vep[(f[0], int(f[1]), f[3], f[4])] = best
            vep.setdefault((f[0], int(f[1])), best)
    os.remove(reg)
log("VEP records pulled: %d for %d small-variant calls" % (len(vep), len(small)))

# ------------------------------------------------------------------ CDS intervals for SV / TR
cds = collections.defaultdict(list)
p = os.path.join(DSV, "gencode_cds.bed")
for ln in (open(p) if os.path.exists(p) else []):
    f = ln.rstrip("\n").split("\t")
    if len(f) >= 4:
        cds[f[0]].append((int(f[1]), int(f[2]), f[3]))
for ch in cds:
    cds[ch].sort()
log("CDS intervals: %d over %d contigs" % (sum(len(v) for v in cds.values()), len(cds)))

def genes_at(chrom, s0, e0):
    out = set()
    for a, b, g in cds.get(chrom, ()):
        if a >= e0:
            break
        if b > s0:
            out.add(g)
    return out

# ------------------------------------------------------------------ the original pipeline's call set, for NEW
orig_small, orig_iv = set(), collections.defaultdict(list)
p = os.path.join(T, "denovo_tiered.tsv")
for r in (csv.DictReader(open(p), delimiter="\t") if os.path.exists(p) else []):
    orig_small.add((r["proband"], r["chrom"], int(r["pos"])))
n_osv = 0
p = os.path.join(DSV, "denovo_sv_PRIORITIZED_v2.tsv")
if not os.path.exists(p):
    p = os.path.join(DSV, "denovo_sv_PRIORITIZED.tsv")
for r in (csv.DictReader(open(p), delimiter="\t") if os.path.exists(p) else []):
    try:
        orig_iv[(r["proband"], r["chrom"])].append((int(r["start"]), int(r["end"])))
        n_osv += 1
    except (KeyError, ValueError):
        pass
p = os.path.join(DSV, "denovo_sv_all.bed")
for ln in (open(p) if os.path.exists(p) else []):
    f = ln.rstrip("\n").split("\t")
    if len(f) >= 5:
        try:
            orig_iv[(f[4], f[0])].append((int(f[1]), int(f[2])))
            n_osv += 1
        except ValueError:
            pass
log("original call set: %d small variants, %d SV intervals" % (len(orig_small), n_osv))

def is_new(c):
    sid, ch = c.get("sample_id", ""), c["chrom"]
    s0 = int(float(c["start"]))
    e0 = int(float(c.get("end") or 0) or s0 + 1)
    if c["_cls"] == "snv_indel":
        return not any((sid, ch, s0 + d) in orig_small for d in (-1, 0, 1))
    for a, b in orig_iv.get((sid, ch), ()):
        if a - A.sv_window < e0 and s0 < b + A.sv_window:
            return False
    return True

# ------------------------------------------------------------------ tier every call
rows = []
for c in calls:
    cls = c["_cls"]
    s0 = int(float(c["start"]))
    e0 = int(float(c.get("end") or 0) or s0 + 1)
    gene = consequence = lof = lof_tier = gaf = ""
    if cls == "snv_indel":
        d = vep.get((c["chrom"], s0, c.get("ref", ""), c.get("alt", ""))) or vep.get((c["chrom"], s0)) or {}
        gene = d.get("SYMBOL", "")
        consequence = d.get("Consequence", "")
        lof = d.get("LoF", "")
        gaf = d.get("gnomAD4.1_joint_AF", "")
        s = shet.get(gene)
        if lof == "HC" and not d.get("LoF_filter"):
            lof_tier = "lof_t1" if (s is not None and s >= SHET_T1) else "lof_t2" if (s is not None and s >= SHET_T2) else "lof_t3"
        genes = [gene] if gene else []
    else:
        allg = sorted(genes_at(c["chrom"], s0, max(e0, s0 + 1)))
        hit = [g for g in allg if g in NDD or g in CONSTRAINED]
        gene = ";".join(hit[:6]) if hit else ";".join(allg[:3])
        consequence = ("%s %s bp" % (c.get("svtype", ""), c.get("svlen", ""))) if cls == "sv" else \
                      ("%s x%s units" % (c.get("motif", ""), c.get("delta_motif_units_nearest_parent", "")))
        genes = hit or allg
    s = max([shet.get(g, 0.0) for g in genes], default=0.0) if genes else 0.0
    fl = ";".join(sorted(set(x for g in genes for x in [gene_flags(g)] if x)))
    rows.append(dict(cls=cls, new="NEW" if is_new(c) else "known", family=c.get("family_id", ""),
                     proband=c.get("sample_id", ""), chrom=c["chrom"], start=s0, end=e0,
                     ref=(c.get("ref") or "")[:12], alt=(c.get("alt") or "")[:12],
                     gene=gene, consequence=consequence, gnomad_af=gaf, s_het=("%.4f" % s) if s else "",
                     lof=lof, lof_tier=lof_tier, gene_flags=fl,
                     dnm_tier=c.get("dnm_tier", ""), decision=c.get("decision_reason", ""),
                     rule_score=c.get("rule_score", ""), rf_prob=c.get("rf_prob", ""),
                     phase_class=c.get("phase_class", ""), poo=c.get("parent_of_origin", ""),
                     svtype=c.get("svtype", ""), svlen=c.get("svlen", ""), trid=c.get("trid", ""),
                     strchive=c.get("strchive_locus", ""), mask=c.get("mask_overlap", ""),
                     flags=(c.get("flags") or "")[:60]))

COLS = ["cls", "new", "family", "proband", "chrom", "start", "end", "ref", "alt", "gene", "consequence", "gnomad_af",
        "s_het", "lof", "lof_tier", "gene_flags", "dnm_tier", "decision", "rule_score", "rf_prob", "phase_class",
        "poo", "svtype", "svlen", "trid", "strchive", "mask", "flags"]

def dump(name, rr):
    p = os.path.join(OUT, name)
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rr:
            w.writerow(r)
    log("wrote %s (%d rows)" % (p, len(rr)))

dump("denovo_tiered.module.tsv", rows)

# ------------------------------------------------------------------ what is reportable
def priority(r):
    if r["cls"] == "snv_indel":
        return bool(r["lof_tier"] in ("lof_t1", "lof_t2") or (r["gene_flags"] and r["lof_tier"]))
    return bool(r["gene_flags"]) or bool(r["strchive"])

pri = [r for r in rows if priority(r)]
pri.sort(key=lambda r: (r["new"] != "NEW", r["cls"] != "sv", r["lof_tier"] != "lof_t1", -float(r["s_het"] or 0)))
dump("denovo_priority.module.tsv", pri)

# ------------------------------------------------------------------ hand the missense calls to the cohort's tierer
# 06_clinical/miss_tier.py counts the four dbNSFP rankscore thresholds of the lab's v2 forward selection. It is driven
# by environment variables rather than reimplemented here, so a missense tier from this module means exactly what a
# missense tier means in the cohort report. It needs a python with duckdb (the `lrtier` env), not the module's.
mt = os.path.join(OUT, "denovo_tiered.module.small.tsv")
with open(mt, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["family", "proband", "chrom", "pos", "ref", "alt", "gene", "consequence"])
    for r in rows:
        if r["cls"] != "snv_indel" or not r["gene"]:
            continue
        w.writerow([r["family"], r["proband"], r["chrom"], r["start"], r["ref"], r["alt"], r["gene"], r["consequence"]])
log("wrote %s -- input for 06_clinical/miss_tier.py (set DENOVO_TSV to it, MISS_OUT to where the tiers land)" % mt)

print("\n== module call set ==")
for cls in ("snv_indel", "sv", "tr"):
    sub = [r for r in rows if r["cls"] == cls]
    print("  %-9s %5d calls  (tier1 %d, tier2 %d)  NEW vs original pipeline: %d"
          % (cls, len(sub), sum(1 for r in sub if r["dnm_tier"] == "1"), sum(1 for r in sub if r["dnm_tier"] == "2"),
             sum(1 for r in sub if r["new"] == "NEW")))
print("\n== reportable: NDD / constrained gene, or a known pathogenic repeat locus ==  (%d)" % len(pri))
hdr = "%-4s %-9s %-12s %-7s %-22s %-8s %-5s %-5s %s"
print(hdr % ("new", "class", "gene", "s_het", "consequence", "lof_tier", "tier", "rule", "decision"))
for r in pri[:80]:
    print(hdr % (r["new"][:4], r["cls"], r["gene"][:12], r["s_het"], str(r["consequence"])[:22], r["lof_tier"],
                 r["dnm_tier"], r["rule_score"], str(r["decision"])[:26]))
