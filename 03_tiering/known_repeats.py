#!/usr/bin/env python3
"""Systematic screen of every known pathogenic repeat locus in a family.

All 73 STRchive v2.4.4 disease loci x all family members, classified against the
published ranges. Four rules learned from the first run on the pilot family, which produced
17 pathogenic-range calls of which essentially all were false:

  1. SCORE THE PATHOGENIC MOTIF, NOT THE LOCUS LENGTH. TRGT reports MC as
     per-motif counts, underscore-separated, aligned to INFO/MOTIFS. At FAME7
     the father shows MC=230_0 over MOTIFS=TTTTA,TTTCA: 230 copies of the benign
     tract and zero of the pathogenic TTTCA. Scoring locus length calls that
     pathogenic. It is not.
  2. QC GATES EVERY VERDICT, ESPECIALLY PATHOGENIC. A pathogenic-range call on
     1-3 spanning reads is not a call. The first version only demoted benign
     verdicts, which is exactly backwards.
  3. FAMILY STRUCTURE IS EVIDENCE. An allele an unaffected parent or sibling
     carries at the same or greater size does not explain an affected child. Such
     loci stay in the output, demoted and labelled, never deleted.
  4. NOTHING IS DROPPED. Every locus x sample appears with a verdict and a reason,
     and the funnel is printed, so a silent zero is impossible to mistake for a
     negative result.

  02_known_repeats.py <family_trgt_vcf> <strchive.json> <out.tsv> [dropout_dir] [ped]
"""
import json, os, subprocess, sys, glob, collections

VCF, STRJSON, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
DROPDIR = sys.argv[4] if len(sys.argv) > 4 else None
PED = sys.argv[5] if len(sys.argv) > 5 else None
MIN_SD = int(os.environ.get("TR_MIN_SPANNING", "5"))
MIN_AP = float(os.environ.get("TR_MIN_PURITY", "0.80"))
SING = os.environ.get("SINGULARITY", "singularity")
CACHE = os.environ.get("SIF_CACHE", "")   # directory holding the pb_wdl_base .sif
os.environ["SINGULARITY_BIND"] = "/expanse:/expanse"

def sif(p):
    m = sorted(glob.glob(os.path.join(CACHE, "*%s*.sif" % p)))
    return m[0] if m else None

def run(argv, timeout=1800):
    p = subprocess.run(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                       universal_newlines=True, timeout=timeout)
    return p.returncode, p.stdout, p.stderr

B = sif("pb_wdl_base")
if not B:
    sys.exit("FATAL: no pb_wdl_base container")

loci = json.load(open(STRJSON))
by_id = {L["id"]: L for L in loci}

rc, sout, err = run([SING, "exec", B, "bcftools", "query", "-l", VCF])
if rc != 0:
    sys.exit("FATAL: bcftools query -l failed: %s" % err.strip())
samples = [s for s in sout.split() if s]

# affected status drives the family filter; without a ped we cannot demote
affected, parents = {}, {}
if PED and os.path.exists(PED):
    for line in open(PED):
        f = line.rstrip().split("\t")
        if len(f) >= 6:
            affected[f[1]] = (f[5] == "2")
            parents[f[1]] = (f[2], f[3])
probands = [s for s in samples if affected.get(s)]
unaffected = [s for s in samples if s in affected and not affected[s]]

TAB, NL = chr(92) + "t", chr(92) + "n"
FMT = ("%INFO/TRID" + TAB + "%CHROM" + TAB + "%POS" + TAB + "%INFO/MOTIFS"
       + "[" + TAB + "%MC|%SD|%AP|%AM]" + NL)
rc, out, err = run([SING, "exec", B, "bcftools", "query", "-f", FMT, VCF])
if rc != 0:
    sys.exit("FATAL: bcftools query failed: %s" % err.strip())

drop = collections.defaultdict(set)
if DROPDIR:
    for f in glob.glob(os.path.join(DROPDIR, "*", "*.txt")):
        s = os.path.basename(f).split(".")[0]
        for line in open(f, errors="replace"):
            p = line.split()
            if p:
                drop[s].add(p[0])

def as_list(v):
    if v is None:
        return []
    return v if isinstance(v, list) else [v]

def path_motif_index(motifs, L):
    """Index into MOTIFS of the disease-causing motif, per STRchive.

    Returns (index, note). A locus whose pathogenic motif is absent from the
    caller MOTIFS cannot be scored by motif; we say so rather than guessing."""
    pm = [m.upper() for m in as_list(L.get("pathogenic_motif_reference_orientation"))]
    if not pm:
        return (None, "no_pathogenic_motif_in_catalog")
    up = [m.upper() for m in motifs]
    for i, m in enumerate(up):
        if m in pm:
            return (i, "")
    return (None, "pathogenic_motif_absent_from_caller_motifs")

def allele_counts(mc_field, idx):
    """Per-allele count of the motif at index idx. MC is 'a_b_c,a_b_c'."""
    outv = []
    for allele in (mc_field or "").split(","):
        parts = allele.split("_")
        if idx is not None and idx < len(parts):
            tok = parts[idx]
        elif len(parts) == 1:
            tok = parts[0]
        else:
            tok = None
        if tok is None:
            continue
        try:
            outv.append(int(float(tok)))
        except ValueError:
            pass
    return outv

def nums(field):
    v = []
    for a in (field or "").split(","):
        try:
            v.append(float(a))
        except ValueError:
            pass
    return v

def classify(counts, L, distinct_path_motif=False):
    if not counts:
        return ("NO_CALL", None)
    if distinct_path_motif and max(counts) == 0:
        # the disease-causing motif is simply absent from this allele
        return ("benign", 0)
    inh = [str(i).upper() for i in as_list(L.get("inheritance"))]
    recessive = any(i in ("AR", "XR") for i in inh)
    allele = min(counts) if (recessive and len(counts) > 1) else max(counts)
    pmin, imin, bmax = (L.get("pathogenic_min"), L.get("intermediate_min"),
                        L.get("benign_max"))
    if pmin is not None and allele >= pmin:
        return ("PATHOGENIC_RANGE", allele)
    if imin is not None and allele >= imin:
        return ("INTERMEDIATE", allele)
    if bmax is not None and allele <= bmax:
        return ("benign", allele)
    return ("unclassified", allele)

# ---- D55 panel-invariance gate. A locus whose modal genotype covers ~all
# unrelated founders in the panel of normals cannot be informative: the four
# 'identical in every member' loci (pre-MIR7-2, TAF1, VWA1, POLG) were a
# locus-definition bug, not four findings. Table from cohort/trgt_invariance.sb.
INV_TABLE = os.environ.get('PON_INVARIANCE',
    '')
INV_FRAC, INV_MIN_CALLED = 0.95, 50
invariant = {}
if os.path.exists(INV_TABLE):
    with open(INV_TABLE) as _fh:
        _hdr = _fh.readline().rstrip(chr(10)).split(chr(9))
        _ix = dict((h, i) for i, h in enumerate(_hdr))
        for _ln in _fh:
            _f = _ln.rstrip(chr(10)).split(chr(9))
            try:
                _nc, _mf = int(_f[_ix['n_called']]), float(_f[_ix['modal_frac']])
            except (ValueError, IndexError):
                continue
            if _nc >= INV_MIN_CALLED and _mf >= INV_FRAC:
                invariant[_f[_ix['trid']]] = (_f[_ix['modal_mc']], _mf, _nc)
    print('  D55 panel-invariance gate ON: %d invariant loci (modal_frac>=%.2f, n>=%d)'
          % (len(invariant), INV_FRAC, INV_MIN_CALLED))
else:
    print('  WARNING: panel invariance table not found -> D55 gate OFF: %s' % INV_TABLE)

raw = collections.defaultdict(dict)     # trid -> sample -> record
meta = {}
for line in out.splitlines():
    f = line.rstrip().split("\t")
    if len(f) < 4 + len(samples):
        continue
    trid, chrom, pos, motifs_s = f[0], f[1], f[2], f[3]
    if trid not in by_id:
        continue
    L = by_id[trid]
    motifs = [m for m in motifs_s.split(",") if m]
    idx, note = path_motif_index(motifs, L)
    meta[trid] = (chrom, pos, motifs, idx, note, L)
    for i, s in enumerate(samples):
        p = (f[4 + i].split("|") + ["", "", "", ""])[:4]
        counts = allele_counts(p[0], idx)
        sd = [int(x) for x in nums(p[1])]
        ap, am = nums(p[2]), nums(p[3])
        ref_m = [m.upper() for m in as_list(
            L.get("reference_motif_reference_orientation"))]
        path_m = [m.upper() for m in as_list(
            L.get("pathogenic_motif_reference_orientation"))]
        distinct = bool(path_m) and bool(ref_m) and not set(path_m) & set(ref_m)
        verdict, allele = classify(counts, L, distinct)
        qc = []
        sd_call = sd
        if allele is not None and isinstance(counts, (list, tuple)) and allele in counts and len(sd) == len(counts):
            sd_call = [sd[counts.index(allele)]]      # the support of the allele being CALLED; max(SD) let the other allele carry it (T9)
        if not sd_call or min(sd_call) < MIN_SD:
            qc.append("low_spanning(SD=%s)" % (",".join(map(str, sd)) or "none"))
        if ap and min(ap) < MIN_AP:
            qc.append("low_purity(AP=%.2f)" % min(ap))
        if trid in drop.get(s, ()):
            qc.append("coverage_dropout")
        if note:
            qc.append(note)
        # QC gates EVERY verdict. A pathogenic call on 2 reads is not a call.
        if qc and verdict != "NO_CALL":
            verdict = "NOT_CALLABLE"
        elif qc:
            verdict = "NOT_CALLABLE"
        # D55 (refined): a PATHOGENIC/INTERMEDIATE call whose allele IS the
        # panel-modal allele is a misconfigured threshold, not a finding. Benign
        # calls at genuinely monomorphic loci (polyalanine tracts) stay benign.
        if trid in invariant and verdict in ('PATHOGENIC_RANGE', 'INTERMEDIATE') and allele is not None:
            _m, _fr, _n = invariant[trid]
            if str(allele) in set(_m.split(',')):
                qc.append('pathogenic_range_is_population_modal(modal=%s,frac=%.2f,n=%d)' % (_m, _fr, _n))
                verdict = 'NOT_CALLABLE'
        raw[trid][s] = {"counts": counts, "allele": allele, "sd": sd, "ap": ap,
                        "am": am, "qc": qc, "verdict": verdict,
                        "all_motif_counts": p[0]}

# ---- family filter: an allele an unaffected relative carries at the same or
# greater size cannot explain an affected child. Demote, label, never delete.
rows = []
for trid in sorted(raw):
    chrom, pos, motifs, idx, note, L = meta[trid]
    for s in samples:
        r = raw[trid][s]
        expl = ""
        if affected.get(s) and r["verdict"] in ("PATHOGENIC_RANGE", "INTERMEDIATE"):
            ua = [raw[trid][u]["allele"] for u in unaffected
                  if raw[trid].get(u) and raw[trid][u]["allele"] is not None]
            if ua and r["allele"] is not None and max(ua) >= r["allele"]:
                expl = "shared_with_unaffected(max=%s)" % max(ua)
                r["verdict"] = r["verdict"] + "_NOT_EXPLANATORY"
        rows.append({
            "trid": trid, "gene": L.get("gene", ""), "disease": L.get("disease", ""),
            "inheritance": ",".join(str(i) for i in as_list(L.get("inheritance"))),
            "mechanism": L.get("mechanism", ""), "chrom": chrom, "pos": pos,
            "sample": s, "affected": ("yes" if affected.get(s) else
                                      ("no" if s in affected else "unknown")),
            "motifs": ",".join(motifs),
            "pathogenic_motif": ",".join(as_list(
                L.get("pathogenic_motif_reference_orientation"))),
            "scored_motif_index": "" if idx is None else idx,
            "all_motif_counts": r["all_motif_counts"],
            "pathogenic_motif_counts": ",".join(map(str, r["counts"])),
            "scored_allele": "" if r["allele"] is None else r["allele"],
            "benign_max": L.get("benign_max", ""),
            "intermediate_min": L.get("intermediate_min", ""),
            "pathogenic_min": L.get("pathogenic_min", ""),
            "spanning_reads": ",".join(map(str, r["sd"])),
            "purity": ",".join("%.2f" % x for x in r["ap"]),
            "methylation": ",".join("%.2f" % x for x in r["am"]),
            "qc_flags": ";".join(r["qc"]), "family_note": expl,
            "verdict": r["verdict"],
        })

cols = ["trid", "gene", "disease", "inheritance", "mechanism", "chrom", "pos",
        "sample", "affected", "motifs", "pathogenic_motif", "scored_motif_index",
        "all_motif_counts", "pathogenic_motif_counts", "scored_allele",
        "benign_max", "intermediate_min", "pathogenic_min", "spanning_reads",
        "purity", "methylation", "qc_flags", "family_note", "verdict"]
with open(OUT, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")

tally = collections.Counter(r["verdict"] for r in rows)
print("  STRchive v2.4.4 disease loci : %d" % len(by_id))
print("  present in this callset      : %d" % len(raw))
print("  samples                      : %d (affected: %s)"
      % (len(samples), ", ".join(probands) or "none declared"))
print("  locus x sample verdicts:")
for v, n in tally.most_common():
    print("    %-34s %d" % (v, n))

act = [r for r in rows if r["verdict"] in ("PATHOGENIC_RANGE", "INTERMEDIATE")
       and r["affected"] == "yes"]
print("")
if act:
    print("  ACTIONABLE in an affected member (motif-scored, QC-passed, not shared):")
    for r in act:
        print("    %-18s %-9s %-30s allele=%s path>=%s SD=%s"
              % (r["trid"], r["gene"], r["disease"][:30], r["scored_allele"],
                 r["pathogenic_min"], r["spanning_reads"]))
else:
    print("  NO locus reaches intermediate or pathogenic range in an affected member")
    print("  after motif-specific scoring, QC gating and family filtering.")

dem = [r for r in rows if r["verdict"].endswith("_NOT_EXPLANATORY")]
if dem:
    print("  demoted -- present in an unaffected relative at the same or greater size:")
    for r in dem:
        print("    %-18s %-9s %s allele=%s  %s"
              % (r["trid"], r["gene"], r["sample"], r["scored_allele"], r["family_note"]))

nc = sorted(set((r["trid"], r["gene"]) for r in rows if r["verdict"] == "NOT_CALLABLE"))
print("  NOT_CALLABLE: %d distinct loci -- a normal-looking call here is unsupported"
      % len(nc))
for t, g in nc[:25]:
    ex = next(r for r in rows if r["trid"] == t and r["verdict"] == "NOT_CALLABLE")
    print("    %-18s %-9s %s" % (t, g, ex["qc_flags"][:70]))
print("  wrote %s" % OUT)
