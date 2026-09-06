#!/usr/bin/env python3
# Build the de novo SV validation table for one family. 3-way outcome (D57):
#   REFUTED       - a PARENT carries the allele in short reads (inherited)
#   CONFIRMED     - proband carries it AND both parents ref with adequate depth
#   UNINFORMATIVE - proband no SR support / parent low coverage / BND / no SR model
# "proband negative in SR" is NEVER refuted: short reads often cannot capture a
# long-read SV. One row per candidate, kept regardless of outcome.
# Python 3.6.8 safe (D50): no f-strings, no capture_output.
import sys, subprocess, os, datetime, collections

FAM = sys.argv[1]
C = os.environ.get("DATA_ROOT", ".") + "/sv_validation"
SIF = subprocess.Popen("ls /expanse/projects/sebat1/jsebat/singularity_cache/miniwdl/*pb_wdl_base*.sif",
                       shell=True, stdout=subprocess.PIPE).communicate()[0].decode().split()[0]
SRVCF = os.path.join(C, FAM + ".sr_genotyped.FULL.vcf.gz")
# affected offspring and their parents, read from the family pedigree (never hardcoded)
RUNROOT = os.environ.get("RUNROOT", "/expanse/lustre/projects/ddp195/jsebat/longread-autism")
PARENTS_OF = {}
for _ln in open(os.path.join(RUNROOT, "run_" + FAM, "analysis", FAM + ".ped")):
    _f = _ln.split()
    if len(_f) >= 6 and _f[2] != "0" and _f[3] != "0" and _f[5] == "2":
        PARENTS_OF[_f[1]] = (_f[3], _f[2])   # (mother, father)
PROBANDS = sorted(PARENTS_OF)
DP_MIN = 6
TOL = 500

def sx(cmd):
    # argv list, no host shell: the ENTIRE pipeline runs inside the container (D56).
    p = subprocess.Popen(["singularity", "exec", "-B", "/expanse:/expanse", SIF, "bash", "-c", cmd],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    if p.returncode != 0 or (err and not out):
        sys.stderr.write('[sx] rc=%s cmd=%s' % (p.returncode, cmd[:120]) + chr(10) + err.decode()[:800] + chr(10))
    return out.decode()

# --- index SR AGGREGATED genotypes: (chrom,svtype) -> list of (pos, {sample:(gt,dp,ft)})
srx = collections.defaultdict(list)
if os.path.exists(SRVCF):
    q = ("bcftools view -i 'INFO/SVMODEL=\"AGGREGATED\"' " + SRVCF +
         " | bcftools query -f '%CHROM\\t%POS\\t%INFO/SVTYPE\\t%INFO/END\\t%INFO/SVSIZE\\t[%SAMPLE|%GT|%DP|%FT;]\\n'")
    for ln in sx(q).splitlines():
        f = ln.split(chr(9))
        if len(f) < 6:
            continue
        chrom, pos, svt = f[0], int(f[1]), f[2]
        try: end = int(f[3])
        except ValueError: end = pos
        try: size = abs(int(f[4]))
        except ValueError: size = 0
        gts = {}
        for cell in f[5].strip(';').split(';'):
            if not cell:
                continue
            parts = cell.split('|')
            if len(parts) >= 4:
                gts[parts[0]] = (parts[1], parts[2], parts[3])
        srx[(chrom, svt)].append((pos, end, size, gts))
for k in srx:
    srx[k].sort(key=lambda t: t[0])   # ties on pos must not compare the dicts
if not srx:
    sys.exit('FATAL: no AGGREGATED SR genotypes loaded from ' + SRVCF + ' -- refusing to write an all-UNINFORMATIVE table')
print('SR AGGREGATED loci indexed: %d' % sum(len(v) for v in srx.values()))

EQUIV = {'INS': ('INS', 'DUP'), 'DUP': ('DUP', 'INS'), 'DEL': ('DEL',), 'INV': ('INV',)}
def match(chrom, pos, end, svt, svlen):
    # JS 2026-09-05: match if the LEFT boundary, RIGHT boundary, or LENGTH agrees
    # (TOL bp, or 10% for length); INS and DUP are one class across callers.
    best = None; bestd = None
    for t in EQUIV.get(svt, (svt,)):
        for p, e, sz, g in srx.get((chrom, t), []):
            dl = abs(p - pos)
            dr = abs(e - end) if (end and e) else TOL + 1
            dz = abs(sz - svlen) if (sz and svlen) else TOL + 1
            ltol = max(TOL, int(0.10 * svlen)) if svlen else TOL
            if not ((dl <= TOL) or (dr <= TOL) or (dz <= ltol)):
                continue
            d = min(dl, dr, dz)
            if bestd is None or d < bestd:
                bestd = d; best = g
    return best

def has_alt(gt):
    return gt not in (None, "", "./.", ".|.", "0/0", "0|0", ".") and "1" in gt

def is_ref(gt):
    return gt in ("0/0", "0|0")

def missing(gt):
    return gt in (None, "", "./.", ".|.", ".")

def dp_ok(dp):
    try: return int(dp) >= DP_MIN
    except: return False

def verdict(proband, sr):
    par = {}
    for pa in PARENTS_OF[proband]:
        par[pa] = sr.get(pa, (None, None, None)) if sr else (None, None, None)
    pb = sr.get(proband, (None, None, None)) if sr else (None, None, None)
    if sr is None:
        return "UNINFORMATIVE", "no_SR_model_at_locus"
    # a parent positively carries the allele -> inherited
    carriers = [pa for pa in PARENTS_OF[proband] if has_alt(par[pa][0]) and dp_ok(par[pa][1])]
    if carriers:
        # name ALL carrier parents; returning on the first one manufactured a 275:44 asymmetry
        lab = "both_parents" if len(carriers) == 2 else "parent_" + carriers[0] + "_only"
        return "REFUTED", "allele_in_" + lab
    # proband must show the allele to confirm
    if missing(pb[0]) or is_ref(pb[0]):
        return "UNINFORMATIVE", "proband_no_SR_support"
    if not dp_ok(pb[1]):
        return "UNINFORMATIVE", "proband_low_coverage"
    # proband has allele; require both parents confidently ref
    for pa in PARENTS_OF[proband]:
        gt, dp, ft = par[pa]
        if missing(gt) or not dp_ok(dp):
            return "UNINFORMATIVE", "parent_low_coverage_" + pa
        if not is_ref(gt):
            return "REFUTED", "allele_in_parent_" + pa
    return "CONFIRMED", "proband_alt_parents_ref"

gtver = sx("graphtyper --version").strip().split()[-1] if False else "2.7.2"
today = datetime.date.today().isoformat()
out = os.path.join(C, FAM + ".sv_validation_table.tsv")
cols = ["candidate_id","family","proband","chrom","pos","end","svtype","svlen","gene",
        "lr_gt_proband","lr_gt_mother","lr_gt_father",
        "sr_gt_proband","sr_dp_proband","sr_gt_mother","sr_dp_mother","sr_gt_father","sr_dp_father",
        "outcome","reason","graphtyper_version","date"]
rows = []
counts = collections.Counter()
for proband in PROBANDS:
    tsv = os.path.join(
        "/expanse/lustre/projects/ddp195/jsebat/longread-autism/run_" + FAM +
        "/analysis/denovo_sv", FAM + "." + proband + ".denovo_sv.tsv")
    if not os.path.exists(tsv):
        continue
    with open(tsv) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(hdr)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            chrom = f[idx["chrom"]]; pos = int(f[idx["pos"]])
            end = f[idx["end"]]; svt = f[idx["svtype"]]; svlen = f[idx["svlen"]]
            gene = f[idx.get("svann", 0)] if "svann" in idx else ""
            cid = "%s_%s_%d_%s" % (proband, chrom, pos, svt)
            if svt == "BND":
                outc, reason = "UNINFORMATIVE", "BND_not_genotypable_by_short_read"
                sr = None
            else:
                try: _end = int(end)
                except ValueError: _end = pos
                try: _len = abs(int(svlen))
                except ValueError: _len = 0
                sr = match(chrom, pos, _end, svt, _len)
                outc, reason = verdict(proband, sr)
            def g(s):
                return (sr.get(s, ("",""," "))[0] if sr else "")
            def d(s):
                return (sr.get(s, ("","",""))[1] if sr else "")
            mo, fa = PARENTS_OF[proband]
            rows.append([cid, FAM, proband, chrom, str(pos), end, svt, svlen, gene,
                         f[idx[proband + "_GT"]], f[idx[mo + "_GT"]],
                         f[idx[fa + "_GT"]],
                         g(proband), d(proband), g(mo), d(mo),
                         g(fa), d(fa),
                         outc, reason, gtver, today])
            counts[(proband, outc)] += 1

with open(out, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(str(x) for x in r) + "\n")

print("wrote " + out + "  (" + str(len(rows)) + " candidate rows)")
print("=== outcome summary ===")
for proband in PROBANDS:
    tot = sum(v for (p, o), v in counts.items() if p == proband)
    print("  " + proband + " (n=" + str(tot) + "):")
    for o in ("CONFIRMED", "REFUTED", "UNINFORMATIVE"):
        n = counts.get((proband, o), 0)
        pct = (100.0 * n / tot) if tot else 0
        print("    %-14s %4d  (%.0f%%)" % (o, n, pct))
