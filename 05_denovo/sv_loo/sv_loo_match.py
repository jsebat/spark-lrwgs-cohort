# LOO SV panel match (D62). JS rule 2026-09-05: same chrom; INS and DUP are one class;
# accept if LEFT boundary <= TOL, or RIGHT boundary <= TOL, or LENGTH within 10% (>= TOL),
# or reciprocal overlap >= 0.5 for interval types. The panel record must lie within a
# len+TOL window of the family record so a length-only match cannot come from anywhere.
# Records WHICH basis matched. Takes the max-AC panel match (conservative for an
# artifact filter). Python 3.6 safe. Runs INSIDE the pb_wdl_base container (bcftools).
import sys, bisect, collections, subprocess
FAM_VCF, PANEL_VCF, OUT, DN_TABLE = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
TOL, PCT, RO_MIN = 500, 0.10, 0.5
TAB, NL = chr(9), chr(10)
EQ = {"INS": "INSDUP", "DUP": "INSDUP", "DEL": "DEL", "INV": "INV"}

def q(vcf, fmt):
    p = subprocess.Popen(["bcftools", "query", "-f", fmt, vcf],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    o, e = p.communicate()
    if p.returncode != 0:
        sys.exit("bcftools query failed on " + vcf + NL + e.decode()[:500])
    return o.decode().splitlines()

def ival(pos, end, svt, ln):
    if svt == "INS":
        return pos, pos
    e = end if (end and end > pos) else pos + abs(ln)
    return pos, e

def toint(x, default=0):
    try:
        return int(x.split(",")[0])
    except (ValueError, AttributeError):
        return default

# ---- panel index: chrom -> sorted (pos, end, cls, len, ac, an, id)
panel = collections.defaultdict(list)
fmt = TAB.join(["%CHROM", "%POS", "%INFO/END", "%INFO/SVTYPE", "%INFO/SVLEN", "%INFO/AC", "%INFO/AN", "%ID"]) + NL
for ln in q(PANEL_VCF, fmt):
    f = ln.split(TAB)
    if len(f) < 8 or f[3] not in EQ:
        continue
    pos = toint(f[1]); end = toint(f[2]); L = abs(toint(f[4])); ac = toint(f[5]); an = toint(f[6])
    if pos <= 0:
        continue
    s, e = ival(pos, end, f[3], L)
    panel[f[0]].append((s, e, EQ[f[3]], L, ac, an, f[7]))
for c in panel:
    panel[c].sort()
ppos = dict((c, [r[0] for r in panel[c]]) for c in panel)
sys.stderr.write("panel SV records indexed: %d" % sum(len(v) for v in panel.values()) + NL)

def best_match(chrom, pos, end, svt, L):
    cls = EQ.get(svt)
    if cls is None or chrom not in panel:
        return None
    lo = bisect.bisect_left(ppos[chrom], pos - L - TOL)
    hi = bisect.bisect_right(ppos[chrom], end + TOL)
    best = None
    for (pp, pe, pcls, pL, ac, an, pid) in panel[chrom][lo:hi]:
        if pcls != cls:
            continue
        basis = []
        if abs(pp - pos) <= TOL:
            basis.append("left")
        if end > pos and pe > pp and abs(pe - end) <= TOL:
            basis.append("right")
        if L and pL and abs(pL - L) <= max(TOL, int(PCT * L)):
            basis.append("length")
        if end > pos and pe > pp:
            ov = min(end, pe) - max(pos, pp)
            if ov > 0 and ov / float(max(end - pos, pe - pp, 1)) >= RO_MIN:
                basis.append("ro")
        if not basis:
            continue
        if best is None or ac > best[0]:
            best = (ac, an, pid, "+".join(basis))
    return best

# ---- family SVs -> panel AC
fmt2 = TAB.join(["%CHROM", "%POS", "%INFO/END", "%INFO/SVTYPE", "%INFO/SVLEN", "%ID"]) + NL
rows = {}
n = m = 0
with open(OUT, "w") as fh:
    fh.write(TAB.join(["chrom", "pos", "end", "svtype", "svlen", "fam_id", "panel_id",
                       "match_basis", "PON_AC", "PON_AN", "PON_AF"]) + NL)
    for ln in q(FAM_VCF, fmt2):
        f = ln.split(TAB)
        if len(f) < 6 or f[3] not in EQ:
            continue
        pos = toint(f[1]); end = toint(f[2]); L = abs(toint(f[4]))
        if pos <= 0:
            continue
        s, e = ival(pos, end, f[3], L)
        n += 1
        b = best_match(f[0], s, e, f[3], L)
        if b:
            m += 1
            ac, an, pid, basis = b
            af = ("%.4f" % (ac / float(an))) if an else ""
        else:
            ac = an = 0; pid = basis = af = ""
        rows[(f[0], pos, f[3])] = (ac, an, basis)
        fh.write(TAB.join([f[0], str(pos), str(e), f[3], str(L), f[5], pid, basis,
                           str(ac), str(an), af]) + NL)
sys.stderr.write("family SV records: %d  matched to panel: %d (%.0f%%)" % (n, m, 100.0 * m / max(n, 1)) + NL)

# ---- payoff: de novo SV candidates (GraphTyper table) x LOO-panel recurrence
hdr = None
cnt = collections.Counter()
ex = []
for ln in open(DN_TABLE):
    f = ln.rstrip(NL).split(TAB)
    if hdr is None:
        hdr = dict((h, i) for i, h in enumerate(f))
        continue
    key = (f[hdr["chrom"]], toint(f[hdr["pos"]]), f[hdr["svtype"]])
    r = rows.get(key)
    if r is None:
        rec = "not_in_family_vcf_or_BND"
    elif r[0] >= 3:
        rec = "PON_AC>=3"
    elif r[0] > 0:
        rec = "PON_AC=1-2"
    else:
        rec = "PON_AC=0"
    cnt[(f[hdr["outcome"]], rec)] += 1
    if r and r[0] >= 3 and len(ex) < 8:
        ex.append("%s:%s %s len=%s GraphTyper=%s PON_AC=%d/%d via %s"
                  % (key[0], key[1], key[2], f[hdr["svlen"]], f[hdr["outcome"]], r[0], r[1], r[2]))
print("=== de novo SV candidates: GraphTyper outcome x LOO-panel recurrence ===")
for (o, rec), v in sorted(cnt.items()):
    print("  %-14s %-26s %d" % (o, rec, v))
print("=== examples PON_AC>=3 ===")
for x in ex:
    print("  " + x)
