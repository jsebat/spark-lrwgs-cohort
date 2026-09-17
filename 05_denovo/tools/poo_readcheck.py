"""Independent, read-backed parent-of-origin for called de novo SNVs.

For every tier-1 de novo SNV in one child, take the child's reads that carry the ALT allele and, on those same reads,
read the alleles at nearby informative sites (sites where the child is heterozygous and one parent lacks the allele).
An allele the father does not carry must have come from the mother, and vice versa. This uses NO haplotag, NO phase
block and NO orientation table: it is the read itself telling us which parent the DNM haplotype came from.

Three assignments are reported per DNM so the source of any disagreement is isolated:
  pipeline   parent_of_origin as written in the final table
  hp_orient  recomputed from the reads' HP tags + the child's orientation table (the pipeline's method, redone here)
  readback   the independent read-backed vote above

REF-carrying reads are the other haplotype and must vote for the OTHER parent: a built-in control per site.

Also emits a genetic-sex table for the trio from chrX (non-PAR) heterozygosity, to catch a parent swap.
"""
import csv
import os
import sys
from collections import Counter, defaultdict

import pysam

TAB, NL = chr(9), chr(10)
PAR = ((10001, 2781479), (155701383, 156030895))
MARGIN = 30000  # same nearest-segment rule as hapmatrix.LabelTables


def base_at(read, pos1):
    seq = read.query_sequence
    if seq is None or read.cigartuples is None:
        return None
    target = pos1 - 1
    ref, q = read.reference_start, 0
    for op, L in read.cigartuples:
        if op in (0, 7, 8):
            if ref <= target < ref + L:
                return seq[q + (target - ref)]
            ref += L; q += L
        elif op in (2, 3):
            if ref <= target < ref + L:
                return None
            ref += L
        elif op in (1, 4):
            q += L
        if ref > target:
            break
    return None


def vcf_index_for(p):
    for c in (p + ".tbi", p + ".csi",
              p.replace("/joint_small_variants_vcf/", "/joint_small_variants_vcf_index/") + ".tbi",
              p.replace("/joint_small_variants_vcf/", "/joint_small_variants_vcf_index/") + ".csi"):
        if os.path.exists(c):
            return c
    return None


def load_orientation(path):
    segs = defaultdict(list)
    if not os.path.exists(path):
        return segs
    for r in csv.DictReader(open(path, newline=""), delimiter=TAB):
        segs[(r["chrom"], int(r["phase_block_id"]))].append((int(r["start"]), int(r["end"]), r["orientation"]))
    for v in segs.values():
        v.sort()
    return segs


def hap1_is(segs, chrom, ps, pos):
    best, bestd = None, None
    for s, e, lab in segs.get((chrom, ps), ()):
        d = 0 if s <= pos <= e else (s - pos if pos < s else pos - e)
        if d <= MARGIN and (bestd is None or d < bestd):
            best, bestd = lab, d
    return {"HAP1_PAT": "P", "HAP1_MAT": "M"}.get(best)


def informative_sites(vf, chrom, a, b, child, father, mother):
    """pos -> {base: 'P'|'M'} for SNV sites where the child is het and a base is absent from one parent.
    Also returns the child's phased GT/PS at each site for the HP cross-check."""
    out, cgt = {}, {}
    for v in vf.fetch(chrom, max(0, a), b):
        if len(v.ref) != 1 or any(len(x) != 1 for x in v.alts or ()):
            continue
        sc, sf, sm = v.samples.get(child), v.samples.get(father), v.samples.get(mother)
        if sc is None or sf is None or sm is None:
            continue
        gc = tuple(x for x in (sc.get("GT") or ()) if x is not None)
        gf = tuple(x for x in (sf.get("GT") or ()) if x is not None)
        gm = tuple(x for x in (sm.get("GT") or ()) if x is not None)
        if len(gc) != 2 or gc[0] == gc[1] or len(gf) != 2 or len(gm) != 2:
            continue
        if (sc.get("GQ") or 0) < 20 or (sf.get("GQ") or 0) < 20 or (sm.get("GQ") or 0) < 20:
            continue
        F, M = set(gf), set(gm)
        alle = (v.ref,) + tuple(v.alts)
        d = {}
        for x in gc:
            if x in F and x not in M:
                d[alle[x]] = "P"
            elif x in M and x not in F:
                d[alle[x]] = "M"
        if d:
            out[v.pos] = d
            cgt[v.pos] = (alle[gc[0]], alle[gc[1]], bool(sc.phased), sc.get("PS"))
    return out, cgt


def genetic_sex(vf, samples, chrom="chrX", a=20_000_000, b=60_000_000):
    """(n_called, n_het, het_frac) per sample from non-PAR chrX SNVs with GQ>=20."""
    n = {s: [0, 0] for s in samples}
    for v in vf.fetch(chrom, a, b):
        if len(v.ref) != 1 or any(len(x) != 1 for x in v.alts or ()):
            continue
        for s in samples:
            sm = v.samples.get(s)
            if sm is None or (sm.get("GQ") or 0) < 20:
                continue
            g = tuple(x for x in (sm.get("GT") or ()) if x is not None)
            if len(g) != 2:
                continue
            n[s][0] += 1
            n[s][1] += int(g[0] != g[1])
    return {s: (c, h, round(h / c, 4) if c else None) for s, (c, h) in n.items()}


def main():
    family, child, manifest, final_dir, phase_dir, out = sys.argv[1:7]
    man = {r["sample_id"]: r for r in csv.DictReader(open(manifest, newline=""), delimiter=TAB)}
    row = man[child]
    father, mother = row["father_id"], row["mother_id"]
    # the manifest's Lustre paths are stale; resolve through the callset layout exactly as m2_review_family.sb does
    import glob
    O = os.path.join(os.environ["CALLSET_ROOT"], family, "out")
    bam_path = sorted(glob.glob(os.path.join(O, "merged_haplotagged_bam", "*", child + ".*.haplotagged.bam")))[0]
    bai_path = bam_path.replace("merged_haplotagged_bam/", "merged_haplotagged_bam_index/") + ".bai"
    bam = pysam.AlignmentFile(bam_path, "rb", index_filename=bai_path)
    vcf_path = sorted(glob.glob(os.path.join(O, "joint_small_variants_vcf", "*.vcf.gz")))[0]
    vf = pysam.VariantFile(vcf_path, index_filename=vcf_index_for(vcf_path))
    print("bam:", bam_path); print("vcf:", vcf_path)
    segs = load_orientation(os.path.join(phase_dir, family, child + ".orientation.tsv"))

    sx = genetic_sex(vf, [child, father, mother])
    with open(out + ".sex.tsv", "w") as fh:
        fh.write(TAB.join(["family", "sample", "role", "manifest_sex", "chrX_called", "chrX_het", "chrX_het_frac"]) + NL)
        for s, role in ((child, "child"), (father, "father"), (mother, "mother")):
            c, h, f = sx[s]
            fh.write(TAB.join(str(x) for x in (family, s, role, man[s]["sex"], c, h, f)) + NL)

    tsv = os.path.join(final_dir, "%s.%s.snv_indel.dnm.tsv" % (family, child))
    dnms = [r for r in csv.DictReader(open(tsv, newline=""), delimiter=TAB)
            if r["variant_class"] == "SNV" and r["dnm_call"] in ("YES", "CANDIDATE")]
    cols = ["family", "child", "chrom", "pos", "ref", "alt", "dnm_call", "dnm_tier", "pipeline_poo", "poo_reason", "flags",
            "n_alt_reads", "n_ref_reads",
            "alt_hp1", "alt_hp2", "alt_hp0", "alt_ps", "orient_hap1_is", "hp_orient_poo",
            "alt_votes_P", "alt_votes_M", "alt_sites_used", "readback_poo",
            "ref_votes_P", "ref_votes_M", "ref_control_ok",
            "hp1_votes_P", "hp1_votes_M", "hp2_votes_P", "hp2_votes_M", "readback_hap1_is",
            "agree_pipeline_readback", "agree_hporient_readback", "agree_orient_table_readback"]
    with open(out, "w") as fh:
        fh.write(TAB.join(cols) + NL)
        for r in dnms:
            chrom, pos, ref, alt = r["chrom"], int(r["start"]), r["ref"], r["alt"]
            reads = []
            for rd in bam.fetch(chrom, max(0, pos - 1), pos):
                if rd.is_unmapped or rd.is_secondary or rd.is_supplementary or rd.is_duplicate or rd.mapping_quality < 20:
                    continue
                b = base_at(rd, pos)
                if b == alt:
                    reads.append(("ALT", rd))
                elif b == ref:
                    reads.append(("REF", rd))
            alt_reads = [rd for s, rd in reads if s == "ALT"]
            ref_reads = [rd for s, rd in reads if s == "REF"]
            if not alt_reads:
                fh.write(TAB.join(str(x) for x in (family, child, chrom, pos, ref, alt, r["dnm_call"], r["dnm_tier"],
                                                   r["parent_of_origin"], r["poo_reason"], r["flags"], 0, len(ref_reads)) + ("",) * 21) + NL)
                continue
            a = min(rd.reference_start for rd in alt_reads + ref_reads)
            b = max(rd.reference_end or 0 for rd in alt_reads + ref_reads)
            info, cgt = informative_sites(vf, chrom, a, b, child, father, mother)

            def votes(rs):
                c, used = Counter(), 0
                for rd in rs:
                    hit = False
                    for p, d in info.items():
                        if p == pos or not (rd.reference_start < p <= (rd.reference_end or 0)):
                            continue
                        bb = base_at(rd, p)
                        if bb in d:
                            c[d[bb]] += 1; hit = True
                    used += hit
                return c, used

            va, used_a = votes(alt_reads)
            vr, _ = votes(ref_reads)
            hp = Counter(); ps = Counter()
            for rd in alt_reads:
                h = rd.get_tag("HP") if rd.has_tag("HP") else 0
                hp[h] += 1
                if rd.has_tag("PS"):
                    ps[rd.get_tag("PS")] += 1
            alt_ps = ps.most_common(1)[0][0] if ps else None
            o = hap1_is(segs, chrom, alt_ps, pos) if alt_ps is not None else None
            # pipeline's method redone: the ALT haplotype is the majority HP among ALT reads, mapped through orientation
            hp_major = 1 if hp[1] > hp[2] else (2 if hp[2] > hp[1] else None)
            hp_orient = None
            if o and hp_major:
                hp_orient = ("paternal" if o == "P" else "maternal") if hp_major == 1 else ("maternal" if o == "P" else "paternal")

            def call(c, minv=2, frac=0.8):
                n = c["P"] + c["M"]
                if n < minv:
                    return "undetermined"
                if c["P"] / n >= frac:
                    return "paternal"
                if c["M"] / n >= frac:
                    return "maternal"
                return "conflict"

            rb = call(va)
            rbr = call(vr)
            ctrl = "" if rb == "undetermined" or rbr == "undetermined" else int({"paternal": "maternal", "maternal": "paternal"}.get(rb) == rbr)
            # what do the reads say hap1 is here?  (all tagged reads, both alleles)
            h1, h2 = Counter(), Counter()
            for rd in alt_reads + ref_reads:
                h = rd.get_tag("HP") if rd.has_tag("HP") else 0
                if h not in (1, 2):
                    continue
                for p, d in info.items():
                    if p == pos or not (rd.reference_start < p <= (rd.reference_end or 0)):
                        continue
                    bb = base_at(rd, p)
                    if bb in d:
                        (h1 if h == 1 else h2)[d[bb]] += 1
            rb_h1 = call(h1)
            rb_hap1_is = {"paternal": "P", "maternal": "M"}.get(rb_h1, "")
            pip = r["parent_of_origin"]
            fh.write(TAB.join(str(x) if x is not None else "" for x in (
                family, child, chrom, pos, ref, alt, r["dnm_call"], r["dnm_tier"], pip, r["poo_reason"], r["flags"],
                len(alt_reads), len(ref_reads),
                hp[1], hp[2], hp[0], alt_ps, o, hp_orient,
                va["P"], va["M"], used_a, rb,
                vr["P"], vr["M"], ctrl,
                h1["P"], h1["M"], h2["P"], h2["M"], rb_hap1_is,
                int(pip == rb) if pip in ("paternal", "maternal") and rb in ("paternal", "maternal") else "",
                int(hp_orient == rb) if hp_orient and rb in ("paternal", "maternal") else "",
                int(o == rb_hap1_is) if o and rb_hap1_is else "")) + NL)
    print("done %s %s: %d SNV rows" % (family, child, len(dnms)))


if __name__ == "__main__":
    main()
