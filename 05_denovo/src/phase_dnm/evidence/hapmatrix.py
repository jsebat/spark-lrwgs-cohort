"""M2 core — the six-haplotype evidence matrix and the rule layer (DESIGN P7, P8, P10, P11).

Pure Python. The read-level adapters (evidence/readers.py, pysam) turn a candidate into a list of ReadObs - one per
spanning read in each of the three samples, carrying the read's HP/PS tags and the class adapter's verdict
(ALT / REF / AMB) - and everything from there on is class-agnostic:

  build_matrix   ReadObs + M1 label tables -> rows keyed by (sample, HP) with per-row counts and read-quality
                 summaries, plus the labels that turn HP into parent of origin (child) or transmitted /
                 untransmitted (parents) at this position.
  features       the registry's D block. Everything here that is rf_safe is a function of the child's two rows
                 (sorted alt-first, never labelled) and of parent rows made symmetric over both parents and both
                 haplotypes. Transmission-dependent quantities are computed separately and marked rf_safe=False.
  classify       the P8 rule table, in order; parent of origin with a reason code; flags.

Identifiers never appear in code; read names arrive already hashed.
"""
from __future__ import annotations

import bisect
import csv
import statistics as st
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple

ROLES = ("C", "F", "M")
CLASSES = ("phase_conflict_artifact", "inherited_missed_in_parent", "parental_mosaic_transmitted",
           "child_postzygotic_mosaic", "germline_DNM_phased", "germline_DNM_unphased", "inconclusive")


@dataclass(slots=True)
class ReadObs:
    role: str                    # 'C' | 'F' | 'M'
    hp: Optional[int]            # 1 | 2 | None (untagged)
    ps: Optional[int]
    support: str                 # 'ALT' | 'REF' | 'AMB'
    mapq: int = 60
    nm_rate: Optional[float] = None
    clipped: bool = False
    supplementary: bool = False
    read_len: int = 0
    al: Optional[int] = None     # TR: per-read allele length
    rid: str = ""
    rq: Optional[float] = None   # HiFi read quality tag


@dataclass
class HapParams:
    k: Tuple[int, ...] = (3, 5)
    error_reads: int = 1
    error_frac: float = 0.05
    min_mapq: int = 20
    amb_flag_frac: float = 0.3      # flag AMBIGUOUS_READS when a haplotype row has > this fraction of unreadable reads


@dataclass
class ClassParams:
    germline_min_hap_frac: float = 0.8
    mosaic_max_hap_frac: float = 0.7
    mosaic_min_dp: int = 15
    mosaic_min_minority_reads: int = 3
    inherited_min_frac_on_T: float = 0.3
    min_alt_reads: int = 3          # evidence floor for any germline / mosaic class: alt reads on A (+ untagged for unphased)
    working_k: int = 5
    del_max_ratio: float = 0.7             # child depth inside / flank at or below this = one haplotype lost (expect ~0.5)
    del_parent_min_ratio: float = 0.85     # a parent at or above this is NOT depleted, so the event is not inherited from them
    del_max_het_persistence: float = 0.35  # heterozygous sites inside / flanking rate; ~0 for a constitutional deletion



# ----------------------------------------------------------------------------------------------
# M1 label tables
# ----------------------------------------------------------------------------------------------
class LabelTables:
    """Child orientation (block+pos -> hap1 is P/M), parent transmission (parent, block, pos -> HAP1/HAP2
    transmitted), change points (for distance flags). All keyed by (chrom, phase_block_id) then position."""

    def __init__(self):
        self.orient: Dict[Tuple[str, int], List[Tuple[int, int, str]]] = defaultdict(list)
        self.trans: Dict[Tuple[str, str, int], List[Tuple[int, int, str]]] = defaultdict(list)
        self.changes: Dict[Tuple[str, str], List[int]] = defaultdict(list)
        self.child_switches: Dict[str, List[int]] = defaultdict(list)

    @classmethod
    def load(cls, orientation_tsv: Optional[str], transmission_tsv: Optional[str], changepoints_tsv: Optional[str]) -> "LabelTables":
        t = cls()
        if orientation_tsv:
            for r in csv.DictReader(open(orientation_tsv, newline=""), delimiter="\t"):
                t.orient[(r["chrom"], int(r["phase_block_id"]))].append((int(r["start"]), int(r["end"]), r["orientation"]))
                if int(r.get("switch_pos", 0) or 0) > 0:
                    t.child_switches[r["chrom"]].append(int(r["switch_pos"]))
        if transmission_tsv:
            for r in csv.DictReader(open(transmission_tsv, newline=""), delimiter="\t"):
                t.trans[(r["parent"], r["chrom"], int(r["parent_phase_block_id"]))].append((int(r["start"]), int(r["end"]), r["transmitted"]))
        if changepoints_tsv:
            for r in csv.DictReader(open(changepoints_tsv, newline=""), delimiter="\t"):
                t.changes[(r["parent"], r["chrom"])].append((int(r["left_pos"]) + int(r["right_pos"])) // 2)
        for d in (t.orient, t.trans):
            for v in d.values():
                v.sort()
        for v in t.changes.values():
            v.sort()
        for v in t.child_switches.values():
            v.sort()
        return t

    # Segments are bounded by the block's heterozygous sites; reads carrying the block's PS extend past the
    # first/last het by up to a read length. A candidate there belongs to the block, so the NEAREST segment of the
    # same phase set applies within `margin` bp (first smoke run: 35% of candidates were left unoriented by an
    # exact-cover lookup against a cohort-wide 3% ambiguity).
    MARGIN = 30000

    @staticmethod
    def _pick(segs, pos: int, margin: int):
        best, bestd = None, None
        for s, e, lab in segs:
            d = 0 if s <= pos <= e else (s - pos if pos < s else pos - e)
            if d <= margin and (bestd is None or d < bestd):
                best, bestd = lab, d
                if d == 0:
                    break
        return best

    def child_hap1_is(self, chrom: str, ps: Optional[int], pos: int) -> Optional[str]:
        if ps is None:
            return None
        o = self._pick(self.orient.get((chrom, ps), ()), pos, self.MARGIN)
        return "P" if o == "HAP1_PAT" else ("M" if o == "HAP1_MAT" else None)

    def parent_transmitted(self, parent: str, chrom: str, ps: Optional[int], pos: int) -> Optional[int]:
        if ps is None:
            return None
        tr = self._pick(self.trans.get((parent, chrom, ps), ()), pos, self.MARGIN)
        return 1 if tr == "HAP1" else (2 if tr == "HAP2" else None)

    def dist_change_point(self, parent: str, chrom: str, pos: int) -> Optional[int]:
        lst = self.changes.get((parent, chrom))
        if not lst:
            return None
        i = bisect.bisect_left(lst, pos)
        c = [abs(lst[j] - pos) for j in (i - 1, i) if 0 <= j < len(lst)]
        return min(c) if c else None

    def dist_child_switch(self, chrom: str, pos: int) -> Optional[int]:
        lst = self.child_switches.get(chrom)
        if not lst:
            return None
        i = bisect.bisect_left(lst, pos)
        c = [abs(lst[j] - pos) for j in (i - 1, i) if 0 <= j < len(lst)]
        return min(c) if c else None


# ----------------------------------------------------------------------------------------------
# matrix
# ----------------------------------------------------------------------------------------------
@dataclass
class Row:
    dp: int = 0
    alt: int = 0
    ref: int = 0
    amb: int = 0
    mapq: List[int] = field(default_factory=list)
    nm_alt: List[float] = field(default_factory=list)
    nm_ref: List[float] = field(default_factory=list)
    clip_alt: int = 0
    supp_alt: int = 0
    al: List[int] = field(default_factory=list)
    len_alt: List[int] = field(default_factory=list)
    rq_alt: List[float] = field(default_factory=list)

    def add(self, r: ReadObs):
        self.dp += 1
        self.mapq.append(r.mapq)
        if r.support == "ALT":
            self.alt += 1
            if r.nm_rate is not None: self.nm_alt.append(r.nm_rate)
            self.clip_alt += r.clipped
            self.supp_alt += r.supplementary
            if r.read_len: self.len_alt.append(r.read_len)
            if r.rq is not None: self.rq_alt.append(r.rq)
        elif r.support == "REF":
            self.ref += 1
            if r.nm_rate is not None: self.nm_ref.append(r.nm_rate)
        else:
            self.amb += 1
        if r.al is not None:
            self.al.append(r.al)

    @property
    def n(self) -> int:
        """Readable reads: those that report REF or ALT at the site. A read that is deleted / soft-clipped / carries a
        third base there (AMB) counts towards depth but does not OBSERVE the allele - observability (k), the error
        allowance and the mosaic depth floors are all defined on readable reads (P7; cohort finding 2026-09-13:
        low-GQ SNVs next to indels had parental haplotypes with 6 of 7 reads unreadable, counted as 'observed')."""
        return self.alt + self.ref

    @property
    def alt_frac(self) -> Optional[float]:
        n = self.alt + self.ref
        return self.alt / n if n else None

    @property
    def amb_frac(self) -> Optional[float]:
        return self.amb / self.dp if self.dp else None

    def summary(self) -> dict:
        return dict(dp=self.dp, alt=self.alt, ref=self.ref, amb=self.amb,
                    mapq_mean=round(st.mean(self.mapq), 1) if self.mapq else None,
                    mapq0_frac=round(sum(1 for m in self.mapq if m == 0) / len(self.mapq), 3) if self.mapq else None,
                    nm_alt_mean=round(st.mean(self.nm_alt), 4) if self.nm_alt else None,
                    nm_ref_mean=round(st.mean(self.nm_ref), 4) if self.nm_ref else None,
                    clip_alt_frac=round(self.clip_alt / self.alt, 3) if self.alt else None,
                    al_mean=round(st.mean(self.al), 1) if self.al else None,
                    al_sd=round(st.pstdev(self.al), 1) if len(self.al) > 1 else None, al_n=len(self.al))


@dataclass
class Matrix:
    rows: Dict[Tuple[str, int], Row]              # (role, hp) for hp in (1, 2)
    untagged: Dict[str, Row]                      # role -> untagged reads
    child_hap1_is: Optional[str]                  # 'P' | 'M' | None (child block unoriented at this position)
    transmitted: Dict[str, Optional[int]]         # 'F'/'M' -> 1 | 2 | None
    n_reads: int

    def row(self, role: str, hp: int) -> Row:
        return self.rows[(role, hp)]

    def child_labelled(self) -> Dict[str, Row]:
        """{'CP': row, 'CM': row} when the child's block is oriented here, else {}."""
        if self.child_hap1_is == "P":
            return {"CP": self.row("C", 1), "CM": self.row("C", 2)}
        if self.child_hap1_is == "M":
            return {"CP": self.row("C", 2), "CM": self.row("C", 1)}
        return {}

    def parent_labelled(self, parent: str) -> Dict[str, Row]:
        """{'T': row, 'U': row} for a parent whose transmitted haplotype is known here, else {}."""
        t = self.transmitted.get(parent)
        if t is None:
            return {}
        return {"T": self.row(parent, t), "U": self.row(parent, 3 - t)}

    @classmethod
    def from_evidence_row(cls, row: Dict[str, object]) -> "Matrix":
        """Rebuild the count part of the matrix from an evidence-table row (review.py columns), so the rule layer can
        be re-run without touching the BAMs. Read-quality lists (MAPQ, NM, allele lengths) are not recoverable and
        stay empty; their summaries are carried over from the original row by the caller."""
        def i(k):
            try:
                return int(float(row.get(k)))  # type: ignore[arg-type]
            except (TypeError, ValueError):
                return 0
        rows = {}
        for (role, hp), pre in ((("C", 1), "C1"), (("C", 2), "C2"), (("F", 1), "F1"), (("F", 2), "F2"), (("M", 1), "M1"), (("M", 2), "M2")):
            rows[(role, hp)] = Row(dp=i(pre + "_dp"), alt=i(pre + "_alt"), ref=i(pre + "_ref"), amb=i(pre + "_amb"))
        untagged = {role: Row(dp=i(role + "_untagged_dp"), alt=i(role + "_untagged_alt")) for role in ROLES}
        h = row.get("child_hap1_is")
        tr = {}
        for par in ("F", "M"):
            v = str(row.get(par + "_transmitted_hap"))
            tr[par] = int(v) if v in ("1", "2") else None
        return cls(rows=rows, untagged=untagged, child_hap1_is=h if h in ("P", "M") else None, transmitted=tr, n_reads=i("n_reads_used"))


def build_matrix(reads: Iterable[ReadObs], labels: LabelTables, chrom: str, pos: int, p: HapParams) -> Matrix:
    rows = {(role, hp): Row() for role in ROLES for hp in (1, 2)}
    untagged = {role: Row() for role in ROLES}
    child_ps: Dict[int, int] = defaultdict(int)
    parent_ps: Dict[str, Dict[int, int]] = {"F": defaultdict(int), "M": defaultdict(int)}
    n = 0
    for r in reads:
        if r.mapq < p.min_mapq:
            continue
        n += 1
        if r.hp in (1, 2):
            rows[(r.role, r.hp)].add(r)
            if r.ps is not None:
                (child_ps if r.role == "C" else parent_ps[r.role])[r.ps] += 1
        else:
            untagged[r.role].add(r)
    # the block at this position is the phase set most of the tagged reads carry (reads from a neighbouring
    # block can overhang the boundary; the majority PS is the block whose labels apply here)
    cps = max(child_ps, key=child_ps.get) if child_ps else None
    transmitted = {}
    for parent in ("F", "M"):
        pps = max(parent_ps[parent], key=parent_ps[parent].get) if parent_ps[parent] else None
        transmitted[parent] = labels.parent_transmitted(parent, chrom, pps, pos)
    return Matrix(rows=rows, untagged=untagged, child_hap1_is=labels.child_hap1_is(chrom, cps, pos),
                  transmitted=transmitted, n_reads=n)


# ----------------------------------------------------------------------------------------------
# features (registry D block) — rf_safe part is symmetric / unlabelled by construction
# ----------------------------------------------------------------------------------------------
def _e(dp: int, p: HapParams) -> float:
    return max(p.error_reads, p.error_frac * dp)


def features(m: Matrix, p: HapParams) -> Dict[str, object]:
    f: Dict[str, object] = {}
    c1, c2 = m.row("C", 1), m.row("C", 2)
    A, O = (c1, c2) if c1.alt >= c2.alt else (c2, c1)           # alt-carrying child haplotype first, never labelled
    f.update(c_alt_hapA=A.alt, c_alt_hapO=O.alt, c_dp_hapA=A.dp, c_dp_hapO=O.dp,
             c_alt_hap_frac=round(A.alt_frac, 4) if A.alt_frac is not None else None,
             c_alt_confined=int(O.alt <= _e(O.n, p)),
             c_alt_tagged_frac=round((c1.alt + c2.alt) / (c1.alt + c2.alt + m.untagged["C"].alt), 3)
             if (c1.alt + c2.alt + m.untagged["C"].alt) else None,
             c_untagged_dp=m.untagged["C"].dp, c_untagged_alt=m.untagged["C"].alt)
    for k in p.k:
        f["c_both_haps_obs_k%d" % k] = int(c1.n >= k and c2.n >= k)
    prow = [m.row(par, hp) for par in ("F", "M") for hp in (1, 2)]
    f.update(p_min_hap_dp=min(r.dp for r in prow), p_max_alt_any_hap=max(r.alt for r in prow),
             p_sum_alt_all_haps=sum(r.alt for r in prow),
             p_n_haps_with_alt=sum(1 for r in prow if r.alt > _e(r.n, p)),
             p_max_alt_hap_frac=round(max((r.alt_frac or 0.0) for r in prow), 4),
             p_untagged_dp_max=max(m.untagged["F"].dp, m.untagged["M"].dp),
             p_untagged_alt_max=max(m.untagged["F"].alt, m.untagged["M"].alt))
    for k in p.k:
        f["p_n_haps_obs_k%d" % k] = sum(1 for r in prow if r.n >= k)
        f["hap_obs_k%d" % k] = f["p_n_haps_obs_k%d" % k] + int(c1.n >= k) + int(c2.n >= k)
    # unreadable reads (deleted / clipped / third base at the site): an indel-context / artefact signal, and the
    # reason a haplotype with depth can still be unobserved
    f.update(c_amb_frac_hapA=round(A.amb_frac, 3) if A.amb_frac is not None else None,
             p_amb_frac_max=round(max((r.amb_frac or 0.0) for r in prow), 3) if any(r.dp for r in prow) else None)
    # read-quality summaries for the child rows (alt vs ref)
    sA = A.summary()
    f.update(c_alt_mapq_mean=sA["mapq_mean"], c_alt_nm_rate=sA["nm_alt_mean"], c_ref_nm_rate=sA["nm_ref_mean"], c_alt_clip_frac=sA["clip_alt_frac"],
             c_alt_mapq0_frac=sA["mapq0_frac"], c_alt_supp_frac=round(A.supp_alt / A.alt, 3) if A.alt else None,
             c_alt_readlen_median=int(st.median(A.len_alt)) if A.len_alt else None,
             c_alt_rq_mean=round(st.mean(A.rq_alt), 4) if A.rq_alt else None)
    # TR: per-haplotype allele lengths
    if A.al or O.al:
        f.update(c_tr_al_hapA_mean=A.summary()["al_mean"], c_tr_al_hapA_sd=A.summary()["al_sd"], c_tr_al_hapO_mean=O.summary()["al_mean"])
    return f


def transmission_features(m: Matrix, p: HapParams) -> Dict[str, object]:
    """rf_safe = FALSE. Parent of origin of the alt-carrying child haplotype, and alt reads on the TRANSMITTED /
    UNTRANSMITTED haplotype of that parent. Only defined in a real trio with M1 labels at this position."""
    out: Dict[str, object] = dict(parent_of_origin="undetermined", poo_reason=None, poo_confidence=None,
                                  t_alt_reads=None, t_dp=None, u_alt_reads=None, u_dp=None, nt_parent_alt_reads=None, t_hap_resolved=0)
    lab = m.child_labelled()
    if not lab:
        out["poo_reason"] = "CHILD_BLOCK_UNORIENTED"
        return out
    cp, cm = lab["CP"], lab["CM"]
    if cp.alt == 0 and cm.alt == 0:
        out["poo_reason"] = "NO_TAGGED_ALT_READS"
        return out
    if cp.alt > _e(cp.n, p) and cm.alt > _e(cm.n, p):
        out["poo_reason"] = "ALT_ON_BOTH_CHILD_HAPLOTYPES"
        return out
    if cp.alt == cm.alt:
        # A TIE carries no information about which parent contributed, and `>=` silently resolved every one of them
        # PATERNAL: genome-wide, 39753 rows sit at exactly equal alt counts and 100% of them were called paternal,
        # 0 maternal. None reach the current call set, so today's paternal fractions are unaffected, but the bias is
        # latent and would appear the moment a threshold loosened or the rows were used in bulk.
        out["poo_reason"] = "TIED_ALT_READS"
        return out
    origin, other = ("F", "M") if cp.alt > cm.alt else ("M", "F")
    alt_on, alt_off = (cp.alt, cm.alt) if origin == "F" else (cm.alt, cp.alt)
    out["parent_of_origin"] = "paternal" if origin == "F" else "maternal"
    # NOT a confidence in the parent assignment, despite the name. The ALT_ON_BOTH_CHILD_HAPLOTYPES guard above has
    # already returned whenever both child haplotypes carry alt above the error expectation, so alt_off is 0 or noise
    # by the time we get here and this is 1.0 in 96.2% of rows (99.6% of SNVs). What it actually measures is the
    # PURITY of the child's read partition. Filtering on it removes almost nothing -- in particular it keeps every one
    # of the 72 change-point calls whose paternal fraction is exactly 0.500. Kept under its existing name so downstream
    # tables do not silently change meaning; redefining it is a separate decision (see DESIGN P34).
    out["poo_confidence"] = round(alt_on / (alt_on + alt_off), 3)
    out["poo_reason"] = "OK"
    pl = m.parent_labelled(origin)
    if pl:
        out.update(t_alt_reads=pl["T"].alt, t_dp=pl["T"].dp, u_alt_reads=pl["U"].alt, u_dp=pl["U"].dp, t_hap_resolved=1)
    else:
        out["poo_reason"] = "OK_TRANSMITTED_HAP_UNRESOLVED"
    nt = m.row(other, 1).alt + m.row(other, 2).alt + m.untagged[other].alt
    out["nt_parent_alt_reads"] = nt
    return out


# ----------------------------------------------------------------------------------------------
# rule layer (P8), in order; first match wins
# ----------------------------------------------------------------------------------------------
def sv_depth_features(raw: Dict[str, object]) -> Dict[str, object]:
    """The registry's SV depth/junction features from the raw per-haplotype counts of readers.sv_interval_evidence.

    Haplotypes are ordered by the evidence itself and never labelled: `_A` is the child haplotype whose depth deviates
    most from the flanking rate (the putatively affected one for a deletion or duplication), `_O` is the other. For a
    heterozygous deletion the expectation is _A near 0 and _O near 1. `p_sv_max_hap_depth_change` is the largest
    deviation over the four parental haplotypes: a parent that also loses depth carries the event, so the candidate is
    inherited rather than de novo — the comparison a single-sample read review cannot make."""
    if not raw:
        return {}
    out: Dict[str, object] = {}
    c = [raw.get("C_sv_ratio_hap1"), raw.get("C_sv_ratio_hap2")]
    c = [x for x in c if x is not None]
    if c:
        dev = sorted(c, key=lambda r: -abs(1.0 - r))
        out["c_sv_hap_depth_change_A"] = round(dev[0], 4)
        if len(dev) > 1:
            out["c_sv_hap_depth_change_O"] = round(dev[1], 4)
    if raw.get("C_sv_ratio_all") is not None:
        out["sv_depth_ratio_inside_flank"] = raw["C_sv_ratio_all"]
    pr = [raw.get("%s_sv_ratio_hap%d" % (role, hp)) for role in ("F", "M") for hp in (1, 2)]
    pr = [x for x in pr if x is not None]
    if pr:
        out["p_sv_max_hap_depth_change"] = round(max(abs(1.0 - x) for x in pr), 4)
    j1, j2 = raw.get("C_sv_junc_hap1") or 0, raw.get("C_sv_junc_hap2") or 0
    if j1 + j2:
        out["c_sv_junction_hap_concentration"] = round(max(j1, j2) / (j1 + j2), 4)
    out["c_sv_junc_untagged_frac"] = round((raw.get("C_sv_junc_hapu") or 0) / max((raw.get("C_sv_junc_hapu") or 0) + j1 + j2, 1), 4)
    out["c_sv_junc_both_ends"] = raw.get("C_sv_junc_both_ends")
    return out


def _rule_score_sv(sv: Dict[str, object], cls: str, cp: ClassParams) -> int:
    """Germline criteria met on the deletion path (0-6), the same scale the alt-read rules use."""
    r = sv.get("sv_depth_ratio_inside_flank")
    het = sv.get("sv_het_snv_persistence")
    p = [sv.get("%s_sv_ratio_all" % role) for role in ("F", "M")]
    p = [x for x in p if x is not None]
    n = 0
    n += int(r is not None and r <= cp.del_max_ratio)                       # child depth halved
    n += int(bool(sv.get("c_sv_junc_both_ends")))                            # junction reads at both breakpoints
    n += int(het is not None and het <= cp.del_max_het_persistence)          # loss of heterozygosity inside
    n += int(bool(p) and all(x >= cp.del_parent_min_ratio for x in p))       # neither parent depleted
    n += int(len(p) == 2)                                                    # both parents measurable
    n += int((sv.get("C_sv_tagged_loss") or 0) > 0.1 or (r is not None and r <= 0.6))
    return n


def classify_sv_interval(sv: Dict[str, object], t: Dict[str, object], cp: ClassParams) -> Optional[Tuple[str, List[str]]]:
    """P8 amendment 2026-09-14 (JS): classify a DELETION from depth across the interval, not from alt-read confinement.

    Alt-read confinement is a point-variant instrument and it cannot work for a large deletion. A read crossing the
    junction is split between a primary and a supplementary alignment; inside the deleted segment there is no
    heterozygous site to phase against, so the read's HP tag is often absent or wrong. Measured on the cohort's 35 kb
    MECP2 deletion: 5 junction reads per breakpoint, 3 untagged and the other 2 assigned to OPPOSITE haplotypes, giving
    an alt fraction of 0.33 on the 'alt' haplotype and a verdict of `inconclusive` for a real pathogenic de novo event.

    The evidence that does work for a deletion, all measured here for all six haplotypes:
      * total depth inside the interval against the flanks — about 0.5 for a constitutional heterozygous deletion;
      * the same ratio in each parent — a parent that is also depleted carries the event, so it is inherited;
      * loss of heterozygosity inside the interval (het-SNV persistence towards 0), which separates a constitutional
        deletion from mosaicism or chimeric reads, where heterozygous sites persist at the flanking rate;
      * junction reads at BOTH breakpoints.
    Returns (class, flags) when the depth evidence decides, otherwise None and the alt-read rules apply."""
    r = sv.get("sv_depth_ratio_inside_flank")
    if r is None:
        return None
    flags: List[str] = ["SV_DEPTH_EVIDENCE"]
    het = sv.get("sv_het_snv_persistence")
    both_ends = sv.get("c_sv_junc_both_ends")
    p_ratios = [sv.get("%s_sv_ratio_all" % role) for role in ("F", "M")]
    p_ratios = [x for x in p_ratios if x is not None]
    parent_depleted = any(x <= cp.del_max_ratio for x in p_ratios)
    if parent_depleted:
        return "inherited_missed_in_parent", flags + ["PARENT_DEPTH_DEPLETED"]
    child_depleted = r <= cp.del_max_ratio
    if child_depleted and p_ratios and all(x >= cp.del_parent_min_ratio for x in p_ratios):
        if not both_ends:
            return "inconclusive", flags + ["NO_JUNCTION_BOTH_ENDS"]
        if het is not None and het > cp.del_max_het_persistence:
            # depth fell but heterozygous sites persist inside: both haplotypes are present, so not constitutional
            return "child_postzygotic_mosaic", flags + ["HET_SNVS_PERSIST_INSIDE"]
        cls = "germline_DNM_phased" if t.get("parent_of_origin") in ("paternal", "maternal") else "germline_DNM_unphased"
        return cls, flags + (["LOSS_OF_HETEROZYGOSITY"] if (het is not None and het <= cp.del_max_het_persistence) else [])
    if not child_depleted and het is not None and het > cp.del_max_het_persistence:
        return "phase_conflict_artifact", flags + ["NO_DEPTH_LOSS_HETS_PERSIST"]
    return None


def classify(m: Matrix, f: Dict[str, object], t: Dict[str, object], hp: HapParams, cp: ClassParams,
             labels: Optional[LabelTables] = None, chrom: str = "", pos: int = 0,
             sv: Optional[Dict[str, object]] = None) -> Dict[str, object]:
    flags: List[str] = []
    k = cp.working_k
    c1, c2 = m.row("C", 1), m.row("C", 2)
    A, O = (c1, c2) if c1.alt >= c2.alt else (c2, c1)
    prow = [m.row(par, hp_) for par in ("F", "M") for hp_ in (1, 2)]
    e = lambda r: _e(r.n, hp)
    parents_het_gt = False  # caller GT knowledge is applied upstream (candidates are child-only alt by construction)
    if f["hap_obs_k%d" % k] < 6:
        flags.append("LOW_HAP_DEPTH")
    if any(r.dp >= 3 and (r.amb_frac or 0.0) > hp.amb_flag_frac for r in (c1, c2, *prow)):
        flags.append("AMBIGUOUS_READS")
    if m.child_hap1_is is None:
        flags.append("CHILD_BLOCK_UNORIENTED")
    if labels is not None and chrom:
        for par in ("F", "M"):
            d = labels.dist_change_point(par, chrom, pos)
            if d is not None and d <= 50000:
                flags.append("NEAR_CHANGE_POINT_%s" % par)
        d = labels.dist_child_switch(chrom, pos)
        if d is not None and d <= 50000:
            flags.append("NEAR_CHILD_SWITCH")
    sv_call = classify_sv_interval(sv, t, cp) if sv else None
    if sv_call is not None:
        cls, sv_flags = sv_call
        flags.extend(sv_flags)
        return dict(phase_class=cls, rule_score=_rule_score_sv(sv, cls, cp), flags=";".join(flags))
    alt_both_child = A.alt > e(A) and O.alt > e(O)
    n_par_haps_alt = sum(1 for r in prow if r.alt > e(r))
    T = m.parent_labelled(t.get("parent_of_origin", "")[:1].upper().replace("P", "F").replace("M", "M")) if t.get("parent_of_origin") in ("paternal", "maternal") else {}
    # 1 conflict
    if alt_both_child or (n_par_haps_alt >= 2 and not parents_het_gt):
        cls = "phase_conflict_artifact"
    # 2 inherited-and-missed
    elif T and T["T"].n >= k and T["T"].alt >= max(3, cp.inherited_min_frac_on_T * T["T"].n):
        cls = "inherited_missed_in_parent"
    # 3 parental mosaic transmitted
    elif T and T["T"].alt > e(T["T"]) and T["T"].alt < cp.inherited_min_frac_on_T * T["T"].n:
        if T["T"].n >= cp.mosaic_min_dp and T["T"].alt >= cp.mosaic_min_minority_reads:
            cls = "parental_mosaic_transmitted"
        else:
            cls = "inconclusive"; flags.append("PARENTAL_ALT_LOW_DEPTH")
    # evidence floor for the remaining classes: enough alt reads to be a claim at all (P10; first smoke run made
    # single-read "germline" calls at ~2,000 per child)
    elif A.alt + m.untagged["C"].alt < cp.min_alt_reads:
        cls = "inconclusive"; flags.append("TOO_FEW_ALT_READS")
    # 4 child postzygotic mosaic
    elif (A.alt_frac is not None and A.alt_frac <= cp.mosaic_max_hap_frac and O.alt <= e(O)
          and all(r.n >= k and r.alt <= e(r) for r in prow)):
        if A.n >= cp.mosaic_min_dp and A.alt >= cp.min_alt_reads and A.ref >= cp.mosaic_min_minority_reads:
            cls = "child_postzygotic_mosaic"
        else:
            cls = "inconclusive"; flags.append("MOSAIC_UNDERPOWERED")
    # 5 germline phased: alt confined to one child haplotype, all FOUR parental haplotypes observed and alt-free,
    #   both child haplotypes observed, transmitted haplotype known and alt-free
    elif (A.alt >= cp.min_alt_reads and A.alt_frac is not None and A.alt_frac >= cp.germline_min_hap_frac and O.alt <= e(O)
          and A.n >= k and O.n >= k and all(r.n >= k and r.alt <= e(r) for r in prow)
          and T and T["T"].alt <= e(T["T"])):
        cls = "germline_DNM_phased"
    # 6 germline unphased: consistent, but some haplotype unobserved or transmission unresolved
    elif (O.alt <= e(O) and all(r.alt <= e(r) for r in prow if r.n > 0)
          and (A.alt_frac is None or A.alt_frac >= cp.germline_min_hap_frac or A.alt + O.alt == 0)):
        cls = "germline_DNM_unphased"
    else:
        cls = "inconclusive"
        if not flags:
            flags.append("NO_RULE_MATCHED")
    # transparent germline score: number of germline criteria met (0-6)
    score = sum([
        A.alt_frac is not None and A.alt_frac >= cp.germline_min_hap_frac,
        O.alt <= e(O),
        bool(T) and T["T"].n >= k and T["T"].alt <= e(T["T"]),
        bool(T) and T["U"].n >= k,
        f["hap_obs_k%d" % k] == 6,
        all(r.alt <= e(r) for r in prow),
    ])
    return dict(phase_class=cls, rule_score=int(score), flags=";".join(flags) if flags else "")


# ----------------------------------------------------------------------------------------------
# re-running the rule layer from an evidence row (no BAM access)
# ----------------------------------------------------------------------------------------------
READ_QUALITY_KEYS = ("c_alt_mapq_mean", "c_alt_nm_rate", "c_ref_nm_rate", "c_alt_clip_frac", "c_alt_mapq0_frac", "c_alt_supp_frac",
                     "c_alt_readlen_median", "c_alt_rq_mean", "c_tr_al_hapA_mean", "c_tr_al_hapA_sd", "c_tr_al_hapO_mean")
POSITIONAL_FLAGS = ("NEAR_CHANGE_POINT_F", "NEAR_CHANGE_POINT_M", "NEAR_CHILD_SWITCH")


def poo_clear_near_switch(t: Dict[str, object], flags: object) -> Dict[str, object]:
    """Blank the parent of origin where the phase geometry says it cannot be believed.

    Within 50 kb of a parental change point or a child haplotype switch, the child's haplotype labels may be the wrong
    way round, so the parent the alt reads point to may be the wrong parent. Measured on the cohort's SNV calls, rows
    flagged NEAR_CHANGE_POINT have a paternal fraction of 0.5000 EXACTLY (36 of 72) against 0.786 for unflagged rows,
    p = 1.1e-8 -- a coin toss, which is what an assignment with no information looks like. NEAR_CHILD_SWITCH gives
    0.623 vs 0.781 (p = 0.0036). Negative controls do not behave this way (LOW_HAP_DEPTH p = 0.59, AMBIGUOUS_READS
    p = 0.43), so this is specific to phase geometry rather than to flags in general.

    Until now those rows were emitted as `poo_reason = OK` with `poo_confidence` 1.0 and the warning visible only in
    free-text flags, so anything reading the parent_of_origin column took a coin toss at face value. They are now
    `undetermined` and say why. Cost: 120 of 3853 calls, 3.1%."""
    fl = {x for x in str(flags or "").split(";") if x}
    hit = sorted(fl & set(POSITIONAL_FLAGS))
    if not hit or t.get("parent_of_origin") in (None, "", "undetermined"):
        return t
    out = dict(t)
    out["parent_of_origin"] = "undetermined"
    out["poo_reason"] = "PHASE_SWITCH_RISK:" + ",".join(hit)
    out["poo_confidence"] = None
    return out


def reclassify_row(row: Dict[str, object], hp: HapParams, cp: ClassParams, thresholds_version: Optional[str] = None) -> Dict[str, object]:
    """Recompute the count-derived features, the transmission block and the P8 class from an evidence row.
    Read-quality summaries and the M1 positional flags (which need the label tables) are carried over."""
    m = Matrix.from_evidence_row(row)
    f = features(m, hp)
    t = transmission_features(m, hp)
    # The SV INTERVAL evidence must be carried forward. It is produced by a BAM pass that this function deliberately
    # avoids, so calling classify() without it silently drops the whole depth branch: the cohort's 35 kb pathogenic
    # deletion went from rule_score 6 and TIER1_SV_DEPTH to rule_score 2 and BELOW_TAU purely by being reclassified.
    # The columns are already in the row (review writes both the raw counts and the derived ratios), so rebuild the
    # dict from them rather than recomputing anything.
    # The row comes from csv.DictReader, so every value is a STRING, and classify_sv_interval / _rule_score_sv
    # compare them with floats ("0.59" <= 0.7 raises TypeError). This went unnoticed because the one cohort
    # reclassify that ran on SV tables (2026-09-16) ran while the sv_ columns were blank, so sv_row was empty and
    # the depth branch was skipped; the "carry the SV evidence forward" intent above had never actually executed.
    def _num(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return v
    sv_row = {k: _num(v) for k, v in row.items()
              if (k.startswith("sv_") or (len(k) > 2 and k[0] in "CFM" and k[1] == "_" and "sv_" in k))
              and v not in (None, "")}
    # classify_sv_interval reads the DERIVED block (c_sv_junc_both_ends, p_sv_max_hap_depth_change, ...), which
    # review.py builds with sv_depth_features(raw) | raw. The filter above keeps only the raw upper-case columns, so
    # without rebuilding the derived block here `both_ends` is None and every reclassified deletion falls to
    # inconclusive / NO_JUNCTION_BOTH_ENDS -- the second cause of the "rule 6 -> 2 after reclassify" episode.
    sv_full = (sv_depth_features(sv_row) | sv_row) if sv_row else None
    c = classify(m, f, t, hp, cp, sv=sv_full)
    old_flags = [x for x in str(row.get("flags") or "").split(";") if x in POSITIONAL_FLAGS]
    c["flags"] = ";".join([x for x in c["flags"].split(";") if x] + old_flags)
    t = poo_clear_near_switch(t, c["flags"])
    out = dict(row)
    out.update({k: v for k, v in f.items() if k not in READ_QUALITY_KEYS})
    out.update(t)
    out.update(c)
    if thresholds_version is not None:
        out["thresholds_version"] = thresholds_version
    return out
