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


@dataclass
class HapParams:
    k: Tuple[int, ...] = (3, 5)
    error_reads: int = 1
    error_frac: float = 0.05
    min_mapq: int = 20


@dataclass
class ClassParams:
    germline_min_hap_frac: float = 0.8
    mosaic_max_hap_frac: float = 0.7
    mosaic_min_dp: int = 15
    mosaic_min_minority_reads: int = 3
    inherited_min_frac_on_T: float = 0.3
    working_k: int = 5


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

    def child_hap1_is(self, chrom: str, ps: Optional[int], pos: int) -> Optional[str]:
        if ps is None:
            return None
        for s, e, o in self.orient.get((chrom, ps), ()):
            if s <= pos <= e:
                return "P" if o == "HAP1_PAT" else ("M" if o == "HAP1_MAT" else None)
        return None

    def parent_transmitted(self, parent: str, chrom: str, ps: Optional[int], pos: int) -> Optional[int]:
        if ps is None:
            return None
        for s, e, tr in self.trans.get((parent, chrom, ps), ()):
            if s <= pos <= e:
                return 1 if tr == "HAP1" else (2 if tr == "HAP2" else None)
        return None

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

    def add(self, r: ReadObs):
        self.dp += 1
        self.mapq.append(r.mapq)
        if r.support == "ALT":
            self.alt += 1
            if r.nm_rate is not None: self.nm_alt.append(r.nm_rate)
            self.clip_alt += r.clipped
            self.supp_alt += r.supplementary
        elif r.support == "REF":
            self.ref += 1
            if r.nm_rate is not None: self.nm_ref.append(r.nm_rate)
        else:
            self.amb += 1
        if r.al is not None:
            self.al.append(r.al)

    @property
    def alt_frac(self) -> Optional[float]:
        n = self.alt + self.ref
        return self.alt / n if n else None

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
             c_alt_confined=int(O.alt <= _e(O.dp, p)),
             c_alt_tagged_frac=round((c1.alt + c2.alt) / (c1.alt + c2.alt + m.untagged["C"].alt), 3)
             if (c1.alt + c2.alt + m.untagged["C"].alt) else None,
             c_untagged_dp=m.untagged["C"].dp, c_untagged_alt=m.untagged["C"].alt)
    for k in p.k:
        f["c_both_haps_obs_k%d" % k] = int(c1.dp >= k and c2.dp >= k)
    prow = [m.row(par, hp) for par in ("F", "M") for hp in (1, 2)]
    f.update(p_min_hap_dp=min(r.dp for r in prow), p_max_alt_any_hap=max(r.alt for r in prow),
             p_sum_alt_all_haps=sum(r.alt for r in prow),
             p_n_haps_with_alt=sum(1 for r in prow if r.alt > _e(r.dp, p)),
             p_max_alt_hap_frac=round(max((r.alt_frac or 0.0) for r in prow), 4),
             p_untagged_dp_max=max(m.untagged["F"].dp, m.untagged["M"].dp),
             p_untagged_alt_max=max(m.untagged["F"].alt, m.untagged["M"].alt))
    for k in p.k:
        f["p_n_haps_obs_k%d" % k] = sum(1 for r in prow if r.dp >= k)
        f["hap_obs_k%d" % k] = f["p_n_haps_obs_k%d" % k] + int(c1.dp >= k) + int(c2.dp >= k)
    # read-quality summaries for the child rows (alt vs ref)
    f.update(c_alt_mapq_mean=A.summary()["mapq_mean"], c_alt_nm_rate=A.summary()["nm_alt_mean"], c_ref_nm_rate=A.summary()["nm_ref_mean"],
             c_alt_clip_frac=A.summary()["clip_alt_frac"])
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
    if cp.alt > _e(cp.dp, p) and cm.alt > _e(cm.dp, p):
        out["poo_reason"] = "ALT_ON_BOTH_CHILD_HAPLOTYPES"
        return out
    origin, other = ("F", "M") if cp.alt >= cm.alt else ("M", "F")
    alt_on, alt_off = (cp.alt, cm.alt) if origin == "F" else (cm.alt, cp.alt)
    out["parent_of_origin"] = "paternal" if origin == "F" else "maternal"
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
def classify(m: Matrix, f: Dict[str, object], t: Dict[str, object], hp: HapParams, cp: ClassParams,
             labels: Optional[LabelTables] = None, chrom: str = "", pos: int = 0) -> Dict[str, object]:
    flags: List[str] = []
    k = cp.working_k
    c1, c2 = m.row("C", 1), m.row("C", 2)
    A, O = (c1, c2) if c1.alt >= c2.alt else (c2, c1)
    prow = [m.row(par, hp_) for par in ("F", "M") for hp_ in (1, 2)]
    e = lambda r: _e(r.dp, hp)
    parents_het_gt = False  # caller GT knowledge is applied upstream (candidates are child-only alt by construction)
    if f["hap_obs_k%d" % k] < 6:
        flags.append("LOW_HAP_DEPTH")
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
    alt_both_child = A.alt > e(A) and O.alt > e(O)
    n_par_haps_alt = sum(1 for r in prow if r.alt > e(r))
    T = m.parent_labelled(t.get("parent_of_origin", "")[:1].upper().replace("P", "F").replace("M", "M")) if t.get("parent_of_origin") in ("paternal", "maternal") else {}
    # 1 conflict
    if alt_both_child or (n_par_haps_alt >= 2 and not parents_het_gt):
        cls = "phase_conflict_artifact"
    # 2 inherited-and-missed
    elif T and T["T"].dp >= k and T["T"].alt >= max(3, cp.inherited_min_frac_on_T * T["T"].dp):
        cls = "inherited_missed_in_parent"
    # 3 parental mosaic transmitted
    elif T and T["T"].alt > e(T["T"]) and T["T"].alt < cp.inherited_min_frac_on_T * T["T"].dp:
        if T["T"].dp >= cp.mosaic_min_dp and T["T"].alt >= cp.mosaic_min_minority_reads:
            cls = "parental_mosaic_transmitted"
        else:
            cls = "inconclusive"; flags.append("PARENTAL_ALT_LOW_DEPTH")
    # 4 child postzygotic mosaic
    elif (A.alt_frac is not None and A.alt_frac <= cp.mosaic_max_hap_frac and O.alt <= e(O)
          and all(r.dp >= k and r.alt <= e(r) for r in prow)):
        if A.dp >= cp.mosaic_min_dp and (A.ref) >= cp.mosaic_min_minority_reads:
            cls = "child_postzygotic_mosaic"
        else:
            cls = "inconclusive"; flags.append("MOSAIC_UNDERPOWERED")
    # 5 germline phased
    elif (A.alt_frac is not None and A.alt_frac >= cp.germline_min_hap_frac and O.alt <= e(O)
          and T and T["T"].dp >= k and T["T"].alt <= e(T["T"]) and T["U"].dp >= k and A.dp >= k and O.dp >= k):
        cls = "germline_DNM_phased"
    # 6 germline unphased
    elif (O.alt <= e(O) and all(r.alt <= e(r) for r in prow if r.dp > 0) and (A.alt > 0 or m.untagged["C"].alt > 0)):
        cls = "germline_DNM_unphased"
    else:
        cls = "inconclusive"
        if not flags:
            flags.append("NO_RULE_MATCHED")
    # transparent germline score: number of germline criteria met (0-6)
    score = sum([
        A.alt_frac is not None and A.alt_frac >= cp.germline_min_hap_frac,
        O.alt <= e(O),
        bool(T) and T["T"].dp >= k and T["T"].alt <= e(T["T"]),
        bool(T) and T["U"].dp >= k,
        f["hap_obs_k%d" % k] == 6,
        all(r.alt <= e(r) for r in prow),
    ])
    return dict(phase_class=cls, rule_score=int(score), flags=";".join(flags) if flags else "")
