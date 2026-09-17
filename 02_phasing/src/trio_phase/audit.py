"""Assert that what this module DECLARES is what its output tables actually CONTAIN.

Same shape and the same reason as `phase_dnm.audit` in 13_phase_dnm: every defect that cost time in
this pipeline was something declared and never produced, with nothing comparing the declaration
against the artefact. A green result here is a measurement, not an absence of complaints, so every
check reports its evidence either way, and the process exits non-zero on any FAIL.

The declarations are not restated by hand. The column lists come from the modules that WRITE the
tables (`orient.ORIENTATION_COLUMNS`, `transmission.SEGMENT_COLUMNS` / `CHANGE_COLUMNS`,
`xo_reads.COLUMNS`), so a writer that gains a column and a reader that never hears about it are
caught here rather than three modules downstream. Only the closed value vocabularies are listed
below, because they exist in the writers as literals; each is annotated with what the cohort of
2026-09-16 actually contained.

Checks
  trios        every complete trio in the manifest has all seven per-child artefacts
  nonempty     every table has at least one data row, and resolved == unresolved change points
  schema       the header on disk is exactly the writer's declared column list
  vocabulary   every categorical column's values lie inside the declared vocabulary
  overlap      orientation and transmission segments do not overlap WITHIN one phase block
  hapdepth     one depth table and summary per manifest sample (only when --hapdepth-dir is given)
"""
from __future__ import annotations

import collections
import csv
import glob
import gzip
import json
import os
from typing import Dict, Iterable, List, Optional, Set, Tuple

from .io.vcf import read_manifest, trio_of
from .phasing import orient as O
from .phasing import transmission as T
from .phasing import xo_reads as X

TAB = "\t"

# ---------------------------------------------------------------------------------------------
# declarations
# ---------------------------------------------------------------------------------------------

# The seven per-child artefacts a complete trio must have. `.tsv` entries are also schema- and
# vocabulary-checked; the summaries are only checked for existence and for being parseable JSON.
CHILD_TABLES = (
    "orientation.tsv",
    "orientation.dissent.tsv",
    "transmission.tsv",
    "changepoints.tsv",
    "changepoints.resolved.tsv",
)
CHILD_SUMMARIES = (
    "orientation.summary.json",
    "transmission.summary.json",
)

# `orientation.dissent.tsv` is the one table that may legitimately be empty: it lists the sites that
# voted against their block's orientation, so nothing in it is the ideal result. It is reported and
# warned about rather than failed (cohort 2026-09-16: 214-693 rows per child, so an empty one on a
# whole genome still deserves a look).
MAY_BE_EMPTY = ("orientation.dissent.tsv",)

DISSENT_COLUMNS = ["chrom", "pos", "phase_block_id", "block_orientation"]

SCHEMA: Dict[str, List[str]] = {
    "orientation.tsv": O.ORIENTATION_COLUMNS,
    "orientation.dissent.tsv": DISSENT_COLUMNS,
    "transmission.tsv": T.SEGMENT_COLUMNS,
    "changepoints.tsv": T.CHANGE_COLUMNS,
    "changepoints.resolved.tsv": X.COLUMNS,
}

# Closed vocabularies. The counts in the comments are the whole cohort on 2026-09-16 (35 complete
# trios), which is what a value outside a vocabulary should be weighed against.
VOCAB: Dict[str, Dict[str, Set[str]]] = {
    "orientation.tsv": {
        # HAP1_PAT 114,765 / HAP1_MAT 113,748 / AMBIGUOUS 263,469
        "orientation": {O.HAP1_PAT, O.HAP1_MAT, O.AMBIGUOUS},
        # OK 183,543 / LOW_SITES 215,595 / NO_INFORMATIVE_SITES 47,627 /
        # OK_SMALL_UNANIMOUS 25,966 / SPLIT_AT_SWITCH 19,004 / MIXED_VOTES 247
        "reason": {"OK", "OK_SMALL_UNANIMOUS", "SPLIT_AT_SWITCH", "LOW_SITES", "MIXED_VOTES",
                   "NO_INFORMATIVE_SITES"},
    },
    "orientation.dissent.tsv": {
        "block_orientation": {O.HAP1_PAT, O.HAP1_MAT, O.AMBIGUOUS},
    },
    "transmission.tsv": {
        "parent": {"F", "M"},
        # HAP1 225,397 / HAP2 207,304 / UNRESOLVED 462,019
        "transmitted": {T.HAP1, T.HAP2, T.UNRESOLVED},
        # OK 432,701 / LOW_SITES 353,583 / NO_INFORMATIVE_SITES 100,488 / MIXED_VOTES 7,948
        "reason": {"OK", "LOW_SITES", "MIXED_VOTES", "NO_INFORMATIVE_SITES"},
        # 18,315 CHANGE_POINT on each side — one per change point, by construction
        "left_boundary": {"BLOCK_EDGE", "CHANGE_POINT"},
        "right_boundary": {"BLOCK_EDGE", "CHANGE_POINT"},
    },
    "changepoints.tsv": {
        "parent": {"F", "M"},
        "left_hap": {T.HAP1, T.HAP2},
        "right_hap": {T.HAP1, T.HAP2},
        # every row, always: the VCF level cannot tell a crossover from a parental switch error
        "status": {"CANDIDATE"},
    },
    "changepoints.resolved.tsv": {
        "parent": {"F", "M"},
        "left_hap": {T.HAP1, T.HAP2},
        "right_hap": {T.HAP1, T.HAP2},
        # AMBIGUOUS 14,711 / CROSSOVER 2,987 / SWITCH_ERROR 617 / UNTESTABLE 0.
        # UNTESTABLE is declared by the classifier (xo_reads.classify's counter) and has never been
        # emitted on this cohort; it stays in the vocabulary because the writer can still produce it.
        "status": {"CROSSOVER", "SWITCH_ERROR", "AMBIGUOUS", "UNTESTABLE"},
        "child_switch_in_interval": {"Y", "N"},
    },
}

# A resolved change-point table with no CROSSOVER at all means the read-level step ran and decided
# nothing: ~25-45 crossovers per meiosis is the biological expectation (cohort: F 27, M 58 median
# after excluding intervals containing a located child switch).
MIN_CROSSOVERS_PER_CHILD = 1


class Finding:
    __slots__ = ("level", "check", "detail")

    def __init__(self, level: str, check: str, detail: str):
        self.level, self.check, self.detail = level, check, detail

    def __str__(self):
        return "%-5s %-22s %s" % (self.level, self.check, self.detail)


def _rows(path: str) -> Tuple[List[str], List[dict]]:
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter=TAB)
        cols = list(rd.fieldnames or [])
        return cols, list(rd)


def complete_trios(manifest: str, family: Optional[str] = None) -> List[Tuple[str, str]]:
    """(family_id, child_id) for every offspring whose BOTH parents are samples of the manifest.

    This is the same rule the workflow uses to decide what to run, so the audit's expectation and the
    pipeline's work list cannot drift apart. Duos are excluded on purpose: nothing is produced for them."""
    man = read_manifest(manifest)
    out = []
    for sid, r in sorted(man.items()):
        if r.get("role") != "offspring":
            continue
        if family and r["family_id"] != family:
            continue
        try:
            trio_of(man, sid)
        except (ValueError, KeyError):
            continue
        out.append((r["family_id"], sid))
    return out


def audit_trios(phase_dir: str, manifest: str, family: Optional[str] = None) -> Tuple[List[Finding], List[Tuple[str, str]]]:
    """Every complete trio the manifest declares must have every artefact this module promises."""
    out: List[Finding] = []
    trios = complete_trios(manifest, family)
    if not trios:
        return [Finding("FAIL", "trios", "the manifest declares no complete trio (%s)" % manifest)], []
    present, missing = [], collections.defaultdict(list)
    for fam, child in trios:
        stem = os.path.join(phase_dir, fam, child)
        gone = [suf for suf in CHILD_TABLES + CHILD_SUMMARIES if not os.path.exists(stem + "." + suf)]
        if gone:
            for suf in gone:
                missing[suf].append(child)
        else:
            present.append((fam, child))
    out.append(Finding("INFO", "trios", "manifest declares %d complete trios; %d have all %d artefacts"
                       % (len(trios), len(present), len(CHILD_TABLES) + len(CHILD_SUMMARIES))))
    for suf, kids in sorted(missing.items()):
        out.append(Finding("FAIL", "trios", "%d complete trio(s) have no .%s: %s"
                           % (len(kids), suf, ", ".join(sorted(kids)[:8]) + (" ..." if len(kids) > 8 else ""))))
    if not missing:
        out.append(Finding("OK", "trios", "every complete trio has every declared artefact"))

    # the converse: tables for a child the manifest does not call a complete trio
    declared = {c for _, c in trios}
    stray = sorted({os.path.basename(p).split(".orientation.tsv")[0]
                    for p in glob.glob(os.path.join(phase_dir, "*", "*.orientation.tsv"))} - declared)
    if stray:
        out.append(Finding("WARN", "trios", "%d orientation table(s) for children the manifest does not make a "
                                            "complete trio: %s" % (len(stray), ", ".join(stray[:8]))))
    return out, present


def audit_nonempty(phase_dir: str, trios: Iterable[Tuple[str, str]]) -> List[Finding]:
    """A table that exists and is empty is the failure mode a `-s` test in a shell script does catch and a
    downstream join does not: the join simply produces nothing and reports success."""
    out: List[Finding] = []
    counts: Dict[str, List[int]] = collections.defaultdict(list)
    empty: Dict[str, List[str]] = collections.defaultdict(list)
    mismatch: List[str] = []
    no_crossover: List[str] = []
    for fam, child in trios:
        stem = os.path.join(phase_dir, fam, child)
        n_by_suf = {}
        for suf in CHILD_TABLES:
            p = stem + "." + suf
            if not os.path.exists(p):
                continue
            _, rows = _rows(p)
            n_by_suf[suf] = len(rows)
            counts[suf].append(len(rows))
            if not rows:
                empty[suf].append(child)
            if suf == "changepoints.resolved.tsv" and rows:
                if sum(1 for r in rows if r.get("status") == "CROSSOVER") < MIN_CROSSOVERS_PER_CHILD:
                    no_crossover.append(child)
        if ("changepoints.tsv" in n_by_suf and "changepoints.resolved.tsv" in n_by_suf
                and n_by_suf["changepoints.tsv"] != n_by_suf["changepoints.resolved.tsv"]):
            mismatch.append("%s (%d -> %d)" % (child, n_by_suf["changepoints.tsv"],
                                               n_by_suf["changepoints.resolved.tsv"]))
    for suf in CHILD_TABLES:
        c = counts.get(suf) or []
        if c:
            out.append(Finding("INFO", "nonempty", "%-26s %d children, rows min %d max %d"
                               % (suf, len(c), min(c), max(c))))
    hard_empty = {s: k for s, k in empty.items() if s not in MAY_BE_EMPTY}
    for suf, kids in sorted(empty.items()):
        out.append(Finding("FAIL" if suf not in MAY_BE_EMPTY else "WARN", "nonempty",
                           "%s has a header and no rows for %d child(ren): %s"
                           % (suf, len(kids), ", ".join(sorted(kids)[:8]))))
    if mismatch:
        out.append(Finding("FAIL", "nonempty", "the read-level step did not classify every change point: %s"
                           % ", ".join(mismatch[:8])))
    if no_crossover:
        out.append(Finding("FAIL", "nonempty", "%d child(ren) have no CROSSOVER at all (a meiosis has ~25-45): %s"
                           % (len(no_crossover), ", ".join(no_crossover[:8]))))
    if not (hard_empty or mismatch or no_crossover):
        out.append(Finding("OK", "nonempty", "every table that must have rows has them, every change point was "
                                             "classified, every child has crossovers"))
    return out


def audit_schema(phase_dir: str, trios: Iterable[Tuple[str, str]]) -> List[Finding]:
    """The header on disk against the column list the writing module declares."""
    out: List[Finding] = []
    bad: Dict[str, List[str]] = collections.defaultdict(list)
    seen: collections.Counter = collections.Counter()
    for fam, child in trios:
        for suf, declared in SCHEMA.items():
            p = os.path.join(phase_dir, fam, child + "." + suf)
            if not os.path.exists(p):
                continue
            cols, _ = _rows(p)
            seen[suf] += 1
            if cols != list(declared):
                extra = [c for c in cols if c not in declared]
                gone = [c for c in declared if c not in cols]
                bad[suf].append("%s (missing %s; unexpected %s; order %s)"
                                % (child, gone or "-", extra or "-",
                                   "differs" if not (gone or extra) else "n/a"))
    for suf in SCHEMA:
        if seen[suf]:
            out.append(Finding("INFO", "schema", "%-26s %d files, %d declared columns"
                               % (suf, seen[suf], len(SCHEMA[suf]))))
    for suf, msgs in sorted(bad.items()):
        out.append(Finding("FAIL", "schema", "%s header differs from the writer's declaration: %s"
                           % (suf, "; ".join(msgs[:4]))))
    if not bad:
        out.append(Finding("OK", "schema", "every header matches its writer's declared columns exactly"))
    return out


def audit_vocabulary(phase_dir: str, trios: Iterable[Tuple[str, str]]) -> List[Finding]:
    """Every categorical column's values against the declared vocabulary, with what was actually seen.

    Reported both ways: a value outside the vocabulary is a FAIL (a writer changed and nothing downstream
    knows), a declared value that never appears is an INFO (it may simply not occur in this cohort)."""
    out: List[Finding] = []
    observed: Dict[Tuple[str, str], collections.Counter] = collections.defaultdict(collections.Counter)
    for fam, child in trios:
        for suf, cols in VOCAB.items():
            p = os.path.join(phase_dir, fam, child + "." + suf)
            if not os.path.exists(p):
                continue
            _, rows = _rows(p)
            for r in rows:
                for col in cols:
                    observed[(suf, col)][r.get(col)] += 1
    fails = 0
    for (suf, col), ctr in sorted(observed.items()):
        allowed = VOCAB[suf][col]
        outside = {v: n for v, n in ctr.items() if v not in allowed}
        unseen = sorted(allowed - set(ctr))
        detail = "%s:%s %s" % (suf.split(".tsv")[0], col,
                               ", ".join("%s=%d" % (v, n) for v, n in sorted(ctr.items(), key=lambda kv: -kv[1])
                                         if v in allowed))
        if unseen:
            detail += " [declared, never seen: %s]" % ", ".join(unseen)
        out.append(Finding("INFO", "vocabulary", detail))
        if outside:
            fails += 1
            out.append(Finding("FAIL", "vocabulary", "%s:%s has values OUTSIDE the declared vocabulary: %s"
                               % (suf, col, ", ".join("%r=%d" % (v, n) for v, n in sorted(outside.items(),
                                                                                          key=lambda kv: -kv[1])[:6]))))
    if not observed:
        out.append(Finding("FAIL", "vocabulary", "no table was read; nothing was checked"))
    elif not fails:
        out.append(Finding("OK", "vocabulary", "every categorical column stays inside its declared vocabulary"))
    return out


def _overlaps(segs: List[Tuple[int, int]]) -> int:
    n = 0
    segs = sorted(segs)
    for a, b in zip(segs, segs[1:]):
        if b[0] <= a[1]:
            n += 1
    return n


def audit_overlap(phase_dir: str, trios: Iterable[Tuple[str, str]], strict_across_block: bool = False) -> List[Finding]:
    """Segments of ONE phase block must tile it, never overlap.

    A block is split only at a located switch, so its segments are disjoint by construction; an overlap
    means the splitter emitted the same interval twice and every position in it has two orientations.

    Overlaps between DIFFERENT blocks are a separate matter and are NOT an error: HiPhase emits small
    blocks nested inside larger ones (cohort 2026-09-16: 193 such pairs across 35 children out of
    491,982 segments, all but 2 involving an AMBIGUOUS segment). They are counted and reported, and
    --strict-across-block promotes them to a FAIL for a caller that wants an exact cover."""
    out: List[Finding] = []
    within_o = within_t = across_o = across_o_oriented = 0
    bad_o: List[str] = []
    bad_t: List[str] = []
    n_seg_o = n_seg_t = 0
    for fam, child in trios:
        p = os.path.join(phase_dir, fam, child + ".orientation.tsv")
        if os.path.exists(p):
            _, rows = _rows(p)
            n_seg_o += len(rows)
            by_block: Dict[Tuple[str, str], List[Tuple[int, int]]] = collections.defaultdict(list)
            by_chrom: Dict[str, List[Tuple[int, int, str, str]]] = collections.defaultdict(list)
            for r in rows:
                s, e = int(r["start"]), int(r["end"])
                by_block[(r["chrom"], r["phase_block_id"])].append((s, e))
                by_chrom[r["chrom"]].append((s, e, r["phase_block_id"], r["orientation"]))
            n = sum(_overlaps(v) for v in by_block.values())
            within_o += n
            if n:
                bad_o.append("%s (%d)" % (child, n))
            for segs in by_chrom.values():
                segs.sort()
                for a, b in zip(segs, segs[1:]):
                    if b[0] <= a[1] and a[2] != b[2]:
                        across_o += 1
                        if O.AMBIGUOUS not in (a[3], b[3]):
                            across_o_oriented += 1
        p = os.path.join(phase_dir, fam, child + ".transmission.tsv")
        if os.path.exists(p):
            _, rows = _rows(p)
            n_seg_t += len(rows)
            by_pblock: Dict[Tuple[str, str, str], List[Tuple[int, int]]] = collections.defaultdict(list)
            for r in rows:
                by_pblock[(r["parent"], r["chrom"], r["parent_phase_block_id"])].append((int(r["start"]), int(r["end"])))
            n = sum(_overlaps(v) for v in by_pblock.values())
            within_t += n
            if n:
                bad_t.append("%s (%d)" % (child, n))
    out.append(Finding("INFO", "overlap", "orientation: %d segments, %d within-block overlaps, %d across-block "
                                          "(%d between two ORIENTED segments)"
                       % (n_seg_o, within_o, across_o, across_o_oriented)))
    out.append(Finding("INFO", "overlap", "transmission: %d segments, %d within-(parent, block) overlaps"
                       % (n_seg_t, within_t)))
    if within_o:
        out.append(Finding("FAIL", "overlap", "orientation segments overlap inside one phase block: %s"
                           % ", ".join(bad_o[:8])))
    if within_t:
        out.append(Finding("FAIL", "overlap", "transmission segments overlap inside one (parent, phase block): %s"
                           % ", ".join(bad_t[:8])))
    if across_o:
        out.append(Finding("FAIL" if strict_across_block else "WARN", "overlap",
                           "%d orientation overlaps between different phase blocks (nested HiPhase blocks; "
                           "%d of them between two oriented segments)" % (across_o, across_o_oriented)))
    if not (within_o or within_t):
        out.append(Finding("OK", "overlap", "no segment overlaps within a phase block, in either table"))
    return out


def audit_summaries(phase_dir: str, trios: Iterable[Tuple[str, str]]) -> List[Finding]:
    """The summaries must parse and must carry the thresholds version that produced them, so a table can
    always be traced back to its rule set."""
    out: List[Finding] = []
    versions: collections.Counter = collections.Counter()
    broken: List[str] = []
    for fam, child in trios:
        for suf in CHILD_SUMMARIES:
            p = os.path.join(phase_dir, fam, child + "." + suf)
            if not os.path.exists(p):
                continue
            try:
                d = json.load(open(p))
            except Exception as e:                                  # noqa: BLE001 - any parse failure is the finding
                broken.append("%s %s (%s)" % (child, suf, e))
                continue
            versions[d.get("thresholds_version")] += 1
    if broken:
        out.append(Finding("FAIL", "summaries", "unparseable: %s" % "; ".join(broken[:6])))
    if None in versions:
        out.append(Finding("FAIL", "summaries", "%d summary file(s) carry no thresholds_version" % versions[None]))
    known = {k: v for k, v in versions.items() if k is not None}
    if known:
        out.append(Finding("INFO", "summaries", "thresholds versions in the summaries: %s"
                           % ", ".join("%s x%d" % kv for kv in sorted(known.items()))))
    if len(known) > 1:
        out.append(Finding("WARN", "summaries", "the cohort was produced by MORE THAN ONE threshold version; "
                                                "the tables are not comparable across children"))
    if not broken and None not in versions and len(known) == 1:
        out.append(Finding("OK", "summaries", "every summary parses and names one threshold version"))
    return out


def audit_hapdepth(hapdepth_dir: str, manifest: str) -> List[Finding]:
    """One depth table and one summary per sample of the manifest (samples, not trios: parents count)."""
    out: List[Finding] = []
    man_all = read_manifest(manifest)
    # the workflow computes depth for samples of COMPLETE-TRIO families only (Snakefile SAMPLES); demanding a table for
    # every manifest sample failed the audit on any duo or singleton, which the README says have no outputs (P4)
    trio_fams = {r["family_id"] for r in man_all.values()
                 if r.get("role") == "offspring" and r.get("father_id") in man_all and r.get("mother_id") in man_all}
    man = {s: r for s, r in man_all.items() if r["family_id"] in trio_fams}
    missing_t = [s for s in sorted(man) if not os.path.exists(os.path.join(hapdepth_dir, s + ".hapdepth.tsv.gz"))]
    missing_s = [s for s in sorted(man) if not os.path.exists(os.path.join(hapdepth_dir, s + ".hapdepth.summary.json"))]
    out.append(Finding("INFO", "hapdepth", "%d manifest samples (%d in complete-trio families), %d depth tables, %d summaries"
                       % (len(man_all), len(man), len(man) - len(missing_t), len(man) - len(missing_s))))
    if missing_t:
        out.append(Finding("FAIL", "hapdepth", "%d sample(s) have no depth table: %s"
                           % (len(missing_t), ", ".join(missing_t[:8]))))
    if missing_s:
        out.append(Finding("FAIL", "hapdepth", "%d sample(s) have no depth summary: %s"
                           % (len(missing_s), ", ".join(missing_s[:8]))))
    if not (missing_t or missing_s):
        # a gzip file that exists and cannot be read is the other half of "the file is there"
        bad = []
        for s in sorted(man):
            p = os.path.join(hapdepth_dir, s + ".hapdepth.tsv.gz")
            try:
                with gzip.open(p, "rt") as fh:
                    if not fh.readline():
                        bad.append(s)
            except Exception as e:                                  # noqa: BLE001
                bad.append("%s (%s)" % (s, e))
        if bad:
            out.append(Finding("FAIL", "hapdepth", "unreadable or empty: %s" % ", ".join(str(b) for b in bad[:8])))
        else:
            out.append(Finding("OK", "hapdepth", "every manifest sample has a readable depth table and a summary"))
    return out


def run(phase_dir: str, manifest: str, hapdepth_dir: Optional[str] = None, family: Optional[str] = None,
        strict_across_block: bool = False, log=print) -> int:
    findings, trios = audit_trios(phase_dir, manifest, family)
    if trios:
        findings += audit_nonempty(phase_dir, trios)
        findings += audit_schema(phase_dir, trios)
        findings += audit_vocabulary(phase_dir, trios)
        findings += audit_overlap(phase_dir, trios, strict_across_block)
        findings += audit_summaries(phase_dir, trios)
    if hapdepth_dir:
        findings += audit_hapdepth(hapdepth_dir, manifest)
    for f in findings:
        log(str(f))
    n_fail = sum(1 for f in findings if f.level == "FAIL")
    n_warn = sum(1 for f in findings if f.level == "WARN")
    log("audit: %d FAIL, %d WARN, %d checks" % (n_fail, n_warn, len(findings)))
    return 1 if n_fail else 0
