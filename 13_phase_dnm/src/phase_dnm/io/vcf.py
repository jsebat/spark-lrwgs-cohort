"""Minimal VCF reading for phased trio genotypes. Pure Python, no htslib.

Why a text parser: Module 1 needs only GT / PS / GQ per site and must run identically on a login
node, inside the container, and on a Windows laptop for the unit tests. The per-family WDL
outputs are per-sample splits of ONE joint VCF (HiPhase is run per sample), so the three files
share contig order and site set minus each sample's uncalled sites. Sites are merged on
(contig order, position) and then matched on (REF, ALT) within the position, so a record that is
missing from one sample's file does not desynchronise the merge.
"""
from __future__ import annotations

import csv
import gzip
from dataclasses import dataclass
from typing import Dict, Iterator, List, Optional, Tuple


@dataclass(slots=True)
class Site:
    chrom: str
    pos: int
    ref: str
    alt: str                       # ALT column as written (comma-joined if multi-allelic)
    alleles: Optional[Tuple[int, ...]]   # allele indices, e.g. (0, 1); None if any allele missing
    phased: bool
    ps: Optional[int]
    gq: Optional[int]
    dp: Optional[int]
    ad: Optional[Tuple[int, ...]]

    @property
    def is_het(self) -> bool:
        return self.alleles is not None and len(self.alleles) == 2 and self.alleles[0] != self.alleles[1]

    @property
    def carried(self) -> Optional[set]:
        return set(self.alleles) if self.alleles else None


@dataclass(slots=True)
class TrioSite:
    chrom: str
    pos: int
    ref: str
    alt: str
    child: Site
    father: Optional[Site]
    mother: Optional[Site]


def open_text(path: str):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def _parse_gt(gt: str) -> Tuple[Optional[Tuple[int, ...]], bool]:
    if gt in (".", "./.", ".|.", ""):
        return None, False
    phased = "|" in gt
    parts = gt.replace("|", "/").split("/")
    if any(p == "." for p in parts):
        return None, False
    try:
        return tuple(int(p) for p in parts), phased
    except ValueError:
        return None, False


def _int_or_none(x: str) -> Optional[int]:
    if x in (".", ""):
        return None
    try:
        return int(x)
    except ValueError:
        try:
            return int(float(x))
        except ValueError:
            return None


class VcfReader:
    """Streams one sample's sites from a single- or multi-sample VCF (.vcf or .vcf.gz)."""

    def __init__(self, path: str):
        self.path = str(path)
        self.contigs: List[str] = []
        self.samples: List[str] = []
        self._fh = open_text(self.path)
        self._first_record: Optional[str] = None
        for line in self._fh:
            if line.startswith("##contig="):
                # ##contig=<ID=chr1,length=...>
                inner = line[len("##contig=<"):].rstrip(">\n")
                for kv in inner.split(","):
                    if kv.startswith("ID="):
                        self.contigs.append(kv[3:])
            elif line.startswith("#CHROM"):
                self.samples = line.rstrip("\n").split("\t")[9:]
                break
            elif not line.startswith("#"):
                raise ValueError("VCF without #CHROM header line: %s" % self.path)
        self.contig_index: Dict[str, int] = {c: i for i, c in enumerate(self.contigs)}

    def close(self):
        self._fh.close()

    def _contig_rank(self, chrom: str) -> int:
        if chrom not in self.contig_index:
            self.contig_index[chrom] = len(self.contig_index)
        return self.contig_index[chrom]

    def sites(self, sample: Optional[str] = None) -> Iterator[Site]:
        if sample is None:
            if len(self.samples) != 1:
                raise ValueError("sample must be named for a multi-sample VCF: %s" % self.path)
            col = 9
        else:
            if sample not in self.samples:
                raise KeyError("sample %s not in %s" % (sample, self.path))
            col = 9 + self.samples.index(sample)
        for line in self._fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) <= col:
                continue
            fmt = f[8].split(":")
            vals = f[col].split(":")
            fields = dict(zip(fmt, vals))
            alleles, phased = _parse_gt(fields.get("GT", "."))
            ad = None
            if "AD" in fields and fields["AD"] not in (".", ""):
                ad = tuple(_int_or_none(x) or 0 for x in fields["AD"].split(","))
            yield Site(chrom=f[0], pos=int(f[1]), ref=f[3], alt=f[4], alleles=alleles, phased=phased,
                       ps=_int_or_none(fields.get("PS", ".")), gq=_int_or_none(fields.get("GQ", ".")),
                       dp=_int_or_none(fields.get("DP", ".")), ad=ad)


def _grouped(reader: VcfReader, sample: Optional[str]):
    """Yield ((contig_rank, pos), {(ref, alt): Site}) groups in file order."""
    key = None
    group: Dict[Tuple[str, str], Site] = {}
    for s in reader.sites(sample):
        k = (reader._contig_rank(s.chrom), s.pos)
        if key is not None and k != key:
            yield key, group
            group = {}
        key = k
        group[(s.ref, s.alt)] = s
    if key is not None:
        yield key, group


def iter_trio(child: VcfReader, father: VcfReader, mother: VcfReader,
              child_sample: Optional[str] = None, father_sample: Optional[str] = None,
              mother_sample: Optional[str] = None) -> Iterator[TrioSite]:
    """Merge three sorted VCF streams; yield every child site with the matching parental sites (or None).

    The three readers must share contig order (they do when split from one joint VCF); a parent
    reader with a different contig list raises rather than silently mis-merging.
    """
    for other in (father, mother):
        if other.contigs and child.contigs and other.contigs != child.contigs:
            raise ValueError("contig order differs between %s and %s" % (child.path, other.path))
        other.contig_index = child.contig_index   # share the rank table (covers header-less contigs)

    def advance(gen):
        try:
            return next(gen)
        except StopIteration:
            return None

    fg, mg = _grouped(father, father_sample), _grouped(mother, mother_sample)
    fcur, mcur = advance(fg), advance(mg)
    for ckey, cgroup in _grouped(child, child_sample):
        while fcur is not None and fcur[0] < ckey:
            fcur = advance(fg)
        while mcur is not None and mcur[0] < ckey:
            mcur = advance(mg)
        fmatch = fcur[1] if (fcur is not None and fcur[0] == ckey) else {}
        mmatch = mcur[1] if (mcur is not None and mcur[0] == ckey) else {}
        for (ref, alt), cs in cgroup.items():
            yield TrioSite(chrom=cs.chrom, pos=cs.pos, ref=ref, alt=alt, child=cs,
                           father=fmatch.get((ref, alt)), mother=mmatch.get((ref, alt)))


# ----------------------------------------------------------------------------------------------
# manifest / pedigree
# ----------------------------------------------------------------------------------------------
MANIFEST_REQUIRED = ("sample_id", "family_id", "father_id", "mother_id", "sex", "role")


def read_manifest(path: str) -> Dict[str, dict]:
    """The cohort manifest: one row per sample, identifiers are DATA here, never embedded in code."""
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise ValueError("empty manifest: %s" % path)
    missing = [c for c in MANIFEST_REQUIRED if c not in rows[0]]
    if missing:
        raise ValueError("manifest %s lacks columns: %s" % (path, ", ".join(missing)))
    return {r["sample_id"]: r for r in rows}


def trio_of(manifest: Dict[str, dict], child: str) -> Tuple[str, str]:
    """(father_id, mother_id) for a child with BOTH parents in the manifest; duos raise (DESIGN P19)."""
    row = manifest.get(child)
    if row is None:
        raise KeyError("child %s not in manifest" % child)
    fa, mo = row.get("father_id", ""), row.get("mother_id", "")
    for pid, label in ((fa, "father"), (mo, "mother")):
        if pid in ("", "0", ".") or pid not in manifest:
            raise ValueError("child %s has no sequenced %s in the manifest (duo) - excluded, see DESIGN P19"
                             % (child, label))
    return fa, mo


def normalise_sex(x: str) -> str:
    """Manifest / PED sex codes to 'M' | 'F' | 'U'."""
    x = (x or "").strip().lower()
    if x in ("1", "m", "male", "xy"):
        return "M"
    if x in ("2", "f", "female", "xx"):
        return "F"
    return "U"
