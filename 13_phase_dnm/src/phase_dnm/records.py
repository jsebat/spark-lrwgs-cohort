"""The shared candidate record (README §2.2): one normalised row per putative de novo event of ANY class.

Everything class-specific lives in `class_payload` (JSON) until the six-haplotype extractor asks the class adapter
the one class-specific question - which reads support the alt allele. Identifiers are values in the row, never
in code.
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass, field, fields
from typing import Dict, Iterable, Iterator, List, Optional

CLASSES = ("SNV", "INDEL", "SV", "TR")
TIERS = ("HIGH", "LOW", "UNFILTERED")


@dataclass
class CandidateRecord:
    family_id: str
    sample_id: str                 # the child
    variant_id: str                # chrom:pos:ref:alt for small variants; caller id for SV; TRID for TR
    chrom: str
    start: int                     # 1-based first affected base (VCF POS for small variants and TR loci)
    end: int
    ref: str
    alt: str                       # the candidate alt allele (one allele per row)
    variant_class: str             # SNV | INDEL | SV | TR
    caller: str                    # deepvariant_glnexus | sawfish | trgt | <list name>
    caller_gt: str = "."
    caller_gq: Optional[int] = None
    caller_dp: Optional[int] = None
    caller_qual: Optional[float] = None
    caller_filter: str = "."
    child_ad: str = "."
    father_gt: str = "."
    mother_gt: str = "."
    father_gq: Optional[int] = None
    mother_gq: Optional[int] = None
    father_dp: Optional[int] = None
    mother_dp: Optional[int] = None
    father_ad: str = "."
    mother_ad: str = "."
    child_pl: str = "."                  # Phred-scaled genotype likelihoods, comma-joined (SynthDNM universal set)
    father_pl: str = "."
    mother_pl: str = "."
    source_tier: str = "UNFILTERED"      # HIGH | LOW | UNFILTERED (P5)
    source_list: str = "joint_vcf"       # which list(s) this row came from, ';'-joined
    mask_overlap: int = 0                # filled by the context step; a FLAG, never a filter (P5)
    class_payload: Dict = field(default_factory=dict)

    def __post_init__(self):
        if self.variant_class not in CLASSES:
            raise ValueError("variant_class must be one of %s, got %r" % (CLASSES, self.variant_class))
        if self.source_tier not in TIERS:
            raise ValueError("source_tier must be one of %s" % (TIERS,))

    @property
    def key(self):
        return (self.sample_id, self.variant_class, self.variant_id)


COLUMNS: List[str] = [f.name for f in fields(CandidateRecord)]


def write_candidates(records: Iterable[CandidateRecord], path: str) -> int:
    n = 0
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COLUMNS)
        for r in records:
            d = asdict(r)
            d["class_payload"] = json.dumps(d["class_payload"], separators=(",", ":"), sort_keys=True)
            w.writerow(["" if d[c] is None else d[c] for c in COLUMNS])
            n += 1
    return n


def _opt_int(x: str) -> Optional[int]:
    return None if x in ("", ".") else int(x)


def read_candidates(path: str) -> Iterator[CandidateRecord]:
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            yield CandidateRecord(
                family_id=r["family_id"], sample_id=r["sample_id"], variant_id=r["variant_id"], chrom=r["chrom"],
                start=int(r["start"]), end=int(r["end"]), ref=r["ref"], alt=r["alt"], variant_class=r["variant_class"],
                caller=r["caller"], caller_gt=r["caller_gt"], caller_gq=_opt_int(r["caller_gq"]), caller_dp=_opt_int(r["caller_dp"]),
                caller_qual=None if r["caller_qual"] in ("", ".") else float(r["caller_qual"]), caller_filter=r["caller_filter"],
                child_ad=r["child_ad"], father_gt=r["father_gt"], mother_gt=r["mother_gt"],
                father_gq=_opt_int(r["father_gq"]), mother_gq=_opt_int(r["mother_gq"]),
                father_dp=_opt_int(r["father_dp"]), mother_dp=_opt_int(r["mother_dp"]),
                father_ad=r["father_ad"], mother_ad=r["mother_ad"],
                child_pl=r.get("child_pl", "."), father_pl=r.get("father_pl", "."), mother_pl=r.get("mother_pl", "."),
                source_tier=r["source_tier"], source_list=r["source_list"],
                mask_overlap=int(r["mask_overlap"] or 0), class_payload=json.loads(r["class_payload"] or "{}"))
