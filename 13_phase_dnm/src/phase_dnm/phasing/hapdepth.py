"""hapdepth — per-haplotype read depth in fixed bins from a HiPhase-haplotagged BAM (M1, one pass per sample).

Feeds `hap_obs` (P7) and the SV adapter's per-haplotype depth inside vs. flanks without re-reading BAMs in M2.
Depth is accumulated over each primary read's REFERENCE SPAN (reference_start..reference_end), not over its aligned
blocks: HiFi alignments carry ~100 small indel gaps per read, so block-wise accumulation is ~100x slower for a
depth error far below one read per kb bin. Reads are split into hap1 / hap2 / untagged (MAPQ >= min_mapq) and a
fourth low-MAPQ category; secondary, supplementary, duplicate and unmapped records are skipped.

pysam is imported lazily so the rest of the package (and its tests) stays pysam-free.
"""
from __future__ import annotations

import gzip
import json
from typing import Dict, Iterable, List, Optional, Tuple

COLUMNS = ["chrom", "start", "end", "dp_hap1", "dp_hap2", "dp_untagged", "dp_lowmapq"]


def _bins_for(contig_len: int, bin_size: int) -> int:
    return contig_len // bin_size + 1


def hapdepth(bam: str, bai: Optional[str], out_tsv_gz: str, out_summary_json: Optional[str] = None,
             bin_size: int = 1000, min_mapq: int = 20, contigs: Optional[Iterable[str]] = None,
             region: Optional[Tuple[str, int, int]] = None) -> Dict[str, dict]:
    import pysam  # lazy: only the read-level steps need it
    try:
        import numpy as np
    except ImportError:  # pragma: no cover - numpy ships with pysam environments
        np = None
    af = pysam.AlignmentFile(bam, "rb", index_filename=bai) if bai else pysam.AlignmentFile(bam, "rb")
    lengths = dict(zip(af.references, af.lengths))
    todo: List[Tuple[str, int, int]]
    if region:
        todo = [region]
    else:
        names = list(contigs) if contigs else [c for c in af.references if not c.startswith(("chrUn", "HLA")) and "_" not in c]
        todo = [(c, 0, lengths[c]) for c in names]
    summary: Dict[str, dict] = {}
    with gzip.open(out_tsv_gz, "wt") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for chrom, rstart, rend in todo:
            nb = _bins_for(lengths[chrom], bin_size)
            first_bin, last_bin = rstart // bin_size, (max(rstart, rend - 1)) // bin_size
            width = last_bin - first_bin + 1
            acc = [[0.0] * width for _ in range(4)] if np is None else np.zeros((4, width), dtype="float64")
            n_reads = 0
            for read in af.fetch(chrom, rstart, rend):
                if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate:
                    continue
                s, e = read.reference_start, read.reference_end
                if e is None or e <= rstart or s >= rend:
                    continue
                n_reads += 1
                if read.mapping_quality < min_mapq:
                    cat = 3
                else:
                    hp = read.get_tag("HP") if read.has_tag("HP") else 0
                    cat = 0 if hp == 1 else (1 if hp == 2 else 2)
                b0, b1 = max(s, rstart) // bin_size, (min(e, rend) - 1) // bin_size
                for b in range(b0, b1 + 1):
                    lo, hi = max(s, b * bin_size, rstart), min(e, (b + 1) * bin_size, rend)
                    if hi > lo:
                        acc[cat][b - first_bin] += hi - lo
            tot = [0.0, 0.0, 0.0, 0.0]
            for i in range(width):
                b = first_bin + i
                start, end = b * bin_size, min((b + 1) * bin_size, lengths[chrom])
                span = max(1, end - start)
                vals = [acc[c][i] / span for c in range(4)]
                for c in range(4):
                    tot[c] += vals[c]
                fh.write("%s\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\n" % (chrom, start, end, *vals))
            summary[chrom] = dict(n_bins=width, n_reads=n_reads,
                                  mean_dp_hap1=round(tot[0] / width, 3), mean_dp_hap2=round(tot[1] / width, 3),
                                  mean_dp_untagged=round(tot[2] / width, 3), mean_dp_lowmapq=round(tot[3] / width, 3))
    af.close()
    if out_summary_json:
        with open(out_summary_json, "w") as fh:
            json.dump({"bin_size": bin_size, "min_mapq": min_mapq, "per_chrom": summary}, fh, indent=1, sort_keys=True)
    return summary
