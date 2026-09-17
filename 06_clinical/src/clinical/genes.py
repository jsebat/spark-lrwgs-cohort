"""Gene model, gene panels and constraint for the targeted clinical arm.

Everything here is a lookup: which genes does an interval touch and how (CDS, UTR, intron), is the gene on the panel,
what is its s_het. Nothing here decides anything about a variant.

Gene model: the Ensembl GFF3 shipped with the HiFi WDL resources. Its hierarchy is gene -> mRNA/transcript -> exon /
CDS / *_UTR, and only the GENE row carries the symbol (Name=); exon rows carry an exon id, UTR rows carry no name at
all. The symbol is therefore resolved through gene_of_tx, not read off the feature (the first version read Name= off
exons, matched no symbol, and dropped every UTR for having none).

Panel: SFARI (high-confidence and all) UNION DDG2P (confident and all). Membership is the clinical criterion (DESIGN
P31/P33); s_het is REPORTED but never gates. A constraint threshold excludes DNMT3A (s_het 0.0063), which is a
haploinsufficiency gene whose apparent LoF tolerance is an artefact of clonal haematopoiesis in blood-derived
population exomes (gnomAD flags its pLI unreliable, as for TET2, ASXL1, PPM1D). A filter with that failure mode would
discard one of the two pathogenic events this cohort is known to contain, systematically, for a whole class of genes.
"""
from __future__ import annotations

import bisect
import collections
import csv
import os
import subprocess
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Set, Tuple

TAB, NL = chr(9), chr(10)
FEATURE_KINDS = {"exon", "CDS", "five_prime_UTR", "three_prime_UTR"}


def _attrs(a: str) -> Dict[str, str]:
    d = {}
    for kv in a.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d


@dataclass
class GeneModel:
    """Per-contig sorted feature lists: (start0, end, symbol, kind, biotype). Query with hits()."""
    feats: Dict[str, List[Tuple[int, int, str, str]]] = field(default_factory=dict)
    starts: Dict[str, List[int]] = field(default_factory=dict)
    span: Dict[str, Tuple[str, int, int]] = field(default_factory=dict)      # symbol -> (chrom, start0, end)
    biotype: Dict[str, str] = field(default_factory=dict)
    max_len: Dict[str, int] = field(default_factory=dict)

    @classmethod
    def from_gff3(cls, path: str, protein_coding_only: bool = False) -> "GeneModel":
        sym_of_gene, gene_of_tx, bt_of_gene, raw = {}, {}, {}, []
        # the shipped .gz is a gzip member followed by plain text; zcat -f reads both, Python's gzip module rejects it
        proc = subprocess.Popen(["zcat", "-f", path], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
        with proc.stdout as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip(NL).split(TAB)
                if len(f) < 9:
                    continue
                kind = f[2]
                if kind == "gene" or kind.endswith("_gene") or kind in ("ncRNA_gene", "pseudogene"):
                    d = _attrs(f[8])
                    gid = d.get("ID", "").replace("gene:", "")
                    if gid and d.get("Name"):
                        sym_of_gene[gid] = d["Name"]
                        bt_of_gene[gid] = d.get("biotype", "")
                elif kind in ("mRNA", "transcript", "lnc_RNA", "ncRNA", "snRNA", "miRNA"):
                    d = _attrs(f[8])
                    tid = d.get("ID", "").replace("transcript:", "")
                    par = d.get("Parent", "").replace("gene:", "")
                    if tid and par:
                        gene_of_tx[tid] = par
                elif kind in FEATURE_KINDS:
                    tid = _attrs(f[8]).get("Parent", "").replace("transcript:", "")
                    if tid:
                        c = f[0] if f[0].startswith("chr") else "chr" + f[0]
                        raw.append((c, int(f[3]) - 1, int(f[4]), tid, kind))
        m = cls()
        per: Dict[str, List[Tuple[int, int, str, str]]] = collections.defaultdict(list)
        for c, a, b, tid, kind in raw:
            gid = gene_of_tx.get(tid, "")
            sym = sym_of_gene.get(gid, "")
            if not sym:
                continue
            bt = bt_of_gene.get(gid, "")
            if protein_coding_only and bt != "protein_coding":
                continue
            per[c].append((a, b, sym, kind))
            m.biotype[sym] = bt
            s = m.span.get(sym)
            m.span[sym] = (c, a if s is None else min(s[1], a), b if s is None else max(s[2], b))
        for c, v in per.items():
            v.sort()
            m.feats[c] = v
            m.starts[c] = [x[0] for x in v]
            m.max_len[c] = max((b - a) for a, b, _, _ in v) if v else 0
        return m

    def hits(self, chrom: str, start0: int, end: int) -> Dict[str, Set[str]]:
        """symbol -> set of feature kinds ('CDS', 'UTR5', 'UTR3', 'exon') overlapping [start0, end)."""
        out: Dict[str, Set[str]] = collections.defaultdict(set)
        v = self.feats.get(chrom)
        if not v:
            return out
        # features are sorted by start; anything starting before (start0 - longest feature) cannot reach start0
        i = bisect.bisect_left(self.starts[chrom], start0 - self.max_len[chrom] - 1)
        for a, b, sym, kind in v[i:]:
            if a >= end:
                break
            if b > start0:
                out[sym].add({"five_prime_UTR": "UTR5", "three_prime_UTR": "UTR3"}.get(kind, kind))
        return out

    def _span_index(self):
        """Per-contig gene spans sorted by start, built once. A linear scan of every span per candidate was
        60k genes x 1.2M candidates and did not finish; this is the same bisect window hits() uses."""
        if not hasattr(self, "_spans_by_chrom"):
            by: Dict[str, List[Tuple[int, int, str]]] = collections.defaultdict(list)
            for s, (c, a, b) in self.span.items():
                by[c].append((a, b, s))
            self._spans_by_chrom = {}
            self._span_starts = {}
            self._span_maxlen = {}
            for c, v in by.items():
                v.sort()
                self._spans_by_chrom[c] = v
                self._span_starts[c] = [x[0] for x in v]
                self._span_maxlen[c] = max(b - a for a, b, _ in v) if v else 0
        return self._spans_by_chrom, self._span_starts, self._span_maxlen

    def genes_spanning(self, chrom: str, pos1: int) -> Set[str]:
        """symbols whose gene span (first to last exon) contains pos1 -- intronic hits included."""
        spans, starts, maxlen = self._span_index()
        v = spans.get(chrom)
        if not v:
            return set()
        i = bisect.bisect_left(starts[chrom], pos1 - maxlen[chrom] - 1)
        out = set()
        for a, b, s in v[i:]:
            if a >= pos1:
                break
            if a < pos1 <= b:
                out.add(s)
        return out


def load_set(path: str) -> Set[str]:
    if not path or not os.path.exists(path):
        return set()
    return {x.strip().replace("\r", "") for x in open(path) if x.strip()}


@dataclass
class Panel:
    sfari_hc: Set[str]
    sfari_all: Set[str]
    ddg2p_conf: Set[str]
    ddg2p_all: Set[str]

    @classmethod
    def load(cls, gene_set_dir: str) -> "Panel":
        g = lambda n: load_set(os.path.join(gene_set_dir, n))
        p = cls(g("sfari_hc.txt"), g("sfari_all.txt"), g("ddg2p_confident.txt"), g("ddg2p_all.txt"))
        if not (p.sfari_all or p.ddg2p_all):
            raise SystemExit("no gene sets found in %s (need sfari_hc.txt sfari_all.txt ddg2p_confident.txt ddg2p_all.txt)" % gene_set_dir)
        return p

    @property
    def union(self) -> Set[str]:
        return self.sfari_hc | self.sfari_all | self.ddg2p_conf | self.ddg2p_all

    def label(self, gene: str) -> str:
        """'SFARI_hc;DDG2P_conf' style membership string, '-' when off-panel."""
        L = []
        if gene in self.sfari_hc:
            L.append("SFARI_hc")
        elif gene in self.sfari_all:
            L.append("SFARI")
        if gene in self.ddg2p_conf:
            L.append("DDG2P_conf")
        elif gene in self.ddg2p_all:
            L.append("DDG2P")
        return ";".join(L) or "-"

    def rank(self, gene: str) -> int:
        """1 = SFARI high-confidence, 2 = DDG2P confident, 3 = SFARI all, 4 = DDG2P all, 9 = off panel (the original
        denovo_sv_priority ordering)."""
        if gene in self.sfari_hc:
            return 1
        if gene in self.ddg2p_conf:
            return 2
        if gene in self.sfari_all:
            return 3
        if gene in self.ddg2p_all:
            return 4
        return 9


def load_shet(genebayes_tsv: str, gff3_path: Optional[str] = None, ensg2sym: Optional[Dict[str, str]] = None) -> Dict[str, float]:
    """s_het (GeneBayes post_mean, column 7) by SYMBOL. GeneBayes is keyed by Ensembl gene id; the symbol map comes
    from the GFF3 gene rows (or a prebuilt ensg->symbol dict)."""
    if ensg2sym is None:
        ensg2sym = {}
        if gff3_path:
            proc = subprocess.Popen(["zcat", "-f", gff3_path], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
            with proc.stdout as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    f = line.split(TAB)
                    if len(f) > 8 and (f[2] == "gene" or f[2].endswith("_gene")):
                        d = _attrs(f[8])
                        gid = d.get("ID", "").replace("gene:", "").split(".")[0]
                        if gid and d.get("Name"):
                            ensg2sym[gid] = d["Name"]
    out: Dict[str, float] = {}
    if not genebayes_tsv or not os.path.exists(genebayes_tsv):
        return out
    with open(genebayes_tsv, newline="") as fh:
        rd = csv.reader(fh, delimiter=TAB)
        hdr = next(rd)
        # a two-column SYMBOL<TAB>s_het file (the lab's shet_by_symbol.tsv, built from GeneBayes + GENCODE) is
        # accepted as-is; it has no header, so the first row is data
        if len(hdr) == 2 and not hdr[0].startswith("ENSG"):
            try:
                out[hdr[0]] = float(hdr[1])
            except ValueError:
                pass
            for r in rd:
                if len(r) >= 2:
                    try:
                        out[r[0]] = float(r[1])
                    except ValueError:
                        pass
            return out
        col = 6
        for i, h in enumerate(hdr):
            if h.strip().lower() in ("post_mean", "s_het", "shet"):
                col = i
                break
        for r in rd:
            if len(r) <= col:
                continue
            e = r[0].split(".")[0]
            sym = ensg2sym.get(e)
            if not sym:
                continue
            try:
                out[sym] = float(r[col])
            except ValueError:
                pass
    return out


def shet_band(v: Optional[float]) -> str:
    if v is None:
        return "unknown"
    return "high" if v >= 0.18 else ("moderate" if v >= 0.075 else "low")
