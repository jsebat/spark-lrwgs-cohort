#!/usr/bin/env python
"""Step 02 - harmonize published episignature tables into hg38 region BEDs + INVENTORY.md.

Inputs
  config/config.yaml            thresholds.region_pad, paths, signatures.layer1
  signatures/sources.yaml       one entry per signature; each lists raw supplementary files and how to parse them
                                (written by hand after inspecting signatures/raw/<ID>/SOURCE.md; never guessed)
  resources/manifests/*.hg38.manifest.tsv.gz   Zhou-lab InfiniumAnnotation (probe -> hg38 CpG)
  resources/liftover/hg19ToHg38.over.chain.gz  for DMR tables published in hg19

Outputs
  signatures/hg38/<ID>.bed          chrom start end signature source direction published_delta_beta n_array_probes
  signatures/hg38/<ID>.tier1.bed    (TBRS_LOF only) regions supported by >=2 sources with concordant direction
  results/02_signatures/<ID>.regions.tsv   per-region detail incl. n_sources, per-source deltas, conflicts
  results/02_signatures/<ID>.probes.tsv    per-probe mapping audit (mapped / unmapped / masked)
  signatures/INVENTORY.md

Rules (CLAUDE.md): no approximation - a signature with no downloadable probe/DMR list is DROPPED with the reason;
probe coordinates come only from the manifest; hg19 DMRs are lifted, never hand-converted.
"""
import argparse
import datetime as dt
import gzip
import pathlib
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
import yaml

BED_COLS = ["chrom", "start", "end", "signature", "source", "direction", "published_delta_beta", "n_array_probes"]
MAIN_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


# ----------------------------------------------------------------------------- manifests
def load_manifest(path):
    """Zhou InfiniumAnnotation hg38 manifest -> DataFrame indexed by probe ID with chrom, pos (1-based), masked."""
    df = pd.read_csv(path, sep="\t", low_memory=False, compression="gzip")
    cols = {c.lower(): c for c in df.columns}
    pid = cols.get("probe_id") or cols.get("probeid") or cols.get("name")
    chrom = cols.get("cpg_chrm") or cols.get("chrm") or cols.get("chr")
    beg = cols.get("cpg_beg") or cols.get("beg") or cols.get("start")
    if not (pid and chrom and beg):
        sys.exit(f"manifest {path}: cannot find probe/chrom/start columns in {list(df.columns)[:15]}")
    ids = df[pid].astype(str).str.replace(r"_[TB][CO]\d+$", "", regex=True)   # EPICv2 suffixes -> base cg ID
    masked = pd.Series(False, index=df.index)
    if cols.get("mask_general"):
        masked = df[cols["mask_general"]].astype(str).str.lower().eq("true")
    else:  # InfiniumAnnotation >= v8: masks live in <plat>.hg38.mask.tsv.gz (Probe_ID, mask, maskUniq, M_general)
        mpath = pathlib.Path(str(path).replace(".manifest.tsv.gz", ".mask.tsv.gz"))
        if mpath.exists():
            mk = pd.read_csv(mpath, sep="\t", compression="gzip")
            bad = set(mk.loc[mk["M_general"].astype(str).str.lower().eq("true"), "Probe_ID"].astype(str)
                      .str.replace(r"_[TB][CO]\d+$", "", regex=True))
            masked = ids.isin(bad)
        else:
            print(f"  WARNING: no mask file for {path.name}; no probes masked")
    out = pd.DataFrame({"chrom": df[chrom].astype(str).to_numpy(), "start0": pd.to_numeric(df[beg], errors="coerce").to_numpy(),
                        "masked": masked.to_numpy()}, index=ids.to_numpy())
    out = out[out["chrom"].isin(MAIN_CHROMS) & out["start0"].notna()]
    out["start0"] = out["start0"].astype(int)          # Zhou CpG_beg is 0-based start of the C
    return out[~out.index.duplicated()]


def load_manifests(res_dir):
    mans = {}
    for f in sorted(pathlib.Path(res_dir, "manifests").glob("*.hg38.manifest.tsv.gz")):
        key = f.name.split(".")[0].upper()      # HM450 / EPIC / EPICv2
        mans[key] = load_manifest(f)
        print(f"manifest {key}: {len(mans[key]):,} probes on main chromosomes")
    if not mans:
        sys.exit("no hg38 manifests found under resources/manifests/")
    return mans


# ----------------------------------------------------------------------------- liftover
def liftover_intervals(df, chain, workdir):
    """hg19 -> hg38 for a DataFrame with chrom,start0,end (0-based half-open). Uses UCSC liftOver if present, else pyliftover."""
    exe = shutil.which("liftOver")
    if exe:
        workdir.mkdir(parents=True, exist_ok=True)
        src, dst, unm = workdir / "in.bed", workdir / "out.bed", workdir / "unmapped.bed"
        df.assign(name=range(len(df)))[["chrom", "start0", "end", "name"]].to_csv(src, sep="\t", header=False, index=False)
        subprocess.run([exe, "-minMatch=0.95", str(src), str(chain), str(dst), str(unm)], check=True)
        out = pd.read_csv(dst, sep="\t", header=None, names=["chrom", "start0", "end", "name"])
        lifted = df.iloc[out["name"].values].copy()
        lifted[["chrom", "start0", "end"]] = out[["chrom", "start0", "end"]].values
        return lifted
    try:
        from pyliftover import LiftOver
    except ImportError:
        sys.exit("need UCSC liftOver on PATH or `pip install pyliftover` to lift hg19 DMRs")
    lo = LiftOver(str(chain))
    rows = []
    for i, r in df.iterrows():
        a = lo.convert_coordinate(r["chrom"], int(r["start0"]))
        b = lo.convert_coordinate(r["chrom"], int(r["end"]) - 1)
        if a and b and a[0][0] == b[0][0] == r["chrom"] and 0 < (b[0][1] + 1 - a[0][1]) <= 2 * (r["end"] - r["start0"]):
            rr = r.copy(); rr["start0"], rr["end"] = a[0][1], b[0][1] + 1
            rows.append(rr)
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- source parsing
PROBE_RE = r"^(cg|ch\.)\d+"


def read_table(src):
    """Read one supplementary table as described in sources.yaml.

    Options: sheet, header (row index or "auto" = row above the first probe-like cell), names, sep, query,
    adapter (arefeshghi2020 | smith2021 | pdf_regions) for tables that need more than column selection.
    """
    p = pathlib.Path(src["file"])
    if not p.exists():
        raise FileNotFoundError(p)
    adapter = src.get("adapter")
    if adapter == "pdf_regions":
        return read_pdf_regions(p, src)
    hdr = src.get("header", 0)
    is_xl = p.suffix.lower() in (".xlsx", ".xls")
    if hdr == "auto":
        raw = pd.read_excel(p, sheet_name=src.get("sheet", 0), header=None) if is_xl else pd.read_csv(p, sep=src.get("sep", "\t"), header=None)
        first = next(i for i in range(len(raw)) if raw.iloc[i].astype(str).str.match(PROBE_RE).any())
        df = raw.iloc[first:].copy(); df.columns = [str(c).strip() for c in raw.iloc[first - 1]]
    elif is_xl:
        df = pd.read_excel(p, sheet_name=src.get("sheet", 0), header=hdr, names=src.get("names"))
    else:
        df = pd.read_csv(p, sep=src.get("sep", "\t"), header=hdr, names=src.get("names"))
    if src.get("probe_col"):   # drop legend / footnote rows
        df = df[df[src["probe_col"]].astype(str).str.match(PROBE_RE)]
    if src.get("query"):
        df = df.query(src["query"])
    if adapter == "arefeshghi2020":
        df = adapt_arefeshghi2020(p, df, src)
    elif adapter == "smith2021":
        df = adapt_smith2021(df, src)
    return df


def adapt_arefeshghi2020(p, s2, src):
    """EpiSign v2 table: S2 boolean membership column -> probes; delta = S3 mean(disorder) - S3 mean(Control)."""
    s3 = pd.read_excel(p, sheet_name=src["means_sheet"], header=src.get("means_header", 1))
    s3 = s3.set_index(s3.columns[0])
    sel = s2[s2[src["member_col"]].astype(str).str.lower().isin(["true", "1", "1.0"])].copy()
    sel["delta_beta"] = (s3[src["means_col"]] - s3[src.get("control_col", "Control")]).reindex(sel[src["probe_col"]].astype(str)).to_numpy()
    return sel


def adapt_smith2021(df, src):
    """Smith 2021 WGBS DMR sheets: per-sample mean-methylation columns; delta = mean(cases) - mean(controls)."""
    cols = [c for c in df.columns if c not in ("dmr.id", "chromosome", "start", "end", "closest_gene", "islands", "shores",
                                                "shelves", "genes", "enhancers", "promoters", "tss", "size")]
    ctrl = [c for c in cols if pd.Series([c]).str.contains(src["control_regex"], regex=True)[0]]
    case = [c for c in cols if any(u in c for u in src["case_ids"])]
    if not ctrl or len(case) != len(src["case_ids"]):
        sys.exit(f"smith2021 adapter: matched {len(ctrl)} control and {len(case)} case columns (expected {len(src['case_ids'])} cases)")
    df = df.copy()
    df["delta_beta"] = df[case].astype(float).mean(axis=1) - df[ctrl].astype(float).mean(axis=1)
    df.attrs["smith_groups"] = dict(n_control=len(ctrl), n_case=len(case))
    return df


def read_pdf_regions(p, src):
    """Regex-extract region rows from a PDF table (e.g. Kernohan 2016 Additional file 1)."""
    import pypdf
    rx = src["regex"]
    rows = []
    for page in pypdf.PdfReader(str(p)).pages:
        for line in page.extract_text().splitlines():
            mm = __import__("re").match(rx, line.strip())
            if mm:
                rows.append(mm.groupdict())
    df = pd.DataFrame(rows)
    for c in df.columns:
        if c != "chrom":
            num = pd.to_numeric(df[c], errors="coerce")
            if num.notna().all():
                df[c] = num
    if src.get("expect_rows") and len(df) != src["expect_rows"]:
        sys.exit(f"pdf_regions {src['name']}: parsed {len(df)} rows, expected {src['expect_rows']}")
    return df


def direction_from(df, src):
    """Return +1 / -1 per row: case-vs-control direction. Uses delta sign, or an explicit direction column."""
    if src.get("delta_col"):
        d = pd.to_numeric(df[src["delta_col"]], errors="coerce")
        if src.get("delta_scale") == "percent":
            d = d / 100.0
        if src.get("delta_sign_flip"):
            d = -d
        return np.sign(d).astype("Int64"), d
    if src.get("direction_col"):
        s = df[src["direction_col"]].astype(str).str.lower()
        pos = s.str.contains(src.get("hyper_token", "hyper|gain|up|\\+"), regex=True)
        neg = s.str.contains(src.get("hypo_token", "hypo|loss|down|-"), regex=True)
        d = pd.Series(pd.NA, index=df.index, dtype="Int64")
        d[pos] = 1; d[neg] = -1
        return d, pd.Series(np.nan, index=df.index)
    if src.get("direction_fixed") in (1, -1):
        return pd.Series(src["direction_fixed"], index=df.index, dtype="Int64"), pd.Series(np.nan, index=df.index)
    raise ValueError(f"source {src['name']}: need delta_col, direction_col or direction_fixed")


def probes_to_cpgs(df, src, mans):
    """Probe table -> per-CpG rows (chrom, start0, end, direction, delta, probe, source)."""
    probes = df[src["probe_col"]].astype(str).str.strip()
    key = src.get("array", "EPIC").upper()
    man = mans[key] if key in mans else next(iter(mans.values()))
    # fall back across arrays: a probe missing from EPIC may be in HM450 and vice versa
    lookup = pd.concat([man] + [m for k, m in mans.items() if m is not man]).pipe(lambda x: x[~x.index.duplicated()])
    dirs, deltas = direction_from(df, src)
    hit = lookup.reindex(probes.values)
    out = pd.DataFrame({"probe": probes.values, "chrom": hit["chrom"].values, "start0": hit["start0"].values,
                        "masked": hit["masked"].fillna(False).values, "direction": dirs.values, "delta": deltas.values,
                        "source": src["name"]})
    out["mapped"] = out["chrom"].notna()
    out["end"] = out["start0"] + 2      # CpG dinucleotide
    return out


def dmrs_to_intervals(df, src, chain, workdir):
    if src.get("coord_col"):   # "chrN:start-end" strings
        parts = df[src["coord_col"]].astype(str).str.extract(r"^(chr\w+):(\d+)-(\d+)$")
        d = pd.DataFrame({"chrom": parts[0], "start0": pd.to_numeric(parts[1]), "end": pd.to_numeric(parts[2])}, index=df.index)
    else:
        d = pd.DataFrame({"chrom": df[src["chrom_col"]].astype(str), "start0": pd.to_numeric(df[src["start_col"]], errors="coerce"),
                          "end": pd.to_numeric(df[src["end_col"]], errors="coerce")}, index=df.index)
    d["chrom"] = d["chrom"].where(d["chrom"].str.startswith("chr"), "chr" + d["chrom"])
    if src.get("one_based", True):
        d["start0"] = d["start0"] - 1
    d = d.dropna().astype({"start0": int, "end": int})
    d["direction"], d["delta"] = direction_from(df.loc[d.index], src)
    d["source"] = src["name"]; d["probe"] = pd.NA; d["masked"] = False; d["mapped"] = True
    if src.get("build", "hg38").lower() == "hg19":
        n0 = len(d); d = liftover_intervals(d, chain, workdir); print(f"  liftOver {src['name']}: {len(d)}/{n0} intervals")
    return d


# ----------------------------------------------------------------------------- region building
def merge_regions(cpgs, pad, sig):
    """Pad each CpG/DMR, merge overlaps, aggregate direction/delta/source support."""
    c = cpgs[cpgs["mapped"] & ~cpgs["masked"]].copy()
    c["rs"] = (c["start0"] - pad).clip(lower=0); c["re"] = c["end"] + pad
    c = c.sort_values(["chrom", "rs"]).reset_index(drop=True)
    regions, cur = [], None
    for r in c.itertuples(index=False):
        if cur and r.chrom == cur["chrom"] and r.rs <= cur["end"]:
            cur["end"] = max(cur["end"], r.re); cur["rows"].append(r)
        else:
            if cur: regions.append(cur)
            cur = dict(chrom=r.chrom, start=int(r.rs), end=int(r.re), rows=[r])
    if cur: regions.append(cur)
    out = []
    for g in regions:
        rows = pd.DataFrame(g["rows"])
        dirs = rows.dropna(subset=["direction"])
        by_src = dirs.groupby("source")["direction"].agg(lambda s: int(np.sign(s.astype(float).sum())) or 0)
        n_src = int((by_src != 0).sum())
        concordant = n_src > 0 and (by_src[by_src != 0].nunique() == 1)
        net = int(np.sign(dirs["direction"].astype(float).sum())) if len(dirs) else 0
        deltas = rows["delta"].dropna()
        out.append(dict(chrom=g["chrom"], start=g["start"], end=g["end"], signature=sig,
                        source=",".join(sorted(set(rows["source"]))), direction={1: "+", -1: "-"}.get(net, "."),
                        published_delta_beta=round(float(deltas.mean()), 4) if len(deltas) else "NA",
                        n_array_probes=int(rows["probe"].notna().sum()), n_sources=n_src, concordant=concordant,
                        within_region_conflict=bool(len(dirs) and dirs["direction"].nunique() > 1)))
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", required=True)
    ap.add_argument("--sources", default="signatures/sources.yaml")
    args = ap.parse_args()
    cfg = yaml.safe_load(open(args.config))
    srcs = yaml.safe_load(open(args.sources))
    pad = int(cfg["thresholds"]["region_pad"]); min_regions = int(cfg["thresholds"]["min_regions_per_signature"])
    res = pathlib.Path(cfg["paths"]["resources_dir"]); out_dir = pathlib.Path(cfg["paths"]["signatures_dir"])
    det_dir = pathlib.Path(cfg["paths"]["output_root"], "02_signatures"); det_dir.mkdir(parents=True, exist_ok=True)
    chain = res / "liftover" / "hg19ToHg38.over.chain.gz"
    mans = load_manifests(res)

    inventory = []
    for sig in cfg["signatures"]["layer1"]:
        spec = srcs.get(sig) or {}
        sources = spec.get("sources") or []
        print(f"\n== {sig}: {len(sources)} source(s)")
        if not sources or spec.get("status") == "DROPPED":
            inventory.append(dict(ID=sig, sources="-", n_probes=0, n_regions="-", direction="-", status="DROPPED", reason=spec.get("reason", "no downloadable probe/DMR list")))
            for f in (out_dir / f"{sig}.bed", out_dir / f"{sig}.tier1.bed"):
                f.write_text("")           # empty output keeps Snakemake happy but downstream treats it as absent
            continue
        parts = []
        for s in sources:
            df = read_table(s)
            if s.get("kind", "probe") == "probe":
                parts.append(probes_to_cpgs(df, s, mans))
            else:
                parts.append(dmrs_to_intervals(df, s, chain, det_dir / "liftover_tmp" / sig / s["name"]))
            p = parts[-1]
            print(f"  {s['name']}: {len(p)} rows, mapped {int(p['mapped'].sum())}, masked {int(p['masked'].sum())}, "
                  f"hyper {int((p['direction'] == 1).sum())} / hypo {int((p['direction'] == -1).sum())}")
        cp = pd.concat(parts, ignore_index=True)
        cp.to_csv(det_dir / f"{sig}.probes.tsv", sep="\t", index=False)
        reg = merge_regions(cp, pad, sig)
        reg.to_csv(det_dir / f"{sig}.regions.tsv", sep="\t", index=False)
        reg[BED_COLS].to_csv(out_dir / f"{sig}.bed", sep="\t", index=False, header=False)
        n_hyper, n_hypo = int((reg["direction"] == "+").sum()), int((reg["direction"] == "-").sum())
        note = []
        if sig == "TBRS_LOF":
            t1 = reg[(reg["n_sources"] >= 2) & reg["concordant"] & ~reg["within_region_conflict"]]
            t1[BED_COLS].to_csv(out_dir / f"{sig}.tier1.bed", sep="\t", index=False, header=False)
            note.append(f"tier1 {len(t1)} regions ({int((t1['direction']=='-').sum())} hypo)")
            n_reg = f"{len(t1)} / {len(reg)}"
        else:
            n_reg = str(len(reg))
        status = "OK" if len(reg) >= min_regions else "DROPPED"
        if status == "DROPPED":
            note.append(f"<{min_regions} regions")
        n_conf = int(reg["within_region_conflict"].sum())
        if n_conf:
            note.append(f"{n_conf} merged regions with mixed direction (direction '.')")
        note.append(f"median region {int(reg.eval('end-start').median())} bp")
        inventory.append(dict(ID=sig, sources="; ".join(dict.fromkeys(s["name"] for s in sources)), n_probes=int(cp["probe"].notna().sum()),
                              n_regions=n_reg, direction=f"{n_hyper} hyper / {n_hypo} hypo", status=status, reason="; ".join(note),
                              notes=[f"{s['name']}: {s.get('note', '')}" for s in sources if s.get("note")]))

    lines = ["# Signature inventory (PLAN step 02)", "",
             f"Generated {dt.datetime.now().isoformat(timespec='seconds')} by workflow/scripts/02_harmonize_signatures.py. "
             f"region_pad = {pad} bp; regions merged within signature; hg38; BED 0-based half-open. "
             "Source provenance: signatures/raw/<ID>/SOURCE.md and signatures/sources.yaml.", "",
             "| ID | Sources | n array probes (mapped) | n regions (tier1 / tier2) | Direction | Status | Notes |", "|---|---|---|---|---|---|---|"]
    for r in inventory:
        lines.append(f"| {r['ID']} | {r['sources']} | {r['n_probes']} | {r['n_regions']} | {r['direction']} | {r['status']} | {r['reason']} |")
    lines += ["", "Direction = sign of case-minus-control delta beta; 'hyper' = more methylated in cases. "
              "TBRS_LOF tier1 = regions supported by >=2 independent sources with concordant direction and no within-region conflict; "
              "tier2 = union of all sources. Probe coordinates come from the Zhou-lab hg38 manifests (masked probes excluded); "
              "hg19 DMRs were lifted with hg19ToHg38.over.chain (see resources/*/SOURCE.md).", "", "## Source notes and caveats", ""]
    for r in inventory:
        lines.append(f"**{r['ID']}**")
        for n in r.get("notes", []):
            lines.append(f"- {n}")
        if r["status"] == "DROPPED":
            lines.append(f"- DROPPED: {r['reason']}")
        lines.append("")
    pathlib.Path("signatures/INVENTORY.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines[5:]))


if __name__ == "__main__":
    main()
