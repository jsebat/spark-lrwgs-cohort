"""Swap-closed, family-grouped folds (DESIGN P12, P14).

Units are FAMILIES (a quad's two children stay together). Outer folds are drawn by seed; within an outer fold the
pedigree-swap pairings are drawn from that fold's families only, so no synthetic trio ever combines a training child
with a test family's parents. Constraints: the single blood-derived family (REACH sample ids; found by the CLI's sample-id
prefix) is never alone in a fold with fewer
than `min_fold_size` families; families are balanced by number of complete-trio children.

Nothing here reads data: the input is the manifest's pedigree (family -> complete-trio children, DNA source), the
output is a deterministic table `family_id, outer_fold, seed` and, from swap.py, the pairing table.
"""
from __future__ import annotations

import random
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple


@dataclass(frozen=True)
class Family:
    family_id: str
    children: Tuple[str, ...]          # complete-trio children only (P19: duos excluded upstream)
    source: str = "saliva"             # saliva | blood (stratification only)


def families_from_manifest(rows: Sequence[dict], sources: Dict[str, str] = None) -> List[Family]:
    """rows: manifest dicts (sample_id family_id father_id mother_id ... role). Complete-trio children = both parents in
    the manifest. `sources` maps family_id -> DNA source (default saliva)."""
    ids = {r["sample_id"] for r in rows}
    by_fam: Dict[str, List[str]] = {}
    for r in rows:
        if r.get("role") != "offspring":
            continue
        if r.get("father_id") in ids and r.get("mother_id") in ids and r.get("father_id") not in ("0", "") and r.get("mother_id") not in ("0", ""):
            by_fam.setdefault(r["family_id"], []).append(r["sample_id"])
    out = [Family(f, tuple(sorted(c)), (sources or {}).get(f, "saliva")) for f, c in sorted(by_fam.items())]
    return out


def outer_folds(fams: Sequence[Family], n_folds: int = 5, seed: int = 0, min_fold_size: int = 3) -> Dict[str, int]:
    """family_id -> fold. Greedy balanced assignment by child count in a seeded random order; a blood family is placed
    first into fold 0 so it always has company. Deterministic for (fams, n_folds, seed)."""
    rng = random.Random(seed)
    order = list(fams)
    rng.shuffle(order)
    blood = [f for f in order if f.source == "blood"]
    rest = [f for f in order if f.source != "blood"]
    # largest families first so the greedy balance works; ties broken by the shuffled order
    rest.sort(key=lambda f: -len(f.children))
    load = [0] * n_folds
    assign: Dict[str, int] = {}
    for f in blood:
        assign[f.family_id] = 0
        load[0] += len(f.children)
    for f in rest:
        k = min(range(n_folds), key=lambda i: (load[i], i))
        assign[f.family_id] = k
        load[k] += len(f.children)
    sizes = [sum(1 for v in assign.values() if v == i) for i in range(n_folds)]
    if min(sizes) < min_fold_size:
        raise ValueError("fold too small: sizes %s (n_folds=%d for %d families)" % (sizes, n_folds, len(fams)))
    return assign


def inner_folds(train_families: Sequence[Family], n_folds: int = 4, seed: int = 0) -> Dict[str, int]:
    """Inner folds from the OUTER-TRAINING families only (P14), same construction."""
    return outer_folds(train_families, n_folds=n_folds, seed=seed + 1000, min_fold_size=1)


def fold_table(assign: Dict[str, int], seed: int) -> List[dict]:
    return [dict(family_id=f, outer_fold=k, seed=seed) for f, k in sorted(assign.items(), key=lambda kv: (kv[1], kv[0]))]
