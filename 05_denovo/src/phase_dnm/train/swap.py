"""Pedigree swap inside folds (DESIGN P12, P23): SynthDNM's pairwise offspring exchange, with the pair drawn from ONE
outer fold. Family A's child is placed with family B's parents and vice versa; every complete-trio child of A gets B's
parents (a quad's two children both go to the same surrogate parents, so the synthetic family keeps its sib structure).
An odd family in a fold is paired with a second family of the same fold in a second within-fold draw (its child is used
once more as a synthetic child; the real trio is unaffected).

Output: one row per synthetic trio - synthetic_id, child, surrogate_father, surrogate_mother, child_family,
parent_family, outer_fold, seed, draw. Synthetic ids never contain sample identifiers (a hash of the pair).
"""
from __future__ import annotations

import hashlib
import random
from typing import Dict, List, Sequence

from .folds import Family


def _sid(child: str, father: str, mother: str, seed: int) -> str:
    return "syn_" + hashlib.blake2b(("%s|%s|%s|%d" % (child, father, mother, seed)).encode(), digest_size=6).hexdigest()


def pairings(fams: Sequence[Family], assign: Dict[str, int], parents: Dict[str, tuple], seed: int = 0) -> List[dict]:
    """parents: family_id -> (father_id, mother_id) of the complete trio(s) in that family (one parent pair per family).
    Returns the synthetic trio table for one seed; every family appears as a child donor at least once."""
    rng = random.Random(seed)
    out: List[dict] = []
    by_fold: Dict[int, List[Family]] = {}
    for f in fams:
        by_fold.setdefault(assign[f.family_id], []).append(f)
    for fold, members in sorted(by_fold.items()):
        members = sorted(members, key=lambda f: f.family_id)
        rng.shuffle(members)
        pairs = [(members[i], members[i + 1]) for i in range(0, len(members) - 1, 2)]
        draw = 0
        for a, b in pairs:
            for donor, host in ((a, b), (b, a)):
                fa, mo = parents[host.family_id]
                for child in donor.children:
                    out.append(dict(synthetic_id=_sid(child, fa, mo, seed), child=child, surrogate_father=fa, surrogate_mother=mo,
                                    child_family=donor.family_id, parent_family=host.family_id, outer_fold=fold, seed=seed, draw=draw))
        if len(members) % 2 == 1:
            odd = members[-1]
            partner = rng.choice(members[:-1]) if len(members) > 1 else None
            if partner is not None:
                fa, mo = parents[partner.family_id]
                for child in odd.children:
                    out.append(dict(synthetic_id=_sid(child, fa, mo, seed), child=child, surrogate_father=fa, surrogate_mother=mo,
                                    child_family=odd.family_id, parent_family=partner.family_id, outer_fold=fold, seed=seed, draw=1))
    return out


def check_closed(table: Sequence[dict], assign: Dict[str, int]) -> None:
    """Every synthetic trio's child family and parent family share the outer fold (P12)."""
    for r in table:
        if assign[r["child_family"]] != assign[r["parent_family"]] or r["child_family"] == r["parent_family"]:
            raise ValueError("swap crosses a fold or pairs a family with itself: %s" % {k: r[k] for k in ("child_family", "parent_family", "outer_fold")})
