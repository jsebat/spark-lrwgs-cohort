"""The feature registry (config/features.yaml) as code (DESIGN P11).

Every feature the module emits is declared with the classes it applies to and an `rf_safe` flag. This module is
the single gate between "everything M2 knows" and "what the classifier may see": `rf_matrix_columns` returns the
rf_safe subset for a class, and `assert_rf_safe` refuses any column set that contains a transmission- or
parent-of-origin-dependent feature. The YAML's hash goes into the training manifest so a model can always be
traced to the exact registry it was trained under.
"""
from __future__ import annotations

import hashlib
import os
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional

import yaml

DEFAULT_PATH = os.path.join(os.path.dirname(__file__), "..", "..", "..", "config", "features.yaml")
SECTION_PREFIXES = ("A_", "B_", "C_", "D_")


@dataclass(frozen=True)
class Feature:
    name: str
    classes: tuple
    source: str
    rf_safe: bool
    status: str
    section: str
    why: str = ""
    symmetric: str = ""
    synthdnm: str = ""


class RfSafetyError(ValueError):
    pass


class Registry:
    def __init__(self, path: Optional[str] = None):
        self.path = os.path.abspath(path or DEFAULT_PATH)
        raw = open(self.path, "rb").read()
        self.sha256 = hashlib.sha256(raw).hexdigest()
        doc = yaml.safe_load(raw.decode("utf-8"))
        self.version = str(doc.get("version"))
        self.k_thresholds = tuple(doc.get("k_thresholds", [3, 5]))
        self.features: Dict[str, Feature] = {}
        self.excluded: List[str] = [e["name"] for e in doc.get("excluded", [])]
        for section, items in doc.items():
            if not section.startswith(SECTION_PREFIXES) or not isinstance(items, list):
                continue
            for it in items:
                if not isinstance(it, dict) or "name" not in it:
                    continue
                f = Feature(name=it["name"], classes=tuple(it.get("classes", [])), source=it.get("source", ""),
                            rf_safe=bool(it.get("rf_safe", False)), status=it.get("status", "add"), section=section,
                            why=str(it.get("why", "")), symmetric=str(it.get("symmetric", "")), synthdnm=str(it.get("synthdnm", "")))
                if f.name in self.features:
                    raise ValueError("duplicate feature %s in %s" % (f.name, self.path))
                if f.status != "drop":
                    self.features[f.name] = f

    # -- queries -------------------------------------------------------------------------------
    def for_class(self, vclass: str) -> List[Feature]:
        return [f for f in self.features.values() if vclass in f.classes]

    def rf_matrix_columns(self, vclass: str) -> List[str]:
        """Columns the classifier may see for this class: rf_safe only, registry order."""
        return [f.name for f in self.for_class(vclass) if f.rf_safe]

    def unsafe_columns(self) -> List[str]:
        return [f.name for f in self.features.values() if not f.rf_safe]

    def assert_rf_safe(self, columns: Iterable[str]) -> None:
        bad = [c for c in columns if c in self.features and not self.features[c].rf_safe]
        unknown = [c for c in columns if c not in self.features]
        if bad:
            raise RfSafetyError("rf_safe violation: %s are transmission/parent-of-origin dependent (DESIGN P11)" % bad)
        if unknown:
            raise RfSafetyError("unregistered feature columns in a classifier matrix: %s" % unknown)

    def manifest(self) -> dict:
        return {"features_yaml": self.path, "features_sha256": self.sha256, "features_version": self.version,
                "n_features": len(self.features), "n_rf_safe": sum(f.rf_safe for f in self.features.values()),
                "unsafe": self.unsafe_columns()}
