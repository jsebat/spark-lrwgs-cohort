"""phase_dnm — phase-aware de novo mutation calling for long-read trios.

Module 13 of spark-lrwgs-cohort. Sub-modules build in order:
  hapmatrix (M2)  six-haplotype evidence matrix + class adapters
  classify  (M2)  rule layer and likelihood layer
  features        registry with rf_safe enforcement
  train     (M4)  swap-closed family-grouped nested CV

The phasing stage this module was built on top of (block orientation, transmission map, resolved
crossovers, per-haplotype depth, phase QC) is now module 01_phasing, package `trio_phase`. Its
tables are read BY PATH from $PHASE_DIR; nothing here imports it.

Every entry point takes family / sample identifiers as arguments and reads paths from the
environment or explicit flags; nothing here embeds an identifier.
"""
__version__ = "0.1.0"
