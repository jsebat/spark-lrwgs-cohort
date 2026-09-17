"""trio_phase — trio-resolved and physical long-read phasing for a family.

Module 02_phasing of spark-lrwgs-cohort. Four steps, in order:
  orient        child phase blocks -> parent of origin, by Mendelian vote
  transmission  per-parent transmitted / untransmitted map + within-block change points
  xo-reads      change points -> CROSSOVER / SWITCH_ERROR / AMBIGUOUS, from the parent's reads
  hapdepth      per-haplotype depth in fixed bins from one haplotagged BAM
  phase-qc      per-child phase_qc.json and the cohort QC table with gates

Everything here was `phase_dnm.phasing` (the "M1" stage of 05_denovo) until it was extracted;
the algorithms, thresholds and output schemas are unchanged. Downstream modules consume the
tables BY PATH under $PHASE_DIR and never import this package.

Every entry point takes family / sample identifiers as arguments and reads paths from the
environment or explicit flags; nothing here embeds an identifier.
"""
__version__ = "0.1.0"
