
## Enumerating and prioritising inherited variants

`inherited_report.py` applies clinical inheritance models and returns only variants meeting reportable
thresholds. When that returns nothing, the aggregate burden tests alone do not tell a reader what was present.
These five scripts enumerate the inherited screen instead, across all three variant classes, and rank what they
find.

| script | what it does |
|---|---|
| `inherited_lof_table.py` | rare loss-of-function small variants transmitted to affected offspring, in SFARI or LoF-constrained (LOEUF < 0.35) genes, filtered on population and cohort frequency (DDG2P enters through `inherited_report.py`, not here) |
| `inherited_sv_table.py` | the same for structural variants overlapping coding sequence; breakends are counted but not interpreted |
| `inherited_tr_table.py` | the same for tandem-repeat length outliers, with rarity defined against the cohort because no population reference exists |
| `prioritised_consequence.py` | assigns a predicted functional consequence per variant, computed rather than inherited from the selection step: coding fraction covered for deletions and duplications, frame for insertions, coding overlap for repeats |
| `inherited_prioritise.py` | three-tier ranking: contributory to NDD/autism, incidental clinical finding, or constrained-gene loss of function with no established disease link |

Three design points are load-bearing.

**Cohort frequency does most of the filtering.** Recurrent artifacts, not recurrent biology, are the dominant
failure mode in all three call sets. A paralogous small-variant locus, one very large inversion, and
homopolymer repeat loci each generated a large share of rows before this filter.

**Quality is judged before the tier.** A tier assignment resting on an unreliable genotype is worse than none,
so each class carries its own gate: depth and genotype quality for small variants, genotype quality for
structural variants, spanning-read count and motif complexity for repeats.

**Transmission from an unaffected parent is not evidence against a contributory role.** For a highly penetrant
Mendelian syndrome it would be. For autism risk alleles, incomplete penetrance and the lower rate of expression
in female carriers make inherited variants in unaffected parents the expected case. It argues only against a
fully penetrant syndromic diagnosis. Co-segregation with an affected parent, where a family has one, is scored
separately and is the strongest single piece of evidence the design can supply.
