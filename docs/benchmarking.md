# DIFLOW validation and benchmarking strategy

DIFLOW is not novel because it uses the joint site-frequency spectrum. The
jSFS is an established population-genetic summary and is intentionally used as
a tested statistical substrate.

The contribution DIFLOW is aiming for is a reproducible evidence architecture
for directional gene-flow inference.

## Validation question

The central question is not:

"Can the optimizer return two different migration parameters?"

It is:

"Under which demographic conditions can DIFLOW correctly recover migration
direction and magnitude, quantify uncertainty, and avoid claiming direction
when the data do not support it?"

## Minimum validation suite

The validated release should include simulations for:

- symmetric migration
- A -> B asymmetric migration
- B -> A asymmetric migration
- near-unidirectional migration
- unequal effective population sizes
- bottlenecks and growth
- secondary contact
- ancient migration
- range expansion
- ghost populations
- uneven sampling and missing data
- linked loci / reduced-representation genomic data
- zero migration

## Primary metrics

For migration magnitude:

- bias
- RMSE
- confidence-interval coverage

For direction:

- direction classification accuracy
- false directional-positive rate under symmetric migration
- sensitivity under weak, moderate, and strong asymmetry
- calibration of directional support probabilities

For model selection:

- correct-model selection frequency
- Akaike-weight calibration
- frequency of confusing continuous migration with secondary contact

For computation:

- runtime
- memory use
- scaling with loci, sample size, candidate pairs, optimization starts, and
  bootstrap replicates

## Comparator philosophy

DIFLOW should eventually be benchmarked against established approaches rather
than evaluated only against itself.

Relevant comparator classes include:

- demographic SFS inference packages
- historical migration estimators
- contemporary migration estimators
- fine-resolution asymmetric spatial migration methods
- descriptive directional statistics

Comparisons must respect differences in what each method actually estimates.
Backward-time lineage migration, recent migrant assignment, historical
coalescent migration, and forward-time demographic migration are not
interchangeable quantities.

## Release criterion

A map that looks plausible is not validation.

Before a validated release, DIFLOW should demonstrate:

1. low false directional-positive rates under symmetry,
2. reliable recovery of both directions under known asymmetry,
3. calibrated uncertainty,
4. robustness limits under major demographic confounders,
5. transparent failure modes,
6. reproducible benchmark datasets and scripts.

Where those conditions are not met, DIFLOW should report uncertainty or
ambiguity rather than force a directional conclusion.
