# DIFLOW statistical design notes

## Statistical substrate

DIFLOW uses the joint site-frequency spectrum (jSFS) as an established summary
of allele-frequency variation across populations.

The jSFS is not claimed as a methodological novelty. DIFLOW's research
contribution is intended to come from the evidence framework built around
demographic inference: competing models, stability checks, uncertainty,
direction classification, simulation calibration, and spatial reporting.

## Target quantities

DIFLOW uses forward-time source-to-recipient migration notation.

For populations i and j:

- m(i -> j): migration from i into j.
- m(j -> i): migration from j into i.

These parameters are estimated separately.

Backend-specific scaled migration parameters are not automatically equivalent
to literal numbers of organisms moving per generation. Conversion and
interpretation depend on the demographic model, reference population size,
mutation scaling, and backend convention.

## Differentiation versus migration

Genetic differentiation statistics such as FST may be useful descriptive
outputs, but they are not directional migration estimators in DIFLOW.

DIFLOW will not use generic transforms such as Nm ≈ (1-FST)/(4FST) to infer
migration direction. Such equilibrium approximations do not supply a robust
source-to-recipient direction and can be badly violated by realistic
demographic histories.

## Model comparison

The current two-population model set includes:

- isolation
- symmetric continuous migration
- asymmetric continuous migration
- asymmetric secondary contact

Additional histories should be introduced when they address identifiable
sources of model misspecification rather than simply increasing model count.

## Optimization

Every demographic model should be fit from multiple starting points.

A directional result should be treated cautiously when similarly good
likelihoods correspond to substantially different migration parameters.

## Directional asymmetry

For a pair of positive or non-zero migration rates,

A_ij = (m_ij - m_ji) / (m_ij + m_ji).

A_ij near zero indicates approximate symmetry.
A_ij near +1 indicates stronger i -> j migration.
A_ij near -1 indicates stronger j -> i migration.

This statistic is descriptive and does not replace uncertainty estimation.

## Uncertainty

The current real-data uncertainty layer can either resample usable loci
independently or resample fixed genomic windows as blocks, reconstruct the
projected jSFS, and refit the asymmetric model.

For A and B, directional support includes:

P(m_A_to_B > m_B_to_A).

When B -> A is the preferred direction, support is evaluated using the reverse
probability rather than incorrectly reusing the A -> B probability.

Linked markers violate the simplest locus-independence assumption. DIFLOW now
supports fixed genomic-window block bootstrap. Simulation calibration of
interval coverage under realistic linkage remains required.

## Evidence classification

A direction should not be mapped as strongly supported merely because one
point estimate is larger than the reverse estimate.

The evidence layer may incorporate:

- asymmetric-model weight
- optimizer stability
- absolute asymmetry
- bootstrap directional support
- uncertainty-interval separation

Thresholds are development defaults until calibrated by simulation.

## Validation targets

For migration magnitude:

1. bias
2. RMSE
3. interval coverage

For migration direction:

1. direction classification accuracy
2. false directional-positive rate under symmetry
3. sensitivity across asymmetry strength
4. support-probability calibration

For model discrimination:

1. correct-model selection frequency
2. confusion between continuous migration and secondary contact
3. robustness to ghost populations and ancestral structure

## Required stress tests

The minimum validation suite includes:

- symmetric migration
- A -> B asymmetry
- B -> A asymmetry
- near-unidirectional migration
- zero migration
- unequal effective population size
- bottlenecks and growth
- secondary contact
- ancient migration
- range expansion
- ghost populations
- missing populations
- uneven sampling
- missing genotypes
- linked loci

## Scientific guardrail

DIFLOW will not infer organismal movement direction from a generic genetic
distance, FST value, ancestry coefficient, or visually suggestive geographic
pattern.

Directional arrows must be traceable to an explicit migration parameter,
model-comparison evidence, and quantified uncertainty.
