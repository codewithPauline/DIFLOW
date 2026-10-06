# Asymmetric demographic inference backend

DIFLOW uses its own forward-time source-to-recipient migration naming while
delegating numerical diffusion calculations to dadi.

## DIFLOW convention

For populations A and B:

- m_A_to_B means migration from A into B
- m_B_to_A means migration from B into A

The first single-time-point model has five scaled parameters:

- nu_A
- nu_B
- T_split
- m_A_to_B
- m_B_to_A

## dadi translation

dadi's two-population asymmetric split model uses the parameter order

(nu1, nu2, T, m12, m21)

where:

- m12 = population 2 -> population 1
- m21 = population 1 -> population 2

Therefore DIFLOW must translate

m_A_to_B -> dadi m21
m_B_to_A -> dadi m12

This translation is centralized and unit-tested because reversing it would
reverse the biological direction reported by DIFLOW.

## Likelihood

The observed pairwise jSFS is compared with the model-predicted spectrum using
dadi's multinomial composite likelihood. The expected spectrum is calculated
on multiple frequency grids and extrapolated to reduce grid-discretization
error.

## Current interpretation

The migration parameters estimated by this backend are scaled demographic
migration parameters in dadi units. They are not yet automatically converted
to unscaled migrants per generation because that conversion requires an
estimate of the ancestral reference population size and mutation rate.

## Current limitations

The first model assumes continuous migration since the population split.
DIFLOW will add model comparison against:

- strict isolation
- symmetric migration
- asymmetric continuous migration
- secondary contact
- ancient migration

Directional interpretation should only be emphasized when the asymmetric
model is supported over simpler alternatives.
