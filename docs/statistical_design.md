# DIFLOW statistical design notes

## Target quantity

DIFLOW uses forward-time source-to-recipient migration notation.

For populations i and j:

- m(i -> j): proportion of population j replaced per generation by migrants from i.
- m(j -> i): reverse-direction migration parameter.

These parameters are not assumed to be equal.

## First validation problem

The first estimator will be developed for two populations under controlled
simulation. The initial validation grid will vary:

- m(A -> B)
- m(B -> A)
- N_e(A)
- N_e(B)
- number of loci
- divergence in initial allele frequencies
- number of generations

The estimator must recover both rates without being told which direction is
stronger.

## Primary validation metrics

For each migration parameter:

1. Bias
2. RMSE
3. Confidence-interval coverage
4. Direction classification accuracy
5. False directional-positive rate under symmetric migration

## Asymmetry statistic

For a pair of migration rates,

A_ij = (m_ij - m_ji) / (m_ij + m_ji)

when at least one rate is positive.

A_ij near zero indicates approximately symmetric migration.
A_ij near +1 indicates predominantly i -> j migration.
A_ij near -1 indicates predominantly j -> i migration.

This statistic is descriptive. It is not itself an estimator of migration and
must not be interpreted independently of uncertainty in m_ij and m_ji.

## Scientific guardrail

DIFLOW will not infer organismal movement direction from a generic genetic
distance or ancestry coefficient. Directional arrows on maps must be tied to
an explicit migration parameter estimated under a stated demographic model.
