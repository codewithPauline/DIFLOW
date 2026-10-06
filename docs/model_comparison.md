# Demographic model comparison

DIFLOW does not treat the best-fitting asymmetric migration model as evidence
of directional gene flow by default.

Instead, the current comparison layer fits four explicit alternatives:

1. isolation
2. continuous symmetric migration
3. continuous asymmetric migration
4. asymmetric secondary contact

The underlying demographic functions are supplied by dadi.

## Information criteria

For each model,

AIC = 2k - 2 ln(L)

where k is the number of free parameters and L is the maximized composite
likelihood.

When requested, DIFLOW also calculates

AICc = AIC + [2k(k+1)] / [n-k-1]

where n is currently defined as the number of non-zero cells in the observed
jSFS. This choice is explicit and may be revisited during validation because
the jSFS composite likelihood does not consist of fully independent
observations.

DIFLOW reports:

- model name
- maximized log-likelihood
- parameter count
- AIC or AICc
- delta information criterion
- Akaike weight

## Interpretation rule

A directional biological statement should not be based only on the fitted
values m_A_to_B and m_B_to_A.

At minimum, DIFLOW should require:

- the asymmetric model to outperform isolation and symmetric migration,
- stable parameter estimates across optimization starts,
- uncertainty that supports the directional contrast,
- simulation validation under the relevant demographic regime.

Secondary contact is included because apparent asymmetry can otherwise absorb
historical model misspecification.

## Direction convention

DIFLOW continues to report migration in forward-time source-to-recipient form.
When calling dadi models, DIFLOW reverses the order required by dadi's m12/m21
labels internally and tests this mapping explicitly.
