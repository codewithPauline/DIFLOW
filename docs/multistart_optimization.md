# Multi-start optimization and convergence diagnostics

Demographic likelihood surfaces can contain local optima and flat ridges.
DIFLOW therefore does not treat a single optimization result as sufficient.

## Multi-start strategy

For each demographic model, DIFLOW:

1. runs the canonical starting point,
2. draws additional starting points log-uniformly within parameter bounds,
3. optimizes each start independently,
4. records every success and failure,
5. ranks successful runs by log-likelihood,
6. evaluates whether near-best solutions agree on parameter values.

## Current stability heuristic

Among near-best runs, DIFLOW evaluates:

- log-likelihood proximity to the best solution,
- relative spread of each fitted parameter,
- fraction of successful optimization starts.

A model is provisionally labeled stable only when multiple near-best
optimizations produce similar parameter values.

These thresholds are development-phase heuristics. They will be calibrated by
simulation before a validated release.

## Why this matters for directional gene flow

A result such as

m_A_to_B = 3.1
m_B_to_A = 0.4

is not biologically trustworthy if equally good optimization runs also find

m_A_to_B = 0.5
m_B_to_A = 4.2.

DIFLOW should report that situation as unstable rather than drawing a strong
directional arrow.

## Reporting

The multistart layer stores:

- run number,
- optimization success,
- log-likelihood,
- fitted parameters,
- best parameter set,
- convergence fraction,
- stability flag.

Future releases will add profile likelihoods and likelihood-surface plots for
problematic parameters.
