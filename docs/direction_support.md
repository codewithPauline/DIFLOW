# Direction-support decision engine

DIFLOW separates parameter estimation from biological interpretation.

A migration estimate is not automatically shown as a supported directional
edge. The decision engine combines several lines of evidence.

## Inputs

For each population pair, DIFLOW can use:

- forward migration estimate
- reverse migration estimate
- directional support from bootstrap or resampling
- Akaike weight of the asymmetric demographic model
- optimizer stability
- optional uncertainty intervals

## Default classification

A direction is provisionally labeled supported when all of the following hold:

- asymmetric-model Akaike weight >= 0.70
- directional support >= 0.95
- absolute asymmetry index >= 0.25
- optimizer is stable
- uncertainty intervals separate when interval separation is required

These are initial transparent thresholds, not universal biological constants.
They must be calibrated during simulation validation.

## Status values

supported
: all configured directional evidence criteria are satisfied.

ambiguous
: there is some directional signal, but one or more evidence criteria fail.

unsupported
: there is no migration signal or the asymmetric model receives insufficient
  support.

## Mapping behavior

By default, only supported edges are passed to the geographic map.

Ambiguous edges can optionally be retained for diagnostic figures, but they
should be visually distinguished from supported routes.

This design prevents a large fitted migration value from becoming a strong map
arrow when the model itself is weakly supported or optimization is unstable.
