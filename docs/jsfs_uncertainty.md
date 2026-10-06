# jSFS bootstrap uncertainty

DIFLOW supports two uncertainty resampling strategies for single-time-point
asymmetric migration inference:

- locus bootstrap for approximately independent markers
- fixed genomic-window block bootstrap for linked markers

## Locus bootstrap

By default, usable loci are sampled independently with replacement.

For each population pair:

1. identify loci usable at the requested projection size,
2. calculate each locus's hypergeometric contribution to the projected jSFS,
3. sample loci with replacement,
4. sum sampled contributions into a bootstrap jSFS,
5. refit the asymmetric demographic model,
6. repeat across bootstrap replicates.

## Genomic block bootstrap

For linked SNPs, users may instead define a fixed genomic window size:

    --bootstrap-block-bp 100000

All usable variants on the same chromosome and within the same fixed window are
summed into one bootstrap unit. Whole blocks are then sampled with replacement.

This preserves local linkage within a resampled block better than naive
SNP-by-SNP resampling.

DIFLOW records:

- the resampling unit
- number of blocks used
- number of loci represented
- successful and attempted bootstrap replicates

A block analysis requires at least two usable blocks.

## Choosing a resampling strategy

Use ordinary locus bootstrap when markers are approximately independent, such
as a dataset intentionally filtered to one SNP per independent locus.

Use genomic block bootstrap when nearby SNPs may be linked and genomic
coordinates are meaningful.

The block size is a biological/statistical choice, not a universal constant.
Users should choose it based on linkage scale, marker design, recombination,
and genome structure where possible.

For reduced-representation data, a one-SNP-per-locus filter remains a
defensible option when physical block definitions are unavailable.

## Directional support

The resulting bootstrap distribution provides:

- mean m_A_to_B
- confidence interval for m_A_to_B
- mean m_B_to_A
- confidence interval for m_B_to_A
- P(m_A_to_B > m_B_to_A)
- support for whichever direction has the larger bootstrap mean

If A -> B is preferred,

support = P(m_A_to_B > m_B_to_A).

If B -> A is preferred,

support = 1 - P(m_A_to_B > m_B_to_A).

## Interpretation

Bootstrap uncertainty addresses sampling variability under the fitted model.

It does not eliminate bias caused by demographic misspecification, ghost
populations, ancestral structure, incorrect population definitions, or severe
identifiability problems.

The block-bootstrap implementation is now available, but simulation-based
coverage calibration under realistic linkage remains part of DIFLOW's
validation roadmap.
