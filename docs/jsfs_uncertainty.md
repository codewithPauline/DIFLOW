# jSFS bootstrap uncertainty

DIFLOW estimates uncertainty for single-time-point asymmetric migration by
resampling loci with replacement.

## Procedure

For each population pair:

1. identify loci usable at the requested projection size,
2. calculate each locus's hypergeometric contribution to the projected jSFS,
3. sample loci with replacement,
4. sum sampled contributions into a bootstrap jSFS,
5. refit the asymmetric demographic model,
6. repeat across bootstrap replicates.

The resulting distribution provides:

- mean m_A_to_B
- confidence interval for m_A_to_B
- mean m_B_to_A
- confidence interval for m_B_to_A
- P(m_A_to_B > m_B_to_A)
- support for whichever direction has the larger bootstrap mean

## Directional support

If A -> B is the preferred direction,

support = P(m_A_to_B > m_B_to_A).

If B -> A is preferred,

support = 1 - P(m_A_to_B > m_B_to_A).

This prevents reverse-direction edges from being penalized simply because the
stored probability is defined in A-to-B order.

## Independence assumption

The current bootstrap samples loci as independent units. For datasets where
multiple SNPs are linked within loci or genomic blocks, users should supply
one approximately independent SNP per locus or wait for DIFLOW's planned block
bootstrap layer.

For reduced-representation datasets, locus/block resampling is preferred over
naive SNP bootstrap when linkage within loci is present.

## Interpretation

Bootstrap uncertainty addresses sampling variability under the fitted model.
It does not eliminate bias caused by an incorrect demographic model, ghost
populations, unsampled migration routes, ancestral structure, or severe
identifiability problems.
