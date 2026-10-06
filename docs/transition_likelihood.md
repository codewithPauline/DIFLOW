# First inference milestone: transition likelihood

DIFLOW's first estimator is intentionally narrower than the eventual method.

## Assumptions

The benchmark model assumes:

1. two populations,
2. allele frequencies observed at the start of a generation,
3. one generation of migration,
4. Wright-Fisher binomial drift,
5. known effective population sizes,
6. independent biallelic loci.

For recipient population B,

p_B' = (1 - m_A_to_B) p_B + m_A_to_B p_A

before drift. The observed post-drift allele count is modeled as

X_B ~ Binomial(2 N_e,B, p_B').

The reverse direction is estimated separately with the corresponding
A-recipient transition.

## Why start here?

This is not yet the intended real-data estimator. It is a unit-testable
statistical benchmark with known truth. It allows the project to verify that:

- forward-time source/recipient notation is implemented correctly,
- asymmetric rates can be recovered independently,
- stronger A -> B migration is not accidentally reported as B -> A,
- unequal effective population sizes are handled correctly.

Only after this benchmark passes do we move to single-time-point genomic
inference using frequency-spectrum or related demographic information.
