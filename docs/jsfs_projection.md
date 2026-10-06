# jSFS projection

Real population-genomic datasets frequently have missing genotypes, so the
number of called chromosomes varies among loci. A site-frequency spectrum
cannot simply combine allele counts observed at different sample sizes.

DIFLOW therefore supports hypergeometric projection.

For a locus with n called chromosomes, x observed alternate alleles, and a
target of k chromosomes, the probability of observing y alternate alleles
after projection is

P(Y = y) = Hypergeometric(N=n, K=x, n=k).

Each locus contributes probability mass across projected frequency bins rather
than being forced into one bin. For a pair of populations, the two projection
vectors are combined with an outer product to contribute to the 2D jSFS.

This allows loci with more called chromosomes than the chosen target to be
retained while preserving a common spectrum size.

Loci with fewer called chromosomes than the projection target are still
excluded.
