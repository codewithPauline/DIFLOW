# Real-data input layer

DIFLOW's first production-facing data layer converts ordinary biallelic VCF
genotypes into population allele counts and then into a pairwise joint
site-frequency spectrum (jSFS).

## Population map

The population map is a delimited text file with two required columns:

sample    population

Sample identifiers must exactly match the VCF sample names.

## VCF scope in the first implementation

The initial parser intentionally accepts a conservative subset:

- plain-text VCF
- biallelic records
- GT genotype field
- diploid or other explicit allele counts such as haploid GT
- missing genotypes are ignored
- multiallelic sites are skipped

This reader is designed for correctness and testing, not yet for very large
compressed datasets. Indexed .vcf.gz / BCF support will use a production
backend later.

## Allele-count table

Each emitted row contains:

- chromosome
- position
- REF allele
- ALT allele
- population
- REF count
- ALT count
- number of called chromosomes

## Pairwise jSFS

For populations A and B, DIFLOW creates a matrix S where

S[i, j]

is the number of loci with i alternate alleles in A and j alternate alleles
in B.

The first implementation only combines loci with a fixed called chromosome
count in both populations. This is intentional: mixing different chromosome
sample sizes changes the support of the SFS and can bias demographic
inference.

A hypergeometric projection layer will be added next so loci with unequal
sample sizes can be projected to common sample sizes instead of discarded.
