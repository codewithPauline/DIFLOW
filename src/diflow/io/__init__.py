"""Input/output helpers for population-genomic data."""

from .popmap import read_popmap
from .vcf import allele_counts_from_vcf

__all__ = ["read_popmap", "allele_counts_from_vcf"]
