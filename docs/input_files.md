# Input files

DIFLOW currently uses three core biological input files.

## VCF

Example: data.vcf

The current development parser accepts plain-text VCF and gzip-compressed `.vcf.gz` files containing biallelic records with a GT genotype field. Multiallelic records are skipped by the conservative parser.

Current requirements:

- samples are columns in the VCF
- genotypes contain allele indices 0 and 1
- a GT field is present
- every sample listed in the population map occurs in the VCF

`.vcf.gz` is read transparently with streaming gzip decompression. BCF and indexed random-access readers remain future production-reader work.

## Population map

The population map connects each VCF sample to one population. It may be tab- or comma-separated.

Example:

    sample    population
    IND001    POP_A
    IND002    POP_A
    IND003    POP_B
    IND004    POP_B

Required column names are sample and population. Each sample may appear only once.

Population labels should represent the biological units to be compared. They may be localities, demes, sampling populations, or another defensible population definition.

## Coordinates

Example:

    population,latitude,longitude
    POP_A,39.510,-84.730
    POP_B,39.850,-84.120
    POP_C,38.990,-83.570

Required columns are population, latitude, and longitude. Each population identifier must be unique.

Every population present in the population map must have coordinates.

Coordinates define and visualize geographic candidate networks. They do not by themselves provide evidence for gene flow.

## Matching identifiers

Population names must match exactly between the population map and coordinate file. Sample names must match exactly between the population map and VCF header.

## Minimal example directory

    project/
    ├── data.vcf
    ├── populations.tsv
    └── coordinates.csv

Then inspect the project with:

    diflow inspect \
      --vcf project/data.vcf \
      --popmap project/populations.tsv \
      --coords project/coordinates.csv

## Recommended preprocessing principles

Before using DIFLOW, apply population-genomic quality control appropriate to the study design, including sample/locus quality, missingness, relatedness where relevant, linkage considerations, and defensible population definitions.

DIFLOW should not be used as a substitute for upstream genomic QC.

## Polarization

By default, DIFLOW treats ordinary REF/ALT coding as unpolarized and folds the
observed jSFS during inference.

Use `--polarized` only when upstream processing has established that ALT is
the derived allele at every included site. DIFLOW does not infer ancestral state
from VCF REF/ALT labels.
