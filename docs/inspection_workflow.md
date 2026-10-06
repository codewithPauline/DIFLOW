# User-friendly inspection workflow

DIFLOW separates data inspection from final demographic inference so that
important assumptions remain visible.

## Step 1: inspect the dataset

A user supplies only the three biological inputs:

    diflow inspect \
      --vcf data.vcf \
      --popmap populations.tsv \
      --coords coordinates.csv

DIFLOW reports:

- number of samples
- number of populations
- number of variant loci
- candidate jSFS projection sizes and retained-locus fractions
- a recommended projection chromosome count
- several geographic k-nearest-neighbor graph options
- the number and distance distribution of candidate population pairs
- an analysis-effort preset
- a ready-to-copy infer command

The projection recommendation uses a transparent retention heuristic. By
default, DIFLOW chooses the largest even projection size that retains at least
80% of loci in the worst-retained population. The user can change this target.

## Step 2: review the geographic graph

The suggested neighbor count is a computational starting point, not a
biological claim.

Users should review whether geographic adjacency is biologically sensible for
their system. Future DIFLOW versions will support watershed, habitat,
resistance, and user-defined adjacency graphs.

## Step 3: run inference

The user may copy the recommended command or use a preset:

    diflow infer \
      --vcf data.vcf \
      --popmap populations.tsv \
      --coords coordinates.csv \
      --projection-chromosomes 8 \
      --neighbors 4 \
      --preset standard \
      --output results/

Available presets:

| preset | model starts | bootstrap replicates | bootstrap starts | max iterations |
| --- | ---: | ---: | ---: | ---: |
| quick | 5 | 20 | 2 | 60 |
| standard | 20 | 100 | 5 | 100 |
| publication | 40 | 500 | 8 | 200 |

Explicit command-line settings override preset values.

The random seed defaults to 42 so repeated runs begin reproducibly unless the
user chooses another seed.

## Why inspect is separate from infer

DIFLOW should not silently choose biologically important graph assumptions.

The inspection command makes recommendations visible before computation begins,
while still allowing experienced users to override every setting.


Projection recommendations are based on **shared usable loci for population
pairs**, because pairwise jSFS inference requires both populations to contribute
data at the same sites. This is more informative than evaluating populations
independently.
