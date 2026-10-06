# Command-line workflow

DIFLOW now provides a first integrated command-line workflow.

## Prepare inputs and spectra only

Use this mode to validate the genomic input layer and candidate graph without
running demographic optimization:

    diflow infer \
      --vcf data.vcf \
      --popmap populations.tsv \
      --coords coordinates.csv \
      --projection-chromosomes 8 \
      --neighbors 4 \
      --output results/ \
      --prepare-only

## Run demographic inference

With the optional dadi backend installed:

    python -m pip install -e ".[demography]"

run:

    diflow infer \
      --vcf data.vcf \
      --popmap populations.tsv \
      --coords coordinates.csv \
      --projection-chromosomes 8 \
      --neighbors 4 \
      --starts 30 \
      --output results/

## Outputs

The workflow writes:

- allele_counts.csv
- candidate_pairs.csv
- pairwise_results.csv
- model_rankings.csv
- run_metadata.csv
- spectra/*.npy

## Important current limitation

The single-time-point demographic workflow does not yet have calibrated
bootstrap or profile-likelihood uncertainty for directional migration.

Therefore pairwise status values produced by the CLI are explicitly
provisional:

- candidate
- ambiguous
- unsupported

They must not yet be interpreted as final supported/unsupported biological
direction calls.

The formal supported/ambiguous/unsupported decision engine will be connected
to this workflow after jSFS uncertainty estimation is implemented.
