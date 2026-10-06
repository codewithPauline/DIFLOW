# Command-line reference

DIFLOW currently provides two main commands:

- diflow inspect
- diflow infer

## diflow inspect

Use inspect before expensive model fitting.

    diflow inspect \
      --vcf data.vcf \
      --popmap populations.tsv \
      --coords coordinates.csv

Optional settings include:

    --preset quick|standard|publication
    --retention-target 0.80
    --output inspection_results/

Inspection reports dataset size, projection-retention tradeoffs, candidate geographic graph options, a recommended starting configuration, and a ready-to-copy inference command.

See [inspection_workflow.md](inspection_workflow.md).

## diflow infer

Example:

    diflow infer \
      --vcf data.vcf \
      --popmap populations.tsv \
      --coords coordinates.csv \
      --projection-chromosomes 8 \
      --neighbors 4 \
      --preset standard \
      --output results/

Important options include projection-chromosomes, neighbors, max-distance-km, preset, starts, maxiter, bootstrap-replicates, bootstrap-starts, bootstrap-block-bp, seed, and prepare-only.

Explicit optimization and bootstrap flags override preset values.

## Presets

| preset | model starts | bootstrap replicates | bootstrap starts | max iterations |
| --- | ---: | ---: | ---: | ---: |
| quick | 5 | 20 | 2 | 60 |
| standard | 20 | 100 | 5 | 100 |
| publication | 40 | 500 | 8 | 200 |

The presets control computational effort. They do not replace model diagnostics, uncertainty evaluation, or scientific validation.

## Prepare-only mode

    diflow infer \
      --vcf data.vcf \
      --popmap populations.tsv \
      --coords coordinates.csv \
      --projection-chromosomes 8 \
      --neighbors 4 \
      --output results/ \
      --prepare-only

## Bootstrap-enabled inference

When bootstrap-replicates is at least 2, DIFLOW attaches locus-bootstrap uncertainty and passes the resulting evidence into the formal directional classifier.

Without bootstrap, directional labels remain provisional.

For linked SNPs, use a fixed genomic-window block bootstrap:

    diflow infer \
      --vcf data.vcf \
      --popmap populations.tsv \
      --coords coordinates.csv \
      --projection-chromosomes 8 \
      --neighbors 4 \
      --bootstrap-replicates 100 \
      --bootstrap-block-bp 100000 \
      --output results/

When block size is omitted, DIFLOW resamples loci independently.

## Reproducibility

The CLI currently defaults to random seed 42. Users may override it with --seed.

## Outputs

See [outputs.md](outputs.md) for the full output reference.