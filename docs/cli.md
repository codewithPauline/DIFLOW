# Command-line reference

DIFLOW currently provides twelve main commands:

- diflow inspect
- diflow infer
- diflow benchmark
- diflow calibrate
- diflow profile
- diflow slurm
- diflow report
- diflow compare
- diflow campaign
- diflow benchmark-export
- diflow release-review
- diflow data-requirements

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

Important options include projection-chromosomes, neighbors, max-distance-km, preset, starts, maxiter, bootstrap-replicates, bootstrap-starts, bootstrap-block-bp, min-model-weight, min-directional-support, min-abs-asymmetry, map-crs, polarized, seed, and prepare-only.

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

## Benchmark suites

DIFLOW provides eight benchmark modes:

    diflow benchmark --suite recovery --output recovery/
    diflow benchmark --suite stress --output stress/
    diflow benchmark --suite forward --output forward/
    diflow benchmark --suite grid --output grid/
    diflow benchmark --suite linked --output linked/
    diflow benchmark --suite mechanistic --output mechanistic/
    diflow benchmark --suite mechanistic-grid --output mechanistic_grid/
    diflow benchmark --suite decision --output decision/

`--suite all` runs recovery, stress, and forward-time stress suites.

The large recovery grid is excluded from `all` because it is substantially
more computationally expensive and must be requested explicitly.


The linked calibration suite has additional controls:

    --linked-blocks
    --snps-per-block
    --linked-block-bp
    --linkage-concentration
    --linked-bootstrap-replicates

These control the number and size of correlated marker blocks and the amount of
bootstrap effort applied to each simulated dataset.


## diflow calibrate

Calibrate final directional thresholds from a known-truth evidence table:

    diflow calibrate \
      --evidence decision/decision_evidence.csv \
      --output calibration/ \
      --max-fpr 0.05

## diflow profile

Profile one scaled directional migration parameter from a saved spectrum:

    diflow profile \
      --spectrum results/spectra/A__B.npy \
      --parameter m_a_to_b \
      --output profile/

## diflow slurm

Generate a Slurm script for an existing DIFLOW command:

    diflow slurm \
      --run-command "diflow infer ..." \
      --script diflow_run.slurm \
      --cpus 8 \
      --mem-gb 64 \
      --hours 72


## diflow report

Create a self-contained searchable HTML report:

    diflow report --results results/

## diflow compare

Compare standardized simulation outputs from multiple methods:

    diflow compare \
      --method DIFLOW=diflow.csv \
      --method OTHER=other.csv \
      --output comparison/


## diflow campaign

Generate a reproducible multi-job Slurm campaign for release-grade validation:

    diflow campaign \
      --output validation_campaign/ \
      --replicates 50 \
      --bootstrap-replicates 100 \
      --cpus 8 \
      --mem-gb 64 \
      --hours 72


## Calibrated thresholds

Apply a calibration result directly to inference:

    diflow infer \
      --vcf data.vcf.gz \
      --popmap populations.tsv \
      --coords coordinates.csv \
      --projection-chromosomes 8 \
      --thresholds-file calibration/selected_thresholds.csv \
      --output results/

Explicit threshold flags override values loaded from the file.

## Comparator export

Convert DIFLOW recovery output into the standardized method-comparison schema:

    diflow benchmark-export \
      --input recovery_grid/recovery_grid_replicates.csv \
      --output diflow_standardized.csv

## Release review

Evaluate completed empirical validation outputs against explicit release
criteria:

    diflow release-review --results validation_campaign/results/


## Mechanistic linkage grid

The release-grade linkage study spans multiple recombination and mutation-rate
regimes:

    diflow benchmark \
      --suite mechanistic-grid \
      --output mechanistic_grid/ \
      --mechanistic-grid-workers 8

The worker count parallelizes independent grid cells using isolated processes.


## Empirical data requirements

After the recovery grid completes, derive the smallest **tested** data regimes
that meet the chosen accuracy, false-direction, and fit-success targets:

    diflow data-requirements \
      --grid-summary recovery_grid/recovery_grid_summary.csv \
      --output data_requirements/

Outputs include regime-level performance, Pareto-minimum passing regimes, and
machine-readable/Markdown guidance. These recommendations apply only to the
validated simulation space and are not universal biological minimums.
