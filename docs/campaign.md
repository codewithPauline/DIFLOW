# Validation campaign

DIFLOW can generate a reproducible Slurm campaign for the expensive empirical
studies required before a validated scientific release.

Generate the campaign:

    diflow campaign \
      --output validation_campaign/ \
      --replicates 50 \
      --starts 10 \
      --bootstrap-replicates 100 \
      --cpus 8 \
      --mem-gb 64 \
      --hours 72 \
      --seed 42

The output directory contains:

- `01_recovery_grid.slurm`
- `02_decision_evidence.slurm`
- `03_linked_calibration.slurm`
- `04_mechanistic_linkage.slurm`
- `campaign_manifest.json`
- `submit_all.sh`
- a `results/` directory for benchmark outputs

Submit all jobs with:

    bash validation_campaign/submit_all.sh

## Why the campaign is separated into jobs

Recovery, threshold calibration, correlated-marker calibration, and mechanistic
recombination validation have different failure modes and computational costs.
Keeping them separate allows a failed or under-resourced study to be rerun
without repeating all other benchmarks.

## After the jobs finish

The manifest includes a post-processing threshold-calibration command, an
automated `release-review` command, and a release-review checklist.

Completing the jobs is not enough by itself. The resulting error rates,
confidence-interval coverage, false directional-positive rates, and sensitivity
must meet predefined scientific acceptance criteria before default thresholds
are frozen.


## Release review

After the benchmark jobs, threshold calibration, and external-method comparison
are complete, run:

    diflow release-review \
      --results validation_campaign/results/

The review writes CSV, JSON, and Markdown summaries and treats missing empirical
studies as release blockers rather than silently ignoring them.


## Archive the final evidence

After external comparison and a passing release review, freeze the completed
validation evidence:

    diflow archive-validation \
      --results validation_campaign/results/

This writes JSON and SHA-256 manifests for every file in the validation results
tree. The campaign manifest includes this archive command for reproducibility.
