# Release readiness

DIFLOW now has a broad software and validation infrastructure, but a validated
1.0 scientific release requires empirical evidence in addition to working code.

## Software infrastructure implemented

- real VCF / VCF.GZ ingestion
- folded-by-default jSFS inference
- explicit polarized mode
- pairwise shared-locus projection inspection
- competing demographic models
- multi-start optimization
- AIC model comparison
- directional migration estimates in both directions
- locus and genomic-block bootstrap
- configurable evidence classifier
- profile-likelihood diagnostics
- geographic candidate graphs
- process-isolated pair-level parallel execution
- projected map output
- run configuration + SHA-256 provenance
- Slurm script generation
- searchable HTML results reports
- recovery, stress, forward-time, linked-marker and msprime validation frameworks
- decision-evidence simulation and threshold calibration
- standardized external-method comparison framework

## Empirical release blockers

These are studies that must actually be run at adequate scale before DIFLOW
should claim validated default behavior:

1. Run the full recovery grid with enough replicates for stable error estimates.
2. Run the decision-evidence benchmark across data regimes and freeze thresholds
   only after demonstrating an acceptable false-direction rate.
3. Run large mechanistic recombination/linkage simulations across recombination
   rates, marker density, and block sizes to evaluate CI coverage.
4. Execute matched benchmarks against established external methods, respecting
   differences in estimands.
5. Run at least one real-data case study as a demonstration of workflow and
   interpretation, not as proof of accuracy. Follow [case_study_protocol.md](case_study_protocol.md).
6. Derive and document recommended minimum data requirements from the completed simulation campaign. Known failure regimes are already documented in [limitations.md](limitations.md).
7. Freeze a release candidate, rerun all tests/benchmarks, archive outputs, and
   only then create a formal software release/DOI.

## Current scientific status

DIFLOW should still be described as active research software with extensive
validation infrastructure, not yet as a universally validated black-box
estimator.

That distinction is a strength: release claims should be earned by the
simulation results rather than inferred from software completeness.


## Reproducible campaign generation

Use `diflow campaign` to generate the major pre-release validation jobs,
Slurm scripts, submission helper, and machine-readable campaign manifest.

This reduces the risk that release claims depend on undocumented one-off
simulation commands.


## Automated release review

Run:

    diflow release-review \
      --results validation_campaign/results/

The default project policy targets are:

- minimum direction accuracy: 0.90
- maximum false-direction rate: 0.05
- minimum directional sensitivity: 0.80
- mechanistic-linkage interval coverage: 0.90 to 0.99
- minimum fitting success rate: 0.95
- at least one standardized external method in the comparison study

These are configurable release-policy targets, not universal biological
constants.

The review fails when a required empirical output is missing. A passing software
test suite therefore cannot be mistaken for a completed scientific validation
campaign.


## Final release sequence

A validated release should follow this order:

1. Run the release-grade validation campaign.
2. Calibrate directional thresholds from known truth.
3. Derive empirical minimum-data guidance from the recovery grid.
4. Execute and standardize the established-method benchmark.
5. Run at least one real-data case study as a workflow demonstration.
6. Run `diflow release-review` and resolve every blocker.
7. Freeze the exact validation evidence with:

       diflow archive-validation --results validation_campaign/results/

8. Rerun the complete CI/test matrix on the release candidate.
9. Tag the release, archive the code/results, and create the formal DOI only
   after the empirical evidence supports the release claims.

The validation archive writes SHA-256 manifests so the evidence used to justify
release defaults can be independently verified.
