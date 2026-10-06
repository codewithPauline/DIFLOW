# Real-data case study protocol

A DIFLOW real-data case study demonstrates workflow, reporting, and biological
interpretation. It does **not** prove estimator accuracy because the true
migration history of an empirical system is unknown.

## Purpose

A release-grade case study should show that a researcher can move from genomic
inputs to an auditable directional gene-flow result while preserving the
distinction between:

- observed genetic structure
- demographic model fit
- directional migration estimates
- uncertainty
- final evidence classification

## Required inputs

Document:

- study system and sampling design
- biological definition of each population
- number of samples per population
- marker-generation method
- genomic QC and filtering
- missing-data thresholds
- linkage treatment
- ancestral-polarization status
- coordinate source and geographic candidate-graph rationale

Do not choose population definitions or geographic edges after viewing the
desired directional result.

## Pre-inference inspection

Run `diflow inspect` and archive:

- sample/population counts
- pairwise shared-locus projection table
- selected projection
- candidate graph options
- final graph choice and rationale

The inspection recommendation is a transparent starting point, not a biological
truth.

## Linkage strategy

State explicitly whether the analysis uses:

- one approximately independent SNP per locus,
- locus bootstrap, or
- genomic block bootstrap.

For block bootstrap, justify the chosen block size and include sensitivity to
at least one smaller and one larger plausible block size when practical.

## Direction thresholds

A release-grade case study should use the empirically calibrated threshold file:

    diflow infer \
      ... \
      --thresholds-file selected_thresholds.csv

Do not replace calibrated thresholds with hand-picked values because they make
the empirical arrows look cleaner.

## Primary inference

Archive:

- resolved configuration
- SHA-256 input provenance
- pairwise jSFS files
- model rankings
- optimizer diagnostics
- migration estimates in both directions
- bootstrap intervals and directional support
- final supported / ambiguous / unsupported classifications
- directional network and map

## Identifiability checks

For biologically important supported edges, run profile likelihood for both
migration directions.

A strong point estimate with a broad or unbounded profile should be described as
weakly identified rather than as a precise directional estimate.

## Sensitivity analyses

At minimum, evaluate whether major conclusions are robust to reasonable changes
in:

- projection size
- candidate geographic graph
- bootstrap block size or SNP-independence treatment
- decision thresholds within their empirically supported range

Sensitivity analysis should test robustness, not search for a configuration
that maximizes support.

## Interpretation

Report model-based migration using DIFLOW's source-to-recipient convention.

Do not translate scaled demographic migration into literal migrants per
generation unless the required reference population-size and scaling
assumptions are available and stated.

Distinguish clearly between:

1. differentiation,
2. descriptive population structure,
3. estimated migration parameters,
4. evidence for directional asymmetry.

## Negative and ambiguous results

The case study must retain ambiguous and unsupported comparisons in the
reported tables.

A method designed to control false directional claims should be allowed to say
that direction is unresolved.

## External biological context

Independent natural-history, geographic, ecological, or previous genomic
evidence can be used to interpret a DIFLOW result, but it should not be treated
as a substitute for the model evidence.

Unexpected directional results deserve additional diagnostics rather than
automatic rejection.

## Reproducibility package

Archive:

- exact DIFLOW version/commit
- commands or resolved configuration
- threshold file
- input-file checksums
- output tables
- maps/figures
- profile-likelihood outputs
- sensitivity-analysis outputs
- a short interpretation README

Where raw genomic data cannot be redistributed, provide accession identifiers
and all derived non-sensitive inputs needed to reproduce the analysis.

## What the case study may claim

Appropriate:

> Under the specified demographic models and calibrated evidence criteria,
> population A showed stronger supported migration into population B than the
> reverse direction.

Not appropriate solely from the case study:

> DIFLOW is 95% accurate on real populations.

Accuracy and false-positive control come from known-truth simulation
validation, not from empirical plausibility.

## Release requirement

The case-study blocker is complete only when at least one real dataset has been
run through the release-candidate workflow, sensitivity checks have been
reviewed, and the reproducibility package has been archived.
