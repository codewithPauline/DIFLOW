# Reproducibility and provenance

Every `diflow infer` run writes machine-readable reproducibility records.

## resolved_config.json

Contains the resolved settings actually used after CLI defaults and presets are
applied, including projection, graph settings, optimization effort, bootstrap
settings, directional thresholds, polarization mode, map CRS, and random seed.

## run_provenance.json

Contains:

- UTC creation timestamp
- DIFLOW and Python version
- versions of core numerical libraries
- dadi, msprime, and tskit versions when installed
- operating-system/platform information
- input paths
- SHA-256 hash of each biological input
- input file sizes
- resolved run settings

The hashes allow a later analysis to verify that the exact VCF, population map,
and coordinates file are the same files used originally.

## run_metadata.csv

A compact tabular summary remains available for quick inspection.

For a manuscript or archived analysis, keep all three files together with the
results directory.
