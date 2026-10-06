# External method comparison

DIFLOW provides a standardized comparison layer for simulation results from
other migration methods.

## Standard input schema

Each method CSV must contain:

- `scenario`
- `replicate`
- `truth_m_a_to_b`
- `truth_m_b_to_a`
- `estimated_m_a_to_b`
- `estimated_m_b_to_a`

Convert DIFLOW recovery output to the standard schema first:

    diflow benchmark-export \
      --input recovery_grid_results/recovery_grid_replicates.csv \
      --output diflow_standardized.csv \
      --estimand dadi_scaled_migration

Compare two or more methods with:

    diflow compare \
      --method DIFLOW=diflow_standardized.csv \
      --method OTHER=other_method_results.csv \
      --output comparison/

Outputs include combined replicate records, method/scenario summaries, and
direction-accuracy / false-direction plots when applicable.

## Scientific requirement

Only compare quantities with genuinely compatible interpretations.

A recent-migrant assignment estimate, a backward-time lineage migration rate,
a historical coalescent migration parameter, and a dadi-scaled demographic
migration parameter are not automatically interchangeable.

The comparison framework is implemented; the final established-method
benchmark study still requires actually running and normalizing the selected
external tools under matched simulated truth.


## Matched-comparison requirement

Release-grade comparisons are strict by default.

Every method must contain the same `(scenario, replicate)` keys, and the known
truth values must match exactly across methods. DIFLOW also rejects a combined
table when multiple non-empty `estimand` labels are present.

This prevents a visually attractive comparison from silently mixing:

- different simulation replicates
- different truth parameters
- incompatible migration estimands

For exploratory work only, `--allow-unmatched` disables the complete-key
requirement. Truth inconsistencies and mixed estimands remain scientific
problems and should not be ignored in publication analyses.


For the release-grade execution protocol, including provenance, estimand
compatibility, failed-run reporting, and direction-only comparisons, see
[External Benchmark Protocol](external_benchmark_protocol.md).
