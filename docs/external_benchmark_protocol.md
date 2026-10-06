# External benchmark protocol

A release-grade DIFLOW benchmark must compare methods on the **same simulated
truth**, not on loosely similar datasets.

## 1. Freeze the simulation campaign

Before running external tools, archive the simulation configuration and random
seeds used for the DIFLOW benchmark.

Every method should receive data generated from the same scenario/replicate
identifiers whenever its input format permits.

The standardized comparison key is:

    scenario + replicate

## 2. Record method provenance

For every comparator, record:

- method name
- software version
- executable or package version
- exact command/configuration
- random seed policy
- convergence/filtering rule
- native migration parameterization
- forward/backward-time direction convention
- any transformation used before comparison

Do not compare an undocumented external run.

## 3. Label the estimand

Every standardized result file should include an `estimand` column.

Examples of different classes of estimand include:

- scaled demographic migration
- per-generation migration probability
- migrants per generation
- recent migrant assignment
- backward-time lineage migration
- spatial effective migration

Different labels are not automatically magnitude-compatible.

## 4. Magnitude-compatible comparison

Use the default comparison mode only after establishing that the compared
numbers have the same interpretation and scale.

    diflow compare \
      --method DIFLOW=diflow_standardized.csv \
      --method METHOD_X=method_x_standardized.csv \
      --output comparison/

The default mode requires:

- identical scenario/replicate keys
- identical known truth values
- a single compatible estimand label

It reports bias, RMSE, direction accuracy, and false directional-positive rate.

## 5. Direction-only comparison

When methods estimate different migration quantities but their directional
ordering can still be meaningfully evaluated, use:

    diflow compare \
      --direction-only \
      --method DIFLOW=diflow_standardized.csv \
      --method METHOD_Y=method_y_standardized.csv \
      --output direction_comparison/

Direction-only mode reports directional recovery and false-direction rates. It
does **not** report magnitude bias/RMSE across incompatible estimands.

## 6. DIFLOW export

Convert a DIFLOW recovery-grid replicate table to the standardized schema with:

    diflow benchmark-export \
      --input recovery_grid/recovery_grid_replicates.csv \
      --output diflow_standardized.csv \
      --estimand dadi_scaled_migration

## 7. External result schema

Each external standardized CSV must contain:

- `scenario`
- `replicate`
- `truth_m_a_to_b`
- `truth_m_b_to_a`
- `estimated_m_a_to_b`
- `estimated_m_b_to_a`

Recommended additional columns:

- `estimand`
- `method_version`
- `run_status`
- `runtime_seconds`
- `notes`

The A -> B and B -> A columns must follow DIFLOW's forward-time
source-to-recipient convention after any required translation.

## 8. Recommended comparator roles

A strong study should include at least one established demographic method that
can be parameterized on matched two-population migration scenarios.

Methods whose migration quantity is not directly scale-compatible can still be
useful as direction-only comparators.

Spatial smoothing or descriptive differentiation methods should be treated as
contextual comparisons rather than forced into a migration-magnitude table.

## 9. Failed runs

Do not silently drop failed external-method replicates.

Retain a run-status table and report:

- attempted replicates
- successful replicates
- failure/convergence rate

The standardized numeric comparison may contain only successful estimates, but
the failure rate must be reported alongside accuracy metrics.

## 10. Release requirement

The external benchmark item is complete only after:

1. at least one established external method has actually been executed,
2. matched scenario/replicate truth is verified,
3. method provenance and estimand are documented,
4. failed-run rates are reported,
5. comparison outputs are archived,
6. `diflow release-review` recognizes the external comparison output.

Building the adapter alone does not satisfy this requirement.
