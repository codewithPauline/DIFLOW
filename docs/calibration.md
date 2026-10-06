# Decision-threshold calibration

DIFLOW's final directional classification combines several evidence layers:

- Akaike weight for the asymmetric model
- optimizer stability
- absolute migration asymmetry
- bootstrap directional support
- separation of directional migration intervals

The numeric thresholds should not be treated as universal constants.

## Generate known-truth evidence

Run:

    diflow benchmark \
      --suite decision \
      --output decision_benchmark/ \
      --replicates 50 \
      --chromosomes 20 \
      --sites 5000 \
      --starts 10 \
      --decision-bootstrap-replicates 100

This writes:

    decision_benchmark/decision_evidence.csv

Each successful row contains the same evidence quantities used by the production
direction classifier, together with the known simulated direction.

## Calibrate thresholds

Use the evidence table to select thresholds subject to a false-direction target:

    diflow calibrate \
      --evidence decision_benchmark/decision_evidence.csv \
      --output calibrated_thresholds/ \
      --max-fpr 0.05

Outputs:

- threshold_scan.csv
- selected_thresholds.csv

The selector chooses the most sensitive tested rule whose false directional
positive rate is at or below the requested target. Ties favor higher direction
accuracy when called, followed by more conservative evidence cutoffs.

## Apply calibrated thresholds

The selected values can then be supplied to inference:

    diflow infer \
      --vcf data.vcf.gz \
      --popmap populations.tsv \
      --coords coordinates.csv \
      --projection-chromosomes 8 \
      --min-model-weight 0.80 \
      --min-directional-support 0.975 \
      --min-abs-asymmetry 0.30 \
      --output results/

The values above are only an example. Use values supported by the calibration
study for the intended data regime.

## Release policy

DIFLOW currently provides the calibration machinery, but final release defaults
should only be frozen after broad simulations across sample size, marker count,
linkage, demographic misspecification, and migration strength.
