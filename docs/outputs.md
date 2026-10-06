# Understanding DIFLOW outputs

The output directory preserves both final results and the evidence behind them.

## allele_counts.csv

Population-level REF and ALT allele counts derived from the VCF. This is the bridge between raw genotypes and the jSFS layer.

## candidate_pairs.csv

Population pairs selected for demographic testing. These rows describe the candidate graph; they are not inferred migration edges.

## spectra/

Projected pairwise jSFS arrays used for demographic inference.

## model_rankings.csv

Model-comparison results for each population pair. Candidate models currently include isolation, symmetric migration, asymmetric migration, and asymmetric secondary contact.

## pairwise_results.csv

The main pairwise inference table. Depending on settings it can include fitted migration in both directions, asymmetry index, asymmetric-model weight, optimizer stability, bootstrap confidence intervals, directional support, preferred direction, final evidence status, and decision reason.

### Scaled migration

Migration parameters from the dadi backend are scaled demographic quantities. They should not automatically be described as literal numbers of individual migrants per generation without the appropriate scaling assumptions.

## Evidence status

supported: sufficient combined evidence under the current decision policy.

ambiguous: one direction may be preferred, but evidence is not strong enough for a supported call.

unsupported: available evidence does not support a directional edge.

Current thresholds are development defaults and still require simulation calibration.

## directional_flows.csv

A standardized table of formally classified directed edges. Provisional candidate results are intentionally excluded.

## network_summary.csv

Population-level summaries of incoming and outgoing inferred migration. Source-like and sink-like are network summaries, not proof of ecological source-sink demography.

## directional_map.png and directional_map.pdf

Geographic visualization of directional flow. The map should always be interpreted together with the statistical tables.

## run_metadata.csv

Records important settings used for the run, supporting reproducibility.

## Recommended interpretation order

1. check data retention and candidate pairs,
2. inspect optimizer stability,
3. inspect model rankings,
4. compare both migration estimates,
5. inspect confidence intervals and directional support,
6. inspect evidence status,
7. then interpret the network and map.