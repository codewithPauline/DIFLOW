# Getting started with DIFLOW

This guide takes a new user from input files to interpretable directional gene-flow results.

## 1. Install DIFLOW

    git clone https://github.com/codewithPauline/DIFLOW.git
    cd DIFLOW
    python -m pip install -e ".[dev,demography]"

Check the command:

    diflow --help

## 2. Prepare the three core inputs

DIFLOW uses:

1. a VCF containing individual genotypes,
2. a sample-to-population mapping file,
3. a population coordinate file.

See [input_files.md](input_files.md) for exact formats.

## 3. Inspect before fitting

    diflow inspect \
      --vcf data.vcf \
      --popmap populations.tsv \
      --coords coordinates.csv

Inspection reports sample count, population count, variant-locus count, projection-retention tradeoffs, a recommended jSFS projection, several candidate geographic graphs, a suggested analysis preset, and a ready-to-copy inference command.

## 4. Review the projection recommendation

The default inspection heuristic recommends the largest even projection that retains at least 80% of variant loci in the worst-retained population. This is a starting rule, not a biological constant.

## 5. Review the candidate population graph

The graph defines which pairs are tested. It does not itself infer gene flow.

## 6. Choose analysis effort

| preset | model starts | bootstrap replicates | bootstrap starts | max iterations |
| --- | ---: | ---: | ---: | ---: |
| quick | 5 | 20 | 2 | 60 |
| standard | 20 | 100 | 5 | 100 |
| publication | 40 | 500 | 8 | 200 |

Use quick for workflow testing, standard for routine analysis, and publication as a high-effort preset rather than a guarantee of publication-quality inference.

## 7. Run inference

    diflow infer \
      --vcf data.vcf \
      --popmap populations.tsv \
      --coords coordinates.csv \
      --projection-chromosomes 8 \
      --neighbors 4 \
      --preset standard \
      --output results/

## 8. Understand the demographic comparison

For each candidate pair, DIFLOW currently compares isolation, symmetric continuous migration, asymmetric continuous migration, and asymmetric secondary contact.

## 9. Understand uncertainty and status

With bootstrap enabled, DIFLOW rebuilds projected jSFS replicates, refits the asymmetric model, and measures how consistently one direction exceeds the reverse direction.

By default, loci are resampled independently. For linked SNPs with meaningful genomic positions, users can request fixed-window block bootstrap with `--bootstrap-block-bp`.

Formal status values are supported, ambiguous, and unsupported.

## 10. Inspect the outputs

Important outputs include allele_counts.csv, candidate_pairs.csv, pairwise_results.csv, model_rankings.csv, directional_flows.csv, network_summary.csv, run_metadata.csv, spectra/*.npy, directional_map.png, and directional_map.pdf.

See [outputs.md](outputs.md) for interpretation.

## 11. Do not overinterpret the map

The map is a visualization of model-based evidence. It is not independent proof of organismal movement.

## 12. Current research-stage limitations

Important remaining work includes explicit ancestral-polarization input, simulation calibration of block-bootstrap coverage under realistic linkage, decision-threshold calibration, large benchmark studies, broader real-data testing, parallel/HPC execution, and indexed .vcf.gz/BCF support.

Do not treat a current development build as a validated black-box estimator.