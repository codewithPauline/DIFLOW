# Pair-level parallel execution

DIFLOW can distribute independent population-pair analyses across isolated
Python processes:

    diflow infer \
      --vcf data.vcf.gz \
      --popmap populations.tsv \
      --coords coordinates.csv \
      --projection-chromosomes 8 \
      --neighbors 4 \
      --workers 4 \
      --output results/

`--workers 1` preserves serial execution and is the default.

## Design

Each candidate pair is represented as an isolated task. A worker receives only
the allele-count rows for the two populations in that pair, performs spectrum
construction and inference without writing shared files, and returns its
results to the parent process.

The parent process:

- restores deterministic candidate-pair order
- writes spectra
- writes pairwise/model-ranking outputs
- builds network and map outputs

This prevents multiple workers from writing the same output files concurrently.

## Resource guidance

Demographic optimization and bootstrap fitting can be memory intensive.
Increasing workers multiplies the number of fits that may be resident at once.

Start conservatively on HPC systems and match Slurm CPU/memory requests to the
chosen worker count. More workers are not automatically faster when memory
bandwidth or per-process memory becomes limiting.
