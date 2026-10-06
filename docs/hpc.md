# HPC and Slurm

DIFLOW demographic fitting and bootstrap inference can be computationally
expensive. The CLI can generate a reproducible Slurm script for an existing
DIFLOW command.

Example:

    diflow slurm \
      --run-command "diflow infer --vcf data.vcf.gz --popmap populations.tsv --coords coordinates.csv --projection-chromosomes 8 --neighbors 4 --preset publication --output results/" \
      --script diflow_run.slurm \
      --job-name DIFLOW_main \
      --cpus 8 \
      --mem-gb 64 \
      --hours 72

Then submit:

    sbatch diflow_run.slurm

The generated script includes strict shell error handling, Slurm resource
directives, job logs, and start/finish timestamps.

## Important distinction

This provides an HPC execution workflow for the current inference command.
DIFLOW does not yet claim internally benchmarked pair-level multiprocessing.
True pair-parallel execution remains a separate optimization task because dadi
fits are memory-intensive and backend/process behavior must be tested carefully.
