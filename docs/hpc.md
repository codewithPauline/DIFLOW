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

## Pair-level process parallelism

DIFLOW can also execute independent population pairs in isolated processes with
`--workers N`.

Example:

    diflow infer \
      ... \
      --workers 4

The default remains `--workers 1`.

Each worker receives only the allele-count rows needed for its population pair,
and the parent process restores deterministic pair ordering and performs all
shared output writes.

Increase workers conservatively because each process may run memory-intensive
dadi optimization/bootstrap fits. Slurm CPU and memory requests should be sized
for the chosen worker count.
