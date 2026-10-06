# DIFLOW validation and benchmarking strategy

DIFLOW is not novel because it uses the joint site-frequency spectrum. The
jSFS is an established population-genetic summary and is intentionally used as
a tested statistical substrate.

The contribution DIFLOW is aiming for is a reproducible evidence architecture
for directional gene-flow inference.

## Validation question

The central question is not:

"Can the optimizer return two different migration parameters?"

It is:

"Under which demographic conditions can DIFLOW correctly recover migration
direction and magnitude, quantify uncertainty, and avoid claiming direction
when the data do not support it?"

## Minimum validation suite

The validated release should include simulations for:

- symmetric migration
- A -> B asymmetric migration
- B -> A asymmetric migration
- near-unidirectional migration
- unequal effective population sizes
- bottlenecks and growth
- secondary contact
- ancient migration
- range expansion
- ghost populations
- uneven sampling and missing data
- linked loci / reduced-representation genomic data
- zero migration

## Primary metrics

For migration magnitude:

- bias
- RMSE
- confidence-interval coverage

For direction:

- direction classification accuracy
- false directional-positive rate under symmetric migration
- sensitivity under weak, moderate, and strong asymmetry
- calibration of directional support probabilities

For model selection:

- correct-model selection frequency
- Akaike-weight calibration
- frequency of confusing continuous migration with secondary contact

For computation:

- runtime
- memory use
- scaling with loci, sample size, candidate pairs, optimization starts, and
  bootstrap replicates

## Comparator philosophy

DIFLOW should eventually be benchmarked against established approaches rather
than evaluated only against itself.

Relevant comparator classes include:

- demographic SFS inference packages
- historical migration estimators
- contemporary migration estimators
- fine-resolution asymmetric spatial migration methods
- descriptive directional statistics

Comparisons must respect differences in what each method actually estimates.
Backward-time lineage migration, recent migrant assignment, historical
coalescent migration, and forward-time demographic migration are not
interchangeable quantities.

## Release criterion

A map that looks plausible is not validation.

Before a validated release, DIFLOW should demonstrate:

1. low false directional-positive rates under symmetry,
2. reliable recovery of both directions under known asymmetry,
3. calibrated uncertainty,
4. robustness limits under major demographic confounders,
5. transparent failure modes,
6. reproducible benchmark datasets and scripts.

Where those conditions are not met, DIFLOW should report uncertainty or
ambiguity rather than force a directional conclusion.


## Current executable benchmark

DIFLOW now includes a first model-consistent parameter-recovery benchmark:

    diflow benchmark \
      --output benchmark_results/ \
      --replicates 10 \
      --chromosomes 20 \
      --sites 5000 \
      --starts 10 \
      --maxiter 100 \
      --seed 42

The benchmark generates finite-SNP jSFS replicates under known asymmetric
continuous-migration parameters and refits those spectra through DIFLOW's
unpolarized/folded inference path.

Current default truths include:

- symmetric migration
- moderate A -> B asymmetry
- moderate B -> A asymmetry
- strong A -> B asymmetry
- strong B -> A asymmetry

Outputs:

- recovery_replicates.csv
- recovery_summary.csv

The summary reports:

- optimization success rate
- bias in m(A -> B)
- RMSE in m(A -> B)
- bias in m(B -> A)
- RMSE in m(B -> A)
- direction accuracy
- false directional-positive rate under symmetric truth
- optimizer-stability rate

### What this benchmark proves

It tests numerical recovery and identifiability when the generating demographic
model matches the fitted model.

### What it does not prove

It does not establish robustness to:

- secondary contact when a continuous-migration model is fit
- population-size changes
- range expansion
- ghost populations
- linked loci
- ascertainment
- incorrect population assignment
- geographic graph misspecification
- other real-data violations

Those stress tests remain required before a validated release.


## Demographic stress benchmark

The second benchmark layer intentionally challenges DIFLOW with histories that
can create misleading directional patterns.

Run all currently implemented benchmark suites with:

    diflow benchmark \
      --suite all \
      --output benchmark_results/ \
      --replicates 10 \
      --chromosomes 20 \
      --sites 5000 \
      --starts 10

Run only the stress suite with:

    diflow benchmark \
      --suite stress \
      --output stress_results/

The current executable stress suite includes:

- zero migration
- symmetric migration with strongly unequal effective population sizes
- symmetric secondary contact
- A -> B asymmetric secondary contact
- B -> A asymmetric secondary contact

Each replicate is fit against the full candidate model set:

- isolation
- symmetric migration
- asymmetric continuous migration
- asymmetric secondary contact

The stress summary reports:

- optimization success rate
- correct generating-model selection rate
- provisional directional-signal rate
- directional recovery rate when direction truly exists
- false directional-signal rate when direction should not be inferred

A provisional directional signal currently requires sufficient asymmetric-model
Akaike weight, optimizer stability, and asymmetry magnitude. It is not the same
as the final bootstrap-supported direction classification.

### Why range expansion and ghost populations are not yet executable here

Those scenarios require explicit generating models beyond the current
two-population candidate-model family. DIFLOW will not approximate them with an
unrelated model merely to check a roadmap box. They will be added when the
simulation model represents the intended biology directly.


## Independent forward-time stress benchmark

DIFLOW also includes a stress generator that does **not** use dadi to generate
the data. Instead, it simulates explicit Wright-Fisher allele-frequency
histories and only uses dadi at the fitting stage.

Run it with:

    diflow benchmark \
      --suite forward \
      --output forward_stress_results/ \
      --replicates 10 \
      --chromosomes 20 \
      --sites 5000 \
      --starts 10

The current scenarios are:

- serial-founder range expansion with no ongoing A/B migration
- ghost-population introgression into B with no direct A/B migration
- a recent bottleneck in B with no migration
- uneven chromosome sampling under no migration

These scenarios are intentionally out-of-model. Their purpose is to measure
whether histories not represented by the fitted candidate set can induce a
false asymmetric-migration signal.

The forward-stress summary reports:

- fitting success rate
- false directional-signal rate
- frequency that asymmetric continuous migration is selected
- frequency that asymmetric secondary contact is selected

This layer is especially important because a method that performs well only
when the generating and fitted models are identical has not demonstrated
robustness to realistic demographic misspecification.


## Large recovery grid

For a broader scaling study, DIFLOW provides an explicit recovery-grid suite:

    diflow benchmark \
      --suite grid \
      --output recovery_grid_results/ \
      --replicates 10 \
      --starts 10 \
      --maxiter 100 \
      --seed 42

The default grid crosses:

- 10, 20, and 40 sampled chromosomes per population
- 1,000, 5,000, and 20,000 segregating sites
- symmetric migration
- weak, moderate, and strong A -> B asymmetry
- weak, moderate, and strong B -> A asymmetry

This produces 63 scenario/data-size cells before replicate expansion.

Outputs include:

- `recovery_grid_replicates.csv`
- `recovery_grid_summary.csv`
- `direction_accuracy_grid.png`
- `direction_accuracy_grid.pdf`
- `false_direction_rate_grid.png`
- `false_direction_rate_grid.pdf`

The grid is intentionally **not** included automatically in `--suite all`
because it can require substantial compute. Users must opt in with
`--suite grid`.

The generated figures are validation summaries, not final manuscript figures;
users should still inspect convergence and replicate counts before publication.


## Linked-marker bootstrap calibration

DIFLOW can directly compare naive locus resampling against genomic block
resampling on the same simulated correlated-marker datasets:

    diflow benchmark \
      --suite linked \
      --output linked_calibration/ \
      --replicates 20 \
      --chromosomes 20 \
      --linked-blocks 50 \
      --snps-per-block 10 \
      --linked-block-bp 100000 \
      --linkage-concentration 25 \
      --linked-bootstrap-replicates 100 \
      --starts 5

Each simulated genomic block shares a latent jSFS distribution. This induces
within-block dependence while preserving known migration truth.

The simulator is intentionally described as a **correlated-block model**. It is
not a mechanistic recombination or haplotype simulator.

Smaller `--linkage-concentration` values create stronger block-to-block
heterogeneity and stronger within-block dependence. Larger values make blocks
more similar to the global demographic expectation.

For every simulated dataset, DIFLOW analyzes the exact same allele counts using:

- independent locus bootstrap
- fixed genomic-window block bootstrap

The calibration reports:

- coverage of the true m(A -> B)
- coverage of the true m(B -> A)
- mean confidence-interval width
- strong directional-support rate
- false strong-direction support under symmetric migration
- number of loci and resampling blocks

Outputs include:

- `linked_bootstrap_replicates.csv`
- `linked_bootstrap_summary.csv`
- coverage comparison PNG/PDF figures
- false-direction comparison PNG/PDF figures

This benchmark tests whether block resampling improves uncertainty calibration
when markers are correlated. DIFLOW now also includes a separate msprime/tskit
mechanistic recombination validation layer.


## Mechanistic recombination validation

DIFLOW includes an optional msprime/tskit validation backend:

    python -m pip install -e ".[dev,demography,validation]"

Run:

    diflow benchmark \
      --suite mechanistic \
      --output mechanistic_linkage/ \
      --replicates 20 \
      --chromosomes 20 \
      --mechanistic-nref 10000 \
      --mechanistic-sequence-length 2000000 \
      --mechanistic-recombination-rate 1e-8 \
      --mechanistic-mutation-rate 1e-8 \
      --mechanistic-block-sizes 50000,100000,250000 \
      --linked-bootstrap-replicates 100

This layer simulates ancestry with recombination and neutral mutations directly.
Forward A -> B migration is translated into msprime's backward-time B -> A
lineage migration convention. Known per-generation rates are converted to dadi
scaled migration truth as 2 * Nref * m before interval coverage is evaluated.

The suite compares locus bootstrap with several candidate genomic block sizes.

## Decision-evidence benchmark and threshold calibration

Run a full classifier-evidence simulation:

    diflow benchmark \
      --suite decision \
      --output decision_benchmark/ \
      --replicates 50 \
      --decision-bootstrap-replicates 100

Then calibrate thresholds:

    diflow calibrate \
      --evidence decision_benchmark/decision_evidence.csv \
      --output calibrated_thresholds/ \
      --max-fpr 0.05

This explicitly controls the tolerated false directional-positive rate under
known symmetric truth while maximizing direction recovery under asymmetric
truth.
