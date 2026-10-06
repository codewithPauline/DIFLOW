# DIFLOW

**Directional Inference and Geographic Visualization of Gene Flow**

DIFLOW is an open-source population-genomics framework for estimating asymmetric migration, comparing alternative demographic explanations, quantifying uncertainty in direction and magnitude, and mapping only those directional relationships that are supported by explicit evidence.

> **Project status:** active research development. DIFLOW is not yet a validated release. Simulation benchmarking and calibration remain required before biological conclusions should be based on its directional classifications.

## What DIFLOW is — and is not

DIFLOW uses the **joint site-frequency spectrum (jSFS)** as an established statistical substrate for demographic inference. The jSFS itself is not the novelty of DIFLOW.

The contribution DIFLOW is being designed around is the full evidence architecture:

```text
VCF genotypes
    ↓
population allele counts
    ↓
projected pairwise jSFS
    ↓
competing demographic models
    ↓
multi-start optimization
    ↓
model comparison
    ↓
m(A→B) and m(B→A)
    ↓
bootstrap uncertainty
    ↓
direction-support decision
    ↓
directed population network
    ↓
geographic migration map
```

DIFLOW does **not** infer migration direction directly from FST, ADMIXTURE coefficients, genetic distance, or visual ancestry patterns.

## Scientific question

For two populations A and B, DIFLOW asks whether the genomic data support:

- isolation,
- approximately symmetric migration,
- asymmetric migration,
- or an alternative history such as secondary contact.

When asymmetry is supported, DIFLOW estimates both directional parameters:

[
m_{Aightarrow B}
]

and

[
m_{Bightarrow A}.
]

Directional asymmetry is summarized as

[
A_{AB} =
rac{m_{Aightarrow B}-m_{Bightarrow A}}
     {m_{Aightarrow B}+m_{Bightarrow A}},
]

when at least one rate is non-zero.

## Current capabilities

DIFLOW currently includes:

- VCF + population-map input
- population allele-count construction
- projected pairwise jSFS construction
- explicit forward-time source→recipient migration convention
- isolation, symmetric-migration, asymmetric-migration, and asymmetric secondary-contact candidate models
- dadi-backed demographic inference
- AIC/AICc utilities and Akaike weights
- multi-start optimization
- convergence and stability diagnostics
- locus-bootstrap and fixed genomic-window block-bootstrap uncertainty for pairwise jSFS inference
- directional support probabilities
- supported / ambiguous / unsupported evidence classification
- sparse geographic candidate-pair construction
- directed migration networks
- source-like / sink-like network summaries
- directional map rendering
- an end-to-end `diflow infer` command
- simulation-validation metrics and canonical stress-test scenarios

## Documentation

New to DIFLOW? Start with the [Getting Started guide](docs/getting_started.md).

Key user documentation:

- [Input Files](docs/input_files.md) — prepare the VCF, population map, and coordinates.
- [Inspection Workflow](docs/inspection_workflow.md) — understand projection and graph recommendations.
- [Command-line Reference](docs/cli.md) — use `diflow inspect` and `diflow infer`.
- [Understanding Outputs](docs/outputs.md) — interpret result tables, network summaries, and maps.
- [Full Documentation Index](docs/README.md) — statistical methods, uncertainty, mapping, and validation.

Typical first step:

```bash
diflow inspect \
  --vcf data.vcf \
  --popmap populations.tsv \
  --coords coordinates.csv
```

DIFLOW will summarize the dataset, show projection-retention tradeoffs, compare candidate geographic graphs, and print a ready-to-copy inference command.

## Command-line workflow

Prepare genomic inputs and spectra without fitting demographic models:

```bash
diflow infer \
  --vcf data.vcf \
  --popmap populations.tsv \
  --coords coordinates.csv \
  --projection-chromosomes 8 \
  --neighbors 4 \
  --output results/ \
  --prepare-only
```

Run demographic inference with locus-bootstrap uncertainty:

```bash
diflow infer \
  --vcf data.vcf \
  --popmap populations.tsv \
  --coords coordinates.csv \
  --projection-chromosomes 8 \
  --neighbors 4 \
  --starts 20 \
  --bootstrap-replicates 100 \
  --bootstrap-starts 5 \
  --output results/
```

Bootstrap demographic inference is computationally expensive because each replicate refits an asymmetric demographic model.

For linked SNPs, users can resample fixed genomic windows instead of individual loci:

```bash
diflow infer \
  --vcf data.vcf \
  --popmap populations.tsv \
  --coords coordinates.csv \
  --projection-chromosomes 8 \
  --neighbors 4 \
  --bootstrap-replicates 100 \
  --bootstrap-block-bp 100000 \
  --output results/
```

DIFLOW records the resampling unit and number of blocks in the output so uncertainty provenance remains explicit.

## Known-truth benchmark

DIFLOW includes an executable first-stage recovery benchmark:

```bash
diflow benchmark \
  --output benchmark_results/ \
  --replicates 10 \
  --chromosomes 20 \
  --sites 5000 \
  --starts 10
```

It simulates jSFS datasets under known symmetric and asymmetric migration
parameters, refits them through the same folded-data inference path, and reports
bias, RMSE, direction accuracy, optimization success, and the false
directional-positive rate under symmetric migration.

This is a **model-consistent recovery benchmark**, not complete biological
validation. Demographic misspecification and other stress tests remain part of
the validation roadmap.

See [docs/benchmarking.md](docs/benchmarking.md).

For deliberately out-of-model histories, including serial-founder range
expansion and unsampled ghost introgression, run:

```bash
diflow benchmark --suite forward --output forward_stress_results/
```

## Linked-marker calibration

To compare locus bootstrap with genomic block bootstrap under correlated marker
blocks:

```bash
diflow benchmark \
  --suite linked \
  --output linked_calibration/ \
  --replicates 20 \
  --linked-blocks 50 \
  --snps-per-block 10 \
  --linked-block-bp 100000 \
  --linked-bootstrap-replicates 100
```

The calibration measures confidence-interval coverage, interval width, and false
strong-direction support under symmetric migration.

## Large recovery grid

For a compute-intensive scaling study across sample sizes, SNP counts, and
asymmetry strengths:

```bash
diflow benchmark \
  --suite grid \
  --output recovery_grid_results/ \
  --replicates 10 \
  --starts 10
```

The grid writes replicate and summary CSV files plus PNG/PDF validation plots
for direction accuracy and false directional-positive rates. It requires
explicit opt-in and is not included in `--suite all`.

See [docs/benchmarking.md](docs/benchmarking.md).

## Directional evidence

A large fitted migration rate is not automatically treated as a supported arrow.

The decision layer can combine:

- support for the asymmetric demographic model,
- optimizer stability,
- directional asymmetry,
- bootstrap support,
- uncertainty-interval separation.

The current decision thresholds are development defaults and must be calibrated by simulation.

## Multi-population scaling

For (n) populations, unrestricted pairwise comparison requires

[
rac{n(n-1)}{2}
]

population pairs.

DIFLOW therefore supports sparse candidate graphs based on geographic distance and k-nearest neighbors. Future versions will add biologically informed adjacency such as watersheds, habitat connectivity, resistance surfaces, and user-supplied graphs.

## Validation is the core product requirement

The most important question for DIFLOW is not whether it can return two migration parameters.

It is whether it can recover direction and magnitude under known truth **without producing false directional conclusions under symmetry or demographic misspecification**.

The validation suite is being built around:

- symmetric migration
- A→B and B→A asymmetry
- near-unidirectional migration
- zero migration
- unequal effective population sizes
- bottlenecks and growth
- secondary contact
- ancient migration
- range expansion
- ghost populations
- uneven sampling and missingness
- linked loci

Primary metrics include:

- parameter bias
- RMSE
- confidence-interval coverage
- direction accuracy
- false directional-positive rate
- model-selection accuracy
- runtime and memory scaling

See [docs/benchmarking.md](docs/benchmarking.md).

## Migration convention

DIFLOW uses forward-time source-to-recipient notation:

```text
m_A_to_B = migration from population A into population B
m_B_to_A = migration from population B into population A
```

Backend-specific parameter ordering is translated internally and tested explicitly.

## Repository layout

```text
DIFLOW/
├── src/diflow/
│   ├── core/          # asymmetry and parameter conventions
│   ├── io/            # VCF and population-map input
│   ├── spectra/       # jSFS construction and projection
│   ├── demography/    # demographic models, fitting, uncertainty
│   ├── decision/      # evidence-based direction classification
│   ├── network/       # sparse pair graphs and directed summaries
│   ├── mapping/       # geographic directional visualization
│   ├── pipeline/      # end-to-end orchestration
│   ├── simulation/    # controlled simulation utilities
│   └── validation/    # benchmark scenarios and metrics
├── tests/
├── docs/
├── pyproject.toml
├── LICENSE
└── README.md
```

## Development roadmap

### Statistical foundation
- [x] Forward-time migration convention
- [x] Asymmetry summaries
- [x] Two-population controlled simulator
- [x] jSFS construction and projection
- [x] dadi demographic backend
- [x] Competing migration models
- [x] Multi-start optimization
- [x] Model comparison
- [x] Locus-bootstrap uncertainty
- [x] Genomic block-bootstrap uncertainty
- [ ] Profile-likelihood diagnostics
- [x] Unpolarized/folded default for ordinary VCF data
- [ ] Explicit ancestral-polarization input workflow

### Spatial inference and reporting
- [x] Geographic candidate-pair graphs
- [x] Directed migration networks
- [x] Source-like / sink-like summaries
- [x] Directional map rendering
- [x] End-to-end CLI
- [ ] Projected publication cartography
- [ ] Parallel/HPC execution
- [ ] Interactive exploration

### Validation
- [x] Validation metric framework
- [x] Canonical benchmark scenario registry
- [x] Large recovery-grid runner and validation plots
- [ ] Final simulation calibration and release thresholds
- [ ] False-positive calibration under symmetry
- [x] Secondary-contact stress benchmark
- [x] Range-expansion forward-time stress test
- [x] Ghost-population forward-time stress test
- [x] Correlated-block locus-vs-block bootstrap calibration framework
- [ ] Mechanistic LD/recombination validation and final coverage calibration
- [ ] Benchmark against established methods

## Scientific guardrails

DIFLOW separates four distinct concepts:

1. **genetic differentiation** — e.g. FST,
2. **descriptive allele-frequency structure**,
3. **model-based migration parameters**,
4. **evidence for directional asymmetry**.

These quantities should not be treated as interchangeable.

A visually compelling directional map is not enough. Every strong arrow must trace back to a fitted demographic parameter, uncertainty estimate, model-comparison result, and explicit evidence classification.

## Installation

Development installation:

```bash
git clone https://github.com/codewithPauline/DIFLOW.git
cd DIFLOW
python -m pip install -e ".[dev]"
```

Install the demographic backend:

```bash
python -m pip install -e ".[dev,demography]"
```

Run tests:

```bash
pytest
```

## License

DIFLOW is released under the MIT License.

## Citation

DIFLOW is under active development and has not yet received a formal software release or DOI. Citation information will be added with the first validated release.
