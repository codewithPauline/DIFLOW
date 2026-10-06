# DIFLOW

**Directional Inference and Geographic Visualization of Gene Flow**

DIFLOW is an open-source population-genomics project being developed to estimate **directional migration between populations**, quantify the **magnitude and asymmetry of gene flow**, and visualize supported migration pathways on geographic maps.

> **Project status:** early research development. The statistical estimator is under active development and must be validated by simulation before biological interpretation.

## Why DIFLOW?

Spatial population-genetic methods are powerful, but direction, magnitude, uncertainty, and geographic visualization are often handled by separate tools or are not estimated in the same framework. DIFLOW is being designed around the quantities researchers usually want to interpret directly:

- \(m_{i\rightarrow j}\): migration from population *i* into population *j*
- \(m_{j\rightarrow i}\): migration in the reverse direction
- uncertainty around each estimate
- strength of directional asymmetry
- spatial location of the inferred movement

The goal is not simply to draw arrows between genetically similar populations. The goal is to build a validated inference framework in which every mapped arrow corresponds to an estimated migration parameter and carries explicit uncertainty.

## Core design

```text
Genomic data + population assignments + coordinates
                         |
                         v
                population summaries
                         |
                         v
              directional inference
                /              \
         m(i -> j)          m(j -> i)
                \              /
                 v            v
                 asymmetry + uncertainty
                         |
                         v
               directed migration graph
                         |
                         v
              geographic gene-flow map
```

## Development roadmap

### Phase 1 — mathematical and simulation foundation
- [x] Repository scaffold
- [x] Directional migration parameter conventions
- [x] Two-population forward migration simulator
- [x] Asymmetry utilities
- [ ] Two-population estimator
- [ ] Block/bootstrap uncertainty
- [ ] Recovery tests against known migration rates

### Phase 2 — genomic inference
- [ ] Allele-count and SFS input layer
- [ ] Joint-SFS likelihood
- [ ] Unequal effective population sizes
- [ ] Divergence and secondary-contact models
- [ ] Model comparison and diagnostics

### Phase 3 — spatial networks
- [ ] Multi-population graph construction
- [ ] Sparse migration parameterization
- [ ] Spatial regularization
- [ ] Source/sink and migration-hub summaries

### Phase 4 — geographic visualization
- [ ] Directional arrows
- [ ] Arrow width scaled by migration magnitude
- [ ] Support/uncertainty encoding
- [ ] Publication-quality PDF/SVG output
- [ ] Interactive exploration

### Phase 5 — validation
- [ ] Symmetric migration
- [ ] Strongly asymmetric migration
- [ ] Stepping-stone systems
- [ ] Range expansion
- [ ] Population-size asymmetry
- [ ] Secondary contact
- [ ] Ghost populations
- [ ] Missing populations and uneven sampling
- [ ] Benchmark bias, RMSE, interval coverage and direction accuracy

## Migration convention

DIFLOW uses **forward-time source-to-recipient notation**:

```text
m_A_to_B = proportion of population B replaced by migrants from A per generation
m_B_to_A = proportion of population A replaced by migrants from B per generation
```

For a pair of populations, directional asymmetry can be summarized as

```text
A_ij = (m_ij - m_ji) / (m_ij + m_ji)
```

which ranges from -1 to +1 when at least one rate is non-zero.

## Repository layout

```text
DIFLOW/
├── src/diflow/
│   ├── core/          # parameter conventions and summary statistics
│   └── simulation/    # validation simulators
├── tests/
├── pyproject.toml
├── LICENSE
└── README.md
```

## Scientific principle

A visually compelling directional map is not enough. DIFLOW will not treat a direction as biologically supported until the estimator can recover known migration direction and magnitude across controlled simulations and report its uncertainty.

## Installation

Development installation:

```bash
git clone https://github.com/codewithPauline/DIFLOW.git
cd DIFLOW
python -m pip install -e ".[dev]"
```

Run tests:

```bash
pytest
```

## License

DIFLOW is released under the MIT License.

## Citation

DIFLOW is under active development and has not yet received a formal software release or DOI. Citation information will be added with the first validated release.
