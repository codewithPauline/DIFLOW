# Known limitations and failure regimes

DIFLOW is active research software. Its strongest scientific claims should be
limited to regimes supported by completed simulation validation.

## Migration estimates are model-dependent

DIFLOW estimates migration parameters under explicit demographic models. A
numerical estimate is not a model-free measurement of organism movement.

Current pairwise candidate histories include isolation, symmetric continuous
migration, asymmetric continuous migration, and asymmetric secondary contact.
Histories outside this set can bias migration estimates or model selection.

## Scaled migration is not automatically migrants per generation

The current dadi-backed migration parameters are scaled demographic quantities.
They should not be described as literal numbers of migrants per generation
without the population-size and scaling information required for conversion.

## Pairwise inference can omit important populations

An unsampled or excluded population can create apparent A/B migration signals.
Ghost-population stress tests are included in validation, but pairwise analysis
cannot reconstruct an unsampled population that is absent from the model.

## Population definitions matter

DIFLOW assumes the supplied population map represents meaningful demographic
units. Pooling distinct demes or splitting a panmictic population can distort
the jSFS and migration inference.

Population labels should therefore come from defensible sampling design and
biological context rather than from the desired direction of a result.

## Range expansion and non-equilibrium history can mimic direction

Founder effects, recent expansion, bottlenecks, and other non-equilibrium
histories can produce asymmetric allele-frequency patterns even without ongoing
direct migration.

DIFLOW includes forward-time and demographic stress tests specifically because
direction should not be inferred from asymmetry alone.

## Linkage can make uncertainty too narrow

Nearby SNPs may not be independent. Naive locus/SNP resampling can therefore
overstate effective information and produce confidence intervals that are too
narrow.

DIFLOW supports genomic block bootstrap and includes correlated-block plus
msprime/tskit recombination validation frameworks. Final release claims still
require the large empirical linkage-coverage campaign.

## Block size is not universal

A fixed genomic block size is a user-controlled approximation to the dependence
scale. The correct value depends on recombination, marker design, genome
structure, and study system.

Sensitivity across multiple plausible block sizes should be examined when the
linkage scale is uncertain.

## Variant ascertainment matters

Variant-only VCF data omit monomorphic sites and reduced-representation methods
can have additional ascertainment effects. DIFLOW masks unobservable fixed
corners and treats ordinary REF/ALT coding as unpolarized by default, but this
does not remove every ascertainment bias.

## Ancestral polarization must be justified

The default folded analysis is appropriate when ancestral state is unknown.
Use unfolded/polarized mode only when ALT has genuinely been established as the
derived allele upstream.

Incorrect polarization can create misleading frequency-spectrum asymmetry.

## Sparse data can be weakly identifiable

A fitted optimum does not guarantee that directional parameters are well
identified. DIFLOW provides multi-start stability diagnostics, bootstrap
uncertainty, and profile likelihoods because migration parameters can remain
weakly identified even when optimization succeeds.

No universal minimum SNP count or sample size is currently claimed. Final
recommended minimum data requirements must be derived from the release-grade
simulation campaign across sample sizes, marker counts, and asymmetry strengths.

## Geographic candidate graphs restrict what is tested

k-nearest-neighbor and distance filters are computational and biological
hypotheses about which population pairs deserve testing. Excluding a pair means
DIFLOW does not evaluate that migration edge.

Maps visualize fitted evidence; geography does not itself establish migration
direction.

## Direction thresholds are provisional until empirical calibration is frozen

The software can calibrate model-weight, directional-support, and asymmetry
thresholds against known truth. Development defaults should not be described as
universal significance cutoffs.

Release defaults will only be frozen after the empirical decision-evidence
campaign meets the declared false-direction and sensitivity targets.

## External methods may estimate different quantities

Historical demographic migration, recent assignment-based migration,
effective-lineage migration, and spatial effective migration are not identical
estimands.

DIFLOW's comparison framework therefore requires matched truth and compatible
estimands for magnitude comparisons. Direction-only comparison is available
when only directional ordering is scientifically comparable.

## Real data cannot prove estimator accuracy

A real-data case study can demonstrate usability and biological interpretation,
but the true migration history is unknown. Accuracy, calibration, and false
positive control must be established using simulation with known truth.

## Current release boundary

Before a validated 1.0 release, DIFLOW still requires:

1. completion of the large recovery and decision-evidence campaigns,
2. frozen empirically justified directional thresholds,
3. large mechanistic linkage-coverage calibration,
4. matched execution against at least one established external method,
5. a real-data demonstration,
6. archival of validation outputs and final release review.

See [release_readiness.md](release_readiness.md) for the release policy.
