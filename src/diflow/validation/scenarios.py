"""Canonical simulation scenarios for validating DIFLOW."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ValidationScenario:
    """One named biological stress test for directional inference."""

    name: str
    purpose: str
    expected_direction: str
    confounders: tuple[str, ...] = ()


def core_validation_scenarios() -> tuple[ValidationScenario, ...]:
    """Return the minimum scenario suite required before validated release."""
    return (
        ValidationScenario(
            name="symmetric_migration",
            purpose="Measure false directional-positive rate under equal migration.",
            expected_direction="symmetric",
        ),
        ValidationScenario(
            name="a_to_b_asymmetry",
            purpose="Test recovery when migration is stronger from A into B.",
            expected_direction="A->B",
        ),
        ValidationScenario(
            name="b_to_a_asymmetry",
            purpose="Test recovery when migration is stronger from B into A.",
            expected_direction="B->A",
        ),
        ValidationScenario(
            name="unidirectional_a_to_b",
            purpose="Test near-one-way migration recovery.",
            expected_direction="A->B",
        ),
        ValidationScenario(
            name="unequal_effective_size",
            purpose="Test whether Ne asymmetry is mistaken for migration asymmetry.",
            expected_direction="scenario-dependent",
            confounders=("unequal_Ne",),
        ),
        ValidationScenario(
            name="secondary_contact",
            purpose="Test whether recent contact is confused with continuous migration.",
            expected_direction="scenario-dependent",
            confounders=("secondary_contact",),
        ),
        ValidationScenario(
            name="range_expansion",
            purpose="Measure directional artifacts caused by non-equilibrium expansion.",
            expected_direction="none-by-default",
            confounders=("range_expansion",),
        ),
        ValidationScenario(
            name="ghost_population",
            purpose="Measure bias from an unsampled donor or recipient population.",
            expected_direction="scenario-dependent",
            confounders=("ghost_population",),
        ),
        ValidationScenario(
            name="uneven_sampling",
            purpose="Test robustness to unequal sample sizes and missing genotypes.",
            expected_direction="scenario-dependent",
            confounders=("sample_size_imbalance", "missing_data"),
        ),
        ValidationScenario(
            name="linked_loci",
            purpose="Evaluate uncertainty calibration when markers are not independent.",
            expected_direction="scenario-dependent",
            confounders=("linkage",),
        ),
    )
