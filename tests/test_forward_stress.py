import numpy as np
import pandas as pd
import pytest

from diflow.validation.forward_stress import (
    default_forward_stress_scenarios,
    simulate_forward_stress_spectrum,
    summarize_forward_stress,
)


def _scenario(kind):
    return next(x for x in default_forward_stress_scenarios() if x.kind == kind)


@pytest.mark.parametrize(
    "kind",
    ["range_expansion", "ghost_introgression", "bottleneck", "uneven_sampling"],
)
def test_forward_stress_generates_variant_spectrum(kind):
    spectrum = simulate_forward_stress_spectrum(
        _scenario(kind),
        loci=300,
        sample_sizes=(8, 10),
        seed=12,
    )
    assert spectrum.shape == (9, 11)
    assert spectrum.sum() > 0
    assert spectrum[0, 0] == 0
    assert spectrum[-1, -1] == 0


def test_forward_stress_is_reproducible():
    scenario = _scenario("range_expansion")
    a = simulate_forward_stress_spectrum(
        scenario, loci=300, sample_sizes=(8, 8), seed=5
    )
    b = simulate_forward_stress_spectrum(
        scenario, loci=300, sample_sizes=(8, 8), seed=5
    )
    assert np.array_equal(a, b)


def test_forward_stress_summary_counts_false_signal():
    results = pd.DataFrame(
        {
            "scenario": ["ghost"] * 4,
            "success": [True] * 4,
            "provisional_directional_signal": [False, True, False, False],
            "best_model": [
                "isolation",
                "asymmetric_migration",
                "secondary_contact_asymmetric",
                "symmetric_migration",
            ],
        }
    )
    summary = summarize_forward_stress(results)
    row = summary.iloc[0]
    assert row["false_direction_signal_rate"] == pytest.approx(0.25)
    assert row["asymmetric_model_selected_rate"] == pytest.approx(0.25)
    assert row["secondary_contact_selected_rate"] == pytest.approx(0.25)
