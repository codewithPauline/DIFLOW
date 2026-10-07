import numpy as np
import pandas as pd
import pytest

import diflow.validation.forward_stress as forward_stress
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



def test_forward_stress_uses_selected_directional_model(monkeypatch):
    params = {
        "isolation": {"nu_a": 1.0, "nu_b": 1.0, "split_time": 0.5},
        "symmetric_migration": {
            "nu_a": 1.0,
            "nu_b": 1.0,
            "split_time": 0.5,
            "migration": 0.4,
        },
        "asymmetric_migration": {
            "nu_a": 1.0,
            "nu_b": 1.0,
            "split_time": 0.5,
            "m_a_to_b": 0.9,
            "m_b_to_a": 0.1,
        },
        "secondary_contact_symmetric": {
            "nu_a": 1.0,
            "nu_b": 1.0,
            "isolation_time": 0.5,
            "contact_time": 0.2,
            "migration": 0.4,
        },
        "secondary_contact_asymmetric": {
            "nu_a": 1.0,
            "nu_b": 1.0,
            "isolation_time": 0.5,
            "contact_time": 0.2,
            "m_a_to_b": 0.2,
            "m_b_to_a": 1.1,
        },
    }

    def fake_fit(observed, model_name, **kwargs):
        return type(
            "Fit",
            (),
            {
                "best_parameters": params[model_name],
                "best_log_likelihood": -10.0,
                "stable": True,
                "converged_fraction": 1.0,
            },
        )()

    monkeypatch.setattr(forward_stress, "fit_multistart", fake_fit)

    monkeypatch.setattr(
        forward_stress,
        "rank_models",
        lambda scores, use_aicc=False: pd.DataFrame(
            {
                "model": [
                    "secondary_contact_asymmetric",
                    "secondary_contact_symmetric",
                    "asymmetric_migration",
                    "symmetric_migration",
                    "isolation",
                ],
                "akaike_weight": [0.82, 0.08, 0.05, 0.03, 0.02],
            }
        ),
    )

    result = forward_stress._fit_full_candidate_set(
        np.ones((5, 5)),
        starts=1,
        maxiter=10,
        seed=4,
    )

    assert result["directional_model"] == "secondary_contact_asymmetric"
    assert result["estimated_m_a_to_b"] == 0.2
    assert result["estimated_m_b_to_a"] == 1.1
    assert result["preferred_direction"] == "B->A"


def test_forward_stress_emits_no_direction_for_symmetric_winner(monkeypatch):
    params = {
        "isolation": {"nu_a": 1.0, "nu_b": 1.0, "split_time": 0.5},
        "symmetric_migration": {
            "nu_a": 1.0,
            "nu_b": 1.0,
            "split_time": 0.5,
            "migration": 0.4,
        },
        "asymmetric_migration": {
            "nu_a": 1.0,
            "nu_b": 1.0,
            "split_time": 0.5,
            "m_a_to_b": 0.9,
            "m_b_to_a": 0.1,
        },
        "secondary_contact_symmetric": {
            "nu_a": 1.0,
            "nu_b": 1.0,
            "isolation_time": 0.5,
            "contact_time": 0.2,
            "migration": 0.4,
        },
        "secondary_contact_asymmetric": {
            "nu_a": 1.0,
            "nu_b": 1.0,
            "isolation_time": 0.5,
            "contact_time": 0.2,
            "m_a_to_b": 0.2,
            "m_b_to_a": 1.1,
        },
    }

    def fake_fit(observed, model_name, **kwargs):
        return type(
            "Fit",
            (),
            {
                "best_parameters": params[model_name],
                "best_log_likelihood": -10.0,
                "stable": True,
                "converged_fraction": 1.0,
            },
        )()

    monkeypatch.setattr(forward_stress, "fit_multistart", fake_fit)

    monkeypatch.setattr(
        forward_stress,
        "rank_models",
        lambda scores, use_aicc=False: pd.DataFrame(
            {
                "model": [
                    "secondary_contact_symmetric",
                    "secondary_contact_asymmetric",
                    "asymmetric_migration",
                    "symmetric_migration",
                    "isolation",
                ],
                "akaike_weight": [0.8, 0.1, 0.05, 0.03, 0.02],
            }
        ),
    )

    result = forward_stress._fit_full_candidate_set(
        np.ones((5, 5)),
        starts=1,
        maxiter=10,
        seed=4,
    )

    assert result["directional_model"] is None
    assert result["preferred_direction"] == "none"
    assert result["asymmetric_model_weight"] == 0.0
