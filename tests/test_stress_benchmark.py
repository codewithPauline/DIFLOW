import pandas as pd
import pytest

import diflow.validation.stress as stress
from diflow.validation.stress import (
    default_stress_scenarios,
    summarize_stress,
)


def test_stress_suite_contains_core_misspecification_cases():
    names = {scenario.name for scenario in default_stress_scenarios()}
    assert "zero_migration" in names
    assert "unequal_ne_symmetric_migration" in names
    assert "secondary_contact_symmetric" in names
    assert "secondary_contact_a_to_b" in names
    assert "secondary_contact_b_to_a" in names


def test_stress_summary_counts_false_direction_signal():
    results = pd.DataFrame(
        {
            "scenario": ["zero"] * 4,
            "success": [True] * 4,
            "expected_direction": ["none"] * 4,
            "correct_model_selected": [True, True, False, True],
            "provisional_directional_signal": [False, False, True, False],
            "preferred_direction": ["symmetric", "symmetric", "A->B", "symmetric"],
        }
    )
    summary = summarize_stress(results)
    row = summary.iloc[0]
    assert row["correct_model_selection_rate"] == pytest.approx(0.75)
    assert row["false_direction_signal_rate"] == pytest.approx(0.25)


def test_stress_summary_scores_direction_recovery():
    results = pd.DataFrame(
        {
            "scenario": ["sc_ab"] * 4,
            "success": [True] * 4,
            "expected_direction": ["A->B"] * 4,
            "correct_model_selected": [True] * 4,
            "provisional_directional_signal": [True, True, False, True],
            "preferred_direction": ["A->B", "A->B", "A->B", "B->A"],
        }
    )
    summary = summarize_stress(results)
    assert summary.iloc[0]["direction_recovery_rate"] == pytest.approx(0.50)



def test_symmetric_secondary_contact_uses_symmetric_generator():
    scenario = next(
        x
        for x in default_stress_scenarios()
        if x.name == "secondary_contact_symmetric"
    )

    assert scenario.generating_model == "secondary_contact_symmetric"
    assert len(scenario.parameters) == 5


def test_stress_uses_selected_secondary_contact_direction(monkeypatch):
    scenario = next(
        x
        for x in default_stress_scenarios()
        if x.name == "secondary_contact_a_to_b"
    )

    monkeypatch.setattr(
        stress,
        "expected_candidate_spectrum",
        lambda *args, **kwargs: __import__("numpy").ones((5, 5)),
    )
    monkeypatch.setattr(
        stress,
        "sample_variant_spectrum",
        lambda *args, **kwargs: __import__("numpy").ones((5, 5)),
    )

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
            "m_a_to_b": 0.2,
            "m_b_to_a": 0.8,
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
            "m_a_to_b": 1.1,
            "m_b_to_a": 0.2,
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

    monkeypatch.setattr(stress, "fit_multistart", fake_fit)
    monkeypatch.setattr(
        stress,
        "rank_models",
        lambda scores, use_aicc=False: pd.DataFrame(
            {
                "model": [
                    "secondary_contact_asymmetric",
                    "asymmetric_migration",
                    "secondary_contact_symmetric",
                    "symmetric_migration",
                    "isolation",
                ],
                "akaike_weight": [0.82, 0.08, 0.05, 0.03, 0.02],
            }
        ),
    )

    result = stress.run_stress_benchmark(
        scenarios=(scenario,),
        replicates=1,
        sample_sizes=(4, 4),
        segregating_sites=10,
        starts=1,
        maxiter=10,
        min_model_weight=0.7,
        min_abs_asymmetry=0.25,
    )

    row = result.iloc[0]
    assert row["directional_model"] == "secondary_contact_asymmetric"
    assert row["estimated_m_a_to_b"] == 1.1
    assert row["estimated_m_b_to_a"] == 0.2
    assert row["preferred_direction"] == "A->B"
    assert bool(row["provisional_directional_signal"]) is True


def test_stress_emits_no_direction_for_symmetric_winner(monkeypatch):
    scenario = next(
        x
        for x in default_stress_scenarios()
        if x.name == "secondary_contact_symmetric"
    )

    monkeypatch.setattr(
        stress,
        "expected_candidate_spectrum",
        lambda *args, **kwargs: __import__("numpy").ones((5, 5)),
    )
    monkeypatch.setattr(
        stress,
        "sample_variant_spectrum",
        lambda *args, **kwargs: __import__("numpy").ones((5, 5)),
    )

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

    monkeypatch.setattr(stress, "fit_multistart", fake_fit)
    monkeypatch.setattr(
        stress,
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

    result = stress.run_stress_benchmark(
        scenarios=(scenario,),
        replicates=1,
        sample_sizes=(4, 4),
        segregating_sites=10,
        starts=1,
        maxiter=10,
        min_model_weight=0.7,
        min_abs_asymmetry=0.25,
    )

    row = result.iloc[0]
    assert row["directional_model"] is None
    assert row["preferred_direction"] == "none"
    assert bool(row["provisional_directional_signal"]) is False
