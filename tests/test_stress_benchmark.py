import pandas as pd
import pytest

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
