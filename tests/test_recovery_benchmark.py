import numpy as np
import pandas as pd
import pytest

from diflow.demography import AsymmetricIMParams
from diflow.validation.recovery import (
    RecoveryScenario,
    default_recovery_scenarios,
    run_recovery_benchmark,
    simulate_model_consistent_spectrum,
    summarize_recovery,
)


class FakeFit:
    def __init__(self):
        self.best_parameters = {"m_a_to_b": 1.0, "m_b_to_a": 0.25}
        self.stable = True
        self.converged_fraction = 1.0
        self.best_log_likelihood = -10.0


def fake_spectrum(params, sample_sizes):
    shape = (sample_sizes[0] + 1, sample_sizes[1] + 1)
    arr = np.ones(shape, dtype=float)
    arr[0, 0] = 0.0
    arr[-1, -1] = 0.0
    return arr


def fake_fit(*args, **kwargs):
    return FakeFit()


def test_default_recovery_suite_has_symmetry_and_both_directions():
    directions = {x.expected_direction for x in default_recovery_scenarios()}
    assert {"symmetric", "A->B", "B->A"}.issubset(directions)


def test_model_consistent_sampler_returns_requested_site_count():
    params = AsymmetricIMParams(1.0, 1.0, 0.5, 1.0, 0.25)
    fs = simulate_model_consistent_spectrum(
        params,
        sample_sizes=(6, 8),
        segregating_sites=250,
        seed=4,
        spectrum_function=fake_spectrum,
    )
    assert fs.shape == (7, 9)
    assert fs.sum() == 250
    assert fs[0, 0] == 0
    assert fs[-1, -1] == 0


def test_recovery_runner_records_truth_and_fit():
    scenario = RecoveryScenario(
        "test_ab",
        AsymmetricIMParams(1.0, 1.0, 0.5, 1.0, 0.25),
        "A->B",
    )
    result = run_recovery_benchmark(
        scenarios=(scenario,),
        replicates=3,
        sample_sizes=(6, 6),
        segregating_sites=100,
        starts=2,
        seed=9,
        spectrum_function=fake_spectrum,
        fit_function=fake_fit,
    )
    assert len(result) == 3
    assert result["success"].all()
    assert (result["true_m_a_to_b"] == 1.0).all()
    assert (result["estimated_m_b_to_a"] == 0.25).all()


def test_summary_reports_false_positive_rate_under_symmetry():
    results = pd.DataFrame(
        {
            "scenario": ["sym"] * 4,
            "success": [True] * 4,
            "true_m_a_to_b": [1.0] * 4,
            "true_m_b_to_a": [1.0] * 4,
            "estimated_m_a_to_b": [1.0, 1.1, 3.0, 0.9],
            "estimated_m_b_to_a": [1.0, 1.0, 1.0, 1.0],
            "optimizer_stable": [True] * 4,
        }
    )
    summary = summarize_recovery(results, asymmetry_threshold=0.25)
    assert summary.iloc[0]["false_directional_positive_rate"] == pytest.approx(0.25)



def test_recovery_summary_rejects_mislabeled_direction():
    results = pd.DataFrame(
        {
            "scenario": ["bad", "bad"],
            "expected_direction": ["B->A", "B->A"],
            "success": [True, True],
            "true_m_a_to_b": [1.0, 1.0],
            "true_m_b_to_a": [0.25, 0.25],
            "estimated_m_a_to_b": [0.9, 1.1],
            "estimated_m_b_to_a": [0.3, 0.2],
            "optimizer_stable": [True, True],
        }
    )
    with pytest.raises(ValueError, match="disagrees"):
        summarize_recovery(results)



def test_recovery_summary_does_not_treat_string_false_as_success():
    results = pd.DataFrame(
        {
            "scenario": ["sym", "sym"],
            "success": ["False", "True"],
            "true_m_a_to_b": [1.0, 1.0],
            "true_m_b_to_a": [1.0, 1.0],
            "estimated_m_a_to_b": [4.0, 1.0],
            "estimated_m_b_to_a": [1.0, 1.0],
            "optimizer_stable": [False, True],
        }
    )

    summary = summarize_recovery(results, asymmetry_threshold=0.25)

    row = summary.iloc[0]
    assert row["attempted_replicates"] == 2
    assert row["successful_replicates"] == 1
    assert row["success_rate"] == pytest.approx(0.5)
    assert row["false_directional_positive_rate"] == pytest.approx(0.0)



@pytest.mark.parametrize(
    "column,value,error_match",
    [
        ("true_m_a_to_b", float("nan"), "true_m_a_to_b"),
        ("true_m_b_to_a", float("inf"), "true_m_b_to_a"),
        ("true_m_a_to_b", -0.1, "true_m_a_to_b"),
    ],
)
def test_recovery_summary_rejects_invalid_truth_migration(
    column,
    value,
    error_match,
):
    results = pd.DataFrame(
        {
            "scenario": ["bad"],
            "success": [True],
            "true_m_a_to_b": [1.0],
            "true_m_b_to_a": [0.25],
            "estimated_m_a_to_b": [0.9],
            "estimated_m_b_to_a": [0.3],
            "optimizer_stable": [True],
        }
    )
    results.loc[0, column] = value

    with pytest.raises(ValueError, match=error_match):
        summarize_recovery(results)



@pytest.mark.parametrize(
    "column,value,error_match",
    [
        ("estimated_m_a_to_b", float("nan"), "estimated_m_a_to_b"),
        ("estimated_m_b_to_a", float("inf"), "estimated_m_b_to_a"),
        ("estimated_m_a_to_b", -0.1, "estimated_m_a_to_b"),
    ],
)
def test_recovery_summary_rejects_invalid_successful_estimate(
    column,
    value,
    error_match,
):
    results = pd.DataFrame(
        {
            "scenario": ["directional"],
            "success": [True],
            "true_m_a_to_b": [1.0],
            "true_m_b_to_a": [0.25],
            "estimated_m_a_to_b": [0.9],
            "estimated_m_b_to_a": [0.3],
            "optimizer_stable": [True],
        }
    )
    results.loc[0, column] = value

    with pytest.raises(ValueError, match=error_match):
        summarize_recovery(results)
