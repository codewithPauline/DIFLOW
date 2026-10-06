from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from diflow.demography import AsymmetricIMParams
from diflow.validation.linked_calibration import (
    LinkedCalibrationScenario,
    simulate_correlated_block_counts,
    run_linked_bootstrap_calibration,
    summarize_linked_bootstrap_calibration,
)


def fake_spectrum(params, sample_sizes):
    shape = (sample_sizes[0] + 1, sample_sizes[1] + 1)
    arr = np.ones(shape, dtype=float)
    arr[0, 0] = 0.0
    arr[-1, -1] = 0.0
    return arr


def fake_fit(spectrum, model_name, *, starts, seed, maxiter):
    total = max(float(spectrum.sum()), 1.0)
    axis_a = sum(i * spectrum[i, :].sum() for i in range(spectrum.shape[0]))
    axis_b = sum(j * spectrum[:, j].sum() for j in range(spectrum.shape[1]))
    return SimpleNamespace(
        best_parameters={
            "m_a_to_b": 0.5 + 0.01 * axis_b / total,
            "m_b_to_a": 0.5 + 0.01 * axis_a / total,
        }
    )


def test_correlated_block_counts_have_expected_structure():
    params = AsymmetricIMParams(1.0, 1.0, 0.5, 1.0, 0.25)
    counts = simulate_correlated_block_counts(
        params,
        sample_sizes=(6, 8),
        blocks=4,
        snps_per_block=3,
        block_size_bp=1000,
        concentration=10,
        seed=3,
        spectrum_function=fake_spectrum,
    )
    assert len(counts) == 4 * 3 * 2
    assert counts["pos"].nunique() == 12
    assert set(counts["population"]) == {"A", "B"}


def test_linked_calibration_compares_locus_and_block_methods():
    scenario = LinkedCalibrationScenario(
        "test",
        AsymmetricIMParams(1.0, 1.0, 0.5, 1.0, 0.25),
        "A->B",
    )
    results = run_linked_bootstrap_calibration(
        scenarios=(scenario,),
        replicates=2,
        sample_sizes=(6, 6),
        blocks=4,
        snps_per_block=3,
        block_size_bp=1000,
        concentration=10,
        bootstrap_replicates=6,
        bootstrap_starts=1,
        seed=7,
        spectrum_function=fake_spectrum,
        fit_function=fake_fit,
    )
    assert len(results) == 4
    assert set(results["method"]) == {"locus", "block"}
    assert results["success"].all()


def test_linked_summary_reports_false_direction_by_method():
    results = pd.DataFrame(
        {
            "scenario": ["sym"] * 4,
            "method": ["locus", "locus", "block", "block"],
            "success": [True] * 4,
            "expected_direction": ["symmetric"] * 4,
            "coverage_a_to_b": [True, False, True, True],
            "coverage_b_to_a": [True, False, True, True],
            "interval_width_a_to_b": [0.2, 0.2, 0.4, 0.4],
            "interval_width_b_to_a": [0.2, 0.2, 0.4, 0.4],
            "strong_directional_support": [True, False, False, False],
        }
    )
    summary = summarize_linked_bootstrap_calibration(results)
    locus = summary[summary["method"] == "locus"].iloc[0]
    block = summary[summary["method"] == "block"].iloc[0]
    assert locus["coverage_a_to_b"] == pytest.approx(0.5)
    assert block["coverage_a_to_b"] == pytest.approx(1.0)
    assert locus["false_directional_support_rate"] == pytest.approx(0.5)
    assert block["false_directional_support_rate"] == pytest.approx(0.0)
