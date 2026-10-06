import json

import pandas as pd
import pytest

from diflow.validation.data_requirements import (
    DataRequirementTargets,
    pareto_minimum_regimes,
    recommend_data_requirements,
    summarize_data_regimes,
    write_data_requirements,
)


def _grid_summary():
    rows = []
    regimes = {
        (10, 1000): (0.80, 0.10, 0.99),
        (10, 5000): (0.92, 0.04, 0.99),
        (20, 1000): (0.91, 0.04, 0.98),
        (20, 5000): (0.96, 0.02, 0.99),
    }
    for (chromosomes, sites), (accuracy, fpr, success) in regimes.items():
        rows.extend(
            [
                {
                    "chromosomes": chromosomes,
                    "segregating_sites": sites,
                    "scenario": "symmetric",
                    "expected_direction": "symmetric",
                    "success_rate": success,
                    "direction_accuracy": 1.0,
                    "false_directional_positive_rate": fpr,
                },
                {
                    "chromosomes": chromosomes,
                    "segregating_sites": sites,
                    "scenario": "moderate_a_to_b",
                    "expected_direction": "A->B",
                    "success_rate": success,
                    "direction_accuracy": accuracy,
                    "false_directional_positive_rate": float("nan"),
                },
                {
                    "chromosomes": chromosomes,
                    "segregating_sites": sites,
                    "scenario": "moderate_b_to_a",
                    "expected_direction": "B->A",
                    "success_rate": success,
                    "direction_accuracy": accuracy,
                    "false_directional_positive_rate": float("nan"),
                },
            ]
        )
    return pd.DataFrame(rows)


def test_data_regime_summary_applies_targets():
    regimes = summarize_data_regimes(_grid_summary())
    lookup = regimes.set_index(["chromosomes", "segregating_sites"])
    assert bool(lookup.loc[(10, 1000), "passes_targets"]) is False
    assert bool(lookup.loc[(10, 5000), "passes_targets"]) is True
    assert bool(lookup.loc[(20, 1000), "passes_targets"]) is True


def test_pareto_minima_keep_tradeoff_regimes():
    regimes = summarize_data_regimes(_grid_summary())
    minima = pareto_minimum_regimes(regimes)
    pairs = set(zip(minima["chromosomes"], minima["segregating_sites"]))
    assert pairs == {(10, 5000), (20, 1000)}


def test_recommendation_reports_available_guidance():
    _, minima, recommendation = recommend_data_requirements(_grid_summary())
    assert recommendation["status"] == "recommendations_available"
    assert recommendation["passing_regimes"] == 3
    assert len(minima) == 2


def test_write_data_requirements_outputs_files(tmp_path):
    source = tmp_path / "grid.csv"
    _grid_summary().to_csv(source, index=False)
    _, _, recommendation = write_data_requirements(
        recovery_grid_summary_csv=source,
        output_dir=tmp_path / "requirements",
    )
    assert recommendation["status"] == "recommendations_available"
    assert (tmp_path / "requirements" / "data_regime_performance.csv").exists()
    assert (tmp_path / "requirements" / "minimum_passing_regimes.csv").exists()
    assert (tmp_path / "requirements" / "data_requirements.json").exists()
    assert (tmp_path / "requirements" / "data_requirements.md").exists()



def test_missing_symmetric_fpr_does_not_become_directional():
    frame = _grid_summary()
    mask = (
        (frame["chromosomes"] == 10)
        & (frame["segregating_sites"] == 5000)
        & (frame["expected_direction"] == "symmetric")
    )
    frame.loc[mask, "false_directional_positive_rate"] = float("nan")

    regimes = summarize_data_regimes(frame)
    row = regimes[
        (regimes["chromosomes"] == 10)
        & (regimes["segregating_sites"] == 5000)
    ].iloc[0]

    assert pd.isna(row["maximum_false_direction_rate"])
    assert bool(row["passes_targets"]) is False
