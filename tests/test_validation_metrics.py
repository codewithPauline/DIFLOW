import pytest

from diflow.validation import (
    core_validation_scenarios,
    direction_accuracy,
    false_directional_positive_rate,
    interval_coverage,
    parameter_bias,
    parameter_rmse,
)


def test_bias_and_rmse():
    truth = [1.0, 1.0, 1.0]
    estimates = [0.8, 1.0, 1.2]
    assert parameter_bias(truth, estimates) == pytest.approx(0.0)
    assert parameter_rmse(truth, estimates) > 0


def test_interval_coverage():
    truth = [1.0, 2.0, 3.0]
    lower = [0.5, 2.1, 2.5]
    upper = [1.5, 2.5, 3.5]
    assert interval_coverage(truth, lower, upper) == pytest.approx(2 / 3)


def test_direction_accuracy():
    accuracy = direction_accuracy(
        [3, 1, 1],
        [1, 3, 1],
        [2.5, 1.2, 1.0],
        [1.0, 2.5, 1.0],
    )
    assert accuracy == pytest.approx(1.0)


def test_false_directional_positive_rate():
    rate = false_directional_positive_rate(
        [1.0, 1.1, 3.0, 0.9],
        [1.0, 1.0, 1.0, 1.0],
        asymmetry_threshold=0.25,
    )
    assert rate == pytest.approx(0.25)


def test_core_scenarios_cover_major_confounders():
    names = {scenario.name for scenario in core_validation_scenarios()}
    assert "symmetric_migration" in names
    assert "secondary_contact" in names
    assert "range_expansion" in names
    assert "ghost_population" in names
