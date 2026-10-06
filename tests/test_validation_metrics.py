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



@pytest.mark.parametrize(
    "truth,lower,upper",
    [
        ([1.0, 2.0], [0.5, float("nan")], [1.5, 2.5]),
        ([1.0, 2.0], [0.5, 1.5], [1.5, float("inf")]),
        ([1.0, float("nan")], [0.5, 1.5], [1.5, 2.5]),
    ],
)
def test_interval_coverage_rejects_nonfinite_bounds(truth, lower, upper):
    with pytest.raises(ValueError, match="finite"):
        interval_coverage(truth, lower, upper)



@pytest.mark.parametrize(
    "true_ab,true_ba,est_ab,est_ba",
    [
        ([1.0], [0.5], [float("nan")], [0.4]),
        ([1.0], [0.5], [0.8], [float("inf")]),
        ([float("nan")], [0.5], [0.8], [0.4]),
    ],
)
def test_direction_accuracy_rejects_nonfinite_values(
    true_ab,
    true_ba,
    est_ab,
    est_ba,
):
    with pytest.raises(ValueError, match="finite"):
        direction_accuracy(true_ab, true_ba, est_ab, est_ba)


@pytest.mark.parametrize("tolerance", [-0.1, float("nan"), float("inf")])
def test_direction_accuracy_rejects_invalid_tolerance(tolerance):
    with pytest.raises(ValueError, match="tolerance"):
        direction_accuracy([1.0], [0.5], [0.9], [0.4], tolerance=tolerance)



@pytest.mark.parametrize(
    "estimated_ab,estimated_ba,error_match",
    [
        ([float("nan")], [1.0], "finite"),
        ([1.0], [float("inf")], "finite"),
        ([-0.1], [1.0], "non-negative"),
        ([1.0], [-0.1], "non-negative"),
    ],
)
def test_false_direction_rate_rejects_invalid_estimates(
    estimated_ab,
    estimated_ba,
    error_match,
):
    with pytest.raises(ValueError, match=error_match):
        false_directional_positive_rate(
            estimated_ab,
            estimated_ba,
            asymmetry_threshold=0.25,
        )
