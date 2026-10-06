import pandas as pd
import pytest

from diflow.validation.thresholds import (
    calibrate_thresholds,
    load_calibrated_thresholds,
    evaluate_thresholds,
    select_thresholds,
)


def _evidence():
    return pd.DataFrame(
        {
            "truth_direction": [
                "symmetric", "symmetric", "A->B", "A->B", "B->A", "B->A"
            ],
            "preferred_direction": [
                "A->B", "B->A", "A->B", "A->B", "B->A", "A->B"
            ],
            "asymmetric_model_weight": [0.8, 0.4, 0.9, 0.75, 0.95, 0.95],
            "directional_support": [0.96, 0.99, 0.99, 0.96, 0.99, 0.99],
            "asymmetry_index": [0.30, -0.40, 0.60, 0.35, -0.70, 0.70],
            "optimizer_stable": [True] * 6,
            "interval_separated": [True, True, True, True, True, True],
        }
    )


def test_threshold_evaluation_counts_false_direction_and_sensitivity():
    metrics = evaluate_thresholds(
        _evidence(),
        min_model_weight=0.7,
        min_directional_support=0.95,
        min_abs_asymmetry=0.25,
    )
    assert metrics["false_directional_positive_rate"] == pytest.approx(0.5)
    assert metrics["directional_sensitivity"] == pytest.approx(0.75)
    assert metrics["direction_accuracy_when_called"] == pytest.approx(0.75)


def test_interval_separation_is_required_for_a_call():
    evidence = _evidence()
    evidence.loc[2, "interval_separated"] = False
    metrics = evaluate_thresholds(
        evidence,
        min_model_weight=0.7,
        min_directional_support=0.95,
        min_abs_asymmetry=0.25,
    )
    assert metrics["directional_sensitivity"] == pytest.approx(0.5)


def test_calibration_selects_rule_meeting_fpr_target():
    selected, scan = calibrate_thresholds(
        _evidence(),
        max_false_directional_positive_rate=0.0,
        model_weights=(0.7, 0.9),
        directional_supports=(0.95,),
        asymmetries=(0.25,),
    )
    assert len(scan) == 2
    assert selected.false_directional_positive_rate == pytest.approx(0.0)
    assert selected.min_model_weight == pytest.approx(0.9)



def test_load_calibrated_thresholds_csv(tmp_path):
    path = tmp_path / "selected_thresholds.csv"
    pd.DataFrame(
        [
            {
                "min_model_weight": 0.8,
                "min_directional_support": 0.975,
                "min_abs_asymmetry": 0.3,
                "false_directional_positive_rate": 0.04,
                "directional_sensitivity": 0.82,
                "direction_accuracy_when_called": 0.96,
                "called_fraction": 0.71,
                "evaluated_rows": 500,
            }
        ]
    ).to_csv(path, index=False)

    selected = load_calibrated_thresholds(path)
    assert selected.min_model_weight == pytest.approx(0.8)
    assert selected.min_directional_support == pytest.approx(0.975)
    assert selected.min_abs_asymmetry == pytest.approx(0.3)
    assert selected.evaluated_rows == 500


def test_load_calibrated_thresholds_rejects_multiple_rows(tmp_path):
    path = tmp_path / "bad.csv"
    pd.DataFrame(
        {
            "min_model_weight": [0.7, 0.8],
            "min_directional_support": [0.95, 0.975],
            "min_abs_asymmetry": [0.25, 0.3],
        }
    ).to_csv(path, index=False)

    with pytest.raises(ValueError, match="exactly one"):
        load_calibrated_thresholds(path)



def test_string_false_is_not_treated_as_true():
    evidence = _evidence()
    evidence["optimizer_stable"] = ["False", "True", "True", "True", "True", "True"]

    metrics = evaluate_thresholds(
        evidence,
        min_model_weight=0.7,
        min_directional_support=0.95,
        min_abs_asymmetry=0.25,
    )

    # The first symmetric row would be a false positive if "False" were
    # incorrectly coerced with bool("False") == True.
    assert metrics["false_directional_positive_rate"] == pytest.approx(0.0)


def test_invalid_boolean_text_is_rejected():
    evidence = _evidence()
    evidence["interval_separated"] = evidence["interval_separated"].astype(object)
    evidence.loc[0, "interval_separated"] = "maybe"

    with pytest.raises(ValueError, match="interval_separated"):
        evaluate_thresholds(
            evidence,
            min_model_weight=0.7,
            min_directional_support=0.95,
            min_abs_asymmetry=0.25,
        )



def test_invalid_preferred_direction_is_rejected():
    evidence = _evidence()
    evidence.loc[0, "preferred_direction"] = "A=>B"

    with pytest.raises(ValueError, match="preferred_direction"):
        evaluate_thresholds(
            evidence,
            min_model_weight=0.7,
            min_directional_support=0.95,
            min_abs_asymmetry=0.25,
        )
