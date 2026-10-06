import json

import pandas as pd

from diflow.validation.release_review import (
    ReleaseCriteria,
    review_validation_campaign,
)


def _write_valid_campaign(root):
    recovery = root / "recovery_grid"
    recovery.mkdir(parents=True)
    pd.DataFrame(
        {
            "scenario": ["symmetric", "moderate_a_to_b", "moderate_b_to_a"],
            "success_rate": [1.0, 0.98, 0.99],
            "direction_accuracy": [1.0, 0.95, 0.96],
            "false_directional_positive_rate": [0.02, float("nan"), float("nan")],
        }
    ).to_csv(recovery / "recovery_grid_summary.csv", index=False)

    threshold = root / "threshold_calibration"
    threshold.mkdir()
    pd.DataFrame(
        [
            {
                "false_directional_positive_rate": 0.03,
                "directional_sensitivity": 0.85,
                "direction_accuracy_when_called": 0.97,
            }
        ]
    ).to_csv(threshold / "selected_thresholds.csv", index=False)

    linkage = root / "mechanistic_linkage"
    linkage.mkdir()
    pd.DataFrame(
        {
            "scenario": ["symmetric", "moderate_a_to_b"],
            "method": ["block_100000", "block_100000"],
            "success_rate": [0.98, 0.97],
            "coverage_a_to_b": [0.95, 0.94],
            "coverage_b_to_a": [0.96, 0.95],
        }
    ).to_csv(linkage / "mechanistic_grid_summary.csv", index=False)

    comparison = root / "external_comparison"
    comparison.mkdir()
    pd.DataFrame(
        {
            "method": ["DIFLOW", "ExternalMethod"],
            "scenario": ["moderate_a_to_b", "moderate_a_to_b"],
            "replicates": [50, 50],
        }
    ).to_csv(comparison / "method_comparison_summary.csv", index=False)


def test_release_review_passes_complete_campaign(tmp_path):
    root = tmp_path / "results"
    _write_valid_campaign(root)

    table, summary = review_validation_campaign(root)

    assert summary["release_ready"] is True
    assert summary["checks_failed_or_missing"] == 0
    assert table["passed"].fillna(False).all()
    assert (root / "release_review" / "release_review.csv").exists()
    assert (root / "release_review" / "release_review.json").exists()
    assert (root / "release_review" / "release_review.md").exists()


def test_release_review_flags_missing_empirical_study(tmp_path):
    root = tmp_path / "results"
    _write_valid_campaign(root)
    (root / "external_comparison" / "method_comparison_summary.csv").unlink()

    _, summary = review_validation_campaign(root)

    assert summary["release_ready"] is False
    assert any(
        blocker["section"] == "external_comparison"
        for blocker in summary["blockers"]
    )


def test_release_review_respects_stricter_criteria(tmp_path):
    root = tmp_path / "results"
    _write_valid_campaign(root)

    _, summary = review_validation_campaign(
        root,
        criteria=ReleaseCriteria(
            min_direction_accuracy=0.99,
            max_false_direction_rate=0.01,
            min_directional_sensitivity=0.95,
            min_ci_coverage=0.95,
            max_ci_coverage=0.99,
            min_success_rate=0.99,
        ),
    )

    assert summary["release_ready"] is False
    assert summary["checks_failed_or_missing"] > 0



def test_release_review_accepts_legacy_linkage_filename(tmp_path):
    root = tmp_path / "results"
    _write_valid_campaign(root)
    grid_path = root / "mechanistic_linkage" / "mechanistic_grid_summary.csv"
    legacy_path = root / "mechanistic_linkage" / "mechanistic_linkage_summary.csv"
    grid_path.rename(legacy_path)

    _, summary = review_validation_campaign(root)
    assert summary["release_ready"] is True
