import json

import pandas as pd
import pytest

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

    data_requirements = root / "data_requirements"
    data_requirements.mkdir()
    (data_requirements / "data_requirements.json").write_text(
        json.dumps(
            {
                "status": "recommendations_available",
                "passing_regimes": 2,
                "pareto_minimum_regimes": [
                    {"chromosomes": 10, "segregating_sites": 5000}
                ],
            }
        ),
        encoding="utf-8",
    )

    comparison = root / "external_comparison"
    comparison.mkdir()
    pd.DataFrame(
        {
            "method": ["DIFLOW", "ExternalMethod"],
            "scenario": ["moderate_a_to_b", "moderate_a_to_b"],
            "replicates": [50, 50],
        }
    ).to_csv(comparison / "method_comparison_summary.csv", index=False)
    (comparison / "comparison_metadata.json").write_text(
        json.dumps(
            {
                "direction_only": False,
                "require_complete_match": True,
                "asymmetry_threshold": 0.25,
                "methods": ["DIFLOW", "ExternalMethod"],
                "estimands": ["dadi_scaled_migration"],
                "method_versions": {
                    "DIFLOW": ["0.0.1"],
                    "ExternalMethod": ["1.2.3"],
                },
                "magnitude_metrics_reported": True,
            }
        ),
        encoding="utf-8",
    )


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



def test_release_review_flags_missing_comparison_provenance(tmp_path):
    root = tmp_path / "results"
    _write_valid_campaign(root)
    (root / "external_comparison" / "comparison_metadata.json").unlink()

    _, summary = review_validation_campaign(root)

    assert summary["release_ready"] is False
    assert any(
        blocker["metric"] == "comparison_provenance"
        for blocker in summary["blockers"]
    )



def test_release_review_flags_missing_external_version(tmp_path):
    root = tmp_path / "results"
    _write_valid_campaign(root)
    metadata_path = root / "external_comparison" / "comparison_metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["method_versions"] = {"DIFLOW": ["0.0.1"]}
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")

    _, summary = review_validation_campaign(root)

    assert summary["release_ready"] is False
    assert any(
        blocker["metric"] == "external_method_versions"
        for blocker in summary["blockers"]
    )



@pytest.mark.parametrize(
    "column,value",
    [
        ("false_directional_positive_rate", float("nan")),
        ("directional_sensitivity", 1.2),
        ("direction_accuracy_when_called", -0.1),
    ],
)
def test_release_review_rejects_invalid_threshold_metrics(
    tmp_path,
    column,
    value,
):
    root = tmp_path / "results"
    _write_valid_campaign(root)
    threshold_path = (
        root / "threshold_calibration" / "selected_thresholds.csv"
    )
    selected = pd.read_csv(threshold_path)
    selected.loc[0, column] = value
    selected.to_csv(threshold_path, index=False)

    with pytest.raises(ValueError, match=column):
        review_validation_campaign(root)
