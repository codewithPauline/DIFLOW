import pandas as pd
import pytest

from diflow.validation.comparators import (
    compare_method_files,
    summarize_method_comparison,
    summarize_direction_only_comparison,
    standardize_diflow_benchmark,
    validate_matched_comparison,
    validate_comparator_table,
)


def _table():
    return pd.DataFrame(
        {
            "scenario": ["sym", "ab"],
            "replicate": [1, 1],
            "truth_m_a_to_b": [0.5, 1.0],
            "truth_m_b_to_a": [0.5, 0.25],
            "estimated_m_a_to_b": [0.55, 0.9],
            "estimated_m_b_to_a": [0.45, 0.3],
        }
    )


def test_comparator_table_is_labeled():
    out = validate_comparator_table(_table(), method="DIFLOW")
    assert set(out["method"]) == {"DIFLOW"}


def test_method_summary_has_direction_and_false_positive_metrics():
    a = validate_comparator_table(_table(), method="DIFLOW")
    b = validate_comparator_table(_table(), method="Other")
    summary = summarize_method_comparison(pd.concat([a, b], ignore_index=True))
    assert set(summary["method"]) == {"DIFLOW", "Other"}
    sym = summary[summary["scenario"] == "sym"]
    assert sym["false_directional_positive_rate"].notna().all()


def test_comparator_requires_standard_columns():
    with pytest.raises(ValueError):
        validate_comparator_table(pd.DataFrame({"scenario": ["x"]}), method="bad")



def test_matched_comparison_rejects_missing_replicates():
    a = validate_comparator_table(_table(), method="DIFLOW")
    b = validate_comparator_table(_table().iloc[[0]], method="Other")
    combined = pd.concat([a, b], ignore_index=True)

    with pytest.raises(ValueError, match="not matched"):
        validate_matched_comparison(combined)


def test_matched_comparison_rejects_inconsistent_truth():
    a = validate_comparator_table(_table(), method="DIFLOW")
    altered = _table().copy()
    altered.loc[0, "truth_m_a_to_b"] = 0.9
    b = validate_comparator_table(altered, method="Other")

    with pytest.raises(ValueError, match="disagree on simulation truth"):
        validate_matched_comparison(pd.concat([a, b], ignore_index=True))


def test_matched_comparison_rejects_mixed_estimands():
    a = validate_comparator_table(_table(), method="DIFLOW")
    b = validate_comparator_table(_table(), method="Other")
    a["estimand"] = "dadi_scaled_migration"
    b["estimand"] = "migrants_per_generation"

    with pytest.raises(ValueError, match="multiple estimands"):
        validate_matched_comparison(pd.concat([a, b], ignore_index=True))



def test_standardize_diflow_benchmark_maps_truth_columns():
    frame = pd.DataFrame(
        {
            "scenario": ["ab"],
            "replicate": [1],
            "true_m_a_to_b": [1.0],
            "true_m_b_to_a": [0.25],
            "estimated_m_a_to_b": [0.9],
            "estimated_m_b_to_a": [0.3],
        }
    )
    out = standardize_diflow_benchmark(frame)
    assert "truth_m_a_to_b" in out.columns
    assert "truth_m_b_to_a" in out.columns
    assert "method_version" in out.columns
    assert str(out.iloc[0]["method_version"]).strip()
    assert out.iloc[0]["estimand"] == "dadi_scaled_migration"



def test_direction_only_summary_allows_different_estimands():
    a = validate_comparator_table(_table(), method="DIFLOW")
    b = validate_comparator_table(_table(), method="Other")
    a["estimand"] = "dadi_scaled_migration"
    b["estimand"] = "effective_lineage_migration"
    combined = pd.concat([a, b], ignore_index=True)

    # Direction-only comparison intentionally ignores scale incompatibility.
    stripped = combined.drop(columns=["estimand"])
    matched = validate_matched_comparison(stripped)
    summary = summarize_direction_only_comparison(matched)
    assert set(summary["method"]) == {"DIFLOW", "Other"}
    assert "bias_m_a_to_b" not in summary.columns
    assert "direction_accuracy" in summary.columns



def test_method_comparison_writes_provenance(tmp_path):
    a_path = tmp_path / "a.csv"
    b_path = tmp_path / "b.csv"
    _table().assign(
        estimand="dadi_scaled_migration",
        method_version="0.0.1",
    ).to_csv(a_path, index=False)
    _table().assign(
        estimand="dadi_scaled_migration",
        method_version="1.2.3",
    ).to_csv(b_path, index=False)

    compare_method_files(
        {"DIFLOW": a_path, "Other": b_path},
        output_dir=tmp_path / "out",
    )

    import json
    metadata = json.loads(
        (tmp_path / "out" / "comparison_metadata.json").read_text(
            encoding="utf-8"
        )
    )
    assert metadata["direction_only"] is False
    assert metadata["require_complete_match"] is True
    assert metadata["methods"] == ["DIFLOW", "Other"]
    assert metadata["estimands"] == ["dadi_scaled_migration"]
    assert metadata["method_versions"]["DIFLOW"] == ["0.0.1"]
    assert metadata["method_versions"]["Other"] == ["1.2.3"]



def test_comparator_rejects_duplicate_scenario_replicate_rows():
    duplicated = pd.concat([_table(), _table().iloc[[0]]], ignore_index=True)

    with pytest.raises(ValueError, match="duplicate .*scenario, replicate"):
        validate_comparator_table(duplicated, method="DIFLOW")
