import pandas as pd
import pytest

from diflow.validation.comparators import (
    summarize_method_comparison,
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
