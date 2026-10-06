import pandas as pd
import pytest

from diflow.validation.mechanistic_linkage import (
    MechanisticLinkageScenario,
    default_mechanistic_linkage_scenarios,
    summarize_mechanistic_linkage,
)


def test_mechanistic_scenarios_cover_symmetry_and_both_directions():
    directions = {
        scenario.expected_direction
        for scenario in default_mechanistic_linkage_scenarios()
    }
    assert directions == {"symmetric", "A->B", "B->A"}


def test_mechanistic_summary_reports_false_direction_under_symmetry():
    results = pd.DataFrame(
        {
            "scenario": ["sym"] * 4,
            "method": ["locus", "locus", "block_100000", "block_100000"],
            "success": [True] * 4,
            "expected_direction": ["symmetric"] * 4,
            "coverage_a_to_b": [False, True, True, True],
            "coverage_b_to_a": [False, True, True, True],
            "interval_width_a_to_b": [0.2, 0.2, 0.4, 0.4],
            "interval_width_b_to_a": [0.2, 0.2, 0.4, 0.4],
            "strong_directional_support": [True, False, False, False],
        }
    )
    summary = summarize_mechanistic_linkage(results)
    locus = summary[summary["method"] == "locus"].iloc[0]
    block = summary[summary["method"] == "block_100000"].iloc[0]
    assert locus["false_directional_support_rate"] == pytest.approx(0.5)
    assert block["false_directional_support_rate"] == pytest.approx(0.0)


def test_msprime_dependency_is_explicit_when_unavailable(monkeypatch):
    import builtins
    from diflow.validation import mechanistic_linkage

    original_import = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name == "msprime":
            raise ImportError("blocked for test")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked)
    with pytest.raises(ImportError, match="validation"):
        mechanistic_linkage._require_msprime()
