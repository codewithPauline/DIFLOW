import pandas as pd

from diflow.inspection import (
    PRESETS,
    recommend_neighbors,
    recommend_projection,
)


def test_presets_are_ordered_by_effort():
    assert PRESETS["quick"]["bootstrap_replicates"] < PRESETS["standard"]["bootstrap_replicates"]
    assert PRESETS["standard"]["bootstrap_replicates"] < PRESETS["publication"]["bootstrap_replicates"]


def test_projection_uses_largest_value_meeting_retention_target():
    summary = pd.DataFrame(
        {
            "projection_chromosomes": [2, 4, 6, 8],
            "minimum_population_retention": [1.0, 0.95, 0.82, 0.60],
        }
    )
    assert recommend_projection(summary, retention_target=0.80) == 6


def test_neighbor_recommendation_scales_with_population_count():
    assert recommend_neighbors(2) == 1
    assert recommend_neighbors(8) == 3
    assert recommend_neighbors(69) == 4
