import pandas as pd

from diflow.inspection import (
    PRESETS,
    recommend_neighbors,
    recommend_projection,
    projection_retention_summary,
)


def test_presets_are_ordered_by_effort():
    assert PRESETS["quick"]["bootstrap_replicates"] < PRESETS["standard"]["bootstrap_replicates"]
    assert PRESETS["standard"]["bootstrap_replicates"] < PRESETS["publication"]["bootstrap_replicates"]


def test_projection_uses_largest_value_meeting_retention_target():
    summary = pd.DataFrame(
        {
            "projection_chromosomes": [2, 4, 6, 8],
            "minimum_pair_retention": [1.0, 0.95, 0.82, 0.60],
        }
    )
    assert recommend_projection(summary, retention_target=0.80) == 6


def test_neighbor_recommendation_scales_with_population_count():
    assert recommend_neighbors(2) == 1
    assert recommend_neighbors(8) == 3
    assert recommend_neighbors(69) == 4



def test_projection_summary_uses_shared_pairwise_loci():
    rows = [
        ["1", 1, "A", "G", "A", 0, 4],
        ["1", 1, "A", "G", "B", 0, 4],
        ["1", 2, "A", "G", "A", 0, 4],
        ["1", 2, "A", "G", "B", 0, 2],
        ["1", 3, "A", "G", "A", 0, 2],
        ["1", 3, "A", "G", "B", 0, 4],
    ]
    counts = pd.DataFrame(
        rows,
        columns=[
            "chrom", "pos", "ref", "alt", "population",
            "alt_count", "called_chromosomes",
        ],
    )
    summary = projection_retention_summary(counts)
    p4 = summary[summary["projection_chromosomes"] == 4].iloc[0]
    assert p4["minimum_pair_loci"] == 1
    assert p4["minimum_pair_retention"] == 1 / 3
