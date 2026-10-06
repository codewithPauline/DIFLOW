from pathlib import Path

import pandas as pd

from diflow.validation.grid import (
    default_recovery_grid,
    plot_recovery_grid,
)


def test_default_recovery_grid_spans_data_sizes_and_directions():
    grid = default_recovery_grid()
    assert len(grid) == 63
    assert {x.chromosomes for x in grid} == {10, 20, 40}
    assert {x.segregating_sites for x in grid} == {1000, 5000, 20000}
    directions = {x.scenario.expected_direction for x in grid}
    assert {"symmetric", "A->B", "B->A"}.issubset(directions)


def test_recovery_grid_plots_are_written(tmp_path):
    summary = pd.DataFrame(
        {
            "chromosomes": [10, 10, 20, 20],
            "segregating_sites": [1000, 5000, 1000, 5000],
            "scenario": ["symmetric", "moderate_a_to_b", "symmetric", "moderate_a_to_b"],
            "direction_accuracy": [1.0, 0.7, 1.0, 0.8],
            "false_directional_positive_rate": [0.2, float("nan"), 0.1, float("nan")],
        }
    )
    paths = plot_recovery_grid(summary, tmp_path)
    names = {Path(p).name for p in paths}
    assert "direction_accuracy_grid.png" in names
    assert "direction_accuracy_grid.pdf" in names
    assert "false_direction_rate_grid.png" in names
    assert "false_direction_rate_grid.pdf" in names
    assert all(Path(p).exists() for p in paths)
