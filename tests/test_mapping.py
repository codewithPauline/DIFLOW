import matplotlib
matplotlib.use("Agg")

import pandas as pd
import pytest

from diflow.mapping import plot_directional_map, validate_coordinates, validate_flows


def _coordinates():
    return pd.DataFrame(
        {
            "population": ["A", "B", "C"],
            "latitude": [39.0, 39.5, 38.8],
            "longitude": [-84.0, -83.2, -82.7],
        }
    )


def test_validates_and_plots_directional_edges():
    flows = pd.DataFrame(
        {
            "source": ["A", "B"],
            "destination": ["B", "A"],
            "migration": [0.03, 0.005],
            "support": [0.99, 0.72],
        }
    )
    ax = plot_directional_map(_coordinates(), flows)
    assert len(ax.patches) == 2
    assert len(ax.collections) == 1


def test_support_filter_removes_weak_edge():
    flows = pd.DataFrame(
        {
            "source": ["A", "B"],
            "destination": ["B", "A"],
            "migration": [0.03, 0.005],
            "support": [0.99, 0.72],
        }
    )
    ax = plot_directional_map(_coordinates(), flows, min_support=0.95)
    assert len(ax.patches) == 1


def test_unknown_population_rejected():
    flows = pd.DataFrame(
        {"source": ["A"], "destination": ["Z"], "migration": [0.02]}
    )
    with pytest.raises(ValueError):
        plot_directional_map(_coordinates(), flows)


def test_coordinate_bounds_checked():
    bad = _coordinates()
    bad.loc[0, "latitude"] = 100
    with pytest.raises(ValueError):
        validate_coordinates(bad)


def test_negative_migration_rejected():
    flows = pd.DataFrame(
        {"source": ["A"], "destination": ["B"], "migration": [-0.1]}
    )
    with pytest.raises(ValueError):
        validate_flows(flows)



def test_projected_map_coordinates():
    flows = pd.DataFrame(
        {
            "source": ["A"],
            "destination": ["B"],
            "migration": [0.03],
            "support": [0.99],
        }
    )
    ax = plot_directional_map(
        _coordinates(),
        flows,
        target_crs="EPSG:5070",
    )
    assert "EPSG:5070" in ax.get_xlabel()
    assert len(ax.patches) == 1
