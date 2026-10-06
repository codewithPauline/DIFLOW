import matplotlib
matplotlib.use("Agg")

import pandas as pd

from diflow.mapping import plot_directional_map


def coordinates():
    return pd.DataFrame(
        {
            "population": ["A", "B", "C"],
            "latitude": [39.0, 39.5, 38.8],
            "longitude": [-84.0, -83.2, -82.7],
        }
    )


def test_supported_edges_draw_by_default():
    flows = pd.DataFrame(
        {
            "source": ["A", "B"],
            "destination": ["B", "C"],
            "migration": [0.03, 0.02],
            "support": [0.99, 0.90],
            "status": ["supported", "ambiguous"],
        }
    )
    ax = plot_directional_map(coordinates(), flows)
    assert len(ax.patches) == 1


def test_ambiguous_edges_can_be_included():
    flows = pd.DataFrame(
        {
            "source": ["A", "B"],
            "destination": ["B", "C"],
            "migration": [0.03, 0.02],
            "support": [0.99, 0.90],
            "status": ["supported", "ambiguous"],
        }
    )
    ax = plot_directional_map(coordinates(), flows, include_ambiguous=True)
    assert len(ax.patches) == 2


def test_unsupported_edges_never_draw():
    flows = pd.DataFrame(
        {
            "source": ["A"],
            "destination": ["B"],
            "migration": [0.03],
            "support": [0.99],
            "status": ["unsupported"],
        }
    )
    ax = plot_directional_map(coordinates(), flows, include_ambiguous=True)
    assert len(ax.patches) == 0
