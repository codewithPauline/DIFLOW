import pandas as pd
import pytest

from diflow.network import build_candidate_pairs, build_flow_network, source_sink_summary


def test_candidate_pairs_respect_distance():
    coordinates = pd.DataFrame(
        {
            "population": ["A", "B", "C"],
            "latitude": [39.0, 39.1, 41.0],
            "longitude": [-84.0, -84.1, -84.0],
        }
    )
    pairs = build_candidate_pairs(coordinates, max_distance_km=30)
    observed = {
        tuple(x)
        for x in pairs[["population_a", "population_b"]].itertuples(index=False, name=None)
    }
    assert observed == {("A", "B")}


def test_k_nearest_keeps_sparse_neighbors():
    coordinates = pd.DataFrame(
        {
            "population": ["A", "B", "C", "D"],
            "latitude": [39.0, 39.1, 39.2, 40.0],
            "longitude": [-84.0, -84.0, -84.0, -84.0],
        }
    )
    pairs = build_candidate_pairs(coordinates, k_nearest=1)
    assert len(pairs) >= 2
    assert len(pairs) < 6


def test_source_sink_summary():
    flows = pd.DataFrame(
        {
            "source": ["A", "A", "B"],
            "destination": ["B", "C", "C"],
            "migration": [0.03, 0.02, 0.01],
            "status": ["supported", "supported", "supported"],
        }
    )
    graph = build_flow_network(flows)
    summary = source_sink_summary(graph).set_index("population")

    assert summary.loc["A", "net_flow"] == pytest.approx(0.05)
    assert summary.loc["C", "net_flow"] == pytest.approx(-0.03)
    assert summary.loc["A", "network_role"] == "source-like"
    assert summary.loc["C", "network_role"] == "sink-like"


def test_ambiguous_edges_excluded_by_default():
    flows = pd.DataFrame(
        {
            "source": ["A", "B"],
            "destination": ["B", "C"],
            "migration": [0.03, 0.02],
            "status": ["supported", "ambiguous"],
        }
    )
    graph = build_flow_network(flows)
    assert graph.number_of_edges() == 1
