"""Directed migration network utilities."""

from __future__ import annotations

import networkx as nx
import pandas as pd


def build_flow_network(
    flows: pd.DataFrame,
    *,
    include_ambiguous: bool = False,
) -> nx.DiGraph:
    """Build a directed graph from a standardized DIFLOW flow table."""
    required = {"source", "destination", "migration"}
    if not required.issubset(flows.columns):
        raise ValueError("flows must contain source, destination, migration.")

    graph = nx.DiGraph()

    for row in flows.itertuples(index=False):
        status = getattr(row, "status", "supported")
        if status == "unsupported":
            continue
        if status == "ambiguous" and not include_ambiguous:
            continue

        attrs = {
            "migration": float(row.migration),
            "status": status,
        }
        for optional in ("support", "asymmetry_index", "model_weight"):
            if hasattr(row, optional):
                attrs[optional] = getattr(row, optional)

        graph.add_edge(str(row.source), str(row.destination), **attrs)

    return graph


def source_sink_summary(
    graph: nx.DiGraph,
) -> pd.DataFrame:
    """Summarize outgoing and incoming migration for each population.

    net_flow = outgoing - incoming

    Positive values indicate a population that is more source-like within the
    fitted network; negative values indicate a more sink-like population.

    This is a network summary, not a claim about ecological source-sink
    demography.
    """
    rows = []

    for node in sorted(graph.nodes):
        outgoing = sum(
            float(data.get("migration", 0.0))
            for _, _, data in graph.out_edges(node, data=True)
        )
        incoming = sum(
            float(data.get("migration", 0.0))
            for _, _, data in graph.in_edges(node, data=True)
        )
        net = outgoing - incoming

        if net > 0:
            role = "source-like"
        elif net < 0:
            role = "sink-like"
        else:
            role = "balanced"

        rows.append(
            {
                "population": node,
                "outgoing_migration": outgoing,
                "incoming_migration": incoming,
                "net_flow": net,
                "network_role": role,
                "out_degree": graph.out_degree(node),
                "in_degree": graph.in_degree(node),
            }
        )

    return pd.DataFrame(rows)
