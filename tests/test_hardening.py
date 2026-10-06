import pandas as pd
import pytest

from diflow.cli import build_parser
from diflow.network import source_sink_summary
import networkx as nx


def test_benchmark_accepts_forward_suite():
    parser = build_parser()
    args = parser.parse_args(
        ["benchmark", "--output", "bench", "--suite", "forward"]
    )
    assert args.suite == "forward"


def test_empty_network_summary_has_stable_schema():
    summary = source_sink_summary(nx.DiGraph())
    assert list(summary.columns) == [
        "population",
        "outgoing_migration",
        "incoming_migration",
        "net_flow",
        "network_role",
        "out_degree",
        "in_degree",
    ]
    assert summary.empty
