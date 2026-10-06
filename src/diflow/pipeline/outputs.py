"""Convert pairwise inference results into network and map outputs."""

from __future__ import annotations

import pandas as pd

from diflow.mapping import plot_directional_map
from diflow.network import build_flow_network, source_sink_summary


def pairwise_to_flow_table(pairwise: pd.DataFrame) -> pd.DataFrame:
    rows = []
    columns = [
        "source",
        "destination",
        "migration",
        "support",
        "status",
        "asymmetry_index",
        "model_weight",
    ]

    if pairwise.empty:
        return pd.DataFrame(columns=columns)

    for row in pairwise.itertuples(index=False):
        status = getattr(row, "status", None)
        preferred = getattr(row, "preferred_direction", None)
        if status not in {"supported", "ambiguous"} or not preferred:
            continue

        pop_a = str(row.population_a)
        pop_b = str(row.population_b)

        if preferred == f"{pop_a}->{pop_b}":
            source, destination = pop_a, pop_b
            migration = float(row.m_a_to_b_scaled)
        elif preferred == f"{pop_b}->{pop_a}":
            source, destination = pop_b, pop_a
            migration = float(row.m_b_to_a_scaled)
        else:
            continue

        rows.append(
            {
                "source": source,
                "destination": destination,
                "migration": migration,
                "support": float(row.directional_support),
                "status": status,
                "asymmetry_index": float(row.asymmetry_index),
                "model_weight": float(row.asymmetric_model_weight),
            }
        )

    return pd.DataFrame(rows, columns=columns)


def write_network_outputs(
    *,
    pairwise: pd.DataFrame,
    coordinates: pd.DataFrame,
    output_dir,
    map_crs: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    output_dir.mkdir(parents=True, exist_ok=True)

    flows = pairwise_to_flow_table(pairwise)
    flows.to_csv(output_dir / "directional_flows.csv", index=False)

    graph = build_flow_network(flows)
    summary = source_sink_summary(graph)
    summary.to_csv(output_dir / "network_summary.csv", index=False)

    if not flows.empty:
        ax = plot_directional_map(
            coordinates,
            flows,
            include_ambiguous=True,
            target_crs=map_crs,
        )
        fig = ax.figure
        fig.tight_layout()
        fig.savefig(output_dir / "directional_map.png", dpi=300, bbox_inches="tight")
        fig.savefig(output_dir / "directional_map.pdf", bbox_inches="tight")

        import matplotlib.pyplot as plt
        plt.close(fig)

    return flows, summary
