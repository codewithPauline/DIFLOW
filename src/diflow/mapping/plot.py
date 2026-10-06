"""Publication-oriented directional migration map renderer."""

from __future__ import annotations

import math

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import pandas as pd
from pyproj import Transformer

from .flows import validate_coordinates, validate_flows


def _scaled(value: float, maximum: float, low: float, high: float) -> float:
    if maximum <= 0:
        return low
    fraction = math.sqrt(max(value, 0.0) / maximum)
    return low + (high - low) * fraction


def plot_directional_map(
    coordinates: pd.DataFrame,
    flows: pd.DataFrame,
    *,
    min_migration: float = 0.0,
    min_support: float | None = None,
    include_ambiguous: bool = False,
    label_populations: bool = True,
    node_size: float = 48.0,
    min_arrow_width: float = 0.7,
    max_arrow_width: float = 5.0,
    target_crs: str | None = None,
    ax=None,
):
    """Plot directed migration edges between geographic population coordinates.

    Supported edges are solid. Ambiguous edges are dashed when explicitly
    included. Unsupported edges are never drawn.
    """
    coords = validate_coordinates(coordinates)
    edge_table = validate_flows(flows)

    if min_migration < 0:
        raise ValueError("min_migration must be non-negative.")
    if min_support is not None and not 0 <= min_support <= 1:
        raise ValueError("min_support must lie within [0, 1].")

    edge_table = edge_table[edge_table["migration"] >= min_migration].copy()

    if "status" in edge_table.columns:
        edge_table = edge_table[edge_table["status"] != "unsupported"].copy()
        if not include_ambiguous:
            edge_table = edge_table[edge_table["status"] == "supported"].copy()

    if min_support is not None:
        if "support" not in edge_table.columns:
            raise ValueError("min_support requires a support column.")
        edge_table = edge_table[edge_table["support"] >= min_support].copy()

    plot_coords = coords.copy()
    x_column = "longitude"
    y_column = "latitude"
    x_label = "Longitude"
    y_label = "Latitude"

    if target_crs is not None:
        transformer = Transformer.from_crs(
            "EPSG:4326",
            target_crs,
            always_xy=True,
        )
        x, y = transformer.transform(
            plot_coords["longitude"].to_numpy(),
            plot_coords["latitude"].to_numpy(),
        )
        plot_coords["x"] = x
        plot_coords["y"] = y
        x_column = "x"
        y_column = "y"
        x_label = f"Easting ({target_crs})"
        y_label = f"Northing ({target_crs})"

    lookup = plot_coords.set_index("population")
    used = set(edge_table["source"]) | set(edge_table["destination"])
    missing = sorted(used - set(lookup.index))
    if missing:
        raise ValueError("flow populations missing coordinates: " + ", ".join(missing))

    if ax is None:
        _, ax = plt.subplots(figsize=(8, 6))

    ax.scatter(
        plot_coords[x_column],
        plot_coords[y_column],
        s=node_size,
        zorder=3,
    )

    if label_populations:
        for row in plot_coords.itertuples(index=False):
            ax.annotate(
                str(row.population),
                (getattr(row, x_column), getattr(row, y_column)),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=8,
                zorder=4,
            )

    max_migration = float(edge_table["migration"].max()) if not edge_table.empty else 0.0

    for row in edge_table.itertuples(index=False):
        source = lookup.loc[row.source]
        destination = lookup.loc[row.destination]
        width = _scaled(
            float(row.migration),
            max_migration,
            min_arrow_width,
            max_arrow_width,
        )
        support = float(getattr(row, "support", 1.0))
        alpha = 0.2 + 0.8 * support
        status = getattr(row, "status", "supported")
        linestyle = "--" if status == "ambiguous" else "-"

        arrow = FancyArrowPatch(
            (source[x_column], source[y_column]),
            (destination[x_column], destination[y_column]),
            arrowstyle="-|>",
            mutation_scale=10 + 2 * width,
            linewidth=width,
            linestyle=linestyle,
            alpha=alpha,
            shrinkA=7,
            shrinkB=7,
            connectionstyle="arc3,rad=0.08",
            zorder=2,
        )
        ax.add_patch(arrow)

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_aspect("equal", adjustable="datalim")
    return ax
