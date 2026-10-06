"""Construct biologically plausible population-pair candidate sets."""

from __future__ import annotations

import math

import pandas as pd


def _haversine_km(lat1, lon1, lat2, lon2):
    radius = 6371.0088
    p1 = math.radians(lat1)
    p2 = math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    a = (
        math.sin(dphi / 2) ** 2
        + math.cos(p1) * math.cos(p2) * math.sin(dlambda / 2) ** 2
    )
    return 2 * radius * math.asin(math.sqrt(a))


def build_candidate_pairs(
    coordinates: pd.DataFrame,
    *,
    max_distance_km: float | None = None,
    k_nearest: int | None = None,
) -> pd.DataFrame:
    """Build a sparse undirected set of population pairs for inference.

    At least one sparsification rule must be supplied unless all pairwise
    comparisons are explicitly desired by passing both as None.

    If both max_distance_km and k_nearest are supplied, a pair is retained when
    it satisfies either criterion.
    """
    required = {"population", "latitude", "longitude"}
    if not required.issubset(coordinates.columns):
        raise ValueError("coordinates must contain population, latitude, longitude.")
    if coordinates["population"].duplicated().any():
        raise ValueError("population identifiers must be unique.")
    if max_distance_km is not None and max_distance_km <= 0:
        raise ValueError("max_distance_km must be positive.")
    if k_nearest is not None and k_nearest < 1:
        raise ValueError("k_nearest must be at least 1.")

    rows = []
    records = list(coordinates.itertuples(index=False))

    for i, a in enumerate(records):
        for j in range(i + 1, len(records)):
            b = records[j]
            distance = _haversine_km(
                float(a.latitude),
                float(a.longitude),
                float(b.latitude),
                float(b.longitude),
            )
            rows.append(
                {
                    "population_a": str(a.population),
                    "population_b": str(b.population),
                    "distance_km": distance,
                }
            )

    pairs = pd.DataFrame(rows)
    if pairs.empty:
        return pairs

    keep = pd.Series(False, index=pairs.index)

    if max_distance_km is not None:
        keep |= pairs["distance_km"] <= max_distance_km

    if k_nearest is not None:
        for population in coordinates["population"].astype(str):
            mask = (pairs["population_a"] == population) | (
                pairs["population_b"] == population
            )
            nearest = pairs.loc[mask].nsmallest(k_nearest, "distance_km").index
            keep.loc[nearest] = True

    if max_distance_km is None and k_nearest is None:
        keep[:] = True

    out = pairs.loc[keep].copy()
    return out.sort_values(
        ["population_a", "population_b"],
        kind="stable",
    ).reset_index(drop=True)
