"""Validation and normalization for directional migration map inputs."""

from __future__ import annotations

import numpy as np
import pandas as pd


COORD_COLUMNS = ("population", "latitude", "longitude")
FLOW_COLUMNS = ("source", "destination", "migration")


def validate_coordinates(coordinates: pd.DataFrame) -> pd.DataFrame:
    """Validate and return a normalized population-coordinate table."""
    missing = [c for c in COORD_COLUMNS if c not in coordinates.columns]
    if missing:
        raise ValueError(f"coordinates missing required columns: {missing}")

    frame = coordinates.loc[:, COORD_COLUMNS].copy()
    if frame["population"].duplicated().any():
        raise ValueError("each population must have exactly one coordinate.")

    for col in ("latitude", "longitude"):
        frame[col] = pd.to_numeric(frame[col], errors="raise")
        if not np.all(np.isfinite(frame[col])):
            raise ValueError(f"{col} must contain only finite values.")

    if ((frame["latitude"] < -90) | (frame["latitude"] > 90)).any():
        raise ValueError("latitude must lie within [-90, 90].")
    if ((frame["longitude"] < -180) | (frame["longitude"] > 180)).any():
        raise ValueError("longitude must lie within [-180, 180].")

    return frame


def validate_flows(flows: pd.DataFrame) -> pd.DataFrame:
    """Validate and normalize source-to-recipient migration edges.

    Required columns:
    source, destination, migration

    Optional columns:
    support: probability/confidence-like score in [0, 1]
    lower: lower uncertainty bound
    upper: upper uncertainty bound
    """
    missing = [c for c in FLOW_COLUMNS if c not in flows.columns]
    if missing:
        raise ValueError(f"flows missing required columns: {missing}")

    keep = list(FLOW_COLUMNS)
    for optional in ("support", "lower", "upper"):
        if optional in flows.columns:
            keep.append(optional)

    frame = flows.loc[:, keep].copy()
    if (frame["source"] == frame["destination"]).any():
        raise ValueError("self-migration edges are not supported in map input.")

    frame["migration"] = pd.to_numeric(frame["migration"], errors="raise")
    if (~np.isfinite(frame["migration"])).any() or (frame["migration"] < 0).any():
        raise ValueError("migration values must be finite and non-negative.")

    if "support" in frame.columns:
        frame["support"] = pd.to_numeric(frame["support"], errors="raise")
        if ((frame["support"] < 0) | (frame["support"] > 1)).any():
            raise ValueError("support must lie within [0, 1].")

    for col in ("lower", "upper"):
        if col in frame.columns:
            frame[col] = pd.to_numeric(frame[col], errors="raise")
            if (~np.isfinite(frame[col])).any() or (frame[col] < 0).any():
                raise ValueError(f"{col} must be finite and non-negative.")

    if {"lower", "upper"}.issubset(frame.columns):
        if (frame["lower"] > frame["upper"]).any():
            raise ValueError("lower uncertainty bound cannot exceed upper.")

    return frame
