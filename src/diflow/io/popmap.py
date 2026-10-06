"""Population assignment file utilities."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def read_popmap(path: str | Path) -> pd.DataFrame:
    """Read a two-column sample-to-population mapping file.

    Accepted files are tab- or comma-separated and must contain columns named
    sample and population. Each sample may occur only once.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    frame = pd.read_csv(path, sep=None, engine="python", dtype=str)
    required = {"sample", "population"}
    if not required.issubset(frame.columns):
        raise ValueError("popmap must contain columns named 'sample' and 'population'.")

    frame = frame.loc[:, ["sample", "population"]].copy()
    if frame.isna().any().any() or (frame == "").any().any():
        raise ValueError("sample and population values must not be missing.")
    if frame["sample"].duplicated().any():
        dupes = frame.loc[frame["sample"].duplicated(), "sample"].tolist()
        raise ValueError(f"duplicate sample IDs in popmap: {dupes}")

    return frame
