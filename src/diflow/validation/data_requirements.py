"""Derive empirical minimum-data recommendations from recovery-grid results."""

from __future__ import annotations

from dataclasses import dataclass, asdict
import json
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class DataRequirementTargets:
    """Release-policy targets for recommending tested data regimes."""

    min_direction_accuracy: float = 0.90
    max_false_direction_rate: float = 0.05
    min_success_rate: float = 0.95


def _validate_target(value: float, name: str) -> float:
    value = float(value)
    if not np.isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f"{name} must lie within [0, 1].")
    return value


def summarize_data_regimes(
    recovery_grid_summary: pd.DataFrame,
    *,
    targets: DataRequirementTargets | None = None,
) -> pd.DataFrame:
    """Aggregate recovery-grid performance by chromosomes and marker count.

    A regime passes only when all tested directional scenarios meet the
    direction-accuracy target, symmetric scenarios meet the false-direction
    target, and every cell meets the fitting-success target.
    """
    targets = targets or DataRequirementTargets()
    _validate_target(targets.min_direction_accuracy, "min_direction_accuracy")
    _validate_target(targets.max_false_direction_rate, "max_false_direction_rate")
    _validate_target(targets.min_success_rate, "min_success_rate")

    required = {
        "chromosomes",
        "segregating_sites",
        "scenario",
        "success_rate",
        "direction_accuracy",
        "false_directional_positive_rate",
    }
    if not required.issubset(recovery_grid_summary.columns):
        missing = sorted(required - set(recovery_grid_summary.columns))
        raise ValueError(
            "recovery grid summary is missing required columns: "
            + ", ".join(missing)
        )

    rows = []
    for (chromosomes, sites), group in recovery_grid_summary.groupby(
        ["chromosomes", "segregating_sites"],
        sort=True,
    ):
        success = pd.to_numeric(group["success_rate"], errors="coerce").dropna()
        min_success = float(success.min()) if not success.empty else np.nan

        directional = group[
            group["false_directional_positive_rate"].isna()
        ].copy()
        accuracy = pd.to_numeric(
            directional["direction_accuracy"], errors="coerce"
        ).dropna()
        min_accuracy = float(accuracy.min()) if not accuracy.empty else np.nan

        symmetric = group[
            group["false_directional_positive_rate"].notna()
        ].copy()
        fpr = pd.to_numeric(
            symmetric["false_directional_positive_rate"], errors="coerce"
        ).dropna()
        max_fpr = float(fpr.max()) if not fpr.empty else np.nan

        complete = (
            np.isfinite(min_success)
            and np.isfinite(min_accuracy)
            and np.isfinite(max_fpr)
        )
        passes = bool(
            complete
            and min_success >= targets.min_success_rate
            and min_accuracy >= targets.min_direction_accuracy
            and max_fpr <= targets.max_false_direction_rate
        )

        rows.append(
            {
                "chromosomes": int(chromosomes),
                "segregating_sites": int(sites),
                "minimum_success_rate": min_success,
                "minimum_direction_accuracy": min_accuracy,
                "maximum_false_direction_rate": max_fpr,
                "passes_targets": passes,
            }
        )

    return pd.DataFrame(rows).sort_values(
        ["segregating_sites", "chromosomes"],
        kind="stable",
    ).reset_index(drop=True)


def pareto_minimum_regimes(regimes: pd.DataFrame) -> pd.DataFrame:
    """Return passing regimes not dominated on both sample size and marker count."""
    required = {"chromosomes", "segregating_sites", "passes_targets"}
    if not required.issubset(regimes.columns):
        raise ValueError("regime table is missing required columns.")

    passing = regimes[regimes["passes_targets"].astype(bool)].copy()
    if passing.empty:
        return passing

    keep = []
    for index, row in passing.iterrows():
        dominated = (
            (passing["chromosomes"] <= row["chromosomes"])
            & (passing["segregating_sites"] <= row["segregating_sites"])
            & (
                (passing["chromosomes"] < row["chromosomes"])
                | (passing["segregating_sites"] < row["segregating_sites"])
            )
        ).any()
        if not dominated:
            keep.append(index)

    return passing.loc[keep].sort_values(
        ["segregating_sites", "chromosomes"],
        kind="stable",
    ).reset_index(drop=True)


def recommend_data_requirements(
    recovery_grid_summary: pd.DataFrame,
    *,
    targets: DataRequirementTargets | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Derive tested passing regimes and Pareto-minimum recommendations."""
    targets = targets or DataRequirementTargets()
    regimes = summarize_data_regimes(
        recovery_grid_summary,
        targets=targets,
    )
    minima = pareto_minimum_regimes(regimes)

    recommendation = {
        "targets": asdict(targets),
        "tested_regimes": int(len(regimes)),
        "passing_regimes": int(regimes["passes_targets"].astype(bool).sum()),
        "pareto_minimum_regimes": minima[
            ["chromosomes", "segregating_sites"]
        ].to_dict(orient="records"),
        "status": (
            "recommendations_available"
            if not minima.empty
            else "no_tested_regime_met_targets"
        ),
        "interpretation": (
            "Recommendations apply only to the simulated regimes and model "
            "conditions represented in the validation campaign."
        ),
    }
    return regimes, minima, recommendation


def write_data_requirements(
    *,
    recovery_grid_summary_csv: str | Path,
    output_dir: str | Path,
    targets: DataRequirementTargets | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Read recovery-grid summary and write empirical data guidance outputs."""
    summary = pd.read_csv(recovery_grid_summary_csv)
    regimes, minima, recommendation = recommend_data_requirements(
        summary,
        targets=targets,
    )

    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    regimes.to_csv(outdir / "data_regime_performance.csv", index=False)
    minima.to_csv(outdir / "minimum_passing_regimes.csv", index=False)
    (outdir / "data_requirements.json").write_text(
        json.dumps(recommendation, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    lines = [
        "# DIFLOW empirical data requirements",
        "",
        f"Status: **{recommendation['status']}**",
        "",
        "## Targets",
        "",
        f"- minimum direction accuracy: {targets.min_direction_accuracy if targets else 0.90}",
        f"- maximum false-direction rate: {targets.max_false_direction_rate if targets else 0.05}",
        f"- minimum fit success rate: {targets.min_success_rate if targets else 0.95}",
        "",
        "## Pareto-minimum tested regimes",
        "",
    ]
    if minima.empty:
        lines.append(
            "No tested sample-size / marker-count regime met all requested targets."
        )
    else:
        for row in minima.itertuples(index=False):
            lines.append(
                f"- {int(row.chromosomes)} chromosomes per population; "
                f"{int(row.segregating_sites)} segregating sites"
            )
    lines.extend(
        [
            "",
            "These are empirical recommendations within the tested simulation "
            "space, not universal population-genetic minimums.",
            "",
        ]
    )
    (outdir / "data_requirements.md").write_text(
        "\n".join(lines),
        encoding="utf-8",
    )
    return regimes, minima, recommendation
