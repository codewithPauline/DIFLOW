"""Standardized comparator benchmarking for directional migration methods."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .metrics import (
    direction_accuracy,
    false_directional_positive_rate,
    parameter_bias,
    parameter_rmse,
)


_REQUIRED = {
    "scenario",
    "replicate",
    "truth_m_a_to_b",
    "truth_m_b_to_a",
    "estimated_m_a_to_b",
    "estimated_m_b_to_a",
}


def validate_comparator_table(
    frame: pd.DataFrame,
    *,
    method: str,
) -> pd.DataFrame:
    """Validate and label one method's standardized benchmark table."""
    missing = sorted(_REQUIRED - set(frame.columns))
    if missing:
        raise ValueError(
            f"{method} comparator table is missing columns: " + ", ".join(missing)
        )

    out = frame.copy()
    out["method"] = str(method)
    numeric = [
        "truth_m_a_to_b",
        "truth_m_b_to_a",
        "estimated_m_a_to_b",
        "estimated_m_b_to_a",
    ]
    for column in numeric:
        out[column] = pd.to_numeric(out[column], errors="coerce")

    out = out.dropna(subset=numeric)
    if out.empty:
        raise ValueError(f"{method} comparator table has no usable rows.")
    if (out[numeric] < 0).any().any():
        raise ValueError("migration truth and estimates must be non-negative.")
    return out


def summarize_method_comparison(
    results: pd.DataFrame,
    *,
    asymmetry_threshold: float = 0.25,
) -> pd.DataFrame:
    """Compute common recovery metrics by method and scenario."""
    required = _REQUIRED | {"method"}
    if not required.issubset(results.columns):
        raise ValueError("comparison results are missing required columns.")
    if not 0 <= asymmetry_threshold <= 1:
        raise ValueError("asymmetry_threshold must lie within [0, 1].")

    rows = []
    for (method, scenario), group in results.groupby(
        ["method", "scenario"],
        sort=False,
    ):
        true_ab = group["truth_m_a_to_b"].to_numpy(float)
        true_ba = group["truth_m_b_to_a"].to_numpy(float)
        est_ab = group["estimated_m_a_to_b"].to_numpy(float)
        est_ba = group["estimated_m_b_to_a"].to_numpy(float)
        symmetric = np.allclose(true_ab, true_ba)

        rows.append(
            {
                "method": method,
                "scenario": scenario,
                "replicates": len(group),
                "bias_m_a_to_b": parameter_bias(true_ab, est_ab),
                "rmse_m_a_to_b": parameter_rmse(true_ab, est_ab),
                "bias_m_b_to_a": parameter_bias(true_ba, est_ba),
                "rmse_m_b_to_a": parameter_rmse(true_ba, est_ba),
                "direction_accuracy": direction_accuracy(
                    true_ab,
                    true_ba,
                    est_ab,
                    est_ba,
                ),
                "false_directional_positive_rate": (
                    false_directional_positive_rate(
                        est_ab,
                        est_ba,
                        true_migration=float(true_ab[0]),
                        asymmetry_threshold=asymmetry_threshold,
                    )
                    if symmetric
                    else np.nan
                ),
            }
        )
    return pd.DataFrame(rows)


def compare_method_files(
    methods: dict[str, str | Path],
    *,
    output_dir: str | Path,
    asymmetry_threshold: float = 0.25,
) -> tuple[pd.DataFrame, pd.DataFrame, list[Path]]:
    """Read standardized method CSVs and write side-by-side benchmark outputs."""
    if len(methods) < 2:
        raise ValueError("at least two methods are required for comparison.")

    frames = []
    for method, path in methods.items():
        frame = pd.read_csv(path)
        frames.append(validate_comparator_table(frame, method=method))

    combined = pd.concat(frames, ignore_index=True)
    summary = summarize_method_comparison(
        combined,
        asymmetry_threshold=asymmetry_threshold,
    )

    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    combined.to_csv(outdir / "method_comparison_replicates.csv", index=False)
    summary.to_csv(outdir / "method_comparison_summary.csv", index=False)

    figures: list[Path] = []

    directional = summary[summary["false_directional_positive_rate"].isna()]
    if not directional.empty:
        pivot = directional.pivot_table(
            index="scenario",
            columns="method",
            values="direction_accuracy",
            aggfunc="mean",
        )
        ax = pivot.plot(kind="bar", figsize=(9, 5))
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("Direction accuracy")
        ax.set_title("Directional recovery by method")
        fig = ax.figure
        fig.tight_layout()
        for suffix in ("png", "pdf"):
            path = outdir / f"method_direction_accuracy.{suffix}"
            fig.savefig(path, dpi=300 if suffix == "png" else None, bbox_inches="tight")
            figures.append(path)
        plt.close(fig)

    symmetric = summary[summary["false_directional_positive_rate"].notna()]
    if not symmetric.empty:
        pivot = symmetric.pivot_table(
            index="scenario",
            columns="method",
            values="false_directional_positive_rate",
            aggfunc="mean",
        )
        ax = pivot.plot(kind="bar", figsize=(9, 5))
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("False directional-positive rate")
        ax.set_title("False directional calls by method")
        fig = ax.figure
        fig.tight_layout()
        for suffix in ("png", "pdf"):
            path = outdir / f"method_false_direction.{suffix}"
            fig.savefig(path, dpi=300 if suffix == "png" else None, bbox_inches="tight")
            figures.append(path)
        plt.close(fig)

    return combined, summary, figures
