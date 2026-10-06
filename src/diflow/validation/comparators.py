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






def standardize_diflow_benchmark(
    frame: pd.DataFrame,
    *,
    estimand: str = "dadi_scaled_migration",
) -> pd.DataFrame:
    """Convert DIFLOW recovery-style outputs to the comparator schema."""
    aliases = {
        "true_m_a_to_b": "truth_m_a_to_b",
        "true_m_b_to_a": "truth_m_b_to_a",
    }
    out = frame.rename(columns=aliases).copy()
    ordered = [
        "scenario",
        "replicate",
        "truth_m_a_to_b",
        "truth_m_b_to_a",
        "estimated_m_a_to_b",
        "estimated_m_b_to_a",
    ]
    required = set(ordered)
    missing = sorted(required - set(out.columns))
    if missing:
        raise ValueError(
            "DIFLOW benchmark table cannot be standardized; missing: "
            + ", ".join(missing)
        )
    out = out.loc[:, ordered].copy()
    out["estimand"] = str(estimand)
    return out


def write_standardized_diflow_benchmark(
    input_csv: str | Path,
    output_csv: str | Path,
    *,
    estimand: str = "dadi_scaled_migration",
) -> Path:
    """Write a DIFLOW benchmark table in the external-comparison schema."""
    frame = pd.read_csv(input_csv)
    standardized = standardize_diflow_benchmark(
        frame,
        estimand=estimand,
    )
    output = Path(output_csv)
    output.parent.mkdir(parents=True, exist_ok=True)
    standardized.to_csv(output, index=False)
    return output


def validate_matched_comparison(
    results: pd.DataFrame,
    *,
    require_complete_match: bool = True,
) -> pd.DataFrame:
    """Verify that methods are compared on identical known-truth replicates.

    For each (scenario, replicate) key, truth values must agree across methods.
    By default every method must contain the same set of keys.
    """
    required = _REQUIRED | {"method"}
    if not required.issubset(results.columns):
        raise ValueError("comparison results are missing required columns.")

    frame = results.copy()
    truth_counts = (
        frame.groupby(["scenario", "replicate"])[
            ["truth_m_a_to_b", "truth_m_b_to_a"]
        ]
        .nunique(dropna=False)
    )
    inconsistent = truth_counts[
        (truth_counts["truth_m_a_to_b"] > 1)
        | (truth_counts["truth_m_b_to_a"] > 1)
    ]
    if not inconsistent.empty:
        first = inconsistent.index[0]
        raise ValueError(
            "methods disagree on simulation truth for "
            f"scenario={first[0]!r}, replicate={first[1]!r}."
        )

    if require_complete_match:
        methods = sorted(set(frame["method"].astype(str)))
        expected = None
        expected_method = None
        for method in methods:
            subset = frame[frame["method"].astype(str) == method]
            keys = set(
                zip(
                    subset["scenario"].astype(str),
                    subset["replicate"].astype(str),
                )
            )
            if expected is None:
                expected = keys
                expected_method = method
                continue
            if keys != expected:
                missing = sorted(expected - keys)[:5]
                extra = sorted(keys - expected)[:5]
                raise ValueError(
                    "method benchmark keys are not matched: "
                    f"{method} differs from {expected_method}; "
                    f"missing examples={missing}, extra examples={extra}."
                )

    if "estimand" in frame.columns:
        estimands = sorted(
            {
                str(value).strip()
                for value in frame["estimand"].dropna()
                if str(value).strip()
            }
        )
        if len(estimands) > 1:
            raise ValueError(
                "comparison contains multiple estimands: "
                + ", ".join(estimands)
                + ". Convert to a compatible quantity or compare separately."
            )

    return frame


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




def summarize_direction_only_comparison(
    results: pd.DataFrame,
    *,
    asymmetry_threshold: float = 0.25,
) -> pd.DataFrame:
    """Compare directional recovery without assuming magnitude-scale equivalence.

    This mode is appropriate when methods estimate different migration
    quantities but still provide two directional scores/rates whose ordering can
    be compared against known truth.
    """
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
    require_complete_match: bool = True,
    direction_only: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, list[Path]]:
    """Read standardized method CSVs and write side-by-side benchmark outputs."""
    if len(methods) < 2:
        raise ValueError("at least two methods are required for comparison.")

    frames = []
    for method, path in methods.items():
        frame = pd.read_csv(path)
        frames.append(validate_comparator_table(frame, method=method))

    combined = pd.concat(frames, ignore_index=True)
    if direction_only:
        # Truth and replicate keys must still match, but estimands may differ
        # because only directional ordering is compared.
        estimand = combined.pop("estimand") if "estimand" in combined.columns else None
        combined = validate_matched_comparison(
            combined,
            require_complete_match=require_complete_match,
        )
        if estimand is not None:
            combined["estimand"] = estimand.to_numpy()
        summary = summarize_direction_only_comparison(
            combined,
            asymmetry_threshold=asymmetry_threshold,
        )
    else:
        combined = validate_matched_comparison(
            combined,
            require_complete_match=require_complete_match,
        )
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
