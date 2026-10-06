"""Simulation-driven calibration of DIFLOW directional decision thresholds."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class CalibratedThresholds:
    """Selected directional evidence thresholds."""

    min_model_weight: float
    min_directional_support: float
    min_abs_asymmetry: float
    false_directional_positive_rate: float
    directional_sensitivity: float
    direction_accuracy_when_called: float
    called_fraction: float
    evaluated_rows: int


def _validate_evidence_table(frame: pd.DataFrame) -> pd.DataFrame:
    required = {
        "truth_direction",
        "preferred_direction",
        "asymmetric_model_weight",
        "directional_support",
        "asymmetry_index",
        "optimizer_stable",
    }
    if not required.issubset(frame.columns):
        missing = sorted(required - set(frame.columns))
        raise ValueError(
            "calibration table is missing required columns: " + ", ".join(missing)
        )

    clean = frame.copy()
    for column in (
        "asymmetric_model_weight",
        "directional_support",
        "asymmetry_index",
    ):
        clean[column] = pd.to_numeric(clean[column], errors="coerce")
    clean = clean.dropna(
        subset=[
            "truth_direction",
            "asymmetric_model_weight",
            "directional_support",
            "asymmetry_index",
        ]
    )
    clean["optimizer_stable"] = clean["optimizer_stable"].astype(bool)

    probability_columns = ("asymmetric_model_weight", "directional_support")
    for column in probability_columns:
        if ((clean[column] < 0) | (clean[column] > 1)).any():
            raise ValueError(f"{column} must lie within [0, 1].")
    if (clean["asymmetry_index"].abs() > 1).any():
        raise ValueError("asymmetry_index must lie within [-1, 1].")

    valid_truth = {"symmetric", "none", "A->B", "B->A"}
    unknown = sorted(set(clean["truth_direction"]) - valid_truth)
    if unknown:
        raise ValueError(
            "truth_direction contains unsupported labels: " + ", ".join(unknown)
        )
    if clean.empty:
        raise ValueError("calibration table contains no usable rows.")
    return clean


def evaluate_thresholds(
    evidence: pd.DataFrame,
    *,
    min_model_weight: float,
    min_directional_support: float,
    min_abs_asymmetry: float,
) -> dict:
    """Evaluate one threshold combination against known simulation truth."""
    frame = _validate_evidence_table(evidence)

    called = (
        frame["optimizer_stable"]
        & (frame["asymmetric_model_weight"] >= min_model_weight)
        & (frame["directional_support"] >= min_directional_support)
        & (frame["asymmetry_index"].abs() >= min_abs_asymmetry)
        & frame["preferred_direction"].isin(["A->B", "B->A"])
    )

    null = frame["truth_direction"].isin(["symmetric", "none"])
    directional = frame["truth_direction"].isin(["A->B", "B->A"])

    false_positive_rate = (
        float(called[null].mean()) if int(null.sum()) > 0 else np.nan
    )
    correct_direction = called & (
        frame["preferred_direction"] == frame["truth_direction"]
    )
    sensitivity = (
        float(correct_direction[directional].mean())
        if int(directional.sum()) > 0
        else np.nan
    )

    called_directional = called & directional
    accuracy_when_called = (
        float(
            (
                frame.loc[called_directional, "preferred_direction"]
                == frame.loc[called_directional, "truth_direction"]
            ).mean()
        )
        if int(called_directional.sum()) > 0
        else np.nan
    )

    return {
        "min_model_weight": float(min_model_weight),
        "min_directional_support": float(min_directional_support),
        "min_abs_asymmetry": float(min_abs_asymmetry),
        "false_directional_positive_rate": false_positive_rate,
        "directional_sensitivity": sensitivity,
        "direction_accuracy_when_called": accuracy_when_called,
        "called_fraction": float(called.mean()),
        "evaluated_rows": int(len(frame)),
        "null_rows": int(null.sum()),
        "directional_rows": int(directional.sum()),
    }


def scan_thresholds(
    evidence: pd.DataFrame,
    *,
    model_weights: tuple[float, ...] = (0.50, 0.60, 0.70, 0.80, 0.90),
    directional_supports: tuple[float, ...] = (0.90, 0.95, 0.975, 0.99),
    asymmetries: tuple[float, ...] = (0.10, 0.20, 0.25, 0.30, 0.40),
) -> pd.DataFrame:
    """Evaluate a grid of candidate directional decision thresholds."""
    rows = [
        evaluate_thresholds(
            evidence,
            min_model_weight=mw,
            min_directional_support=ds,
            min_abs_asymmetry=aa,
        )
        for mw, ds, aa in product(
            model_weights,
            directional_supports,
            asymmetries,
        )
    ]
    return pd.DataFrame(rows)


def select_thresholds(
    scan: pd.DataFrame,
    *,
    max_false_directional_positive_rate: float = 0.05,
) -> CalibratedThresholds:
    """Select the most sensitive threshold set meeting the false-positive target.

    Ties are broken toward higher direction accuracy when called, then toward a
    more conservative support/asymmetry/model-weight rule.
    """
    if not 0 <= max_false_directional_positive_rate <= 1:
        raise ValueError(
            "max_false_directional_positive_rate must lie within [0, 1]."
        )

    required = {
        "min_model_weight",
        "min_directional_support",
        "min_abs_asymmetry",
        "false_directional_positive_rate",
        "directional_sensitivity",
        "direction_accuracy_when_called",
        "called_fraction",
        "evaluated_rows",
    }
    if not required.issubset(scan.columns):
        raise ValueError("threshold scan is missing required columns.")

    eligible = scan[
        scan["false_directional_positive_rate"]
        <= max_false_directional_positive_rate
    ].copy()
    eligible = eligible.dropna(
        subset=["false_directional_positive_rate", "directional_sensitivity"]
    )

    if eligible.empty:
        raise RuntimeError(
            "no tested threshold combination met the requested false-direction "
            "target; expand the threshold grid or relax the target."
        )

    eligible["_accuracy"] = eligible["direction_accuracy_when_called"].fillna(-1)
    eligible = eligible.sort_values(
        [
            "directional_sensitivity",
            "_accuracy",
            "min_directional_support",
            "min_abs_asymmetry",
            "min_model_weight",
        ],
        ascending=[False, False, False, False, False],
        kind="stable",
    )
    best = eligible.iloc[0]

    return CalibratedThresholds(
        min_model_weight=float(best["min_model_weight"]),
        min_directional_support=float(best["min_directional_support"]),
        min_abs_asymmetry=float(best["min_abs_asymmetry"]),
        false_directional_positive_rate=float(
            best["false_directional_positive_rate"]
        ),
        directional_sensitivity=float(best["directional_sensitivity"]),
        direction_accuracy_when_called=float(
            best["direction_accuracy_when_called"]
        ),
        called_fraction=float(best["called_fraction"]),
        evaluated_rows=int(best["evaluated_rows"]),
    )


def calibrate_thresholds(
    evidence: pd.DataFrame,
    *,
    max_false_directional_positive_rate: float = 0.05,
    model_weights: tuple[float, ...] = (0.50, 0.60, 0.70, 0.80, 0.90),
    directional_supports: tuple[float, ...] = (0.90, 0.95, 0.975, 0.99),
    asymmetries: tuple[float, ...] = (0.10, 0.20, 0.25, 0.30, 0.40),
) -> tuple[CalibratedThresholds, pd.DataFrame]:
    """Scan and select directional thresholds from known-truth evidence."""
    scan = scan_thresholds(
        evidence,
        model_weights=model_weights,
        directional_supports=directional_supports,
        asymmetries=asymmetries,
    )
    selected = select_thresholds(
        scan,
        max_false_directional_positive_rate=max_false_directional_positive_rate,
    )
    return selected, scan


def write_threshold_calibration(
    *,
    evidence_csv: str | Path,
    output_dir: str | Path,
    max_false_directional_positive_rate: float = 0.05,
) -> tuple[CalibratedThresholds, pd.DataFrame]:
    """Read known-truth evidence, calibrate thresholds, and write outputs."""
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    evidence = pd.read_csv(evidence_csv)
    selected, scan = calibrate_thresholds(
        evidence,
        max_false_directional_positive_rate=max_false_directional_positive_rate,
    )
    scan.to_csv(outdir / "threshold_scan.csv", index=False)
    pd.DataFrame([selected.__dict__]).to_csv(
        outdir / "selected_thresholds.csv",
        index=False,
    )
    return selected, scan
