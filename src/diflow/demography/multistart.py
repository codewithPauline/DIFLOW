"""Multi-start optimization and convergence diagnostics for demographic models."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
import pandas as pd

from .fit_models import MODEL_SPECS, _fit_candidate


@dataclass(frozen=True)
class MultiStartResult:
    """Collection of optimization replicates for one demographic model."""

    model: str
    runs: pd.DataFrame
    best_parameters: dict[str, float]
    best_log_likelihood: float
    converged_fraction: float
    stable: bool


def _random_log_uniform(
    rng: np.random.Generator,
    lower: float,
    upper: float,
) -> float:
    if lower <= 0 or upper <= 0 or lower >= upper:
        raise ValueError("log-uniform bounds must satisfy 0 < lower < upper.")
    return float(np.exp(rng.uniform(np.log(lower), np.log(upper))))


def generate_starting_points(
    model_name: str,
    *,
    starts: int,
    seed: int | None = None,
) -> list[list[float]]:
    """Generate reproducible log-uniform starting points within model bounds."""
    if model_name not in MODEL_SPECS:
        raise ValueError(f"unknown demographic model: {model_name}")
    if starts < 1:
        raise ValueError("starts must be at least 1.")

    spec = MODEL_SPECS[model_name]
    rng = np.random.default_rng(seed)
    points: list[list[float]] = []

    # Always include the canonical starting point first.
    points.append([float(x) for x in spec["initial"]])

    for _ in range(starts - 1):
        point = [
            _random_log_uniform(rng, lo, hi)
            for lo, hi in zip(spec["lower"], spec["upper"])
        ]
        points.append(point)
    return points


def _relative_spread(values: np.ndarray) -> float:
    finite = values[np.isfinite(values)]
    if finite.size < 2:
        return math.inf
    median = float(np.median(finite))
    scale = max(abs(median), 1e-12)
    return float((np.max(finite) - np.min(finite)) / scale)


def fit_multistart(
    observed_spectrum,
    model_name: str,
    *,
    starts: int = 20,
    seed: int | None = None,
    grid_points: tuple[int, int, int] | None = None,
    maxiter: int = 100,
    top_fraction: float = 0.25,
    ll_tolerance: float = 2.0,
    parameter_spread_tolerance: float = 0.5,
) -> MultiStartResult:
    """Fit one demographic model from multiple starting points.

    Stability is assessed among near-best runs:
    - keep the best top_fraction of successful runs,
    - require their log-likelihoods to lie within ll_tolerance of the best,
    - require each fitted parameter's relative spread to be below
      parameter_spread_tolerance.

    These thresholds are transparent heuristics for the development phase and
    will be calibrated through simulation.
    """
    if model_name not in MODEL_SPECS:
        raise ValueError(f"unknown demographic model: {model_name}")
    if not 0 < top_fraction <= 1:
        raise ValueError("top_fraction must lie in (0, 1].")
    if ll_tolerance < 0:
        raise ValueError("ll_tolerance must be non-negative.")
    if parameter_spread_tolerance < 0:
        raise ValueError("parameter_spread_tolerance must be non-negative.")

    spec = MODEL_SPECS[model_name]
    points = generate_starting_points(model_name, starts=starts, seed=seed)
    original_initial = list(spec["initial"])

    rows: list[dict] = []

    try:
        for run_id, point in enumerate(points, start=1):
            spec["initial"] = point
            try:
                fit = _fit_candidate(
                    observed_spectrum,
                    model_name,
                    grid_points=grid_points,
                    maxiter=maxiter,
                )
                row = {
                    "run": run_id,
                    "success": True,
                    "log_likelihood": fit.log_likelihood,
                    **fit.parameters,
                }
            except Exception:
                row = {
                    "run": run_id,
                    "success": False,
                    "log_likelihood": -math.inf,
                }
                for name in spec["names"]:
                    row[name] = np.nan
            rows.append(row)
    finally:
        spec["initial"] = original_initial

    frame = pd.DataFrame(rows)
    successful = frame[frame["success"]].copy()
    converged_fraction = float(len(successful) / len(frame))

    if successful.empty:
        raise RuntimeError(f"all optimization starts failed for model {model_name}.")

    successful = successful.sort_values("log_likelihood", ascending=False).reset_index(drop=True)
    best = successful.iloc[0]
    best_ll = float(best["log_likelihood"])

    near_best = successful[
        successful["log_likelihood"] >= best_ll - ll_tolerance
    ].copy()

    keep_n = max(2, int(math.ceil(len(successful) * top_fraction)))
    near_best = near_best.head(keep_n)

    spreads = {
        name: _relative_spread(near_best[name].to_numpy(dtype=float))
        for name in spec["names"]
    }
    stable = (
        len(near_best) >= 2
        and all(spread <= parameter_spread_tolerance for spread in spreads.values())
    )

    best_parameters = {
        name: float(best[name])
        for name in spec["names"]
    }

    return MultiStartResult(
        model=model_name,
        runs=frame,
        best_parameters=best_parameters,
        best_log_likelihood=best_ll,
        converged_fraction=converged_fraction,
        stable=bool(stable),
    )
