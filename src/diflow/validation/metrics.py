"""Validation metrics for directional migration inference."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


@dataclass(frozen=True)
class ValidationMetrics:
    """Summary metrics for one migration parameter or validation scenario."""

    bias: float
    rmse: float
    interval_coverage: float | None = None
    direction_accuracy: float | None = None
    false_directional_positive_rate: float | None = None


def _arrays(truth, estimates) -> tuple[np.ndarray, np.ndarray]:
    truth_arr = np.asarray(truth, dtype=float)
    estimate_arr = np.asarray(estimates, dtype=float)
    if truth_arr.shape != estimate_arr.shape:
        raise ValueError("truth and estimates must have identical shapes.")
    if truth_arr.size == 0:
        raise ValueError("validation arrays cannot be empty.")
    if not np.all(np.isfinite(truth_arr)) or not np.all(np.isfinite(estimate_arr)):
        raise ValueError("validation arrays must contain only finite values.")
    return truth_arr, estimate_arr


def parameter_bias(truth, estimates) -> float:
    """Mean signed estimation error."""
    truth_arr, estimate_arr = _arrays(truth, estimates)
    return float(np.mean(estimate_arr - truth_arr))


def parameter_rmse(truth, estimates) -> float:
    """Root mean squared estimation error."""
    truth_arr, estimate_arr = _arrays(truth, estimates)
    return float(np.sqrt(np.mean((estimate_arr - truth_arr) ** 2)))


def interval_coverage(truth, lower, upper) -> float:
    """Fraction of true values contained in reported intervals."""
    truth_arr = np.asarray(truth, dtype=float)
    lower_arr = np.asarray(lower, dtype=float)
    upper_arr = np.asarray(upper, dtype=float)

    if not (truth_arr.shape == lower_arr.shape == upper_arr.shape):
        raise ValueError("truth, lower, and upper must have identical shapes.")
    if truth_arr.size == 0:
        raise ValueError("coverage arrays cannot be empty.")
    if np.any(lower_arr > upper_arr):
        raise ValueError("lower interval bounds cannot exceed upper bounds.")

    covered = (truth_arr >= lower_arr) & (truth_arr <= upper_arr)
    return float(np.mean(covered))


def _direction(m_ab: float, m_ba: float, tolerance: float = 0.0) -> int:
    difference = float(m_ab) - float(m_ba)
    if abs(difference) <= tolerance:
        return 0
    return 1 if difference > 0 else -1


def direction_accuracy(
    true_m_ab,
    true_m_ba,
    estimated_m_ab,
    estimated_m_ba,
    *,
    tolerance: float = 0.0,
) -> float:
    """Fraction of replicates with correctly recovered migration direction."""
    arrays = [
        np.asarray(x, dtype=float)
        for x in (true_m_ab, true_m_ba, estimated_m_ab, estimated_m_ba)
    ]
    if len({arr.shape for arr in arrays}) != 1:
        raise ValueError("all direction arrays must have identical shapes.")
    if arrays[0].size == 0:
        raise ValueError("direction arrays cannot be empty.")

    correct = []
    for t_ab, t_ba, e_ab, e_ba in zip(*arrays):
        correct.append(
            _direction(t_ab, t_ba, tolerance)
            == _direction(e_ab, e_ba, tolerance)
        )
    return float(np.mean(correct))


def false_directional_positive_rate(
    estimated_m_ab,
    estimated_m_ba,
    *,
    true_migration: float | None = None,
    asymmetry_threshold: float = 0.25,
) -> float:
    """Directional false-positive rate under symmetric migration truth.

    true_migration is retained for reporting/API clarity; symmetry is the key
    null condition. A directional positive occurs when the estimated absolute
    asymmetry index exceeds asymmetry_threshold.
    """
    ab = np.asarray(estimated_m_ab, dtype=float)
    ba = np.asarray(estimated_m_ba, dtype=float)
    if ab.shape != ba.shape:
        raise ValueError("estimated migration arrays must have identical shapes.")
    if ab.size == 0:
        raise ValueError("estimated migration arrays cannot be empty.")
    if not 0 <= asymmetry_threshold <= 1:
        raise ValueError("asymmetry_threshold must lie within [0, 1].")
    if true_migration is not None and (not math.isfinite(true_migration) or true_migration < 0):
        raise ValueError("true_migration must be finite and non-negative.")

    total = ab + ba
    asymmetry = np.zeros_like(total, dtype=float)
    nonzero = total > 0
    asymmetry[nonzero] = (ab[nonzero] - ba[nonzero]) / total[nonzero]
    return float(np.mean(np.abs(asymmetry) >= asymmetry_threshold))
