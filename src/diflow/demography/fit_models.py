"""Fit and compare candidate two-population demographic models."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .candidate_models import (
    asymmetric_migration_model,
    no_migration_model,
    secondary_contact_asymmetric_model,
    symmetric_migration_model,
)
from .comparison import ModelScore, rank_models
from .dadi_backend import _require_dadi


@dataclass(frozen=True)
class CandidateFit:
    """One optimized candidate demographic model."""

    model: str
    parameters: dict[str, float]
    log_likelihood: float
    theta: float
    parameter_count: int


MODEL_SPECS = {
    "isolation": {
        "builder": no_migration_model,
        "names": ["nu_a", "nu_b", "split_time"],
        "initial": [1.0, 1.0, 0.5],
        "lower": [1e-3, 1e-3, 1e-4],
        "upper": [100.0, 100.0, 20.0],
    },
    "symmetric_migration": {
        "builder": symmetric_migration_model,
        "names": ["nu_a", "nu_b", "split_time", "migration"],
        "initial": [1.0, 1.0, 0.5, 0.5],
        "lower": [1e-3, 1e-3, 1e-4, 1e-5],
        "upper": [100.0, 100.0, 20.0, 50.0],
    },
    "asymmetric_migration": {
        "builder": asymmetric_migration_model,
        "names": ["nu_a", "nu_b", "split_time", "m_a_to_b", "m_b_to_a"],
        "initial": [1.0, 1.0, 0.5, 0.5, 0.5],
        "lower": [1e-3, 1e-3, 1e-4, 1e-5, 1e-5],
        "upper": [100.0, 100.0, 20.0, 50.0, 50.0],
    },
    "secondary_contact_asymmetric": {
        "builder": secondary_contact_asymmetric_model,
        "names": [
            "nu_a",
            "nu_b",
            "isolation_time",
            "contact_time",
            "m_a_to_b",
            "m_b_to_a",
        ],
        "initial": [1.0, 1.0, 0.5, 0.2, 0.5, 0.5],
        "lower": [1e-3, 1e-3, 1e-4, 1e-4, 1e-5, 1e-5],
        "upper": [100.0, 100.0, 20.0, 20.0, 50.0, 50.0],
    },
}


def _fit_candidate(
    observed_spectrum,
    model_name: str,
    *,
    grid_points: tuple[int, int, int] | None = None,
    maxiter: int = 100,
) -> CandidateFit:
    dadi = _require_dadi()
    data_array = np.asarray(observed_spectrum, dtype=float)
    if data_array.ndim != 2 or data_array.sum() <= 0:
        raise ValueError("observed_spectrum must be a non-empty 2D spectrum.")

    spec = MODEL_SPECS[model_name]
    ns = (data_array.shape[0] - 1, data_array.shape[1] - 1)
    data = dadi.Spectrum(data_array)

    if grid_points is None:
        largest = max(ns)
        grid_points = (largest + 10, largest + 20, largest + 30)

    model = spec["builder"](dadi)
    extrapolated = dadi.Numerics.make_extrap_log_func(model)

    optimized = dadi.Inference.optimize_log_lbfgsb(
        spec["initial"],
        data,
        extrapolated,
        list(grid_points),
        lower_bound=spec["lower"],
        upper_bound=spec["upper"],
        multinom=True,
        maxiter=maxiter,
        full_output=False,
    )

    model_sfs = extrapolated(
        optimized,
        ns,
        list(grid_points),
    )
    ll = float(dadi.Inference.ll_multinom(model_sfs, data))
    theta = float(dadi.Inference.optimal_sfs_scaling(model_sfs, data))

    params = {
        name: float(value)
        for name, value in zip(spec["names"], optimized)
    }

    return CandidateFit(
        model=model_name,
        parameters=params,
        log_likelihood=ll,
        theta=theta,
        parameter_count=len(spec["names"]),
    )


def compare_candidate_models(
    observed_spectrum,
    *,
    models: tuple[str, ...] = (
        "isolation",
        "symmetric_migration",
        "asymmetric_migration",
        "secondary_contact_asymmetric",
    ),
    grid_points: tuple[int, int, int] | None = None,
    maxiter: int = 100,
    use_aicc: bool = True,
):
    """Fit candidate models and rank them by AIC or AICc.

    The effective observation count is the number of non-zero jSFS cells.
    This is a pragmatic information-criterion sample-size definition for the
    current implementation and is reported transparently in the output.
    """
    unknown = [name for name in models if name not in MODEL_SPECS]
    if unknown:
        raise ValueError("unknown demographic models: " + ", ".join(unknown))

    data_array = np.asarray(observed_spectrum, dtype=float)
    observations = int(np.count_nonzero(data_array))
    if observations < 1:
        raise ValueError("observed_spectrum must contain at least one non-zero cell.")

    fits = [
        _fit_candidate(
            data_array,
            name,
            grid_points=grid_points,
            maxiter=maxiter,
        )
        for name in models
    ]

    ranking = rank_models(
        [
            ModelScore(
                name=fit.model,
                log_likelihood=fit.log_likelihood,
                parameters=fit.parameter_count,
                observations=observations,
            )
            for fit in fits
        ],
        use_aicc=use_aicc,
    )
    return fits, ranking
