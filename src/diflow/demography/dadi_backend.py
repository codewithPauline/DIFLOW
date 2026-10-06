"""dadi numerical backend for asymmetric two-population jSFS inference.

DIFLOW owns the population labeling and migration-direction convention.
dadi supplies the diffusion-equation solver, spectrum extrapolation,
multinomial likelihood, and numerical optimizer.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .models import AsymmetricIMParams, to_dadi_split_asym_mig


@dataclass(frozen=True)
class DadiFitResult:
    """Result of fitting DIFLOW's asymmetric isolation-with-migration model."""

    params: AsymmetricIMParams
    log_likelihood: float
    theta: float
    model_spectrum: np.ndarray


def _require_dadi():
    try:
        import dadi
    except ImportError as exc:
        raise ImportError(
            "The dadi backend is optional. Install it with: "
            "python -m pip install -e '.[demography]'"
        ) from exc
    return dadi


def _diflow_model(dadi):
    """Return a dadi-compatible model function using DIFLOW parameter order."""

    def model(raw_params, ns, pts):
        nu_a, nu_b, split_time, m_a_to_b, m_b_to_a = raw_params
        dadi_params = (
            nu_a,
            nu_b,
            split_time,
            m_b_to_a,
            m_a_to_b,
        )
        return dadi.Demographics2D.split_asym_mig(dadi_params, ns, pts)

    return model


def expected_spectrum(
    params: AsymmetricIMParams,
    sample_sizes: tuple[int, int],
    *,
    grid_points: tuple[int, int, int] | None = None,
) -> np.ndarray:
    """Calculate the expected two-population jSFS with dadi.

    The returned spectrum is unscaled; the multinomial fitting step estimates
    the optimal theta separately.
    """
    dadi = _require_dadi()
    ns = tuple(int(x) for x in sample_sizes)
    if len(ns) != 2 or min(ns) < 1:
        raise ValueError("sample_sizes must contain two positive chromosome counts.")

    if grid_points is None:
        largest = max(ns)
        grid_points = (largest + 10, largest + 20, largest + 30)

    model = _diflow_model(dadi)
    extrapolated = dadi.Numerics.make_extrap_log_func(model)
    spectrum = extrapolated(
        [
            params.nu_a,
            params.nu_b,
            params.split_time,
            params.m_a_to_b,
            params.m_b_to_a,
        ],
        ns,
        list(grid_points),
    )
    return np.asarray(spectrum, dtype=float)


def fit_asymmetric_im(
    observed_spectrum,
    *,
    initial: AsymmetricIMParams | None = None,
    lower: AsymmetricIMParams | None = None,
    upper: AsymmetricIMParams | None = None,
    grid_points: tuple[int, int, int] | None = None,
    maxiter: int = 100,
) -> DadiFitResult:
    """Fit a two-population asymmetric isolation-with-migration model.

    This is the first single-time-point demographic inference backend in
    DIFLOW. It fits relative population sizes, split time, and the two
    directional scaled migration rates from an observed 2D jSFS.

    Multiple-start optimization and model comparison are intentionally handled
    in a separate layer so this function remains a transparent single fit.
    """
    dadi = _require_dadi()
    data_array = np.asarray(observed_spectrum, dtype=float)
    if data_array.ndim != 2:
        raise ValueError("observed_spectrum must be a two-dimensional array.")
    if np.any(~np.isfinite(data_array)) or np.any(data_array < 0):
        raise ValueError("observed_spectrum must contain finite non-negative values.")
    if data_array.sum() <= 0:
        raise ValueError("observed_spectrum must contain positive total mass.")

    ns = (data_array.shape[0] - 1, data_array.shape[1] - 1)
    data = dadi.Spectrum(data_array)

    initial = initial or AsymmetricIMParams(1.0, 1.0, 0.5, 0.5, 0.5)
    lower = lower or AsymmetricIMParams(1e-3, 1e-3, 1e-4, 1e-5, 1e-5)
    upper = upper or AsymmetricIMParams(100.0, 100.0, 20.0, 50.0, 50.0)

    p0 = [
        initial.nu_a,
        initial.nu_b,
        initial.split_time,
        initial.m_a_to_b,
        initial.m_b_to_a,
    ]
    lower_bound = [
        lower.nu_a,
        lower.nu_b,
        lower.split_time,
        lower.m_a_to_b,
        lower.m_b_to_a,
    ]
    upper_bound = [
        upper.nu_a,
        upper.nu_b,
        upper.split_time,
        upper.m_a_to_b,
        upper.m_b_to_a,
    ]

    if grid_points is None:
        largest = max(ns)
        grid_points = (largest + 10, largest + 20, largest + 30)

    model = _diflow_model(dadi)
    extrapolated = dadi.Numerics.make_extrap_log_func(model)

    optimized = dadi.Inference.optimize_log_lbfgsb(
        p0,
        data,
        extrapolated,
        list(grid_points),
        lower_bound=lower_bound,
        upper_bound=upper_bound,
        multinom=True,
        maxiter=maxiter,
        full_output=False,
    )

    fitted = AsymmetricIMParams(
        nu_a=float(optimized[0]),
        nu_b=float(optimized[1]),
        split_time=float(optimized[2]),
        m_a_to_b=float(optimized[3]),
        m_b_to_a=float(optimized[4]),
    )

    model_sfs = extrapolated(
        [
            fitted.nu_a,
            fitted.nu_b,
            fitted.split_time,
            fitted.m_a_to_b,
            fitted.m_b_to_a,
        ],
        ns,
        list(grid_points),
    )
    ll = float(dadi.Inference.ll_multinom(model_sfs, data))
    theta = float(dadi.Inference.optimal_sfs_scaling(model_sfs, data))

    return DadiFitResult(
        params=fitted,
        log_likelihood=ll,
        theta=theta,
        model_spectrum=np.asarray(model_sfs, dtype=float),
    )
