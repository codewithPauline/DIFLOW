"""Profile-likelihood diagnostics for asymmetric migration parameters."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2

from .candidate_models import asymmetric_migration_model
from .dadi_backend import _require_dadi
from .fit_models import MODEL_SPECS
from .multistart import fit_multistart
from .spectrum import prepare_observed_spectrum


@dataclass(frozen=True)
class ProfileLikelihoodResult:
    """One-dimensional profile likelihood for a migration parameter."""

    parameter: str
    table: pd.DataFrame
    mle: float
    maximum_log_likelihood: float
    confidence: float
    lower: float | None
    upper: float | None


def profile_grid(
    mle: float,
    *,
    lower: float,
    upper: float,
    points: int = 15,
    fold_range: float = 10.0,
) -> np.ndarray:
    """Generate a positive log-spaced profile grid centered on the MLE."""
    if not lower > 0 or not upper > lower:
        raise ValueError("profile bounds must satisfy 0 < lower < upper.")
    if points < 5:
        raise ValueError("points must be at least 5.")
    if mle <= 0:
        raise ValueError("mle must be positive.")
    if fold_range <= 1:
        raise ValueError("fold_range must be greater than 1.")

    lo = max(lower, mle / fold_range)
    hi = min(upper, mle * fold_range)
    values = np.geomspace(lo, hi, points)
    values = np.unique(np.concatenate([values, [mle]]))
    return np.sort(values)


def confidence_interval_from_profile(
    table: pd.DataFrame,
    *,
    confidence: float = 0.95,
) -> tuple[float | None, float | None]:
    """Approximate a profile-likelihood interval using a chi-square cutoff."""
    if not 0 < confidence < 1:
        raise ValueError("confidence must lie within (0, 1).")
    required = {"fixed_value", "log_likelihood"}
    if not required.issubset(table.columns):
        raise ValueError("profile table is missing required columns.")

    usable = table.replace([np.inf, -np.inf], np.nan).dropna(
        subset=["fixed_value", "log_likelihood"]
    )
    if usable.empty:
        return None, None

    best = float(usable["log_likelihood"].max())
    cutoff = 0.5 * float(chi2.ppf(confidence, df=1))
    inside = usable[usable["log_likelihood"] >= best - cutoff]
    if inside.empty:
        return None, None
    return (
        float(inside["fixed_value"].min()),
        float(inside["fixed_value"].max()),
    )


def profile_asymmetric_migration(
    observed_spectrum,
    *,
    parameter: str,
    values: np.ndarray | None = None,
    points: int = 15,
    starts: int = 10,
    maxiter: int = 100,
    confidence: float = 0.95,
    seed: int = 42,
    polarized: bool = False,
) -> ProfileLikelihoodResult:
    """Profile one scaled directional migration parameter.

    At each fixed migration value, all remaining asymmetric-model parameters are
    reoptimized. This is an identifiability diagnostic, not a replacement for
    bootstrap uncertainty.
    """
    if parameter not in {"m_a_to_b", "m_b_to_a"}:
        raise ValueError("parameter must be m_a_to_b or m_b_to_a.")

    data_array = np.asarray(observed_spectrum, dtype=float)
    if data_array.ndim != 2 or data_array.sum() <= 0:
        raise ValueError("observed_spectrum must be a non-empty 2D spectrum.")

    best = fit_multistart(
        data_array,
        "asymmetric_migration",
        starts=starts,
        maxiter=maxiter,
        seed=seed,
        polarized=polarized,
    )

    names = list(MODEL_SPECS["asymmetric_migration"]["names"])
    spec = MODEL_SPECS["asymmetric_migration"]
    target_index = names.index(parameter)
    mle = float(best.best_parameters[parameter])

    if values is None:
        values = profile_grid(
            mle,
            lower=float(spec["lower"][target_index]),
            upper=float(spec["upper"][target_index]),
            points=points,
        )
    else:
        values = np.asarray(values, dtype=float)
        if values.ndim != 1 or values.size < 2:
            raise ValueError("values must be a one-dimensional array with >=2 points.")
        if np.any(~np.isfinite(values)) or np.any(values <= 0):
            raise ValueError("profile values must be finite and positive.")

    dadi = _require_dadi()
    data = prepare_observed_spectrum(data_array, polarized=polarized)
    ns = (data_array.shape[0] - 1, data_array.shape[1] - 1)
    largest = max(ns)
    pts = [largest + 10, largest + 20, largest + 30]
    base_builder = asymmetric_migration_model(dadi)

    free_names = [name for name in names if name != parameter]
    free_indices = [names.index(name) for name in free_names]
    free_lower = [spec["lower"][i] for i in free_indices]
    free_upper = [spec["upper"][i] for i in free_indices]

    rows: list[dict] = []
    for grid_index, fixed in enumerate(values):
        fixed = float(fixed)
        if fixed < spec["lower"][target_index] or fixed > spec["upper"][target_index]:
            rows.append(
                {
                    "fixed_value": fixed,
                    "log_likelihood": np.nan,
                    "success": False,
                    "error": "fixed value outside model bounds",
                }
            )
            continue

        def reduced_model(free_params, ns_inner, pts_inner):
            full = []
            free_iter = iter(free_params)
            for name in names:
                if name == parameter:
                    full.append(fixed)
                else:
                    full.append(next(free_iter))
            return base_builder(full, ns_inner, pts_inner)

        extrapolated = dadi.Numerics.make_extrap_log_func(reduced_model)
        start = [best.best_parameters[name] for name in free_names]

        try:
            optimized = dadi.Inference.optimize_log_lbfgsb(
                start,
                data,
                extrapolated,
                pts,
                lower_bound=free_lower,
                upper_bound=free_upper,
                multinom=True,
                maxiter=maxiter,
                full_output=False,
            )
            model_sfs = extrapolated(optimized, ns, pts)
            ll = float(dadi.Inference.ll_multinom(model_sfs, data))
            row = {
                "fixed_value": fixed,
                "log_likelihood": ll,
                "success": True,
            }
            row.update(
                {
                    name: float(value)
                    for name, value in zip(free_names, optimized)
                }
            )
        except Exception as exc:
            row = {
                "fixed_value": fixed,
                "log_likelihood": np.nan,
                "success": False,
                "error": str(exc),
            }
        rows.append(row)

    table = pd.DataFrame(rows).sort_values("fixed_value").reset_index(drop=True)
    finite = table["log_likelihood"].replace([np.inf, -np.inf], np.nan).dropna()
    if finite.empty:
        raise RuntimeError("all profile-likelihood evaluations failed.")

    maximum = float(finite.max())
    table["delta_log_likelihood"] = maximum - table["log_likelihood"]
    lower, upper = confidence_interval_from_profile(
        table,
        confidence=confidence,
    )

    return ProfileLikelihoodResult(
        parameter=parameter,
        table=table,
        mle=mle,
        maximum_log_likelihood=maximum,
        confidence=confidence,
        lower=lower,
        upper=upper,
    )


def write_profile_likelihood(
    spectrum_path: str | Path,
    *,
    parameter: str,
    output_dir: str | Path,
    points: int = 15,
    starts: int = 10,
    maxiter: int = 100,
    confidence: float = 0.95,
    seed: int = 42,
    polarized: bool = False,
) -> ProfileLikelihoodResult:
    """Profile a saved DIFLOW .npy spectrum and write tables/figures."""
    spectrum = np.load(spectrum_path)
    result = profile_asymmetric_migration(
        spectrum,
        parameter=parameter,
        points=points,
        starts=starts,
        maxiter=maxiter,
        confidence=confidence,
        seed=seed,
        polarized=polarized,
    )

    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    result.table.to_csv(
        outdir / f"profile_{parameter}.csv",
        index=False,
    )

    successful = result.table[result.table["success"].astype(bool)]
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(
        successful["fixed_value"],
        successful["delta_log_likelihood"],
        marker="o",
    )
    cutoff = 0.5 * float(chi2.ppf(confidence, df=1))
    ax.axhline(cutoff, linestyle="--", linewidth=1)
    ax.axvline(result.mle, linestyle=":", linewidth=1)
    ax.set_xscale("log")
    ax.set_xlabel(f"Fixed scaled {parameter}")
    ax.set_ylabel("Delta log-likelihood")
    ax.set_title(f"Profile likelihood: {parameter}")
    fig.tight_layout()
    fig.savefig(
        outdir / f"profile_{parameter}.png",
        dpi=300,
        bbox_inches="tight",
    )
    fig.savefig(
        outdir / f"profile_{parameter}.pdf",
        bbox_inches="tight",
    )
    plt.close(fig)

    return result
