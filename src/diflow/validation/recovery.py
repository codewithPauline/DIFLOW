"""Known-truth recovery benchmarks for DIFLOW's jSFS inference layer."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from diflow.demography import AsymmetricIMParams, expected_spectrum, fit_multistart
from .metrics import (
    direction_accuracy,
    false_directional_positive_rate,
    parameter_bias,
    parameter_rmse,
)


@dataclass(frozen=True)
class RecoveryScenario:
    """Model-consistent truth used for parameter-recovery benchmarking."""

    name: str
    params: AsymmetricIMParams
    expected_direction: str


def default_recovery_scenarios() -> tuple[RecoveryScenario, ...]:
    """Core symmetric and asymmetric truths for the first benchmark suite."""
    base = dict(nu_a=1.0, nu_b=1.0, split_time=0.75)
    return (
        RecoveryScenario(
            "symmetric_migration",
            AsymmetricIMParams(**base, m_a_to_b=0.5, m_b_to_a=0.5),
            "symmetric",
        ),
        RecoveryScenario(
            "moderate_a_to_b",
            AsymmetricIMParams(**base, m_a_to_b=1.0, m_b_to_a=0.25),
            "A->B",
        ),
        RecoveryScenario(
            "moderate_b_to_a",
            AsymmetricIMParams(**base, m_a_to_b=0.25, m_b_to_a=1.0),
            "B->A",
        ),
        RecoveryScenario(
            "strong_a_to_b",
            AsymmetricIMParams(**base, m_a_to_b=2.0, m_b_to_a=0.10),
            "A->B",
        ),
        RecoveryScenario(
            "strong_b_to_a",
            AsymmetricIMParams(**base, m_a_to_b=0.10, m_b_to_a=2.0),
            "B->A",
        ),
    )


def simulate_model_consistent_spectrum(
    params: AsymmetricIMParams,
    *,
    sample_sizes: tuple[int, int] = (20, 20),
    segregating_sites: int = 5000,
    seed: int | None = None,
    spectrum_function: Callable | None = None,
) -> np.ndarray:
    """Draw a finite-SNP jSFS from a dadi demographic expectation.

    The absent and globally fixed corners are excluded before multinomial
    sampling because ordinary variant-only VCF data do not contain those
    categories. The downstream fitter masks the same corners.
    """
    if len(sample_sizes) != 2 or min(sample_sizes) < 2:
        raise ValueError("sample_sizes must contain two values >= 2.")
    if segregating_sites < 1:
        raise ValueError("segregating_sites must be positive.")

    generator = expected_spectrum if spectrum_function is None else spectrum_function
    expectation = np.asarray(generator(params, sample_sizes), dtype=float)
    expected_shape = (sample_sizes[0] + 1, sample_sizes[1] + 1)
    if expectation.shape != expected_shape:
        raise ValueError(
            f"expected spectrum shape {expected_shape}, got {expectation.shape}."
        )
    if np.any(~np.isfinite(expectation)) or np.any(expectation < 0):
        raise ValueError("expected spectrum must be finite and non-negative.")

    probabilities = expectation.copy()
    probabilities[0, 0] = 0.0
    probabilities[-1, -1] = 0.0
    total = float(probabilities.sum())
    if total <= 0:
        raise ValueError("expected spectrum has no segregating-site mass.")
    probabilities /= total

    rng = np.random.default_rng(seed)
    sampled = rng.multinomial(segregating_sites, probabilities.ravel())
    return sampled.reshape(expected_shape).astype(float)


def run_recovery_benchmark(
    *,
    scenarios: tuple[RecoveryScenario, ...] | None = None,
    replicates: int = 10,
    sample_sizes: tuple[int, int] = (20, 20),
    segregating_sites: int = 5000,
    starts: int = 10,
    maxiter: int = 100,
    seed: int = 42,
    spectrum_function: Callable | None = None,
    fit_function: Callable | None = None,
) -> pd.DataFrame:
    """Run model-consistent known-truth recovery replicates.

    This benchmark tests numerical/identifiability recovery when the generating
    model matches the fitted asymmetric continuous-migration model. It is a
    necessary benchmark, not a sufficient validation of real biological data.
    """
    if replicates < 1:
        raise ValueError("replicates must be at least 1.")
    if starts < 1:
        raise ValueError("starts must be at least 1.")

    scenarios = scenarios or default_recovery_scenarios()
    fitter = fit_multistart if fit_function is None else fit_function
    rows: list[dict] = []

    for scenario_index, scenario in enumerate(scenarios):
        for replicate in range(1, replicates + 1):
            run_seed = seed + scenario_index * 100000 + replicate
            spectrum = simulate_model_consistent_spectrum(
                scenario.params,
                sample_sizes=sample_sizes,
                segregating_sites=segregating_sites,
                seed=run_seed,
                spectrum_function=spectrum_function,
            )
            row = {
                "scenario": scenario.name,
                "replicate": replicate,
                "expected_direction": scenario.expected_direction,
                "true_m_a_to_b": scenario.params.m_a_to_b,
                "true_m_b_to_a": scenario.params.m_b_to_a,
                "sample_chromosomes_a": sample_sizes[0],
                "sample_chromosomes_b": sample_sizes[1],
                "segregating_sites": segregating_sites,
                "seed": run_seed,
            }
            try:
                fit = fitter(
                    spectrum,
                    "asymmetric_migration",
                    starts=starts,
                    maxiter=maxiter,
                    seed=run_seed,
                    polarized=False,
                )
                row.update(
                    {
                        "success": True,
                        "estimated_m_a_to_b": float(
                            fit.best_parameters["m_a_to_b"]
                        ),
                        "estimated_m_b_to_a": float(
                            fit.best_parameters["m_b_to_a"]
                        ),
                        "optimizer_stable": bool(fit.stable),
                        "optimizer_success_fraction": float(
                            fit.converged_fraction
                        ),
                        "log_likelihood": float(fit.best_log_likelihood),
                    }
                )
            except Exception as exc:
                row.update(
                    {
                        "success": False,
                        "estimated_m_a_to_b": np.nan,
                        "estimated_m_b_to_a": np.nan,
                        "optimizer_stable": False,
                        "optimizer_success_fraction": 0.0,
                        "log_likelihood": np.nan,
                        "error": str(exc),
                    }
                )
            rows.append(row)

    return pd.DataFrame(rows)


def summarize_recovery(
    results: pd.DataFrame,
    *,
    asymmetry_threshold: float = 0.25,
) -> pd.DataFrame:
    """Summarize recovery accuracy by known-truth scenario."""
    required = {
        "scenario",
        "success",
        "true_m_a_to_b",
        "true_m_b_to_a",
        "estimated_m_a_to_b",
        "estimated_m_b_to_a",
    }
    if not required.issubset(results.columns):
        raise ValueError("benchmark results are missing required columns.")

    summaries: list[dict] = []
    for scenario, group in results.groupby("scenario", sort=False):
        successful = group[group["success"].astype(bool)].copy()
        row = {
            "scenario": scenario,
            "attempted_replicates": len(group),
            "successful_replicates": len(successful),
            "success_rate": float(len(successful) / len(group)),
        }

        if successful.empty:
            row.update(
                {
                    "bias_m_a_to_b": np.nan,
                    "rmse_m_a_to_b": np.nan,
                    "bias_m_b_to_a": np.nan,
                    "rmse_m_b_to_a": np.nan,
                    "direction_accuracy": np.nan,
                    "false_directional_positive_rate": np.nan,
                    "optimizer_stable_rate": np.nan,
                }
            )
            summaries.append(row)
            continue

        true_ab = successful["true_m_a_to_b"].to_numpy(float)
        true_ba = successful["true_m_b_to_a"].to_numpy(float)
        est_ab = successful["estimated_m_a_to_b"].to_numpy(float)
        est_ba = successful["estimated_m_b_to_a"].to_numpy(float)

        row.update(
            {
                "bias_m_a_to_b": parameter_bias(true_ab, est_ab),
                "rmse_m_a_to_b": parameter_rmse(true_ab, est_ab),
                "bias_m_b_to_a": parameter_bias(true_ba, est_ba),
                "rmse_m_b_to_a": parameter_rmse(true_ba, est_ba),
                "direction_accuracy": direction_accuracy(
                    true_ab, true_ba, est_ab, est_ba
                ),
                "optimizer_stable_rate": float(
                    successful["optimizer_stable"].astype(bool).mean()
                ),
            }
        )

        symmetric_truth = np.allclose(true_ab, true_ba)
        row["false_directional_positive_rate"] = (
            false_directional_positive_rate(
                est_ab,
                est_ba,
                true_migration=float(true_ab[0]),
                asymmetry_threshold=asymmetry_threshold,
            )
            if symmetric_truth
            else np.nan
        )
        summaries.append(row)

    return pd.DataFrame(summaries)


def write_recovery_benchmark(
    *,
    output_dir: str | Path,
    **kwargs,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run the benchmark and save raw and summarized CSV outputs."""
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    raw = run_recovery_benchmark(**kwargs)
    summary = summarize_recovery(raw)
    raw.to_csv(outdir / "recovery_replicates.csv", index=False)
    summary.to_csv(outdir / "recovery_summary.csv", index=False)
    return raw, summary
