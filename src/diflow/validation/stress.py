"""Demographic misspecification stress tests for DIFLOW."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from diflow.demography import ModelScore, fit_multistart, rank_models
from diflow.demography.candidate_models import (
    asymmetric_migration_model,
    no_migration_model,
    secondary_contact_asymmetric_model,
    secondary_contact_symmetric_model,
    symmetric_migration_model,
)
from diflow.demography.dadi_backend import _require_dadi


@dataclass(frozen=True)
class StressScenario:
    """Known generating history used to challenge the inference workflow."""

    name: str
    generating_model: str
    parameters: tuple[float, ...]
    expected_direction: str
    purpose: str


def default_stress_scenarios() -> tuple[StressScenario, ...]:
    """First executable demographic stress-test suite."""
    return (
        StressScenario(
            name="zero_migration",
            generating_model="isolation",
            parameters=(1.0, 1.0, 0.75),
            expected_direction="none",
            purpose="Measure false direction when true migration is absent.",
        ),
        StressScenario(
            name="unequal_ne_symmetric_migration",
            generating_model="symmetric_migration",
            parameters=(0.25, 2.0, 0.75, 0.5),
            expected_direction="symmetric",
            purpose="Test whether strong population-size asymmetry is mistaken for migration asymmetry.",
        ),
        StressScenario(
            name="secondary_contact_symmetric",
            generating_model="secondary_contact_symmetric",
            parameters=(1.0, 1.0, 0.75, 0.15, 0.5),
            expected_direction="symmetric",
            purpose="Test whether recent symmetric contact creates false directional evidence.",
        ),
        StressScenario(
            name="secondary_contact_a_to_b",
            generating_model="secondary_contact_asymmetric",
            parameters=(1.0, 1.0, 0.75, 0.15, 1.0, 0.15),
            expected_direction="A->B",
            purpose="Test direction recovery when gene flow occurs only after secondary contact.",
        ),
        StressScenario(
            name="secondary_contact_b_to_a",
            generating_model="secondary_contact_asymmetric",
            parameters=(1.0, 1.0, 0.75, 0.15, 0.15, 1.0),
            expected_direction="B->A",
            purpose="Reverse-direction secondary-contact recovery test.",
        ),
    )


_BUILDERS = {
    "isolation": no_migration_model,
    "symmetric_migration": symmetric_migration_model,
    "asymmetric_migration": asymmetric_migration_model,
    "secondary_contact_symmetric": secondary_contact_symmetric_model,
    "secondary_contact_asymmetric": secondary_contact_asymmetric_model,
}

_PARAMETER_COUNTS = {
    "isolation": 3,
    "symmetric_migration": 4,
    "asymmetric_migration": 5,
    "secondary_contact_symmetric": 5,
    "secondary_contact_asymmetric": 6,
}


def expected_candidate_spectrum(
    scenario: StressScenario,
    sample_sizes: tuple[int, int],
    *,
    grid_points: tuple[int, int, int] | None = None,
) -> np.ndarray:
    """Generate the expected jSFS under a named candidate demographic history."""
    if scenario.generating_model not in _BUILDERS:
        raise ValueError(f"unsupported generating model: {scenario.generating_model}")

    dadi = _require_dadi()
    ns = tuple(int(x) for x in sample_sizes)
    if len(ns) != 2 or min(ns) < 2:
        raise ValueError("sample_sizes must contain two values >= 2.")

    if grid_points is None:
        largest = max(ns)
        grid_points = (largest + 10, largest + 20, largest + 30)

    model = _BUILDERS[scenario.generating_model](dadi)
    extrapolated = dadi.Numerics.make_extrap_log_func(model)
    spectrum = extrapolated(list(scenario.parameters), ns, list(grid_points))
    return np.asarray(spectrum, dtype=float)


def sample_variant_spectrum(
    expectation: np.ndarray,
    *,
    segregating_sites: int,
    seed: int,
) -> np.ndarray:
    """Sample a finite variant-only jSFS from an expected spectrum."""
    if segregating_sites < 1:
        raise ValueError("segregating_sites must be positive.")
    probs = np.asarray(expectation, dtype=float).copy()
    if probs.ndim != 2 or np.any(~np.isfinite(probs)) or np.any(probs < 0):
        raise ValueError("expectation must be a finite non-negative 2D spectrum.")

    probs[0, 0] = 0.0
    probs[-1, -1] = 0.0
    total = float(probs.sum())
    if total <= 0:
        raise ValueError("expected spectrum has no segregating-site mass.")
    probs /= total

    rng = np.random.default_rng(seed)
    draws = rng.multinomial(segregating_sites, probs.ravel())
    return draws.reshape(probs.shape).astype(float)


def _preferred_direction(m_ab: float, m_ba: float) -> str:
    if m_ab > m_ba:
        return "A->B"
    if m_ba > m_ab:
        return "B->A"
    return "symmetric"


def run_stress_benchmark(
    *,
    scenarios: tuple[StressScenario, ...] | None = None,
    replicates: int = 10,
    sample_sizes: tuple[int, int] = (20, 20),
    segregating_sites: int = 5000,
    starts: int = 10,
    maxiter: int = 100,
    seed: int = 42,
    min_model_weight: float = 0.70,
    min_abs_asymmetry: float = 0.25,
) -> pd.DataFrame:
    """Fit the full candidate set to data from challenging known histories.

    Directional signal here is deliberately provisional: no bootstrap support is
    used. A spurious-direction flag requires asymmetric-model weight, optimizer
    stability, and asymmetry magnitude to pass the supplied thresholds.
    """
    if replicates < 1:
        raise ValueError("replicates must be at least 1.")
    if starts < 1:
        raise ValueError("starts must be at least 1.")
    if not 0 <= min_model_weight <= 1:
        raise ValueError("min_model_weight must lie in [0, 1].")
    if not 0 <= min_abs_asymmetry <= 1:
        raise ValueError("min_abs_asymmetry must lie in [0, 1].")

    scenarios = scenarios or default_stress_scenarios()
    model_names = tuple(_BUILDERS)
    rows: list[dict] = []

    for scenario_index, scenario in enumerate(scenarios):
        expectation = expected_candidate_spectrum(scenario, sample_sizes)

        for replicate in range(1, replicates + 1):
            run_seed = seed + scenario_index * 100000 + replicate
            observed = sample_variant_spectrum(
                expectation,
                segregating_sites=segregating_sites,
                seed=run_seed,
            )

            row = {
                "scenario": scenario.name,
                "replicate": replicate,
                "generating_model": scenario.generating_model,
                "expected_direction": scenario.expected_direction,
                "purpose": scenario.purpose,
                "seed": run_seed,
            }

            fits = {}
            scores = []
            observations = int(np.count_nonzero(observed))

            try:
                for model_offset, model_name in enumerate(model_names):
                    fit = fit_multistart(
                        observed,
                        model_name,
                        starts=starts,
                        maxiter=maxiter,
                        seed=run_seed + model_offset * 1000,
                        polarized=False,
                    )
                    fits[model_name] = fit
                    scores.append(
                        ModelScore(
                            name=model_name,
                            log_likelihood=fit.best_log_likelihood,
                            parameters=_PARAMETER_COUNTS[model_name],
                            observations=observations,
                        )
                    )

                ranking = rank_models(scores, use_aicc=False)
                best_model = str(ranking.iloc[0]["model"])
                generating_row = ranking[
                    ranking["model"] == scenario.generating_model
                ].iloc[0]
                directional_models = {
                    "asymmetric_migration",
                    "secondary_contact_asymmetric",
                }
                directional_model = (
                    best_model if best_model in directional_models else None
                )

                if directional_model is None:
                    m_ab = 0.0
                    m_ba = 0.0
                    asymmetry = 0.0
                    preferred = "none"
                    asymmetric_weight = 0.0
                    optimizer_stable = False
                    directional_signal = False
                else:
                    directional_row = ranking[
                        ranking["model"] == directional_model
                    ].iloc[0]
                    directional_fit = fits[directional_model]
                    m_ab = float(
                        directional_fit.best_parameters["m_a_to_b"]
                    )
                    m_ba = float(
                        directional_fit.best_parameters["m_b_to_a"]
                    )
                    total = m_ab + m_ba
                    asymmetry = (
                        0.0 if total == 0 else (m_ab - m_ba) / total
                    )
                    preferred = _preferred_direction(m_ab, m_ba)
                    asymmetric_weight = float(
                        directional_row["akaike_weight"]
                    )
                    optimizer_stable = bool(directional_fit.stable)
                    directional_signal = bool(
                        asymmetric_weight >= min_model_weight
                        and optimizer_stable
                        and abs(asymmetry) >= min_abs_asymmetry
                    )

                row.update(
                    {
                        "success": True,
                        "best_model": best_model,
                        "directional_model": directional_model,
                        "correct_model_selected": best_model
                        == scenario.generating_model,
                        "generating_model_weight": float(
                            generating_row["akaike_weight"]
                        ),
                        "asymmetric_model_weight": asymmetric_weight,
                        "estimated_m_a_to_b": m_ab,
                        "estimated_m_b_to_a": m_ba,
                        "estimated_asymmetry": asymmetry,
                        "preferred_direction": preferred,
                        "asymmetric_optimizer_stable": optimizer_stable,
                        "provisional_directional_signal": directional_signal,
                    }
                )
            except Exception as exc:
                row.update(
                    {
                        "success": False,
                        "best_model": None,
                        "correct_model_selected": False,
                        "provisional_directional_signal": False,
                        "error": str(exc),
                    }
                )

            rows.append(row)

    return pd.DataFrame(rows)


def summarize_stress(results: pd.DataFrame) -> pd.DataFrame:
    """Summarize model selection and false-direction behavior by scenario."""
    required = {
        "scenario",
        "success",
        "expected_direction",
        "correct_model_selected",
        "provisional_directional_signal",
    }
    if not required.issubset(results.columns):
        raise ValueError("stress benchmark results are missing required columns.")

    rows: list[dict] = []
    for scenario, group in results.groupby("scenario", sort=False):
        successful = group[group["success"].astype(bool)].copy()
        expected = str(group.iloc[0]["expected_direction"])
        summary = {
            "scenario": scenario,
            "attempted_replicates": len(group),
            "successful_replicates": len(successful),
            "success_rate": float(len(successful) / len(group)),
        }

        if successful.empty:
            summary.update(
                {
                    "correct_model_selection_rate": np.nan,
                    "provisional_directional_signal_rate": np.nan,
                    "direction_recovery_rate": np.nan,
                    "false_direction_signal_rate": np.nan,
                }
            )
        else:
            signal = successful["provisional_directional_signal"].astype(bool)
            summary["correct_model_selection_rate"] = float(
                successful["correct_model_selected"].astype(bool).mean()
            )
            summary["provisional_directional_signal_rate"] = float(signal.mean())

            if expected in {"A->B", "B->A"}:
                recovered = (
                    signal
                    & (successful["preferred_direction"].astype(str) == expected)
                )
                summary["direction_recovery_rate"] = float(recovered.mean())
                summary["false_direction_signal_rate"] = np.nan
            else:
                summary["direction_recovery_rate"] = np.nan
                summary["false_direction_signal_rate"] = float(signal.mean())

        rows.append(summary)

    return pd.DataFrame(rows)


def write_stress_benchmark(
    *,
    output_dir: str | Path,
    **kwargs,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    raw = run_stress_benchmark(**kwargs)
    summary = summarize_stress(raw)
    raw.to_csv(outdir / "stress_replicates.csv", index=False)
    summary.to_csv(outdir / "stress_summary.csv", index=False)
    return raw, summary
