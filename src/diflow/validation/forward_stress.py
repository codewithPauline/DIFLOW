"""Out-of-model forward-time stress simulations for DIFLOW.

These simulations are deliberately independent of dadi's diffusion models.
They create sampled two-population jSFS data from explicit Wright-Fisher
allele-frequency histories, then challenge DIFLOW's demographic model set.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from diflow.demography import ModelScore, fit_multistart, rank_models


@dataclass(frozen=True)
class ForwardStressScenario:
    """A forward-time history that may violate fitted demographic assumptions."""

    name: str
    kind: str
    expected_direction: str
    purpose: str


def default_forward_stress_scenarios() -> tuple[ForwardStressScenario, ...]:
    return (
        ForwardStressScenario(
            "serial_founder_range_expansion",
            "range_expansion",
            "none",
            "Test whether founder-driven allele-frequency shifts are mistaken for ongoing directional migration.",
        ),
        ForwardStressScenario(
            "ghost_introgression_into_b",
            "ghost_introgression",
            "none",
            "Test whether migration from an unsampled third population creates a false A/B directional edge.",
        ),
        ForwardStressScenario(
            "bottleneck_in_b",
            "bottleneck",
            "none",
            "Test whether a recent population-size collapse creates false directional migration.",
        ),
        ForwardStressScenario(
            "uneven_sampling_no_migration",
            "uneven_sampling",
            "none",
            "Test whether unequal chromosome sampling alone creates a directional signal.",
        ),
    )


def _drift(
    rng: np.random.Generator,
    frequencies: np.ndarray,
    ne: int,
) -> np.ndarray:
    if ne < 2:
        raise ValueError("effective population size must be at least 2.")
    return rng.binomial(2 * ne, np.clip(frequencies, 0, 1)) / (2 * ne)


def _sample_jsfs(
    rng: np.random.Generator,
    p_a: np.ndarray,
    p_b: np.ndarray,
    sample_sizes: tuple[int, int],
) -> np.ndarray:
    n_a, n_b = map(int, sample_sizes)
    if min(n_a, n_b) < 2:
        raise ValueError("sample chromosome counts must be at least 2.")

    x_a = rng.binomial(n_a, np.clip(p_a, 0, 1))
    x_b = rng.binomial(n_b, np.clip(p_b, 0, 1))

    spectrum = np.zeros((n_a + 1, n_b + 1), dtype=float)
    for a, b in zip(x_a, x_b):
        # Mimic a variant-only input: loci fixed reference or fixed alternate
        # across both sampled populations are absent.
        if (a == 0 and b == 0) or (a == n_a and b == n_b):
            continue
        spectrum[int(a), int(b)] += 1.0
    return spectrum


def simulate_forward_stress_spectrum(
    scenario: ForwardStressScenario,
    *,
    loci: int = 5000,
    sample_sizes: tuple[int, int] = (20, 20),
    seed: int = 42,
) -> np.ndarray:
    """Generate a sampled jSFS from an explicit forward-time stress history."""
    if loci < 100:
        raise ValueError("loci must be at least 100 for stress simulation.")

    rng = np.random.default_rng(seed)
    ancestral = rng.beta(0.7, 0.7, size=loci)
    p_a = ancestral.copy()
    p_b = ancestral.copy()

    if scenario.kind == "range_expansion":
        # Source population A persists, while B is founded by a small sample
        # from A and then expands without ongoing migration.
        for _ in range(20):
            p_a = _drift(rng, p_a, 2000)

        founder_gene_copies = 20
        p_b = rng.binomial(
            founder_gene_copies,
            np.clip(p_a, 0, 1),
        ) / founder_gene_copies

        for _ in range(20):
            p_a = _drift(rng, p_a, 2000)
            p_b = _drift(rng, p_b, 2000)

    elif scenario.kind == "ghost_introgression":
        p_g = ancestral.copy()
        for _ in range(35):
            p_a = _drift(rng, p_a, 1500)
            p_b = _drift(rng, p_b, 1500)
            p_g = _drift(rng, p_g, 800)

        # Only G -> B occurs. A and B have no direct migration.
        for _ in range(10):
            p_g = _drift(rng, p_g, 800)
            p_a = _drift(rng, p_a, 1500)
            migrated_b = 0.85 * p_b + 0.15 * p_g
            p_b = _drift(rng, migrated_b, 1500)

    elif scenario.kind == "bottleneck":
        for _ in range(25):
            p_a = _drift(rng, p_a, 2000)
            p_b = _drift(rng, p_b, 2000)
        for _ in range(8):
            p_a = _drift(rng, p_a, 2000)
            p_b = _drift(rng, p_b, 60)
        for _ in range(12):
            p_a = _drift(rng, p_a, 2000)
            p_b = _drift(rng, p_b, 2000)

    elif scenario.kind == "uneven_sampling":
        for _ in range(35):
            p_a = _drift(rng, p_a, 1500)
            p_b = _drift(rng, p_b, 1500)

    else:
        raise ValueError(f"unknown forward stress kind: {scenario.kind}")

    return _sample_jsfs(rng, p_a, p_b, sample_sizes)


def _fit_full_candidate_set(
    observed: np.ndarray,
    *,
    starts: int,
    maxiter: int,
    seed: int,
) -> dict:
    model_names = (
        "isolation",
        "symmetric_migration",
        "asymmetric_migration",
        "secondary_contact_symmetric",
        "secondary_contact_asymmetric",
    )
    fits = {}
    scores = []
    observations = int(np.count_nonzero(observed))
    if observations < 1:
        raise ValueError("simulated spectrum contains no variable cells.")

    for offset, model_name in enumerate(model_names):
        fit = fit_multistart(
            observed,
            model_name,
            starts=starts,
            maxiter=maxiter,
            seed=seed + offset * 1000,
            polarized=False,
        )
        fits[model_name] = fit
        scores.append(
            ModelScore(
                name=model_name,
                log_likelihood=fit.best_log_likelihood,
                parameters=len(fit.best_parameters),
                observations=observations,
            )
        )

    ranking = rank_models(scores, use_aicc=False)
    asym_row = ranking[ranking["model"] == "asymmetric_migration"].iloc[0]
    asym_fit = fits["asymmetric_migration"]
    m_ab = float(asym_fit.best_parameters["m_a_to_b"])
    m_ba = float(asym_fit.best_parameters["m_b_to_a"])
    total = m_ab + m_ba
    asymmetry = 0.0 if total == 0 else (m_ab - m_ba) / total

    return {
        "best_model": str(ranking.iloc[0]["model"]),
        "asymmetric_model_weight": float(asym_row["akaike_weight"]),
        "estimated_m_a_to_b": m_ab,
        "estimated_m_b_to_a": m_ba,
        "estimated_asymmetry": asymmetry,
        "preferred_direction": (
            "A->B" if m_ab > m_ba else "B->A" if m_ba > m_ab else "symmetric"
        ),
        "asymmetric_optimizer_stable": bool(asym_fit.stable),
        "optimizer_success_fraction": float(asym_fit.converged_fraction),
    }


def run_forward_stress_benchmark(
    *,
    scenarios: tuple[ForwardStressScenario, ...] | None = None,
    replicates: int = 10,
    loci: int = 5000,
    sample_sizes: tuple[int, int] = (20, 20),
    starts: int = 10,
    maxiter: int = 100,
    seed: int = 42,
    min_model_weight: float = 0.70,
    min_abs_asymmetry: float = 0.25,
) -> pd.DataFrame:
    """Challenge DIFLOW using histories generated outside the dadi model family."""
    if replicates < 1:
        raise ValueError("replicates must be at least 1.")
    scenarios = scenarios or default_forward_stress_scenarios()
    rows: list[dict] = []

    for scenario_index, scenario in enumerate(scenarios):
        for replicate in range(1, replicates + 1):
            run_seed = seed + scenario_index * 100000 + replicate
            current_sample_sizes = sample_sizes
            if scenario.kind == "uneven_sampling":
                current_sample_sizes = (
                    max(4, sample_sizes[0] // 2),
                    sample_sizes[1],
                )

            row = {
                "scenario": scenario.name,
                "kind": scenario.kind,
                "replicate": replicate,
                "expected_direction": scenario.expected_direction,
                "purpose": scenario.purpose,
                "sample_chromosomes_a": current_sample_sizes[0],
                "sample_chromosomes_b": current_sample_sizes[1],
                "seed": run_seed,
            }

            try:
                observed = simulate_forward_stress_spectrum(
                    scenario,
                    loci=loci,
                    sample_sizes=current_sample_sizes,
                    seed=run_seed,
                )
                fit = _fit_full_candidate_set(
                    observed,
                    starts=starts,
                    maxiter=maxiter,
                    seed=run_seed,
                )
                directional_signal = bool(
                    fit["asymmetric_model_weight"] >= min_model_weight
                    and fit["asymmetric_optimizer_stable"]
                    and abs(fit["estimated_asymmetry"]) >= min_abs_asymmetry
                )
                row.update(
                    {
                        "success": True,
                        **fit,
                        "provisional_directional_signal": directional_signal,
                    }
                )
            except Exception as exc:
                row.update(
                    {
                        "success": False,
                        "provisional_directional_signal": False,
                        "error": str(exc),
                    }
                )
            rows.append(row)

    return pd.DataFrame(rows)


def summarize_forward_stress(results: pd.DataFrame) -> pd.DataFrame:
    required = {
        "scenario",
        "success",
        "provisional_directional_signal",
        "best_model",
    }
    if not required.issubset(results.columns):
        raise ValueError("forward stress results are missing required columns.")

    rows = []
    for scenario, group in results.groupby("scenario", sort=False):
        successful = group[group["success"].astype(bool)]
        row = {
            "scenario": scenario,
            "attempted_replicates": len(group),
            "successful_replicates": len(successful),
            "success_rate": float(len(successful) / len(group)),
        }
        if successful.empty:
            row.update(
                {
                    "false_direction_signal_rate": np.nan,
                    "asymmetric_model_selected_rate": np.nan,
                    "secondary_contact_selected_rate": np.nan,
                }
            )
        else:
            row.update(
                {
                    "false_direction_signal_rate": float(
                        successful["provisional_directional_signal"]
                        .astype(bool)
                        .mean()
                    ),
                    "asymmetric_model_selected_rate": float(
                        (successful["best_model"] == "asymmetric_migration").mean()
                    ),
                    "secondary_contact_selected_rate": float(
                        (
                            successful["best_model"]
                            == "secondary_contact_asymmetric"
                        ).mean()
                    ),
                }
            )
        rows.append(row)
    return pd.DataFrame(rows)


def write_forward_stress_benchmark(
    *,
    output_dir: str | Path,
    **kwargs,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    raw = run_forward_stress_benchmark(**kwargs)
    summary = summarize_forward_stress(raw)
    raw.to_csv(outdir / "forward_stress_replicates.csv", index=False)
    summary.to_csv(outdir / "forward_stress_summary.csv", index=False)
    return raw, summary
