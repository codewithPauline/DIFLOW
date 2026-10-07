"""End-to-end known-truth evidence benchmark for threshold calibration."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from diflow.demography import (
    AsymmetricIMParams,
    ModelScore,
    bootstrap_asymmetric_jsfs,
    directional_support_for_estimate,
    expected_spectrum,
    fit_multistart,
    rank_models,
)


@dataclass(frozen=True)
class DecisionEvidenceScenario:
    """Known truth used to calibrate the final directional classifier."""

    name: str
    params: AsymmetricIMParams
    truth_direction: str


def default_decision_evidence_scenarios() -> tuple[DecisionEvidenceScenario, ...]:
    base = dict(nu_a=1.0, nu_b=1.0, split_time=0.75)
    return (
        DecisionEvidenceScenario(
            "symmetric",
            AsymmetricIMParams(**base, m_a_to_b=0.5, m_b_to_a=0.5),
            "symmetric",
        ),
        DecisionEvidenceScenario(
            "weak_a_to_b",
            AsymmetricIMParams(**base, m_a_to_b=0.75, m_b_to_a=0.5),
            "A->B",
        ),
        DecisionEvidenceScenario(
            "moderate_a_to_b",
            AsymmetricIMParams(**base, m_a_to_b=1.0, m_b_to_a=0.25),
            "A->B",
        ),
        DecisionEvidenceScenario(
            "strong_a_to_b",
            AsymmetricIMParams(**base, m_a_to_b=2.0, m_b_to_a=0.10),
            "A->B",
        ),
        DecisionEvidenceScenario(
            "weak_b_to_a",
            AsymmetricIMParams(**base, m_a_to_b=0.5, m_b_to_a=0.75),
            "B->A",
        ),
        DecisionEvidenceScenario(
            "moderate_b_to_a",
            AsymmetricIMParams(**base, m_a_to_b=0.25, m_b_to_a=1.0),
            "B->A",
        ),
        DecisionEvidenceScenario(
            "strong_b_to_a",
            AsymmetricIMParams(**base, m_a_to_b=0.10, m_b_to_a=2.0),
            "B->A",
        ),
    )


def simulate_counts_from_expected_spectrum(
    params: AsymmetricIMParams,
    *,
    chromosomes: int = 20,
    segregating_sites: int = 5000,
    seed: int = 42,
) -> pd.DataFrame:
    """Generate a long-form allele-count table from a known dadi expectation."""
    if chromosomes < 2:
        raise ValueError("chromosomes must be at least 2.")
    if segregating_sites < 2:
        raise ValueError("segregating_sites must be at least 2.")

    expected = np.asarray(
        expected_spectrum(params, (chromosomes, chromosomes)),
        dtype=float,
    )
    probs = expected.copy()
    probs[0, 0] = 0.0
    probs[-1, -1] = 0.0
    total = float(probs.sum())
    if total <= 0:
        raise ValueError("expected spectrum has no variable-site mass.")
    probs /= total

    rng = np.random.default_rng(seed)
    draws = rng.choice(
        probs.size,
        size=segregating_sites,
        replace=True,
        p=probs.ravel(),
    )

    rows: list[list] = []
    for site_index, flat_index in enumerate(draws, start=1):
        i, j = np.unravel_index(int(flat_index), probs.shape)
        rows.append(
            ["1", site_index, "0", "1", "A", chromosomes - i, i, chromosomes]
        )
        rows.append(
            ["1", site_index, "0", "1", "B", chromosomes - j, j, chromosomes]
        )

    return pd.DataFrame(
        rows,
        columns=[
            "chrom",
            "pos",
            "ref",
            "alt",
            "population",
            "ref_count",
            "alt_count",
            "called_chromosomes",
        ],
    )


def _spectrum_from_counts(
    counts: pd.DataFrame,
    *,
    chromosomes: int,
) -> np.ndarray:
    spectrum = np.zeros((chromosomes + 1, chromosomes + 1), dtype=float)
    a = counts[counts["population"] == "A"].set_index("pos")
    b = counts[counts["population"] == "B"].set_index("pos")
    for pos in a.index.intersection(b.index):
        i = int(a.at[pos, "alt_count"])
        j = int(b.at[pos, "alt_count"])
        spectrum[i, j] += 1.0
    return spectrum


def run_decision_evidence_benchmark(
    *,
    scenarios: tuple[DecisionEvidenceScenario, ...] | None = None,
    replicates: int = 10,
    chromosomes: int = 20,
    segregating_sites: int = 5000,
    starts: int = 10,
    bootstrap_replicates: int = 100,
    bootstrap_starts: int = 5,
    maxiter: int = 100,
    seed: int = 42,
) -> pd.DataFrame:
    """Generate full evidence columns needed for threshold calibration."""
    if replicates < 1:
        raise ValueError("replicates must be at least 1.")
    if bootstrap_replicates < 2:
        raise ValueError("bootstrap_replicates must be at least 2.")

    scenarios = scenarios or default_decision_evidence_scenarios()
    model_names = (
        "isolation",
        "symmetric_migration",
        "asymmetric_migration",
        "secondary_contact_symmetric",
        "secondary_contact_asymmetric",
    )
    parameter_counts = {
        "isolation": 3,
        "symmetric_migration": 4,
        "asymmetric_migration": 5,
        "secondary_contact_symmetric": 5,
        "secondary_contact_asymmetric": 6,
    }

    rows: list[dict] = []

    for scenario_index, scenario in enumerate(scenarios):
        for replicate in range(1, replicates + 1):
            run_seed = seed + scenario_index * 100000 + replicate
            counts = simulate_counts_from_expected_spectrum(
                scenario.params,
                chromosomes=chromosomes,
                segregating_sites=segregating_sites,
                seed=run_seed,
            )
            spectrum = _spectrum_from_counts(
                counts,
                chromosomes=chromosomes,
            )
            observations = int(np.count_nonzero(spectrum))
            row = {
                "scenario": scenario.name,
                "truth_direction": scenario.truth_direction,
                "replicate": replicate,
                "true_m_a_to_b": scenario.params.m_a_to_b,
                "true_m_b_to_a": scenario.params.m_b_to_a,
                "seed": run_seed,
            }

            try:
                fits = {}
                scores = []
                for offset, model_name in enumerate(model_names):
                    fit = fit_multistart(
                        spectrum,
                        model_name,
                        starts=starts,
                        maxiter=maxiter,
                        seed=run_seed + offset * 1000,
                        polarized=False,
                    )
                    fits[model_name] = fit
                    scores.append(
                        ModelScore(
                            name=model_name,
                            log_likelihood=fit.best_log_likelihood,
                            parameters=parameter_counts[model_name],
                            observations=observations,
                        )
                    )

                ranking = rank_models(scores, use_aicc=False)
                asym_row = ranking[
                    ranking["model"] == "asymmetric_migration"
                ].iloc[0]
                asym_fit = fits["asymmetric_migration"]

                m_ab = float(asym_fit.best_parameters["m_a_to_b"])
                m_ba = float(asym_fit.best_parameters["m_b_to_a"])
                total = m_ab + m_ba
                asymmetry = 0.0 if total == 0 else (m_ab - m_ba) / total
                preferred = (
                    "A->B" if m_ab > m_ba else
                    "B->A" if m_ba > m_ab else
                    "symmetric"
                )

                boot = bootstrap_asymmetric_jsfs(
                    counts,
                    "A",
                    "B",
                    chromosomes_a=chromosomes,
                    chromosomes_b=chromosomes,
                    replicates=bootstrap_replicates,
                    starts=bootstrap_starts,
                    maxiter=maxiter,
                    seed=run_seed + 500000,
                )

                interval_separated = bool(
                    boot.m_a_to_b_lower > boot.m_b_to_a_upper
                    or boot.m_b_to_a_lower > boot.m_a_to_b_upper
                )
                point_direction_support = directional_support_for_estimate(
                    boot.probability_a_to_b_stronger,
                    m_ab,
                    m_ba,
                )

                row.update(
                    {
                        "success": True,
                        "best_model": str(ranking.iloc[0]["model"]),
                        "asymmetric_model_weight": float(
                            asym_row["akaike_weight"]
                        ),
                        "optimizer_stable": bool(asym_fit.stable),
                        "estimated_m_a_to_b": m_ab,
                        "estimated_m_b_to_a": m_ba,
                        "asymmetry_index": asymmetry,
                        "preferred_direction": preferred,
                        "directional_support": float(point_direction_support),
                        "bootstrap_preferred_direction_support": float(
                            boot.preferred_direction_support
                        ),
                        "interval_separated": interval_separated,
                        "m_a_to_b_lower": boot.m_a_to_b_lower,
                        "m_a_to_b_upper": boot.m_a_to_b_upper,
                        "m_b_to_a_lower": boot.m_b_to_a_lower,
                        "m_b_to_a_upper": boot.m_b_to_a_upper,
                    }
                )
            except Exception as exc:
                row.update({"success": False, "error": str(exc)})
            rows.append(row)

    return pd.DataFrame(rows)


def write_decision_evidence_benchmark(
    *,
    output_dir: str | Path,
    **kwargs,
) -> pd.DataFrame:
    """Run and save the end-to-end evidence table used for calibration."""
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    result = run_decision_evidence_benchmark(**kwargs)
    result.to_csv(outdir / "decision_evidence.csv", index=False)
    return result
