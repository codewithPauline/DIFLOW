"""Isolated pairwise inference worker for serial or process-parallel execution."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from diflow.decision import DirectionEvidence, classify_direction
from diflow.demography import (
    ModelScore,
    bootstrap_asymmetric_jsfs,
    directional_support_for_estimate,
    fit_multistart,
    rank_models,
)
from diflow.spectra import pairwise_projected_jsfs


@dataclass(frozen=True)
class PairInferenceTask:
    pair_index: int
    population_a: str
    population_b: str
    distance_km: float
    counts: pd.DataFrame
    projection_chromosomes: int
    starts: int
    maxiter: int
    bootstrap_replicates: int
    bootstrap_starts: int
    bootstrap_block_bp: int | None
    min_model_weight: float
    min_directional_support: float
    min_abs_asymmetry: float
    polarized: bool
    prepare_only: bool
    seed: int | None


@dataclass(frozen=True)
class PairInferenceOutput:
    pair_index: int
    row: dict
    rankings: list[dict]
    spectrum: np.ndarray


def _provisional_direction(
    *,
    population_a: str,
    population_b: str,
    m_a_to_b: float,
    m_b_to_a: float,
    asymmetric_weight: float,
    stable: bool,
    min_model_weight: float,
    min_abs_asymmetry: float,
) -> tuple[str, str | None, float]:
    total = m_a_to_b + m_b_to_a
    asymmetry = 0.0 if total == 0 else (m_a_to_b - m_b_to_a) / total

    if asymmetric_weight < min_model_weight:
        return "unsupported", None, asymmetry

    if m_a_to_b > m_b_to_a:
        preferred = f"{population_a}->{population_b}"
    elif m_b_to_a > m_a_to_b:
        preferred = f"{population_b}->{population_a}"
    else:
        preferred = None

    if stable and abs(asymmetry) >= min_abs_asymmetry and preferred is not None:
        return "candidate", preferred, asymmetry
    return "ambiguous", preferred, asymmetry


def infer_pair_task(task: PairInferenceTask) -> PairInferenceOutput:
    """Fit one candidate pair without touching shared filesystem state."""
    pop_a = task.population_a
    pop_b = task.population_b

    projected = pairwise_projected_jsfs(
        task.counts,
        pop_a,
        pop_b,
        chromosomes_a=task.projection_chromosomes,
        chromosomes_b=task.projection_chromosomes,
    )

    spectrum_file = f"spectra/{pop_a}__{pop_b}.npy"
    base = {
        "population_a": pop_a,
        "population_b": pop_b,
        "distance_km": float(task.distance_km),
        "loci_used": int(projected.loci_used),
        "loci_skipped": int(projected.loci_skipped),
        "projection_chromosomes": int(task.projection_chromosomes),
        "spectrum_file": spectrum_file,
    }

    if projected.loci_used == 0:
        return PairInferenceOutput(
            pair_index=task.pair_index,
            row={**base, "status": "no_usable_loci"},
            rankings=[],
            spectrum=projected.spectrum,
        )

    if task.prepare_only:
        return PairInferenceOutput(
            pair_index=task.pair_index,
            row={**base, "status": "prepared"},
            rankings=[],
            spectrum=projected.spectrum,
        )

    model_names = (
        "isolation",
        "symmetric_migration",
        "asymmetric_migration",
        "secondary_contact_symmetric",
        "secondary_contact_asymmetric",
    )
    scores = []
    fits = {}
    observations = int(np.count_nonzero(projected.spectrum))

    for model_offset, model_name in enumerate(model_names):
        result = fit_multistart(
            projected.spectrum,
            model_name,
            starts=task.starts,
            seed=(
                None
                if task.seed is None
                else task.seed + task.pair_index * 100 + model_offset
            ),
            maxiter=task.maxiter,
            polarized=task.polarized,
        )
        fits[model_name] = result
        scores.append(
            ModelScore(
                name=model_name,
                log_likelihood=result.best_log_likelihood,
                parameters=len(result.best_parameters),
                observations=observations,
            )
        )

    ranking = rank_models(scores, use_aicc=False)
    best_model = str(ranking.iloc[0]["model"])
    directional_models = {
        "asymmetric_migration",
        "secondary_contact_asymmetric",
    }
    directional_model = (
        best_model if best_model in directional_models else None
    )

    if directional_model is None:
        row = {
            **base,
            "status": "unsupported",
            "best_model": best_model,
            "directional_model": None,
            "preferred_direction": None,
            "m_a_to_b_scaled": np.nan,
            "m_b_to_a_scaled": np.nan,
            "asymmetry_index": np.nan,
            "asymmetric_model_weight": 0.0,
            "optimizer_stable": False,
            "optimizer_success_fraction": np.nan,
            "decision_reason": (
                "best-supported demographic model is not asymmetric"
            ),
        }
        ranked = ranking.copy()
        ranked.insert(0, "population_b", pop_b)
        ranked.insert(0, "population_a", pop_a)
        return PairInferenceOutput(
            pair_index=task.pair_index,
            row=row,
            rankings=ranked.to_dict(orient="records"),
            spectrum=projected.spectrum,
        )

    directional_row = ranking[
        ranking["model"] == directional_model
    ].iloc[0]
    asym_weight = float(directional_row["akaike_weight"])
    asym_fit = fits[directional_model]
    m_a_to_b = float(asym_fit.best_parameters["m_a_to_b"])
    m_b_to_a = float(asym_fit.best_parameters["m_b_to_a"])

    status, preferred, asymmetry = _provisional_direction(
        population_a=pop_a,
        population_b=pop_b,
        m_a_to_b=m_a_to_b,
        m_b_to_a=m_b_to_a,
        asymmetric_weight=asym_weight,
        stable=asym_fit.stable,
        min_model_weight=task.min_model_weight,
        min_abs_asymmetry=task.min_abs_asymmetry,
    )

    row = {
        **base,
        "status": status,
        "best_model": best_model,
        "directional_model": directional_model,
        "preferred_direction": preferred,
        "m_a_to_b_scaled": m_a_to_b,
        "m_b_to_a_scaled": m_b_to_a,
        "asymmetry_index": asymmetry,
        "asymmetric_model_weight": asym_weight,
        "optimizer_stable": bool(asym_fit.stable),
        "optimizer_success_fraction": float(asym_fit.converged_fraction),
    }

    if task.bootstrap_replicates >= 2:
        boot = bootstrap_asymmetric_jsfs(
            task.counts,
            pop_a,
            pop_b,
            chromosomes_a=task.projection_chromosomes,
            chromosomes_b=task.projection_chromosomes,
            replicates=task.bootstrap_replicates,
            starts=task.bootstrap_starts,
            maxiter=task.maxiter,
            seed=(
                None
                if task.seed is None
                else task.seed + task.pair_index * 10000
            ),
            block_size_bp=task.bootstrap_block_bp,
            polarized=task.polarized,
            model_name=directional_model,
        )

        point_direction_support = directional_support_for_estimate(
            boot.probability_a_to_b_stronger,
            m_a_to_b,
            m_b_to_a,
        )

        evidence = DirectionEvidence(
            source=pop_a,
            destination=pop_b,
            migration_forward=m_a_to_b,
            migration_reverse=m_b_to_a,
            directional_support=point_direction_support,
            asymmetric_model_weight=asym_weight,
            optimizer_stable=bool(asym_fit.stable),
            forward_lower=boot.m_a_to_b_lower,
            forward_upper=boot.m_a_to_b_upper,
            reverse_lower=boot.m_b_to_a_lower,
            reverse_upper=boot.m_b_to_a_upper,
        )
        decision = classify_direction(
            evidence,
            min_model_weight=task.min_model_weight,
            min_directional_support=task.min_directional_support,
            min_abs_asymmetry=task.min_abs_asymmetry,
            require_interval_separation=True,
        )

        row.update(
            {
                "status": decision.status,
                "preferred_direction": decision.preferred_direction,
                "directional_support": point_direction_support,
                "bootstrap_preferred_direction_support": boot.preferred_direction_support,
                "probability_a_to_b_stronger": boot.probability_a_to_b_stronger,
                "m_a_to_b_lower": boot.m_a_to_b_lower,
                "m_a_to_b_upper": boot.m_a_to_b_upper,
                "m_b_to_a_lower": boot.m_b_to_a_lower,
                "m_b_to_a_upper": boot.m_b_to_a_upper,
                "bootstrap_successful": boot.successful_replicates,
                "bootstrap_attempted": boot.attempted_replicates,
                "bootstrap_resampling_unit": boot.resampling_unit,
                "bootstrap_blocks_used": boot.blocks_used,
                "bootstrap_loci_used": boot.loci_used,
                "decision_reason": decision.reason,
            }
        )

    ranked = ranking.copy()
    ranked.insert(0, "population_b", pop_b)
    ranked.insert(0, "population_a", pop_a)

    return PairInferenceOutput(
        pair_index=task.pair_index,
        row=row,
        rankings=ranked.to_dict(orient="records"),
        spectrum=projected.spectrum,
    )
