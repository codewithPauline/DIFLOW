"""End-to-end orchestration for pairwise directional gene-flow inference."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from diflow.decision import DirectionEvidence, classify_direction
from diflow.demography import (
    ModelScore,
    bootstrap_asymmetric_jsfs,
    fit_multistart,
    rank_models,
)
from diflow.io import allele_counts_from_vcf, read_popmap
from diflow.network import build_candidate_pairs
from diflow.spectra import pairwise_projected_jsfs
from .outputs import write_network_outputs


@dataclass(frozen=True)
class PipelineResult:
    candidate_pairs: pd.DataFrame
    pairwise_results: pd.DataFrame
    model_rankings: pd.DataFrame
    output_dir: Path


def _read_coordinates(path: str | Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    required = {"population", "latitude", "longitude"}
    if not required.issubset(frame.columns):
        raise ValueError(
            "coordinates file must contain population, latitude, longitude columns."
        )
    return frame.loc[:, ["population", "latitude", "longitude"]].copy()


def _provisional_direction(
    *,
    population_a: str,
    population_b: str,
    m_a_to_b: float,
    m_b_to_a: float,
    asymmetric_weight: float,
    stable: bool,
    min_model_weight: float = 0.70,
    min_abs_asymmetry: float = 0.25,
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


def run_infer_pipeline(
    *,
    vcf_path: str | Path,
    popmap_path: str | Path,
    coordinates_path: str | Path,
    output_dir: str | Path,
    projection_chromosomes: int,
    k_nearest: int | None = None,
    max_distance_km: float | None = None,
    starts: int = 10,
    maxiter: int = 100,
    bootstrap_replicates: int = 0,
    bootstrap_starts: int = 5,
    bootstrap_block_bp: int | None = None,
    prepare_only: bool = False,
    seed: int | None = None,
) -> PipelineResult:
    if projection_chromosomes < 2:
        raise ValueError("projection_chromosomes must be at least 2.")
    if starts < 1:
        raise ValueError("starts must be at least 1.")
    if bootstrap_replicates != 0 and bootstrap_replicates < 2:
        raise ValueError("bootstrap_replicates must be 0 or at least 2.")
    if bootstrap_starts < 1:
        raise ValueError("bootstrap_starts must be at least 1.")
    if bootstrap_block_bp is not None and bootstrap_block_bp < 1:
        raise ValueError("bootstrap_block_bp must be a positive integer.")

    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    spectra_dir = outdir / "spectra"
    spectra_dir.mkdir(exist_ok=True)

    popmap = read_popmap(popmap_path)
    coordinates = _read_coordinates(coordinates_path)

    mapped_pops = set(popmap["population"].astype(str))
    coord_pops = set(coordinates["population"].astype(str))
    missing_coords = sorted(mapped_pops - coord_pops)
    if missing_coords:
        raise ValueError(
            "populations in popmap missing coordinates: " + ", ".join(missing_coords)
        )

    counts = allele_counts_from_vcf(vcf_path, popmap)
    counts.to_csv(outdir / "allele_counts.csv", index=False)

    pairs = build_candidate_pairs(
        coordinates[coordinates["population"].astype(str).isin(mapped_pops)],
        max_distance_km=max_distance_km,
        k_nearest=k_nearest,
    )
    pairs.to_csv(outdir / "candidate_pairs.csv", index=False)

    pair_rows = []
    ranking_rows = []
    model_names = (
        "isolation",
        "symmetric_migration",
        "asymmetric_migration",
        "secondary_contact_asymmetric",
    )

    for pair_index, pair in enumerate(pairs.itertuples(index=False), start=1):
        pop_a = str(pair.population_a)
        pop_b = str(pair.population_b)

        projected = pairwise_projected_jsfs(
            counts,
            pop_a,
            pop_b,
            chromosomes_a=projection_chromosomes,
            chromosomes_b=projection_chromosomes,
        )

        spectrum_path = spectra_dir / f"{pop_a}__{pop_b}.npy"
        np.save(spectrum_path, projected.spectrum)

        base = {
            "population_a": pop_a,
            "population_b": pop_b,
            "distance_km": float(pair.distance_km),
            "loci_used": int(projected.loci_used),
            "loci_skipped": int(projected.loci_skipped),
            "projection_chromosomes": int(projection_chromosomes),
            "spectrum_file": str(spectrum_path.relative_to(outdir)),
        }

        if projected.loci_used == 0:
            pair_rows.append({**base, "status": "no_usable_loci"})
            continue

        if prepare_only:
            pair_rows.append({**base, "status": "prepared"})
            continue

        scores = []
        fits = {}
        observations = int(np.count_nonzero(projected.spectrum))

        for model_offset, model_name in enumerate(model_names):
            result = fit_multistart(
                projected.spectrum,
                model_name,
                starts=starts,
                seed=None if seed is None else seed + pair_index * 100 + model_offset,
                maxiter=maxiter,
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
        asym_row = ranking[ranking["model"] == "asymmetric_migration"].iloc[0]
        asym_weight = float(asym_row["akaike_weight"])
        asym_fit = fits["asymmetric_migration"]
        m_a_to_b = float(asym_fit.best_parameters["m_a_to_b"])
        m_b_to_a = float(asym_fit.best_parameters["m_b_to_a"])

        status, preferred, asymmetry = _provisional_direction(
            population_a=pop_a,
            population_b=pop_b,
            m_a_to_b=m_a_to_b,
            m_b_to_a=m_b_to_a,
            asymmetric_weight=asym_weight,
            stable=asym_fit.stable,
        )

        row = {
            **base,
            "status": status,
            "best_model": best_model,
            "preferred_direction": preferred,
            "m_a_to_b_scaled": m_a_to_b,
            "m_b_to_a_scaled": m_b_to_a,
            "asymmetry_index": asymmetry,
            "asymmetric_model_weight": asym_weight,
            "optimizer_stable": bool(asym_fit.stable),
            "optimizer_success_fraction": float(asym_fit.converged_fraction),
        }

        if bootstrap_replicates >= 2:
            boot = bootstrap_asymmetric_jsfs(
                counts,
                pop_a,
                pop_b,
                chromosomes_a=projection_chromosomes,
                chromosomes_b=projection_chromosomes,
                replicates=bootstrap_replicates,
                starts=bootstrap_starts,
                maxiter=maxiter,
                seed=None if seed is None else seed + pair_index * 10000,
                block_size_bp=bootstrap_block_bp,
            )

            evidence = DirectionEvidence(
                source=pop_a,
                destination=pop_b,
                migration_forward=m_a_to_b,
                migration_reverse=m_b_to_a,
                directional_support=boot.preferred_direction_support,
                asymmetric_model_weight=asym_weight,
                optimizer_stable=bool(asym_fit.stable),
                forward_lower=boot.m_a_to_b_lower,
                forward_upper=boot.m_a_to_b_upper,
                reverse_lower=boot.m_b_to_a_lower,
                reverse_upper=boot.m_b_to_a_upper,
            )
            decision = classify_direction(
                evidence,
                require_interval_separation=True,
            )

            row.update(
                {
                    "status": decision.status,
                    "preferred_direction": decision.preferred_direction,
                    "directional_support": boot.preferred_direction_support,
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

        pair_rows.append(row)

        ranked = ranking.copy()
        ranked.insert(0, "population_b", pop_b)
        ranked.insert(0, "population_a", pop_a)
        ranking_rows.extend(ranked.to_dict(orient="records"))

    pairwise = pd.DataFrame(pair_rows)
    rankings = pd.DataFrame(ranking_rows)

    pairwise.to_csv(outdir / "pairwise_results.csv", index=False)
    rankings.to_csv(outdir / "model_rankings.csv", index=False)

    if bootstrap_replicates >= 2 and not prepare_only:
        write_network_outputs(
            pairwise=pairwise,
            coordinates=coordinates,
            output_dir=outdir,
        )

    metadata = pd.DataFrame(
        [
            {
                "projection_chromosomes": projection_chromosomes,
                "k_nearest": k_nearest,
                "max_distance_km": max_distance_km,
                "starts": starts,
                "maxiter": maxiter,
                "bootstrap_replicates": bootstrap_replicates,
                "bootstrap_starts": bootstrap_starts,
                "bootstrap_block_bp": bootstrap_block_bp,
                "bootstrap_resampling_unit": (
                    "locus" if bootstrap_block_bp is None
                    else f"{bootstrap_block_bp}-bp genomic block"
                ),
                "prepare_only": prepare_only,
                "direction_status_note": (
                    "supported/ambiguous/unsupported uses bootstrap uncertainty "
                    "only when bootstrap_replicates >= 2; otherwise candidate "
                    "labels remain provisional"
                ),
            }
        ]
    )
    metadata.to_csv(outdir / "run_metadata.csv", index=False)

    return PipelineResult(
        candidate_pairs=pairs,
        pairwise_results=pairwise,
        model_rankings=rankings,
        output_dir=outdir,
    )
