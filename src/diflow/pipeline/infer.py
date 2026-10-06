"""End-to-end orchestration for pairwise directional gene-flow inference."""

from __future__ import annotations

from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json

import numpy as np
import pandas as pd

from diflow.io import allele_counts_from_vcf, read_popmap
from diflow.network import build_candidate_pairs
from diflow.provenance import write_provenance
from .outputs import write_network_outputs
from .pair_worker import PairInferenceTask, infer_pair_task


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
    min_model_weight: float = 0.70,
    min_directional_support: float = 0.95,
    min_abs_asymmetry: float = 0.25,
    polarized: bool = False,
    map_crs: str | None = None,
    workers: int = 1,
    prepare_only: bool = False,
    seed: int | None = None,
) -> PipelineResult:
    if projection_chromosomes < 2:
        raise ValueError("projection_chromosomes must be at least 2.")
    if starts < 1:
        raise ValueError("starts must be at least 1.")
    if workers < 1:
        raise ValueError("workers must be at least 1.")
    if bootstrap_replicates != 0 and bootstrap_replicates < 2:
        raise ValueError("bootstrap_replicates must be 0 or at least 2.")
    if bootstrap_starts < 1:
        raise ValueError("bootstrap_starts must be at least 1.")
    if bootstrap_block_bp is not None and bootstrap_block_bp < 1:
        raise ValueError("bootstrap_block_bp must be a positive integer.")
    for name, value in (
        ("min_model_weight", min_model_weight),
        ("min_directional_support", min_directional_support),
        ("min_abs_asymmetry", min_abs_asymmetry),
    ):
        if not 0 <= value <= 1:
            raise ValueError(f"{name} must lie within [0, 1].")

    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    spectra_dir = outdir / "spectra"
    spectra_dir.mkdir(exist_ok=True)

    resolved_settings = {
        "projection_chromosomes": projection_chromosomes,
        "k_nearest": k_nearest,
        "max_distance_km": max_distance_km,
        "starts": starts,
        "maxiter": maxiter,
        "bootstrap_replicates": bootstrap_replicates,
        "bootstrap_starts": bootstrap_starts,
        "bootstrap_block_bp": bootstrap_block_bp,
        "min_model_weight": min_model_weight,
        "min_directional_support": min_directional_support,
        "min_abs_asymmetry": min_abs_asymmetry,
        "polarized": polarized,
        "map_crs": map_crs,
        "workers": workers,
        "prepare_only": prepare_only,
        "seed": seed,
    }
    (outdir / "resolved_config.json").write_text(
        json.dumps(resolved_settings, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_provenance(
        outdir / "run_provenance.json",
        inputs={
            "vcf": vcf_path,
            "popmap": popmap_path,
            "coordinates": coordinates_path,
        },
        settings=resolved_settings,
    )

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

    tasks: list[PairInferenceTask] = []
    for pair_index, pair in enumerate(pairs.itertuples(index=False), start=1):
        pop_a = str(pair.population_a)
        pop_b = str(pair.population_b)
        pair_counts = counts[
            counts["population"].astype(str).isin([pop_a, pop_b])
        ].copy()
        tasks.append(
            PairInferenceTask(
                pair_index=pair_index,
                population_a=pop_a,
                population_b=pop_b,
                distance_km=float(pair.distance_km),
                counts=pair_counts,
                projection_chromosomes=projection_chromosomes,
                starts=starts,
                maxiter=maxiter,
                bootstrap_replicates=bootstrap_replicates,
                bootstrap_starts=bootstrap_starts,
                bootstrap_block_bp=bootstrap_block_bp,
                min_model_weight=min_model_weight,
                min_directional_support=min_directional_support,
                min_abs_asymmetry=min_abs_asymmetry,
                polarized=polarized,
                prepare_only=prepare_only,
                seed=seed,
            )
        )

    if workers == 1:
        outputs = [infer_pair_task(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            outputs = list(executor.map(infer_pair_task, tasks, chunksize=1))

    outputs = sorted(outputs, key=lambda result: result.pair_index)

    for result in outputs:
        spectrum_path = outdir / result.row["spectrum_file"]
        spectrum_path.parent.mkdir(parents=True, exist_ok=True)
        np.save(spectrum_path, result.spectrum)
        pair_rows.append(result.row)
        ranking_rows.extend(result.rankings)

    pairwise = pd.DataFrame(pair_rows)
    rankings = pd.DataFrame(ranking_rows)

    pairwise.to_csv(outdir / "pairwise_results.csv", index=False)
    rankings.to_csv(outdir / "model_rankings.csv", index=False)

    if bootstrap_replicates >= 2 and not prepare_only:
        write_network_outputs(
            pairwise=pairwise,
            coordinates=coordinates,
            output_dir=outdir,
            map_crs=map_crs,
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
                "min_model_weight": min_model_weight,
                "min_directional_support": min_directional_support,
                "min_abs_asymmetry": min_abs_asymmetry,
                "polarized": polarized,
                "map_crs": map_crs,
                "workers": workers,
                "spectrum_orientation": (
                    "polarized/unfolded (ALT asserted derived)"
                    if polarized
                    else "unpolarized/folded"
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
