"""Dataset inspection and user-facing analysis recommendations."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shlex

import numpy as np
import pandas as pd

from diflow.io import allele_counts_from_vcf, read_popmap
from diflow.network import build_candidate_pairs


PRESETS = {
    "quick": {
        "starts": 5,
        "bootstrap_replicates": 20,
        "bootstrap_starts": 2,
        "maxiter": 60,
    },
    "standard": {
        "starts": 20,
        "bootstrap_replicates": 100,
        "bootstrap_starts": 5,
        "maxiter": 100,
    },
    "publication": {
        "starts": 40,
        "bootstrap_replicates": 500,
        "bootstrap_starts": 8,
        "maxiter": 200,
    },
}


@dataclass(frozen=True)
class InspectionResult:
    samples: int
    populations: int
    loci: int
    recommended_projection: int
    projection_summary: pd.DataFrame
    graph_summary: pd.DataFrame
    recommended_neighbors: int
    recommended_preset: str
    command: str


def read_coordinates(path: str | Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    required = {"population", "latitude", "longitude"}
    if not required.issubset(frame.columns):
        raise ValueError(
            "coordinates file must contain population, latitude, longitude columns."
        )
    frame = frame.loc[:, ["population", "latitude", "longitude"]].copy()
    if frame["population"].duplicated().any():
        raise ValueError("coordinate populations must be unique.")
    return frame


def _candidate_projection_values(max_chromosomes: int) -> list[int]:
    if max_chromosomes < 2:
        return []
    values = list(range(2, max_chromosomes + 1, 2))
    if not values:
        values = [2]
    return values


def projection_retention_summary(counts: pd.DataFrame) -> pd.DataFrame:
    """Summarize shared pairwise locus retention across projection sizes.

    DIFLOW fits pairwise jSFS models, so the relevant quantity is not merely how
    many loci each population retains independently. It is how many loci are
    jointly usable for both populations in each pair.
    """
    if counts.empty:
        raise ValueError("no allele-count records were produced from the VCF.")

    max_chromosomes = int(counts["called_chromosomes"].max())
    projections = _candidate_projection_values(max_chromosomes)
    if not projections:
        raise ValueError(
            "fewer than two called chromosomes are available; a pairwise jSFS "
            "projection cannot be recommended."
        )

    populations = sorted(counts["population"].astype(str).unique())
    if len(populations) < 2:
        raise ValueError("at least two populations are required for pairwise inference.")

    key = ["chrom", "pos", "ref", "alt"]
    total_loci = counts[key].drop_duplicates().shape[0]
    called = counts.pivot_table(
        index=key,
        columns="population",
        values="called_chromosomes",
        aggfunc="first",
    )

    rows = []
    for projection in projections:
        pair_counts = []
        for i, pop_a in enumerate(populations):
            for pop_b in populations[i + 1 :]:
                if pop_a not in called.columns or pop_b not in called.columns:
                    pair_counts.append(0)
                    continue
                shared = (
                    called[[pop_a, pop_b]].fillna(-1).ge(projection).all(axis=1)
                )
                pair_counts.append(int(shared.sum()))

        minimum_pair = min(pair_counts) if pair_counts else 0
        median_pair = float(np.median(pair_counts)) if pair_counts else 0.0

        rows.append(
            {
                "projection_chromosomes": projection,
                "total_loci": total_loci,
                "population_pairs": len(pair_counts),
                "minimum_pair_loci": minimum_pair,
                "median_pair_loci": median_pair,
                "minimum_pair_retention": (
                    minimum_pair / total_loci if total_loci else 0.0
                ),
                "median_pair_retention": (
                    median_pair / total_loci if total_loci else 0.0
                ),
            }
        )

    return pd.DataFrame(rows)


def recommend_projection(summary: pd.DataFrame, *, retention_target: float = 0.80) -> int:
    """Choose the largest even projection retaining enough loci across populations.

    The default heuristic seeks at least 80% shared-locus retention in the
    worst-retained population pair. If no projection meets that target, the
    smallest available projection is returned and the inspection table exposes
    the tradeoff.
    """
    if summary.empty:
        raise ValueError("projection summary cannot be empty.")
    eligible = summary[
        summary["minimum_pair_retention"] >= retention_target
    ]
    if not eligible.empty:
        return int(eligible["projection_chromosomes"].max())
    return int(summary["projection_chromosomes"].min())


def graph_option_summary(coordinates: pd.DataFrame) -> pd.DataFrame:
    """Show the computational consequences of several k-nearest graph choices."""
    n = len(coordinates)
    if n < 2:
        return pd.DataFrame(
            columns=["neighbors", "candidate_pairs", "median_distance_km", "max_distance_km"]
        )

    candidate_ks = sorted({k for k in (2, 3, 4, 5) if k < n})
    rows = []
    for k in candidate_ks:
        pairs = build_candidate_pairs(coordinates, k_nearest=k)
        rows.append(
            {
                "neighbors": k,
                "candidate_pairs": len(pairs),
                "median_distance_km": float(pairs["distance_km"].median()),
                "max_distance_km": float(pairs["distance_km"].max()),
            }
        )
    return pd.DataFrame(rows)


def recommend_neighbors(populations: int) -> int:
    """Provide a transparent starting point, not a biological truth."""
    if populations <= 3:
        return max(1, populations - 1)
    if populations <= 10:
        return 3
    return 4


def inspect_dataset(
    *,
    vcf_path: str | Path,
    popmap_path: str | Path,
    coordinates_path: str | Path,
    output_dir: str | Path | None = None,
    preset: str = "standard",
    retention_target: float = 0.80,
) -> InspectionResult:
    if preset not in PRESETS:
        raise ValueError(f"unknown preset: {preset}")

    popmap = read_popmap(popmap_path)
    coordinates = read_coordinates(coordinates_path)

    mapped = set(popmap["population"].astype(str))
    coordinate_populations = set(coordinates["population"].astype(str))
    missing = sorted(mapped - coordinate_populations)
    if missing:
        raise ValueError(
            "populations in popmap missing coordinates: " + ", ".join(missing)
        )

    counts = allele_counts_from_vcf(vcf_path, popmap)
    projection_summary = projection_retention_summary(counts)
    projection = recommend_projection(
        projection_summary,
        retention_target=retention_target,
    )

    coordinates = coordinates[
        coordinates["population"].astype(str).isin(mapped)
    ].copy()
    graph_summary = graph_option_summary(coordinates)
    neighbors = recommend_neighbors(len(mapped))

    preset_values = PRESETS[preset]
    destination = str(output_dir or "diflow_results")
    command = (
        "diflow infer "
        f"--vcf {shlex.quote(str(vcf_path))} "
        f"--popmap {shlex.quote(str(popmap_path))} "
        f"--coords {shlex.quote(str(coordinates_path))} "
        f"--projection-chromosomes {projection} "
        f"--neighbors {neighbors} "
        f"--starts {preset_values['starts']} "
        f"--bootstrap-replicates {preset_values['bootstrap_replicates']} "
        f"--bootstrap-starts {preset_values['bootstrap_starts']} "
        f"--maxiter {preset_values['maxiter']} "
        f"--seed 42 "
        f"--output {shlex.quote(destination)}"
    )

    if output_dir is not None:
        outdir = Path(output_dir)
        outdir.mkdir(parents=True, exist_ok=True)
        projection_summary.to_csv(
            outdir / "projection_recommendations.csv",
            index=False,
        )
        graph_summary.to_csv(
            outdir / "graph_options.csv",
            index=False,
        )

    loci = counts[["chrom", "pos", "ref", "alt"]].drop_duplicates().shape[0]

    return InspectionResult(
        samples=len(popmap),
        populations=len(mapped),
        loci=loci,
        recommended_projection=projection,
        projection_summary=projection_summary,
        graph_summary=graph_summary,
        recommended_neighbors=neighbors,
        recommended_preset=preset,
        command=command,
    )
