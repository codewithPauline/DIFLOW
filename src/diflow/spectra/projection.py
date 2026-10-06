"""Hypergeometric projection for pairwise joint site-frequency spectra."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import hypergeom


@dataclass(frozen=True)
class ProjectedPairwiseJSFS:
    """Projected two-population joint site-frequency spectrum."""

    population_a: str
    population_b: str
    chromosomes_a: int
    chromosomes_b: int
    spectrum: np.ndarray
    loci_used: int
    loci_skipped: int


def _projection_probabilities(
    alt_count: int,
    called_chromosomes: int,
    target_chromosomes: int,
) -> np.ndarray:
    if target_chromosomes > called_chromosomes:
        raise ValueError("target chromosome count cannot exceed called chromosomes.")

    k = np.arange(target_chromosomes + 1)
    return hypergeom.pmf(
        k,
        called_chromosomes,
        alt_count,
        target_chromosomes,
    )


def pairwise_projected_jsfs(
    counts: pd.DataFrame,
    population_a: str,
    population_b: str,
    *,
    chromosomes_a: int,
    chromosomes_b: int,
) -> ProjectedPairwiseJSFS:
    """Project loci with unequal sample sizes to a common two-population jSFS.

    Each locus contributes probability mass rather than a single integer count.
    Projection uses the hypergeometric distribution, matching random sampling
    without replacement from the observed chromosomes at that locus.
    """
    required = {
        "chrom",
        "pos",
        "ref",
        "alt",
        "population",
        "alt_count",
        "called_chromosomes",
    }
    if not required.issubset(counts.columns):
        raise ValueError("counts table is missing required allele-count columns.")
    if population_a == population_b:
        raise ValueError("population_a and population_b must be different.")
    if chromosomes_a < 1 or chromosomes_b < 1:
        raise ValueError("projection chromosome counts must be positive.")

    sub = counts[counts["population"].isin([population_a, population_b])].copy()
    if sub.empty:
        raise ValueError("requested populations are absent from the counts table.")

    key = ["chrom", "pos", "ref", "alt"]
    alt = sub.pivot_table(index=key, columns="population", values="alt_count", aggfunc="first")
    called = sub.pivot_table(
        index=key,
        columns="population",
        values="called_chromosomes",
        aggfunc="first",
    )

    if population_a not in alt.columns or population_b not in alt.columns:
        raise ValueError("both requested populations must be represented.")

    spectrum = np.zeros((chromosomes_a + 1, chromosomes_b + 1), dtype=float)
    loci_used = 0

    for idx in alt.index:
        if (
            pd.isna(alt.at[idx, population_a])
            or pd.isna(alt.at[idx, population_b])
            or pd.isna(called.at[idx, population_a])
            or pd.isna(called.at[idx, population_b])
        ):
            continue

        n_a = int(called.at[idx, population_a])
        n_b = int(called.at[idx, population_b])
        if n_a < chromosomes_a or n_b < chromosomes_b:
            continue

        x_a = int(alt.at[idx, population_a])
        x_b = int(alt.at[idx, population_b])

        p_a = _projection_probabilities(x_a, n_a, chromosomes_a)
        p_b = _projection_probabilities(x_b, n_b, chromosomes_b)
        spectrum += np.outer(p_a, p_b)
        loci_used += 1

    return ProjectedPairwiseJSFS(
        population_a=population_a,
        population_b=population_b,
        chromosomes_a=chromosomes_a,
        chromosomes_b=chromosomes_b,
        spectrum=spectrum,
        loci_used=loci_used,
        loci_skipped=int(len(alt.index) - loci_used),
    )
