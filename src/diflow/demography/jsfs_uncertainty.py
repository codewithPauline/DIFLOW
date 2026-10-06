"""Bootstrap uncertainty for single-time-point pairwise jSFS inference."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Callable

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from .multistart import fit_multistart


@dataclass(frozen=True)
class JSFSBootstrapResult:
    """Bootstrap uncertainty for asymmetric migration parameters."""

    m_a_to_b_mean: float
    m_a_to_b_lower: float
    m_a_to_b_upper: float
    m_b_to_a_mean: float
    m_b_to_a_lower: float
    m_b_to_a_upper: float
    probability_a_to_b_stronger: float
    preferred_direction_support: float
    successful_replicates: int
    attempted_replicates: int
    resampling_unit: str
    blocks_used: int
    loci_used: int


def _locus_contribution_records(
    counts: pd.DataFrame,
    population_a: str,
    population_b: str,
    chromosomes_a: int,
    chromosomes_b: int,
) -> list[tuple[str, int, np.ndarray]]:
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

    sub = counts[counts["population"].isin([population_a, population_b])].copy()
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

    records: list[tuple[str, int, np.ndarray]] = []

    for idx in alt.index:
        values = (
            alt.at[idx, population_a],
            alt.at[idx, population_b],
            called.at[idx, population_a],
            called.at[idx, population_b],
        )
        if any(pd.isna(v) for v in values):
            continue

        x_a, x_b, n_a, n_b = map(int, values)
        if n_a < chromosomes_a or n_b < chromosomes_b:
            continue

        k_a = np.arange(chromosomes_a + 1)
        k_b = np.arange(chromosomes_b + 1)
        p_a = hypergeom.pmf(k_a, n_a, x_a, chromosomes_a)
        p_b = hypergeom.pmf(k_b, n_b, x_b, chromosomes_b)
        chrom = str(idx[0])
        pos = int(idx[1])
        records.append((chrom, pos, np.outer(p_a, p_b)))

    return records


def _resampling_blocks(
    records: list[tuple[str, int, np.ndarray]],
    *,
    block_size_bp: int | None,
) -> tuple[list[np.ndarray], str]:
    """Collapse locus contributions into resampling units.

    Without block_size_bp, every usable locus is sampled independently.
    With block_size_bp, all usable variants in the same chromosome/window are
    summed and resampled together.
    """
    if block_size_bp is None:
        return [record[2] for record in records], "locus"

    if block_size_bp < 1:
        raise ValueError("block_size_bp must be a positive integer.")

    grouped: dict[tuple[str, int], np.ndarray] = {}
    for chrom, pos, contribution in records:
        # 1-based VCF positions are assigned to zero-based fixed windows.
        block_index = (pos - 1) // block_size_bp
        key = (chrom, block_index)
        if key not in grouped:
            grouped[key] = np.zeros_like(contribution, dtype=float)
        grouped[key] += contribution

    return list(grouped.values()), f"{block_size_bp}-bp genomic block"


def bootstrap_asymmetric_jsfs(
    counts: pd.DataFrame,
    population_a: str,
    population_b: str,
    *,
    chromosomes_a: int,
    chromosomes_b: int,
    replicates: int = 100,
    starts: int = 5,
    maxiter: int = 100,
    confidence: float = 0.95,
    seed: int | None = None,
    block_size_bp: int | None = None,
    polarized: bool = False,
    fit_function: Callable | None = None,
) -> JSFSBootstrapResult:
    """Bootstrap asymmetric migration estimates from projected jSFS data.

    By default, usable loci are resampled independently. When block_size_bp is
    supplied, all variants within each fixed chromosome/window are resampled as
    one unit, preserving local linkage within bootstrap blocks.

    The default fitter is DIFLOW's multi-start dadi optimizer. fit_function is
    injectable to enable lightweight unit testing and alternative backends.
    """
    if replicates < 2:
        raise ValueError("replicates must be at least 2.")
    if starts < 1:
        raise ValueError("starts must be at least 1.")
    if not 0.0 < confidence < 1.0:
        raise ValueError("confidence must lie strictly between 0 and 1.")
    if chromosomes_a < 1 or chromosomes_b < 1:
        raise ValueError("projection chromosome counts must be positive.")
    if block_size_bp is not None and block_size_bp < 1:
        raise ValueError("block_size_bp must be a positive integer.")

    records = _locus_contribution_records(
        counts,
        population_a,
        population_b,
        chromosomes_a,
        chromosomes_b,
    )
    if len(records) < 2:
        raise ValueError("at least two usable loci are required for bootstrap.")

    blocks, resampling_unit = _resampling_blocks(
        records,
        block_size_bp=block_size_bp,
    )
    if len(blocks) < 2:
        raise ValueError(
            "at least two resampling blocks are required; reduce block_size_bp "
            "or provide more genomic data."
        )

    rng = np.random.default_rng(seed)
    fitter = fit_multistart if fit_function is None else fit_function
    estimates: list[tuple[float, float]] = []

    for replicate in range(replicates):
        sampled = rng.integers(0, len(blocks), size=len(blocks))
        spectrum = np.zeros(
            (chromosomes_a + 1, chromosomes_b + 1),
            dtype=float,
        )
        for index in sampled:
            spectrum += blocks[int(index)]

        try:
            fit_kwargs = {
                "starts": starts,
                "seed": None if seed is None else seed + replicate + 1,
                "maxiter": maxiter,
            }
            if fit_function is None:
                fit_kwargs["polarized"] = polarized
            fit = fitter(
                spectrum,
                "asymmetric_migration",
                **fit_kwargs,
            )
            m_ab = float(fit.best_parameters["m_a_to_b"])
            m_ba = float(fit.best_parameters["m_b_to_a"])
            if math.isfinite(m_ab) and math.isfinite(m_ba):
                estimates.append((m_ab, m_ba))
        except Exception:
            continue

    if len(estimates) < 2:
        raise RuntimeError("fewer than two bootstrap replicates completed successfully.")

    values = np.asarray(estimates, dtype=float)
    alpha = 1.0 - confidence
    low_q = alpha / 2.0
    high_q = 1.0 - low_q
    p_ab = float(np.mean(values[:, 0] > values[:, 1]))

    mean_ab = float(np.mean(values[:, 0]))
    mean_ba = float(np.mean(values[:, 1]))
    preferred_support = p_ab if mean_ab >= mean_ba else 1.0 - p_ab

    return JSFSBootstrapResult(
        m_a_to_b_mean=mean_ab,
        m_a_to_b_lower=float(np.quantile(values[:, 0], low_q)),
        m_a_to_b_upper=float(np.quantile(values[:, 0], high_q)),
        m_b_to_a_mean=mean_ba,
        m_b_to_a_lower=float(np.quantile(values[:, 1], low_q)),
        m_b_to_a_upper=float(np.quantile(values[:, 1], high_q)),
        probability_a_to_b_stronger=p_ab,
        preferred_direction_support=float(preferred_support),
        successful_replicates=len(estimates),
        attempted_replicates=replicates,
        resampling_unit=resampling_unit,
        blocks_used=len(blocks),
        loci_used=len(records),
    )
