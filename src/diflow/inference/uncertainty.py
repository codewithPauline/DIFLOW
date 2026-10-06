"""Bootstrap uncertainty for directional migration estimates."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .transition import estimate_one_generation


@dataclass(frozen=True)
class BootstrapSummary:
    """Locus-bootstrap uncertainty summary for two directional migration rates."""

    m_a_to_b_mean: float
    m_a_to_b_lower: float
    m_a_to_b_upper: float
    m_b_to_a_mean: float
    m_b_to_a_lower: float
    m_b_to_a_upper: float
    probability_a_to_b_stronger: float
    replicates: int


def bootstrap_one_generation(
    *,
    p_a_initial,
    p_b_initial,
    p_a_final,
    p_b_final,
    ne_a: int,
    ne_b: int,
    replicates: int = 200,
    confidence: float = 0.95,
    seed: int | None = None,
) -> BootstrapSummary:
    """Estimate locus-bootstrap uncertainty for the transition benchmark.

    Loci are resampled with replacement. The returned directional probability
    is the fraction of bootstrap replicates in which m(A->B) exceeds m(B->A).

    This procedure assumes loci are independent. A future genomic-data layer
    will support block bootstrap resampling for linked markers.
    """
    if replicates < 2:
        raise ValueError("replicates must be at least 2.")
    if not 0.0 < confidence < 1.0:
        raise ValueError("confidence must lie strictly between 0 and 1.")

    arrays = [
        np.asarray(p_a_initial, dtype=float),
        np.asarray(p_b_initial, dtype=float),
        np.asarray(p_a_final, dtype=float),
        np.asarray(p_b_final, dtype=float),
    ]
    n = arrays[0].size
    if n == 0 or any(arr.ndim != 1 or arr.size != n for arr in arrays):
        raise ValueError("all allele-frequency arrays must be non-empty and equal length.")

    rng = np.random.default_rng(seed)
    estimates = np.empty((replicates, 2), dtype=float)

    for b in range(replicates):
        idx = rng.integers(0, n, size=n)
        fit = estimate_one_generation(
            p_a_initial=arrays[0][idx],
            p_b_initial=arrays[1][idx],
            p_a_final=arrays[2][idx],
            p_b_final=arrays[3][idx],
            ne_a=ne_a,
            ne_b=ne_b,
        )
        if not fit.success:
            raise RuntimeError("migration optimization failed during bootstrap.")
        estimates[b, 0] = fit.m_a_to_b
        estimates[b, 1] = fit.m_b_to_a

    alpha = 1.0 - confidence
    lower_q = alpha / 2.0
    upper_q = 1.0 - lower_q

    return BootstrapSummary(
        m_a_to_b_mean=float(np.mean(estimates[:, 0])),
        m_a_to_b_lower=float(np.quantile(estimates[:, 0], lower_q)),
        m_a_to_b_upper=float(np.quantile(estimates[:, 0], upper_q)),
        m_b_to_a_mean=float(np.mean(estimates[:, 1])),
        m_b_to_a_lower=float(np.quantile(estimates[:, 1], lower_q)),
        m_b_to_a_upper=float(np.quantile(estimates[:, 1], upper_q)),
        probability_a_to_b_stronger=float(np.mean(estimates[:, 0] > estimates[:, 1])),
        replicates=replicates,
    )
