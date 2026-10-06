"""One-generation transition-likelihood estimator.

This module is a benchmark estimator for controlled validation experiments.
It assumes allele frequencies are observed before migration and after exactly
one generation of migration plus Wright-Fisher drift, and that effective
population sizes are known.

It is NOT yet an estimator for ordinary single-time-point population-genomic
datasets. Its role is to verify DIFLOW's parameter conventions and recovery
logic before more realistic demographic inference is added.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import xlogy, xlog1py


@dataclass(frozen=True)
class TransitionEstimate:
    """Estimated forward-time migration rates for two populations."""

    m_a_to_b: float
    m_b_to_a: float
    log_likelihood: float
    success: bool


def _as_frequency_array(values, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional.")
    if arr.size == 0:
        raise ValueError(f"{name} must not be empty.")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values.")
    if np.any((arr < 0.0) | (arr > 1.0)):
        raise ValueError(f"{name} must contain allele frequencies in [0, 1].")
    return arr


def _binomial_loglikelihood(
    m: float,
    source_initial: np.ndarray,
    recipient_initial: np.ndarray,
    recipient_final: np.ndarray,
    ne_recipient: int,
) -> float:
    """Binomial transition log-likelihood up to an additive constant."""
    expected = (1.0 - m) * recipient_initial + m * source_initial
    expected = np.clip(expected, 1e-12, 1.0 - 1e-12)

    chromosomes = 2 * ne_recipient
    counts = np.rint(recipient_final * chromosomes).astype(int)
    counts = np.clip(counts, 0, chromosomes)

    ll = xlogy(counts, expected) + xlog1py(chromosomes - counts, -expected)
    return float(np.sum(ll))


def _estimate_direction(
    source_initial: np.ndarray,
    recipient_initial: np.ndarray,
    recipient_final: np.ndarray,
    ne_recipient: int,
) -> tuple[float, float, bool]:
    if ne_recipient <= 0:
        raise ValueError("effective population sizes must be positive.")

    objective = lambda m: -_binomial_loglikelihood(
        m,
        source_initial,
        recipient_initial,
        recipient_final,
        ne_recipient,
    )

    result = minimize_scalar(
        objective,
        bounds=(0.0, 1.0),
        method="bounded",
        options={"xatol": 1e-10},
    )
    return float(result.x), float(-result.fun), bool(result.success)


def estimate_one_generation(
    *,
    p_a_initial,
    p_b_initial,
    p_a_final,
    p_b_final,
    ne_a: int,
    ne_b: int,
) -> TransitionEstimate:
    """Estimate asymmetric migration after one migration-drift generation.

    m_a_to_b is the fraction of population B replaced by migrants from A.
    m_b_to_a is the fraction of population A replaced by migrants from B.
    All frequency vectors must refer to the same loci in the same order.
    """
    p_a_initial = _as_frequency_array(p_a_initial, "p_a_initial")
    p_b_initial = _as_frequency_array(p_b_initial, "p_b_initial")
    p_a_final = _as_frequency_array(p_a_final, "p_a_final")
    p_b_final = _as_frequency_array(p_b_final, "p_b_final")

    n = p_a_initial.size
    if not all(arr.size == n for arr in (p_b_initial, p_a_final, p_b_final)):
        raise ValueError("all allele-frequency arrays must have equal length.")

    m_a_to_b, ll_b, ok_b = _estimate_direction(
        p_a_initial,
        p_b_initial,
        p_b_final,
        ne_b,
    )
    m_b_to_a, ll_a, ok_a = _estimate_direction(
        p_b_initial,
        p_a_initial,
        p_a_final,
        ne_a,
    )

    return TransitionEstimate(
        m_a_to_b=m_a_to_b,
        m_b_to_a=m_b_to_a,
        log_likelihood=ll_a + ll_b,
        success=ok_a and ok_b,
    )
