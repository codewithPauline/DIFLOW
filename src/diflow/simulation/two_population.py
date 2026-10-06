"""Forward-time two-population migration-and-drift simulator.

This module provides controlled data with known directional migration rates
for testing DIFLOW estimators.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def _validate_migration_rate(value: float, name: str) -> float:
    value = float(value)
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must lie between 0 and 1.")
    return value


def simulate_two_population(
    *,
    generations: int,
    loci: int,
    ne_a: int,
    ne_b: int,
    m_a_to_b: float,
    m_b_to_a: float,
    seed: int | None = None,
) -> pd.DataFrame:
    """Simulate allele-frequency change under asymmetric migration and drift.

    Migration is forward-time and source-to-recipient. m_a_to_b is the
    fraction of population B replaced each generation by migrants from A,
    and m_b_to_a is defined analogously.

    The model tracks independent biallelic loci. Each generation applies
    migration to allele frequencies, followed by Wright-Fisher binomial drift.
    """
    if generations < 0:
        raise ValueError("generations must be non-negative.")
    if loci <= 0:
        raise ValueError("loci must be positive.")
    if ne_a <= 0 or ne_b <= 0:
        raise ValueError("effective population sizes must be positive.")

    m_a_to_b = _validate_migration_rate(m_a_to_b, "m_a_to_b")
    m_b_to_a = _validate_migration_rate(m_b_to_a, "m_b_to_a")

    rng = np.random.default_rng(seed)
    p_a0 = rng.beta(0.8, 0.8, size=loci)
    p_b0 = rng.beta(0.8, 0.8, size=loci)
    p_a = p_a0.copy()
    p_b = p_b0.copy()

    for _ in range(generations):
        migrated_a = (1.0 - m_b_to_a) * p_a + m_b_to_a * p_b
        migrated_b = (1.0 - m_a_to_b) * p_b + m_a_to_b * p_a
        p_a = rng.binomial(2 * ne_a, migrated_a) / (2 * ne_a)
        p_b = rng.binomial(2 * ne_b, migrated_b) / (2 * ne_b)

    return pd.DataFrame(
        {
            "locus": np.arange(loci),
            "p_a_initial": p_a0,
            "p_b_initial": p_b0,
            "p_a_final": p_a,
            "p_b_final": p_b,
        }
    )
