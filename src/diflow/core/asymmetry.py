"""Directional gene-flow summary statistics."""

from __future__ import annotations

import math


def _validate_rate(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite.")
    if value < 0:
        raise ValueError(f"{name} must be non-negative.")
    return value


def asymmetry_index(m_ij: float, m_ji: float) -> float:
    """Return normalized directional asymmetry in [-1, 1]."""
    m_ij = _validate_rate(m_ij, "m_ij")
    m_ji = _validate_rate(m_ji, "m_ji")
    total = m_ij + m_ji
    if total == 0:
        return 0.0
    return (m_ij - m_ji) / total


def directional_ratio(m_ij: float, m_ji: float) -> float:
    """Return the directional migration-rate ratio m_ij / m_ji."""
    m_ij = _validate_rate(m_ij, "m_ij")
    m_ji = _validate_rate(m_ji, "m_ji")
    if m_ji == 0:
        return math.inf if m_ij > 0 else 1.0
    return m_ij / m_ji
