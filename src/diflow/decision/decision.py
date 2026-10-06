"""Evidence-based classification of directional migration edges."""

from __future__ import annotations

from dataclasses import dataclass

import math


@dataclass(frozen=True)
class DirectionEvidence:
    """Inputs used to decide whether a directional edge is supported."""

    source: str
    destination: str
    migration_forward: float
    migration_reverse: float
    directional_support: float
    asymmetric_model_weight: float
    optimizer_stable: bool
    forward_lower: float | None = None
    forward_upper: float | None = None
    reverse_lower: float | None = None
    reverse_upper: float | None = None


@dataclass(frozen=True)
class DirectionDecision:
    """Final classification for one source-to-recipient contrast."""

    source: str
    destination: str
    status: str
    preferred_direction: str | None
    migration: float
    reverse_migration: float
    asymmetry_index: float
    directional_support: float
    asymmetric_model_weight: float
    optimizer_stable: bool
    reason: str


def _validate_probability(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must lie within [0, 1].")
    return value


def _validate_rate(value: float, name: str) -> float:
    value = float(value)
    if not math.isfinite(value) or value < 0.0:
        raise ValueError(f"{name} must be finite and non-negative.")
    return value


def _asymmetry(m_forward: float, m_reverse: float) -> float:
    total = m_forward + m_reverse
    if total == 0:
        return 0.0
    return (m_forward - m_reverse) / total


def classify_direction(
    evidence: DirectionEvidence,
    *,
    min_model_weight: float = 0.70,
    min_directional_support: float = 0.95,
    min_abs_asymmetry: float = 0.25,
    require_interval_separation: bool = False,
) -> DirectionDecision:
    """Classify a migration contrast as supported, ambiguous, or unsupported.

    Rules are deliberately conservative and transparent.

    supported:
      asymmetric model has sufficient Akaike weight,
      optimizer is stable,
      directional support is high,
      asymmetry magnitude exceeds threshold,
      and optional uncertainty intervals separate.

    ambiguous:
      there is some directional signal, but one or more support criteria fail.

    unsupported:
      migration is effectively absent or the asymmetric model has weak support.
    """
    mf = _validate_rate(evidence.migration_forward, "migration_forward")
    mr = _validate_rate(evidence.migration_reverse, "migration_reverse")
    ds = _validate_probability(evidence.directional_support, "directional_support")
    mw = _validate_probability(evidence.asymmetric_model_weight, "asymmetric_model_weight")

    if not 0 <= min_model_weight <= 1:
        raise ValueError("min_model_weight must lie within [0, 1].")
    if not 0 <= min_directional_support <= 1:
        raise ValueError("min_directional_support must lie within [0, 1].")
    if not 0 <= min_abs_asymmetry <= 1:
        raise ValueError("min_abs_asymmetry must lie within [0, 1].")

    a = _asymmetry(mf, mr)
    preferred = None
    if mf > mr:
        preferred = f"{evidence.source}->{evidence.destination}"
    elif mr > mf:
        preferred = f"{evidence.destination}->{evidence.source}"

    interval_separated = True
    if require_interval_separation:
        bounds = (
            evidence.forward_lower,
            evidence.forward_upper,
            evidence.reverse_lower,
            evidence.reverse_upper,
        )
        if any(x is None for x in bounds):
            interval_separated = False
        else:
            fl, fu, rl, ru = map(float, bounds)
            if not all(math.isfinite(x) and x >= 0 for x in (fl, fu, rl, ru)):
                raise ValueError("uncertainty bounds must be finite and non-negative.")
            if fl > fu or rl > ru:
                raise ValueError("lower uncertainty bounds cannot exceed upper bounds.")
            interval_separated = fl > ru or rl > fu

    if mf == 0 and mr == 0:
        return DirectionDecision(
            source=evidence.source,
            destination=evidence.destination,
            status="unsupported",
            preferred_direction=None,
            migration=mf,
            reverse_migration=mr,
            asymmetry_index=a,
            directional_support=ds,
            asymmetric_model_weight=mw,
            optimizer_stable=evidence.optimizer_stable,
            reason="no migration signal",
        )

    if mw < min_model_weight:
        return DirectionDecision(
            source=evidence.source,
            destination=evidence.destination,
            status="unsupported",
            preferred_direction=preferred,
            migration=mf,
            reverse_migration=mr,
            asymmetry_index=a,
            directional_support=ds,
            asymmetric_model_weight=mw,
            optimizer_stable=evidence.optimizer_stable,
            reason="asymmetric model support below threshold",
        )

    criteria = [
        evidence.optimizer_stable,
        ds >= min_directional_support,
        abs(a) >= min_abs_asymmetry,
        interval_separated,
    ]

    if all(criteria):
        return DirectionDecision(
            source=evidence.source,
            destination=evidence.destination,
            status="supported",
            preferred_direction=preferred,
            migration=mf,
            reverse_migration=mr,
            asymmetry_index=a,
            directional_support=ds,
            asymmetric_model_weight=mw,
            optimizer_stable=evidence.optimizer_stable,
            reason="all directional evidence criteria satisfied",
        )

    failed = []
    if not evidence.optimizer_stable:
        failed.append("optimizer instability")
    if ds < min_directional_support:
        failed.append("directional support below threshold")
    if abs(a) < min_abs_asymmetry:
        failed.append("asymmetry below threshold")
    if not interval_separated:
        failed.append("uncertainty intervals overlap or are unavailable")

    return DirectionDecision(
        source=evidence.source,
        destination=evidence.destination,
        status="ambiguous",
        preferred_direction=preferred,
        migration=mf,
        reverse_migration=mr,
        asymmetry_index=a,
        directional_support=ds,
        asymmetric_model_weight=mw,
        optimizer_stable=evidence.optimizer_stable,
        reason="; ".join(failed),
    )
