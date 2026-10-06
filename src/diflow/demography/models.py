"""Demographic parameter definitions for directional migration."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class AsymmetricIMParams:
    """Two-population split model with continuous asymmetric migration.

    All quantities are in the scaled units used by dadi.

    Parameters
    ----------
    nu_a, nu_b
        Population sizes after the split relative to the reference size.
    split_time
        Time since split in units of 2*N_ref generations.
    m_a_to_b
        Scaled migration from population A into B.
    m_b_to_a
        Scaled migration from population B into A.

    Notes
    -----
    DIFLOW always names migration in forward-time source-to-recipient order.
    This differs from dadi's split_asym_mig argument labels:

    dadi m12 = population 2 -> population 1 = B -> A
    dadi m21 = population 1 -> population 2 = A -> B
    """

    nu_a: float
    nu_b: float
    split_time: float
    m_a_to_b: float
    m_b_to_a: float

    def __post_init__(self) -> None:
        for name in ("nu_a", "nu_b", "split_time", "m_a_to_b", "m_b_to_a"):
            value = float(getattr(self, name))
            if value <= 0:
                raise ValueError(f"{name} must be > 0 for log-scale optimization.")


def to_dadi_split_asym_mig(params: AsymmetricIMParams) -> tuple[float, ...]:
    """Translate DIFLOW parameters into dadi split_asym_mig order.

    dadi expects (nu1, nu2, T, m12, m21), where m12 is migration
    from population 2 to population 1 and m21 is the reverse.
    """
    return (
        float(params.nu_a),
        float(params.nu_b),
        float(params.split_time),
        float(params.m_b_to_a),
        float(params.m_a_to_b),
    )
