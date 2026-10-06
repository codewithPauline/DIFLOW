"""Convert directional decisions into standardized map-edge tables."""

from __future__ import annotations

import pandas as pd

from .decision import DirectionDecision


def decisions_to_flow_table(
    decisions: list[DirectionDecision],
    *,
    include_ambiguous: bool = False,
) -> pd.DataFrame:
    """Convert decisions into a standardized source-to-recipient flow table."""
    rows = []

    for decision in decisions:
        if decision.status == "unsupported":
            continue
        if decision.status == "ambiguous" and not include_ambiguous:
            continue
        if decision.preferred_direction is None:
            continue

        if decision.migration >= decision.reverse_migration:
            source = decision.source
            destination = decision.destination
            migration = decision.migration
        else:
            source = decision.destination
            destination = decision.source
            migration = decision.reverse_migration

        rows.append(
            {
                "source": source,
                "destination": destination,
                "migration": float(migration),
                "support": float(decision.directional_support),
                "status": decision.status,
                "asymmetry_index": float(decision.asymmetry_index),
                "model_weight": float(decision.asymmetric_model_weight),
                "reason": decision.reason,
            }
        )

    return pd.DataFrame(
        rows,
        columns=[
            "source",
            "destination",
            "migration",
            "support",
            "status",
            "asymmetry_index",
            "model_weight",
            "reason",
        ],
    )
