"""Information-criterion model comparison for demographic fits."""

from __future__ import annotations

from dataclasses import dataclass

import math
import pandas as pd


@dataclass(frozen=True)
class ModelScore:
    """Likelihood summary for one demographic model."""

    name: str
    log_likelihood: float
    parameters: int
    observations: int | None = None


def aic(log_likelihood: float, parameters: int) -> float:
    """Akaike Information Criterion."""
    if parameters < 0:
        raise ValueError("parameters must be non-negative.")
    return 2.0 * parameters - 2.0 * float(log_likelihood)


def aicc(log_likelihood: float, parameters: int, observations: int) -> float:
    """Small-sample corrected AIC."""
    if observations <= parameters + 1:
        return math.inf
    base = aic(log_likelihood, parameters)
    correction = (2.0 * parameters * (parameters + 1)) / (
        observations - parameters - 1
    )
    return base + correction


def rank_models(scores: list[ModelScore], *, use_aicc: bool = True) -> pd.DataFrame:
    """Rank fitted demographic models and calculate Akaike weights."""
    if not scores:
        raise ValueError("at least one model score is required.")

    rows = []
    for score in scores:
        if use_aicc:
            if score.observations is None:
                raise ValueError("observations are required for AICc ranking.")
            criterion = aicc(
                score.log_likelihood,
                score.parameters,
                score.observations,
            )
            criterion_name = "aicc"
        else:
            criterion = aic(score.log_likelihood, score.parameters)
            criterion_name = "aic"

        rows.append(
            {
                "model": score.name,
                "log_likelihood": float(score.log_likelihood),
                "parameters": int(score.parameters),
                criterion_name: float(criterion),
            }
        )

    frame = pd.DataFrame(rows)
    best = frame[criterion_name].min()
    frame["delta"] = frame[criterion_name] - best

    finite = frame["delta"].map(math.isfinite)
    raw = frame["delta"].copy()
    frame["akaike_weight"] = 0.0
    if finite.any():
        weights = (-0.5 * raw[finite]).map(math.exp)
        frame.loc[finite, "akaike_weight"] = weights / weights.sum()

    return frame.sort_values(
        [criterion_name, "model"],
        kind="stable",
    ).reset_index(drop=True)
