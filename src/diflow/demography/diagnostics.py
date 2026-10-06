"""Convergence summaries for DIFLOW optimization replicates."""

from __future__ import annotations

import pandas as pd


def summarize_multistart(result) -> pd.DataFrame:
    """Return a one-row summary suitable for reports and QC tables."""
    return pd.DataFrame(
        [
            {
                "model": result.model,
                "runs": int(len(result.runs)),
                "successful_runs": int(result.runs["success"].sum()),
                "converged_fraction": float(result.converged_fraction),
                "best_log_likelihood": float(result.best_log_likelihood),
                "stable": bool(result.stable),
            }
        ]
    )
