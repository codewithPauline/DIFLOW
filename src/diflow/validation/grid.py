"""Large simulation-grid recovery benchmarks for DIFLOW."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from diflow.demography import AsymmetricIMParams
from .recovery import RecoveryScenario, run_recovery_benchmark, summarize_recovery


@dataclass(frozen=True)
class RecoveryGridConfig:
    """One cell in the simulation-recovery grid."""

    chromosomes: int
    segregating_sites: int
    scenario: RecoveryScenario


def default_recovery_grid() -> tuple[RecoveryGridConfig, ...]:
    """Return a compact but informative default validation grid."""
    truths = (
        RecoveryScenario(
            "symmetric",
            AsymmetricIMParams(1.0, 1.0, 0.75, 0.5, 0.5),
            "symmetric",
        ),
        RecoveryScenario(
            "weak_a_to_b",
            AsymmetricIMParams(1.0, 1.0, 0.75, 0.75, 0.5),
            "A->B",
        ),
        RecoveryScenario(
            "moderate_a_to_b",
            AsymmetricIMParams(1.0, 1.0, 0.75, 1.0, 0.25),
            "A->B",
        ),
        RecoveryScenario(
            "strong_a_to_b",
            AsymmetricIMParams(1.0, 1.0, 0.75, 2.0, 0.10),
            "A->B",
        ),
        RecoveryScenario(
            "weak_b_to_a",
            AsymmetricIMParams(1.0, 1.0, 0.75, 0.5, 0.75),
            "B->A",
        ),
        RecoveryScenario(
            "moderate_b_to_a",
            AsymmetricIMParams(1.0, 1.0, 0.75, 0.25, 1.0),
            "B->A",
        ),
        RecoveryScenario(
            "strong_b_to_a",
            AsymmetricIMParams(1.0, 1.0, 0.75, 0.10, 2.0),
            "B->A",
        ),
    )
    chromosomes = (10, 20, 40)
    sites = (1000, 5000, 20000)

    return tuple(
        RecoveryGridConfig(n, s, truth)
        for n in chromosomes
        for s in sites
        for truth in truths
    )


def run_recovery_grid(
    *,
    configs: tuple[RecoveryGridConfig, ...] | None = None,
    replicates: int = 10,
    starts: int = 10,
    maxiter: int = 100,
    seed: int = 42,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run a factorial recovery study across sample size, SNP count, and truth."""
    configs = configs or default_recovery_grid()
    raw_frames: list[pd.DataFrame] = []
    summary_frames: list[pd.DataFrame] = []

    for index, config in enumerate(configs):
        raw = run_recovery_benchmark(
            scenarios=(config.scenario,),
            replicates=replicates,
            sample_sizes=(config.chromosomes, config.chromosomes),
            segregating_sites=config.segregating_sites,
            starts=starts,
            maxiter=maxiter,
            seed=seed + index * 1000000,
        )
        raw.insert(0, "grid_chromosomes", config.chromosomes)
        raw.insert(1, "grid_segregating_sites", config.segregating_sites)
        raw_frames.append(raw)

        summary = summarize_recovery(raw)
        summary.insert(0, "chromosomes", config.chromosomes)
        summary.insert(1, "segregating_sites", config.segregating_sites)
        summary_frames.append(summary)

    return (
        pd.concat(raw_frames, ignore_index=True),
        pd.concat(summary_frames, ignore_index=True),
    )


def plot_recovery_grid(
    summary: pd.DataFrame,
    output_dir: str | Path,
) -> list[Path]:
    """Create compact validation figures from a recovery-grid summary."""
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []

    directional = summary[summary["scenario"] != "symmetric"].copy()
    if not directional.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        for chromosomes, group in directional.groupby("chromosomes"):
            collapsed = (
                group.groupby("segregating_sites", as_index=False)["direction_accuracy"]
                .mean()
                .sort_values("segregating_sites")
            )
            ax.plot(
                collapsed["segregating_sites"],
                collapsed["direction_accuracy"],
                marker="o",
                label=f"{chromosomes} chromosomes/population",
            )
        ax.set_xscale("log")
        ax.set_ylim(0, 1.05)
        ax.set_xlabel("Segregating sites")
        ax.set_ylabel("Mean direction accuracy")
        ax.set_title("DIFLOW directional recovery across data sizes")
        ax.legend(frameon=False)
        fig.tight_layout()
        for suffix in ("png", "pdf"):
            path = outdir / f"direction_accuracy_grid.{suffix}"
            fig.savefig(path, dpi=300 if suffix == "png" else None, bbox_inches="tight")
            paths.append(path)
        plt.close(fig)

    symmetric = summary[summary["scenario"] == "symmetric"].copy()
    if not symmetric.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        for chromosomes, group in symmetric.groupby("chromosomes"):
            group = group.sort_values("segregating_sites")
            ax.plot(
                group["segregating_sites"],
                group["false_directional_positive_rate"],
                marker="o",
                label=f"{chromosomes} chromosomes/population",
            )
        ax.set_xscale("log")
        ax.set_ylim(0, 1.05)
        ax.set_xlabel("Segregating sites")
        ax.set_ylabel("False directional-positive rate")
        ax.set_title("False direction under symmetric migration")
        ax.legend(frameon=False)
        fig.tight_layout()
        for suffix in ("png", "pdf"):
            path = outdir / f"false_direction_rate_grid.{suffix}"
            fig.savefig(path, dpi=300 if suffix == "png" else None, bbox_inches="tight")
            paths.append(path)
        plt.close(fig)

    return paths


def write_recovery_grid(
    *,
    output_dir: str | Path,
    **kwargs,
) -> tuple[pd.DataFrame, pd.DataFrame, list[Path]]:
    """Run the grid, save raw/summary tables, and render validation figures."""
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    raw, summary = run_recovery_grid(**kwargs)
    raw.to_csv(outdir / "recovery_grid_replicates.csv", index=False)
    summary.to_csv(outdir / "recovery_grid_summary.csv", index=False)
    figures = plot_recovery_grid(summary, outdir)
    return raw, summary, figures
