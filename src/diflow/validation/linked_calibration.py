"""Calibration of locus versus block bootstrap under correlated marker blocks.

This module creates finite SNP datasets with known demographic truth and
within-block dependence in jSFS cell outcomes. It is designed to test whether
naive locus resampling becomes overconfident and whether genomic-block
resampling improves uncertainty calibration.

The correlation model is statistical rather than a mechanistic recombination
or haplotype simulator. It should therefore be described as correlated-block
validation, not as a full population-genetic LD simulation.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from diflow.demography import AsymmetricIMParams, bootstrap_asymmetric_jsfs, expected_spectrum


@dataclass(frozen=True)
class LinkedCalibrationScenario:
    """Known migration truth for bootstrap calibration."""

    name: str
    params: AsymmetricIMParams
    expected_direction: str


def default_linked_calibration_scenarios() -> tuple[LinkedCalibrationScenario, ...]:
    base = dict(nu_a=1.0, nu_b=1.0, split_time=0.75)
    return (
        LinkedCalibrationScenario(
            "symmetric",
            AsymmetricIMParams(**base, m_a_to_b=0.5, m_b_to_a=0.5),
            "symmetric",
        ),
        LinkedCalibrationScenario(
            "moderate_a_to_b",
            AsymmetricIMParams(**base, m_a_to_b=1.0, m_b_to_a=0.25),
            "A->B",
        ),
        LinkedCalibrationScenario(
            "moderate_b_to_a",
            AsymmetricIMParams(**base, m_a_to_b=0.25, m_b_to_a=1.0),
            "B->A",
        ),
    )


def _variant_probabilities(
    params: AsymmetricIMParams,
    sample_sizes: tuple[int, int],
    *,
    spectrum_function: Callable | None = None,
) -> np.ndarray:
    generator = expected_spectrum if spectrum_function is None else spectrum_function
    expectation = np.asarray(generator(params, sample_sizes), dtype=float)
    expected_shape = (sample_sizes[0] + 1, sample_sizes[1] + 1)
    if expectation.shape != expected_shape:
        raise ValueError(
            f"expected spectrum shape {expected_shape}, got {expectation.shape}."
        )
    if np.any(~np.isfinite(expectation)) or np.any(expectation < 0):
        raise ValueError("expected spectrum must be finite and non-negative.")

    probs = expectation.copy()
    probs[0, 0] = 0.0
    probs[-1, -1] = 0.0
    total = float(probs.sum())
    if total <= 0:
        raise ValueError("expected spectrum has no variable-site mass.")
    return probs / total


def simulate_correlated_block_counts(
    params: AsymmetricIMParams,
    *,
    sample_sizes: tuple[int, int] = (20, 20),
    blocks: int = 50,
    snps_per_block: int = 10,
    block_size_bp: int = 100000,
    concentration: float = 25.0,
    seed: int = 42,
    spectrum_function: Callable | None = None,
) -> pd.DataFrame:
    """Generate allele-count records with correlated jSFS outcomes within blocks.

    For each genomic block, a block-specific categorical distribution is drawn
    from a Dirichlet distribution centered on the demographic expected jSFS.
    SNPs in that block are then sampled from the shared latent distribution.

    Smaller concentration values create stronger block-to-block heterogeneity
    and therefore stronger within-block dependence.
    """
    if len(sample_sizes) != 2 or min(sample_sizes) < 2:
        raise ValueError("sample_sizes must contain two values >= 2.")
    if blocks < 2:
        raise ValueError("blocks must be at least 2.")
    if snps_per_block < 1:
        raise ValueError("snps_per_block must be positive.")
    if block_size_bp < snps_per_block:
        raise ValueError("block_size_bp must be at least snps_per_block.")
    if concentration <= 0:
        raise ValueError("concentration must be positive.")

    probs = _variant_probabilities(
        params,
        sample_sizes,
        spectrum_function=spectrum_function,
    )
    flat = probs.ravel()
    nonzero = flat > 0
    cell_indices = np.flatnonzero(nonzero)
    base = flat[nonzero]
    base /= base.sum()

    rng = np.random.default_rng(seed)
    rows: list[list] = []
    n_a, n_b = sample_sizes

    for block in range(blocks):
        alpha = np.maximum(base * concentration, 1e-9)
        latent = rng.dirichlet(alpha)
        draws = rng.choice(cell_indices, size=snps_per_block, p=latent)

        chrom = str(block // 10 + 1)
        block_within_chrom = block % 10
        start = block_within_chrom * block_size_bp + 1
        spacing = max(1, block_size_bp // snps_per_block)

        for snp_index, flat_index in enumerate(draws):
            i, j = np.unravel_index(int(flat_index), probs.shape)
            pos = start + snp_index * spacing
            rows.append(["A", chrom, pos, "A", "G", n_a - i, i, n_a])
            rows.append(["B", chrom, pos, "A", "G", n_b - j, j, n_b])

    frame = pd.DataFrame(
        rows,
        columns=[
            "population",
            "chrom",
            "pos",
            "ref",
            "alt",
            "ref_count",
            "alt_count",
            "called_chromosomes",
        ],
    )
    return frame[
        [
            "chrom",
            "pos",
            "ref",
            "alt",
            "population",
            "ref_count",
            "alt_count",
            "called_chromosomes",
        ]
    ]


def _asymmetry(m_ab: float, m_ba: float) -> float:
    total = m_ab + m_ba
    return 0.0 if total == 0 else (m_ab - m_ba) / total


def run_linked_bootstrap_calibration(
    *,
    scenarios: tuple[LinkedCalibrationScenario, ...] | None = None,
    replicates: int = 10,
    sample_sizes: tuple[int, int] = (20, 20),
    blocks: int = 50,
    snps_per_block: int = 10,
    block_size_bp: int = 100000,
    concentration: float = 25.0,
    bootstrap_replicates: int = 100,
    bootstrap_starts: int = 5,
    maxiter: int = 100,
    confidence: float = 0.95,
    directional_support_threshold: float = 0.95,
    asymmetry_threshold: float = 0.25,
    seed: int = 42,
    spectrum_function: Callable | None = None,
    fit_function: Callable | None = None,
) -> pd.DataFrame:
    """Compare locus and block bootstrap on the same correlated-marker datasets."""
    if replicates < 1:
        raise ValueError("replicates must be at least 1.")
    if bootstrap_replicates < 2:
        raise ValueError("bootstrap_replicates must be at least 2.")
    if not 0 < confidence < 1:
        raise ValueError("confidence must lie within (0, 1).")

    scenarios = scenarios or default_linked_calibration_scenarios()
    rows: list[dict] = []

    for scenario_index, scenario in enumerate(scenarios):
        for replicate in range(1, replicates + 1):
            run_seed = seed + scenario_index * 100000 + replicate
            counts = simulate_correlated_block_counts(
                scenario.params,
                sample_sizes=sample_sizes,
                blocks=blocks,
                snps_per_block=snps_per_block,
                block_size_bp=block_size_bp,
                concentration=concentration,
                seed=run_seed,
                spectrum_function=spectrum_function,
            )

            for method, block_bp in (
                ("locus", None),
                ("block", block_size_bp),
            ):
                row = {
                    "scenario": scenario.name,
                    "expected_direction": scenario.expected_direction,
                    "replicate": replicate,
                    "method": method,
                    "true_m_a_to_b": scenario.params.m_a_to_b,
                    "true_m_b_to_a": scenario.params.m_b_to_a,
                    "blocks_simulated": blocks,
                    "snps_per_block": snps_per_block,
                    "block_size_bp": block_size_bp,
                    "concentration": concentration,
                    "seed": run_seed,
                }

                try:
                    boot = bootstrap_asymmetric_jsfs(
                        counts,
                        "A",
                        "B",
                        chromosomes_a=sample_sizes[0],
                        chromosomes_b=sample_sizes[1],
                        replicates=bootstrap_replicates,
                        starts=bootstrap_starts,
                        maxiter=maxiter,
                        confidence=confidence,
                        seed=run_seed + (0 if method == "locus" else 500000),
                        block_size_bp=block_bp,
                        fit_function=fit_function,
                    )
                    asymmetry = _asymmetry(
                        boot.m_a_to_b_mean,
                        boot.m_b_to_a_mean,
                    )
                    strong_directional_support = bool(
                        boot.preferred_direction_support
                        >= directional_support_threshold
                        and abs(asymmetry) >= asymmetry_threshold
                    )
                    row.update(
                        {
                            "success": True,
                            "m_a_to_b_mean": boot.m_a_to_b_mean,
                            "m_a_to_b_lower": boot.m_a_to_b_lower,
                            "m_a_to_b_upper": boot.m_a_to_b_upper,
                            "m_b_to_a_mean": boot.m_b_to_a_mean,
                            "m_b_to_a_lower": boot.m_b_to_a_lower,
                            "m_b_to_a_upper": boot.m_b_to_a_upper,
                            "coverage_a_to_b": (
                                boot.m_a_to_b_lower
                                <= scenario.params.m_a_to_b
                                <= boot.m_a_to_b_upper
                            ),
                            "coverage_b_to_a": (
                                boot.m_b_to_a_lower
                                <= scenario.params.m_b_to_a
                                <= boot.m_b_to_a_upper
                            ),
                            "interval_width_a_to_b": (
                                boot.m_a_to_b_upper - boot.m_a_to_b_lower
                            ),
                            "interval_width_b_to_a": (
                                boot.m_b_to_a_upper - boot.m_b_to_a_lower
                            ),
                            "preferred_direction_support": (
                                boot.preferred_direction_support
                            ),
                            "estimated_asymmetry": asymmetry,
                            "strong_directional_support": strong_directional_support,
                            "blocks_used": boot.blocks_used,
                            "loci_used": boot.loci_used,
                        }
                    )
                except Exception as exc:
                    row.update(
                        {
                            "success": False,
                            "error": str(exc),
                        }
                    )

                rows.append(row)

    return pd.DataFrame(rows)


def summarize_linked_bootstrap_calibration(results: pd.DataFrame) -> pd.DataFrame:
    """Summarize coverage and false-direction behavior by scenario and method."""
    required = {
        "scenario",
        "method",
        "success",
        "expected_direction",
        "coverage_a_to_b",
        "coverage_b_to_a",
        "interval_width_a_to_b",
        "interval_width_b_to_a",
        "strong_directional_support",
    }
    if not required.issubset(results.columns):
        raise ValueError("linked calibration results are missing required columns.")

    rows: list[dict] = []
    for (scenario, method), group in results.groupby(
        ["scenario", "method"], sort=False
    ):
        successful = group[group["success"].astype(bool)]
        expected = str(group.iloc[0]["expected_direction"])
        row = {
            "scenario": scenario,
            "method": method,
            "attempted_replicates": len(group),
            "successful_replicates": len(successful),
            "success_rate": float(len(successful) / len(group)),
        }

        if successful.empty:
            row.update(
                {
                    "coverage_a_to_b": np.nan,
                    "coverage_b_to_a": np.nan,
                    "mean_interval_width_a_to_b": np.nan,
                    "mean_interval_width_b_to_a": np.nan,
                    "strong_directional_support_rate": np.nan,
                    "false_directional_support_rate": np.nan,
                }
            )
        else:
            row.update(
                {
                    "coverage_a_to_b": float(
                        successful["coverage_a_to_b"].astype(bool).mean()
                    ),
                    "coverage_b_to_a": float(
                        successful["coverage_b_to_a"].astype(bool).mean()
                    ),
                    "mean_interval_width_a_to_b": float(
                        successful["interval_width_a_to_b"].mean()
                    ),
                    "mean_interval_width_b_to_a": float(
                        successful["interval_width_b_to_a"].mean()
                    ),
                    "strong_directional_support_rate": float(
                        successful["strong_directional_support"]
                        .astype(bool)
                        .mean()
                    ),
                    "false_directional_support_rate": (
                        float(
                            successful["strong_directional_support"]
                            .astype(bool)
                            .mean()
                        )
                        if expected == "symmetric"
                        else np.nan
                    ),
                }
            )
        rows.append(row)

    return pd.DataFrame(rows)


def plot_linked_bootstrap_calibration(
    summary: pd.DataFrame,
    output_dir: str | Path,
) -> list[Path]:
    """Create comparison figures for locus versus block bootstrap."""
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []

    if summary.empty:
        return paths

    for metric, ylabel, stem in (
        ("coverage_a_to_b", "Coverage of true m(A→B)", "coverage_a_to_b"),
        ("coverage_b_to_a", "Coverage of true m(B→A)", "coverage_b_to_a"),
    ):
        pivot = summary.pivot(index="scenario", columns="method", values=metric)
        ax = pivot.plot(kind="bar", figsize=(8, 5))
        ax.axhline(0.95, linestyle="--", linewidth=1)
        ax.set_ylim(0, 1.05)
        ax.set_ylabel(ylabel)
        ax.set_title("Bootstrap interval calibration under correlated blocks")
        ax.legend(title="Resampling")
        fig = ax.figure
        fig.tight_layout()
        for suffix in ("png", "pdf"):
            path = outdir / f"{stem}.{suffix}"
            fig.savefig(path, dpi=300 if suffix == "png" else None, bbox_inches="tight")
            paths.append(path)
        plt.close(fig)

    symmetric = summary[summary["scenario"] == "symmetric"]
    if not symmetric.empty:
        ax = symmetric.set_index("method")[
            "false_directional_support_rate"
        ].plot(kind="bar", figsize=(6, 5))
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("False strong-direction support rate")
        ax.set_title("False direction under correlated symmetric migration")
        fig = ax.figure
        fig.tight_layout()
        for suffix in ("png", "pdf"):
            path = outdir / f"false_direction_linkage.{suffix}"
            fig.savefig(path, dpi=300 if suffix == "png" else None, bbox_inches="tight")
            paths.append(path)
        plt.close(fig)

    return paths


def write_linked_bootstrap_calibration(
    *,
    output_dir: str | Path,
    **kwargs,
) -> tuple[pd.DataFrame, pd.DataFrame, list[Path]]:
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    raw = run_linked_bootstrap_calibration(**kwargs)
    summary = summarize_linked_bootstrap_calibration(raw)
    raw.to_csv(outdir / "linked_bootstrap_replicates.csv", index=False)
    summary.to_csv(outdir / "linked_bootstrap_summary.csv", index=False)
    figures = plot_linked_bootstrap_calibration(summary, outdir)
    return raw, summary, figures
