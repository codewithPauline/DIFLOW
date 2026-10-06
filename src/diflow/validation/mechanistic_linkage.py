"""Mechanistic linkage validation using msprime/tskit ancestry simulations.

This layer complements the correlated-block simulator by generating recombining
ancestral genealogies and mutations explicitly. It therefore provides a stronger
test of locus versus block bootstrap calibration under realistic genealogical
linkage.

Migration convention
--------------------
DIFLOW reports forward-time source-to-recipient migration.

msprime specifies continuous migration backwards in time. Therefore forward
A -> B corresponds to a backwards-time lineage transition B -> A.

dadi reports scaled migration M = 2 * N_ref * m, so known msprime per-generation
migration rates are converted to dadi-scaled truth before coverage is assessed.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from diflow.demography import bootstrap_asymmetric_jsfs


@dataclass(frozen=True)
class MechanisticLinkageScenario:
    """Known demographic truth for msprime linkage validation."""

    name: str
    m_a_to_b: float
    m_b_to_a: float
    expected_direction: str


def default_mechanistic_linkage_scenarios() -> tuple[MechanisticLinkageScenario, ...]:
    return (
        MechanisticLinkageScenario("symmetric", 2.5e-5, 2.5e-5, "symmetric"),
        MechanisticLinkageScenario("moderate_a_to_b", 5e-5, 1.25e-5, "A->B"),
        MechanisticLinkageScenario("moderate_b_to_a", 1.25e-5, 5e-5, "B->A"),
    )


def _require_msprime():
    try:
        import msprime
    except ImportError as exc:
        raise ImportError(
            "Mechanistic linkage validation requires msprime. Install with: "
            "python -m pip install -e '.[validation]'"
        ) from exc
    return msprime


def simulate_msprime_counts(
    scenario: MechanisticLinkageScenario,
    *,
    chromosomes_per_population: int = 20,
    nref: int = 10000,
    split_time_scaled: float = 0.75,
    sequence_length: int = 2_000_000,
    recombination_rate: float = 1e-8,
    mutation_rate: float = 1e-8,
    seed: int = 42,
) -> pd.DataFrame:
    """Simulate recombining biallelic markers for two populations.

    The two daughter populations have size nref and split from an ancestor of
    size nref at split_time_scaled * 2 * nref generations in the past.
    """
    if chromosomes_per_population < 2:
        raise ValueError("chromosomes_per_population must be at least 2.")
    if nref < 2:
        raise ValueError("nref must be at least 2.")
    if split_time_scaled <= 0:
        raise ValueError("split_time_scaled must be positive.")
    if sequence_length < 1000:
        raise ValueError("sequence_length must be at least 1000.")
    if recombination_rate < 0 or mutation_rate <= 0:
        raise ValueError("recombination_rate must be >= 0 and mutation_rate > 0.")

    msprime = _require_msprime()

    demography = msprime.Demography()
    demography.add_population(name="ANC", initial_size=nref)
    demography.add_population(name="A", initial_size=nref)
    demography.add_population(name="B", initial_size=nref)

    # msprime migration is backwards in time:
    # forward A -> B == backward lineage B -> A
    demography.set_migration_rate(
        source="B", dest="A", rate=float(scenario.m_a_to_b)
    )
    demography.set_migration_rate(
        source="A", dest="B", rate=float(scenario.m_b_to_a)
    )

    split_generations = float(split_time_scaled) * 2.0 * float(nref)
    demography.add_population_split(
        time=split_generations,
        derived=["A", "B"],
        ancestral="ANC",
    )

    samples = [
        msprime.SampleSet(
            num_samples=chromosomes_per_population,
            population="A",
            ploidy=1,
        ),
        msprime.SampleSet(
            num_samples=chromosomes_per_population,
            population="B",
            ploidy=1,
        ),
    ]

    ts = msprime.sim_ancestry(
        samples=samples,
        demography=demography,
        sequence_length=sequence_length,
        recombination_rate=recombination_rate,
        ploidy=1,
        random_seed=seed,
    )
    mts = msprime.sim_mutations(
        ts,
        rate=mutation_rate,
        model=msprime.BinaryMutationModel(),
        random_seed=seed + 1,
    )

    sample_nodes = np.asarray(mts.samples(), dtype=int)
    sample_populations = np.asarray(
        [mts.node(int(node)).population for node in sample_nodes],
        dtype=int,
    )
    a_id = demography["A"].id
    b_id = demography["B"].id
    a_indices = np.flatnonzero(sample_populations == a_id)
    b_indices = np.flatnonzero(sample_populations == b_id)

    rows: list[list] = []
    seen_positions: set[int] = set()

    for variant in mts.variants():
        # BinaryMutationModel should be biallelic, but keep this guard explicit.
        if len(variant.alleles) != 2:
            continue

        pos = int(variant.site.position) + 1
        if pos in seen_positions:
            continue
        seen_positions.add(pos)

        genotypes = np.asarray(variant.genotypes, dtype=int)
        if np.any(genotypes < 0) or np.any(genotypes > 1):
            continue

        a_alt = int(genotypes[a_indices].sum())
        b_alt = int(genotypes[b_indices].sum())
        n_a = len(a_indices)
        n_b = len(b_indices)

        # Match variant-only VCF behavior by excluding globally fixed sites.
        if (a_alt == 0 and b_alt == 0) or (a_alt == n_a and b_alt == n_b):
            continue

        rows.append(["1", pos, "0", "1", "A", n_a - a_alt, a_alt, n_a])
        rows.append(["1", pos, "0", "1", "B", n_b - b_alt, b_alt, n_b])

    if not rows:
        raise RuntimeError(
            "msprime simulation produced no usable variable sites; increase "
            "sequence_length or mutation_rate."
        )

    return pd.DataFrame(
        rows,
        columns=[
            "chrom",
            "pos",
            "ref",
            "alt",
            "population",
            "ref_count",
            "alt_count",
            "called_chromosomes",
        ],
    )


def _asymmetry(m_ab: float, m_ba: float) -> float:
    total = m_ab + m_ba
    return 0.0 if total == 0 else (m_ab - m_ba) / total


def run_mechanistic_linkage_calibration(
    *,
    scenarios: tuple[MechanisticLinkageScenario, ...] | None = None,
    replicates: int = 10,
    chromosomes_per_population: int = 20,
    nref: int = 10000,
    split_time_scaled: float = 0.75,
    sequence_length: int = 2_000_000,
    recombination_rate: float = 1e-8,
    mutation_rate: float = 1e-8,
    block_sizes_bp: tuple[int, ...] = (50_000, 100_000, 250_000),
    bootstrap_replicates: int = 100,
    bootstrap_starts: int = 5,
    maxiter: int = 100,
    confidence: float = 0.95,
    directional_support_threshold: float = 0.95,
    asymmetry_threshold: float = 0.25,
    seed: int = 42,
) -> pd.DataFrame:
    """Calibrate locus and block bootstrap using recombining genealogies."""
    if replicates < 1:
        raise ValueError("replicates must be at least 1.")
    if bootstrap_replicates < 2:
        raise ValueError("bootstrap_replicates must be at least 2.")
    if not block_sizes_bp or any(int(x) < 1 for x in block_sizes_bp):
        raise ValueError("block_sizes_bp must contain positive integers.")

    scenarios = scenarios or default_mechanistic_linkage_scenarios()
    rows: list[dict] = []

    for scenario_index, scenario in enumerate(scenarios):
        true_scaled_ab = 2.0 * nref * scenario.m_a_to_b
        true_scaled_ba = 2.0 * nref * scenario.m_b_to_a

        for replicate in range(1, replicates + 1):
            run_seed = seed + scenario_index * 100000 + replicate
            counts = simulate_msprime_counts(
                scenario,
                chromosomes_per_population=chromosomes_per_population,
                nref=nref,
                split_time_scaled=split_time_scaled,
                sequence_length=sequence_length,
                recombination_rate=recombination_rate,
                mutation_rate=mutation_rate,
                seed=run_seed,
            )

            methods: list[tuple[str, int | None]] = [("locus", None)]
            methods.extend((f"block_{bp}", bp) for bp in block_sizes_bp)

            for method, block_bp in methods:
                row = {
                    "scenario": scenario.name,
                    "expected_direction": scenario.expected_direction,
                    "replicate": replicate,
                    "method": method,
                    "block_size_bp": block_bp,
                    "true_m_a_to_b_scaled": true_scaled_ab,
                    "true_m_b_to_a_scaled": true_scaled_ba,
                    "nref": nref,
                    "sequence_length": sequence_length,
                    "recombination_rate": recombination_rate,
                    "mutation_rate": mutation_rate,
                    "variant_loci": int(counts["pos"].nunique()),
                    "seed": run_seed,
                }

                try:
                    boot = bootstrap_asymmetric_jsfs(
                        counts,
                        "A",
                        "B",
                        chromosomes_a=chromosomes_per_population,
                        chromosomes_b=chromosomes_per_population,
                        replicates=bootstrap_replicates,
                        starts=bootstrap_starts,
                        maxiter=maxiter,
                        confidence=confidence,
                        seed=run_seed
                        + (0 if block_bp is None else int(block_bp)),
                        block_size_bp=block_bp,
                    )
                    asym = _asymmetry(
                        boot.m_a_to_b_mean,
                        boot.m_b_to_a_mean,
                    )
                    strong = bool(
                        boot.preferred_direction_support
                        >= directional_support_threshold
                        and abs(asym) >= asymmetry_threshold
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
                                <= true_scaled_ab
                                <= boot.m_a_to_b_upper
                            ),
                            "coverage_b_to_a": (
                                boot.m_b_to_a_lower
                                <= true_scaled_ba
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
                            "estimated_asymmetry": asym,
                            "strong_directional_support": strong,
                            "blocks_used": boot.blocks_used,
                            "loci_used": boot.loci_used,
                        }
                    )
                except Exception as exc:
                    row.update({"success": False, "error": str(exc)})
                rows.append(row)

    return pd.DataFrame(rows)


def summarize_mechanistic_linkage(results: pd.DataFrame) -> pd.DataFrame:
    """Summarize coverage, interval width, and false direction by method."""
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
        raise ValueError("mechanistic linkage results are missing required columns.")

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
            strong_rate = float(
                successful["strong_directional_support"].astype(bool).mean()
            )
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
                    "strong_directional_support_rate": strong_rate,
                    "false_directional_support_rate": (
                        strong_rate if expected == "symmetric" else np.nan
                    ),
                }
            )
        rows.append(row)

    return pd.DataFrame(rows)


def plot_mechanistic_linkage(
    summary: pd.DataFrame,
    output_dir: str | Path,
) -> list[Path]:
    """Render coverage and false-direction comparisons across block sizes."""
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []

    if summary.empty:
        return paths

    for metric, ylabel, stem in (
        ("coverage_a_to_b", "Coverage of true scaled m(A→B)", "mechanistic_coverage_ab"),
        ("coverage_b_to_a", "Coverage of true scaled m(B→A)", "mechanistic_coverage_ba"),
    ):
        pivot = summary.pivot(index="scenario", columns="method", values=metric)
        ax = pivot.plot(kind="bar", figsize=(9, 5))
        ax.axhline(0.95, linestyle="--", linewidth=1)
        ax.set_ylim(0, 1.05)
        ax.set_ylabel(ylabel)
        ax.set_title("Bootstrap calibration with recombining genealogies")
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
        ].plot(kind="bar", figsize=(7, 5))
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("False strong-direction support rate")
        ax.set_title("False direction with recombining genealogies")
        fig = ax.figure
        fig.tight_layout()
        for suffix in ("png", "pdf"):
            path = outdir / f"mechanistic_false_direction.{suffix}"
            fig.savefig(path, dpi=300 if suffix == "png" else None, bbox_inches="tight")
            paths.append(path)
        plt.close(fig)

    return paths


def write_mechanistic_linkage_calibration(
    *,
    output_dir: str | Path,
    **kwargs,
) -> tuple[pd.DataFrame, pd.DataFrame, list[Path]]:
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    raw = run_mechanistic_linkage_calibration(**kwargs)
    summary = summarize_mechanistic_linkage(raw)
    raw.to_csv(outdir / "mechanistic_linkage_replicates.csv", index=False)
    summary.to_csv(outdir / "mechanistic_linkage_summary.csv", index=False)
    figures = plot_mechanistic_linkage(summary, outdir)
    return raw, summary, figures
