"""Command-line interface for DIFLOW."""

from __future__ import annotations

import argparse
from pathlib import Path

from diflow.inspection import PRESETS, inspect_dataset
from diflow.pipeline import run_infer_pipeline
from diflow.validation import (
    write_forward_stress_benchmark,
    write_recovery_benchmark,
    write_recovery_grid,
    write_stress_benchmark,
    write_linked_bootstrap_calibration,
    write_mechanistic_linkage_calibration,
    write_decision_evidence_benchmark,
    write_threshold_calibration,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="diflow",
        description="Directional inference and geographic visualization of gene flow.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    inspect = subparsers.add_parser(
        "inspect",
        help="Inspect input data and recommend transparent analysis settings.",
    )
    inspect.add_argument("--vcf", required=True, help="Input plain-text VCF.")
    inspect.add_argument("--popmap", required=True, help="Sample-to-population table.")
    inspect.add_argument("--coords", required=True, help="Population coordinate CSV.")
    inspect.add_argument(
        "--preset",
        choices=sorted(PRESETS),
        default="standard",
        help="Analysis-effort preset used for the recommended infer command.",
    )
    inspect.add_argument(
        "--retention-target",
        type=float,
        default=0.80,
        help="Minimum worst-population locus retention used for projection recommendation.",
    )
    inspect.add_argument(
        "--output",
        default=None,
        help="Optional directory for inspection tables and recommended run output.",
    )


    benchmark = subparsers.add_parser(
        "benchmark",
        help="Run model-consistent known-truth migration recovery benchmarks.",
    )
    benchmark.add_argument("--output", required=True, help="Benchmark output directory.")
    benchmark.add_argument(
        "--suite",
        choices=("recovery", "stress", "forward", "grid", "linked", "mechanistic", "decision", "all"),
        default="all",
        help="Benchmark suite to run.",
    )
    benchmark.add_argument("--replicates", type=int, default=10)
    benchmark.add_argument(
        "--chromosomes",
        type=int,
        default=20,
        help="Sampled chromosomes per population in simulated spectra.",
    )
    benchmark.add_argument(
        "--sites",
        type=int,
        default=5000,
        help="Segregating sites sampled per replicate.",
    )
    benchmark.add_argument("--starts", type=int, default=10)
    benchmark.add_argument("--maxiter", type=int, default=100)
    benchmark.add_argument("--seed", type=int, default=42)
    benchmark.add_argument("--linked-blocks", type=int, default=50)
    benchmark.add_argument(
        "--linked-bootstrap-replicates",
        type=int,
        default=100,
        help="Bootstrap replicates per simulated linked-marker dataset.",
    )
    benchmark.add_argument("--snps-per-block", type=int, default=10)
    benchmark.add_argument("--linked-block-bp", type=int, default=100000)
    benchmark.add_argument(
        "--linkage-concentration",
        type=float,
        default=25.0,
        help="Smaller values create stronger within-block dependence.",
    )
    benchmark.add_argument("--mechanistic-nref", type=int, default=10000)
    benchmark.add_argument("--mechanistic-sequence-length", type=int, default=2000000)
    benchmark.add_argument("--mechanistic-recombination-rate", type=float, default=1e-8)
    benchmark.add_argument("--mechanistic-mutation-rate", type=float, default=1e-8)
    benchmark.add_argument(
        "--decision-bootstrap-replicates",
        type=int,
        default=100,
        help="Bootstrap replicates per dataset in decision-evidence calibration.",
    )
    benchmark.add_argument(
        "--mechanistic-block-sizes",
        default="50000,100000,250000",
        help="Comma-separated block sizes in bp for msprime linkage calibration.",
    )

    calibrate = subparsers.add_parser(
        "calibrate",
        help="Calibrate directional decision thresholds from known-truth evidence.",
    )
    calibrate.add_argument(
        "--evidence",
        required=True,
        help="CSV from 'diflow benchmark --suite decision'.",
    )
    calibrate.add_argument("--output", required=True)
    calibrate.add_argument(
        "--max-fpr",
        type=float,
        default=0.05,
        help="Maximum tolerated false directional-positive rate.",
    )

    infer = subparsers.add_parser(
        "infer",
        help="Prepare and fit pairwise directional demographic models.",
    )
    infer.add_argument("--vcf", required=True, help="Input plain-text VCF.")
    infer.add_argument("--popmap", required=True, help="Sample-to-population table.")
    infer.add_argument("--coords", required=True, help="Population coordinate CSV.")
    infer.add_argument("--output", required=True, help="Output directory.")
    infer.add_argument(
        "--projection-chromosomes",
        required=True,
        type=int,
        help="Chromosome count used for both populations in each projected jSFS.",
    )
    infer.add_argument("--neighbors", type=int, default=None, help="k nearest neighbors.")
    infer.add_argument("--max-distance-km", type=float, default=None)
    infer.add_argument(
        "--preset",
        choices=sorted(PRESETS),
        default=None,
        help="Set optimization/bootstrap effort; explicit flags override the preset.",
    )
    infer.add_argument("--starts", type=int, default=None, help="Optimization starts/model.")
    infer.add_argument("--maxiter", type=int, default=None)
    infer.add_argument(
        "--bootstrap-replicates",
        type=int,
        default=None,
        help="Locus-bootstrap replicates; 0 disables bootstrap.",
    )
    infer.add_argument(
        "--bootstrap-starts",
        type=int,
        default=None,
        help="Optimization starts within each bootstrap replicate.",
    )
    infer.add_argument(
        "--bootstrap-block-bp",
        type=int,
        default=None,
        help=(
            "Resample fixed genomic windows of this size in bp instead of "
            "individual loci. Recommended when nearby SNPs are linked."
        ),
    )
    infer.add_argument(
        "--min-model-weight",
        type=float,
        default=0.70,
        help="Minimum asymmetric-model Akaike weight for directional support.",
    )
    infer.add_argument(
        "--min-directional-support",
        type=float,
        default=0.95,
        help="Minimum bootstrap directional support.",
    )
    infer.add_argument(
        "--min-abs-asymmetry",
        type=float,
        default=0.25,
        help="Minimum absolute migration asymmetry index.",
    )
    infer.add_argument(
        "--polarized",
        action="store_true",
        help=(
            "Treat ALT as the derived allele and use an unfolded spectrum. "
            "Use only when ancestral state was established upstream."
        ),
    )
    infer.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed. Defaults to 42 for reproducibility.",
    )
    infer.add_argument("--prepare-only", action="store_true")

    return parser


def _resolve_effort(args):
    defaults = {
        "starts": 10,
        "bootstrap_replicates": 0,
        "bootstrap_starts": 5,
        "maxiter": 100,
    }
    if args.preset is not None:
        defaults.update(PRESETS[args.preset])

    return {
        key: getattr(args, key) if getattr(args, key) is not None else value
        for key, value in defaults.items()
    }


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "inspect":
        if not 0 < args.retention_target <= 1:
            parser.error("--retention-target must be within (0, 1].")

        result = inspect_dataset(
            vcf_path=args.vcf,
            popmap_path=args.popmap,
            coordinates_path=args.coords,
            output_dir=args.output,
            preset=args.preset,
            retention_target=args.retention_target,
        )

        print("DIFLOW dataset inspection")
        print(f"Samples: {result.samples}")
        print(f"Populations: {result.populations}")
        print(f"Variant loci: {result.loci}")
        print(
            "Recommended projection: "
            f"{result.recommended_projection} chromosomes"
        )
        print(
            "Recommended starting graph: "
            f"{result.recommended_neighbors}-nearest neighbors"
        )
        print(f"Recommended effort preset: {result.recommended_preset}")
        print("")
        print("Projection retention options:")
        print(result.projection_summary.to_string(index=False))
        if not result.graph_summary.empty:
            print("")
            print("Geographic graph options:")
            print(result.graph_summary.to_string(index=False))
        print("")
        print("Recommended command:")
        print(result.command)
        print("")
        print(
            "Graph recommendations are computational starting points, not "
            "biological truths; review them before final inference."
        )
        return 0


    if args.command == "benchmark":
        if args.replicates < 1:
            parser.error("--replicates must be at least 1.")
        if args.chromosomes < 2:
            parser.error("--chromosomes must be at least 2.")
        if args.sites < 1:
            parser.error("--sites must be at least 1.")
        if args.starts < 1:
            parser.error("--starts must be at least 1.")
        if args.linked_blocks < 2:
            parser.error("--linked-blocks must be at least 2.")
        if args.linked_bootstrap_replicates < 2:
            parser.error("--linked-bootstrap-replicates must be at least 2.")
        if args.decision_bootstrap_replicates < 2:
            parser.error("--decision-bootstrap-replicates must be at least 2.")
        if args.snps_per_block < 1:
            parser.error("--snps-per-block must be at least 1.")
        if args.linked_block_bp < args.snps_per_block:
            parser.error("--linked-block-bp must be at least --snps-per-block.")
        if args.linkage_concentration <= 0:
            parser.error("--linkage-concentration must be positive.")
        if args.mechanistic_nref < 2:
            parser.error("--mechanistic-nref must be at least 2.")
        if args.mechanistic_sequence_length < 1000:
            parser.error("--mechanistic-sequence-length must be at least 1000.")
        if args.mechanistic_recombination_rate < 0:
            parser.error("--mechanistic-recombination-rate must be non-negative.")
        if args.mechanistic_mutation_rate <= 0:
            parser.error("--mechanistic-mutation-rate must be positive.")
        try:
            mechanistic_block_sizes = tuple(
                int(value.strip())
                for value in args.mechanistic_block_sizes.split(",")
                if value.strip()
            )
        except ValueError:
            parser.error("--mechanistic-block-sizes must be comma-separated integers.")
        if not mechanistic_block_sizes or any(value < 1 for value in mechanistic_block_sizes):
            parser.error("--mechanistic-block-sizes must contain positive integers.")

        common = dict(
            replicates=args.replicates,
            sample_sizes=(args.chromosomes, args.chromosomes),
            segregating_sites=args.sites,
            starts=args.starts,
            maxiter=args.maxiter,
            seed=args.seed,
        )

        if args.suite in {"recovery", "all"}:
            recovery_dir = (
                args.output if args.suite == "recovery" else str(Path(args.output) / "recovery")
            )
            raw, summary = write_recovery_benchmark(
                output_dir=recovery_dir,
                **common,
            )
            print("DIFLOW known-truth recovery benchmark")
            print(f"Replicate rows: {len(raw)}")
            print(f"Scenarios: {len(summary)}")
            print(summary.to_string(index=False))
            print("")

        if args.suite in {"stress", "all"}:
            stress_dir = (
                args.output if args.suite == "stress" else str(Path(args.output) / "stress")
            )
            raw, summary = write_stress_benchmark(
                output_dir=stress_dir,
                **common,
            )
            print("DIFLOW demographic stress benchmark")
            print(f"Replicate rows: {len(raw)}")
            print(f"Scenarios: {len(summary)}")
            print(summary.to_string(index=False))
            print("")

        if args.suite in {"forward", "all"}:
            forward_dir = (
                args.output if args.suite == "forward" else str(Path(args.output) / "forward")
            )
            raw, summary = write_forward_stress_benchmark(
                output_dir=forward_dir,
                replicates=args.replicates,
                loci=args.sites,
                sample_sizes=(args.chromosomes, args.chromosomes),
                starts=args.starts,
                maxiter=args.maxiter,
                seed=args.seed,
            )
            print("DIFLOW independent forward-time stress benchmark")
            print(f"Replicate rows: {len(raw)}")
            print(f"Scenarios: {len(summary)}")
            print(summary.to_string(index=False))
            print("")

        if args.suite == "grid":
            grid_dir = (
                args.output if args.suite == "grid" else str(Path(args.output) / "grid")
            )
            raw, summary, figures = write_recovery_grid(
                output_dir=grid_dir,
                replicates=args.replicates,
                starts=args.starts,
                maxiter=args.maxiter,
                seed=args.seed,
            )
            print("DIFLOW large recovery-grid benchmark")
            print(f"Replicate rows: {len(raw)}")
            print(f"Grid cells: {len(summary)}")
            print(f"Validation figures: {len(figures)}")
            print("")

        if args.suite == "linked":
            raw, summary, figures = write_linked_bootstrap_calibration(
                output_dir=args.output,
                replicates=args.replicates,
                sample_sizes=(args.chromosomes, args.chromosomes),
                blocks=args.linked_blocks,
                snps_per_block=args.snps_per_block,
                block_size_bp=args.linked_block_bp,
                concentration=args.linkage_concentration,
                bootstrap_replicates=args.linked_bootstrap_replicates,
                bootstrap_starts=max(1, min(args.starts, 5)),
                maxiter=args.maxiter,
                seed=args.seed,
            )
            print("DIFLOW linked-marker bootstrap calibration")
            print(f"Replicate rows: {len(raw)}")
            print(f"Summary rows: {len(summary)}")
            print(f"Validation figures: {len(figures)}")
            print("")

        if args.suite == "mechanistic":
            raw, summary, figures = write_mechanistic_linkage_calibration(
                output_dir=args.output,
                replicates=args.replicates,
                chromosomes_per_population=args.chromosomes,
                nref=args.mechanistic_nref,
                sequence_length=args.mechanistic_sequence_length,
                recombination_rate=args.mechanistic_recombination_rate,
                mutation_rate=args.mechanistic_mutation_rate,
                block_sizes_bp=mechanistic_block_sizes,
                bootstrap_replicates=args.linked_bootstrap_replicates,
                bootstrap_starts=max(1, min(args.starts, 5)),
                maxiter=args.maxiter,
                seed=args.seed,
            )
            print("DIFLOW mechanistic linkage calibration")
            print(f"Replicate rows: {len(raw)}")
            print(f"Summary rows: {len(summary)}")
            print(f"Validation figures: {len(figures)}")
            print("")

        if args.suite == "decision":
            evidence = write_decision_evidence_benchmark(
                output_dir=args.output,
                replicates=args.replicates,
                chromosomes=args.chromosomes,
                segregating_sites=args.sites,
                starts=args.starts,
                bootstrap_replicates=args.decision_bootstrap_replicates,
                bootstrap_starts=max(1, min(args.starts, 5)),
                maxiter=args.maxiter,
                seed=args.seed,
            )
            print("DIFLOW decision-evidence benchmark")
            print(f"Evidence rows: {len(evidence)}")
            print(f"Successful rows: {int(evidence['success'].astype(bool).sum())}")
            print("")

        print(f"Results: {args.output}")
        print(
            "Recovery tests model-consistent identifiability; stress and forward "
            "suites challenge misspecification; grid evaluates scaling across data sizes."
        )
        return 0

    if args.command == "calibrate":
        if not 0 <= args.max_fpr <= 1:
            parser.error("--max-fpr must lie within [0, 1].")
        selected, scan = write_threshold_calibration(
            evidence_csv=args.evidence,
            output_dir=args.output,
            max_false_directional_positive_rate=args.max_fpr,
        )
        print("DIFLOW directional threshold calibration")
        print(f"Threshold combinations evaluated: {len(scan)}")
        print(f"Maximum false-direction rate: {args.max_fpr:.3f}")
        print(f"Selected minimum model weight: {selected.min_model_weight:.3f}")
        print(
            "Selected minimum directional support: "
            f"{selected.min_directional_support:.3f}"
        )
        print(
            "Selected minimum absolute asymmetry: "
            f"{selected.min_abs_asymmetry:.3f}"
        )
        print(
            "Observed false-direction rate: "
            f"{selected.false_directional_positive_rate:.3f}"
        )
        print(
            "Directional sensitivity: "
            f"{selected.directional_sensitivity:.3f}"
        )
        print(f"Results: {args.output}")
        return 0

    if args.command == "infer":
        effort = _resolve_effort(args)
        result = run_infer_pipeline(
            vcf_path=args.vcf,
            popmap_path=args.popmap,
            coordinates_path=args.coords,
            output_dir=args.output,
            projection_chromosomes=args.projection_chromosomes,
            k_nearest=args.neighbors,
            max_distance_km=args.max_distance_km,
            starts=effort["starts"],
            maxiter=effort["maxiter"],
            bootstrap_replicates=effort["bootstrap_replicates"],
            bootstrap_starts=effort["bootstrap_starts"],
            bootstrap_block_bp=args.bootstrap_block_bp,
            min_model_weight=args.min_model_weight,
            min_directional_support=args.min_directional_support,
            min_abs_asymmetry=args.min_abs_asymmetry,
            polarized=args.polarized,
            prepare_only=args.prepare_only,
            seed=args.seed,
        )
        print(f"DIFLOW results: {result.output_dir}")
        print(f"Candidate pairs: {len(result.candidate_pairs)}")
        print(f"Pairwise rows: {len(result.pairwise_results)}")
        if args.prepare_only:
            print("Inference skipped (--prepare-only).")
        elif effort["bootstrap_replicates"] >= 2:
            if args.bootstrap_block_bp is None:
                print("Locus-bootstrap directional evidence classification completed.")
            else:
                print(
                    "Block-bootstrap directional evidence classification completed "
                    f"({args.bootstrap_block_bp} bp windows)."
                )
        else:
            print("Direction labels are provisional; bootstrap is disabled.")
        return 0

    parser.error("unknown command")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
