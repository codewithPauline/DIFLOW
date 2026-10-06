"""Command-line interface for DIFLOW."""

from __future__ import annotations

import argparse

from diflow.inspection import PRESETS, inspect_dataset
from diflow.pipeline import run_infer_pipeline
from diflow.validation import write_recovery_benchmark


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

        raw, summary = write_recovery_benchmark(
            output_dir=args.output,
            replicates=args.replicates,
            sample_sizes=(args.chromosomes, args.chromosomes),
            segregating_sites=args.sites,
            starts=args.starts,
            maxiter=args.maxiter,
            seed=args.seed,
        )
        print("DIFLOW known-truth recovery benchmark")
        print(f"Replicate rows: {len(raw)}")
        print(f"Scenarios: {len(summary)}")
        print(f"Results: {args.output}")
        print("")
        print(summary.to_string(index=False))
        print("")
        print(
            "This is a model-consistent recovery benchmark. Passing it is necessary "
            "but does not establish robustness to demographic misspecification."
        )
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
            prepare_only=args.prepare_only,
            seed=args.seed,
        )
        print(f"DIFLOW results: {result.output_dir}")
        print(f"Candidate pairs: {len(result.candidate_pairs)}")
        print(f"Pairwise rows: {len(result.pairwise_results)}")
        if args.prepare_only:
            print("Inference skipped (--prepare-only).")
        elif effort["bootstrap_replicates"] >= 2:
            print("Bootstrap directional evidence classification completed.")
        else:
            print("Direction labels are provisional; bootstrap is disabled.")
        return 0

    parser.error("unknown command")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
