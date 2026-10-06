"""Command-line interface for DIFLOW."""

from __future__ import annotations

import argparse

from diflow.pipeline import run_infer_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="diflow",
        description="Directional inference and geographic visualization of gene flow.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

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
    infer.add_argument(
        "--max-distance-km",
        type=float,
        default=None,
        help="Maximum geographic distance for candidate pairs.",
    )
    infer.add_argument("--starts", type=int, default=10, help="Optimization starts/model.")
    infer.add_argument("--maxiter", type=int, default=100, help="Optimizer iterations/start.")
    infer.add_argument("--seed", type=int, default=None, help="Random seed.")
    infer.add_argument(
        "--prepare-only",
        action="store_true",
        help="Stop after candidate pairs and projected jSFS construction.",
    )

    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "infer":
        result = run_infer_pipeline(
            vcf_path=args.vcf,
            popmap_path=args.popmap,
            coordinates_path=args.coords,
            output_dir=args.output,
            projection_chromosomes=args.projection_chromosomes,
            k_nearest=args.neighbors,
            max_distance_km=args.max_distance_km,
            starts=args.starts,
            maxiter=args.maxiter,
            prepare_only=args.prepare_only,
            seed=args.seed,
        )
        print(f"DIFLOW results: {result.output_dir}")
        print(f"Candidate pairs: {len(result.candidate_pairs)}")
        print(f"Pairwise rows: {len(result.pairwise_results)}")
        if args.prepare_only:
            print("Inference skipped (--prepare-only).")
        else:
            print("Direction labels are provisional pending jSFS uncertainty calibration.")
        return 0

    parser.error("unknown command")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
