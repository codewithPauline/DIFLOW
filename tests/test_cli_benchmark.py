from diflow.cli import build_parser


def test_benchmark_defaults_to_all_suites():
    parser = build_parser()
    args = parser.parse_args(
        [
            "benchmark",
            "--output", "bench",
        ]
    )
    assert args.suite == "all"
    assert args.replicates == 10
    assert args.chromosomes == 20
    assert args.sites == 5000


def test_benchmark_accepts_stress_suite():
    parser = build_parser()
    args = parser.parse_args(
        [
            "benchmark",
            "--output", "bench",
            "--suite", "stress",
            "--replicates", "3",
        ]
    )
    assert args.suite == "stress"
    assert args.replicates == 3
