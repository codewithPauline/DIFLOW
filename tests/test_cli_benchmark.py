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



def test_benchmark_accepts_grid_suite():
    parser = build_parser()
    args = parser.parse_args(
        [
            "benchmark",
            "--output", "bench",
            "--suite", "grid",
            "--replicates", "2",
        ]
    )
    assert args.suite == "grid"
    assert args.replicates == 2



def test_benchmark_accepts_linked_suite():
    parser = build_parser()
    args = parser.parse_args(
        [
            "benchmark",
            "--output", "linked",
            "--suite", "linked",
            "--replicates", "5",
            "--linked-blocks", "20",
            "--linked-bootstrap-replicates", "50",
            "--snps-per-block", "8",
            "--linked-block-bp", "50000",
            "--linkage-concentration", "15",
        ]
    )
    assert args.suite == "linked"
    assert args.linked_blocks == 20
    assert args.linked_bootstrap_replicates == 50
    assert args.snps_per_block == 8
    assert args.linked_block_bp == 50000
    assert args.linkage_concentration == 15



def test_benchmark_accepts_mechanistic_grid_suite():
    parser = build_parser()
    args = parser.parse_args(
        [
            "benchmark",
            "--output", "mechanistic_grid",
            "--suite", "mechanistic-grid",
            "--replicates", "2",
        ]
    )
    assert args.suite == "mechanistic-grid"
