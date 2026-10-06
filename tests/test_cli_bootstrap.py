from diflow.cli import build_parser


def test_infer_cli_parses_bootstrap_arguments():
    parser = build_parser()
    args = parser.parse_args(
        [
            "infer",
            "--vcf", "data.vcf",
            "--popmap", "popmap.tsv",
            "--coords", "coords.csv",
            "--output", "results",
            "--projection-chromosomes", "8",
            "--neighbors", "4",
            "--starts", "20",
            "--bootstrap-replicates", "100",
            "--bootstrap-starts", "4",
        ]
    )

    assert args.command == "infer"
    assert args.bootstrap_replicates == 100
    assert args.bootstrap_starts == 4
