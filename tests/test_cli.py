from diflow.cli import build_parser


def test_infer_cli_parses_required_arguments():
    parser = build_parser()
    args = parser.parse_args(
        [
            "infer",
            "--vcf",
            "data.vcf",
            "--popmap",
            "popmap.tsv",
            "--coords",
            "coords.csv",
            "--output",
            "results",
            "--projection-chromosomes",
            "8",
            "--neighbors",
            "4",
            "--starts",
            "20",
        ]
    )

    assert args.command == "infer"
    assert args.projection_chromosomes == 8
    assert args.neighbors == 4
    assert args.starts == 20
