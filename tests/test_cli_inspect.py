from diflow.cli import _resolve_effort, build_parser


def test_inspect_cli_parses_minimal_inputs():
    parser = build_parser()
    args = parser.parse_args(
        [
            "inspect",
            "--vcf", "data.vcf",
            "--popmap", "popmap.tsv",
            "--coords", "coords.csv",
        ]
    )
    assert args.command == "inspect"
    assert args.preset == "standard"
    assert args.retention_target == 0.80


def test_infer_preset_can_be_overridden():
    parser = build_parser()
    args = parser.parse_args(
        [
            "infer",
            "--vcf", "data.vcf",
            "--popmap", "popmap.tsv",
            "--coords", "coords.csv",
            "--output", "results",
            "--projection-chromosomes", "8",
            "--preset", "publication",
            "--starts", "25",
        ]
    )
    effort = _resolve_effort(args)
    assert effort["starts"] == 25
    assert effort["bootstrap_replicates"] == 500
    assert effort["bootstrap_starts"] == 8
    assert effort["maxiter"] == 200


def test_infer_defaults_are_reproducible():
    parser = build_parser()
    args = parser.parse_args(
        [
            "infer",
            "--vcf", "data.vcf",
            "--popmap", "popmap.tsv",
            "--coords", "coords.csv",
            "--output", "results",
            "--projection-chromosomes", "8",
        ]
    )
    assert args.seed == 42
