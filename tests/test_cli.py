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



def test_infer_accepts_polarized_mode():
    parser = build_parser()
    args = parser.parse_args(
        [
            "infer",
            "--vcf", "data.vcf",
            "--popmap", "popmap.tsv",
            "--coords", "coords.csv",
            "--output", "results",
            "--projection-chromosomes", "8",
            "--polarized",
        ]
    )
    assert args.polarized is True



def test_profile_cli_parses():
    parser = build_parser()
    args = parser.parse_args(
        [
            "profile",
            "--spectrum", "pair.npy",
            "--parameter", "m_a_to_b",
            "--output", "profile_out",
        ]
    )
    assert args.command == "profile"
    assert args.parameter == "m_a_to_b"
    assert args.points == 15



def test_infer_accepts_calibrated_threshold_file():
    parser = build_parser()
    args = parser.parse_args(
        [
            "infer",
            "--vcf", "data.vcf",
            "--popmap", "popmap.tsv",
            "--coords", "coords.csv",
            "--output", "results",
            "--projection-chromosomes", "8",
            "--thresholds-file", "selected_thresholds.csv",
            "--min-model-weight", "0.85",
        ]
    )
    assert args.thresholds_file == "selected_thresholds.csv"
    assert args.min_model_weight == 0.85
    assert args.min_directional_support is None
    assert args.min_abs_asymmetry is None



def test_benchmark_export_command_parses():
    parser = build_parser()
    args = parser.parse_args(
        [
            "benchmark-export",
            "--input", "recovery.csv",
            "--output", "diflow_standardized.csv",
            "--estimand", "dadi_scaled_migration",
        ]
    )
    assert args.command == "benchmark-export"
    assert args.input == "recovery.csv"
    assert args.output == "diflow_standardized.csv"



def test_data_requirements_command_parses():
    parser = build_parser()
    args = parser.parse_args(
        [
            "data-requirements",
            "--grid-summary", "recovery_grid_summary.csv",
            "--output", "requirements",
            "--min-direction-accuracy", "0.92",
            "--max-false-direction-rate", "0.04",
            "--min-success-rate", "0.96",
        ]
    )
    assert args.command == "data-requirements"
    assert args.grid_summary == "recovery_grid_summary.csv"
    assert args.min_direction_accuracy == 0.92
    assert args.max_false_direction_rate == 0.04
    assert args.min_success_rate == 0.96



def test_archive_validation_command_parses():
    parser = build_parser()
    args = parser.parse_args(
        [
            "archive-validation",
            "--results", "validation_campaign/results",
            "--output", "archive/manifest.json",
        ]
    )
    assert args.command == "archive-validation"
    assert args.results == "validation_campaign/results"
    assert args.output == "archive/manifest.json"
