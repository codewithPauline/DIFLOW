"""Simulation validation and benchmark utilities for DIFLOW."""

from .metrics import (
    ValidationMetrics,
    direction_accuracy,
    false_directional_positive_rate,
    interval_coverage,
    parameter_bias,
    parameter_rmse,
)
from .scenarios import ValidationScenario, core_validation_scenarios

__all__ = [
    "ValidationMetrics",
    "direction_accuracy",
    "false_directional_positive_rate",
    "interval_coverage",
    "parameter_bias",
    "parameter_rmse",
    "ValidationScenario",
    "core_validation_scenarios",
]


from .recovery import (
    RecoveryScenario,
    default_recovery_scenarios,
    run_recovery_benchmark,
    simulate_model_consistent_spectrum,
    summarize_recovery,
    write_recovery_benchmark,
)

__all__ += [
    "RecoveryScenario",
    "default_recovery_scenarios",
    "run_recovery_benchmark",
    "simulate_model_consistent_spectrum",
    "summarize_recovery",
    "write_recovery_benchmark",
]


from .stress import (
    StressScenario,
    default_stress_scenarios,
    run_stress_benchmark,
    summarize_stress,
    write_stress_benchmark,
)

__all__ += [
    "StressScenario",
    "default_stress_scenarios",
    "run_stress_benchmark",
    "summarize_stress",
    "write_stress_benchmark",
]


from .forward_stress import (
    ForwardStressScenario,
    default_forward_stress_scenarios,
    run_forward_stress_benchmark,
    simulate_forward_stress_spectrum,
    summarize_forward_stress,
    write_forward_stress_benchmark,
)

__all__ += [
    "ForwardStressScenario",
    "default_forward_stress_scenarios",
    "run_forward_stress_benchmark",
    "simulate_forward_stress_spectrum",
    "summarize_forward_stress",
    "write_forward_stress_benchmark",
]


from .grid import (
    RecoveryGridConfig,
    default_recovery_grid,
    plot_recovery_grid,
    run_recovery_grid,
    write_recovery_grid,
)

__all__ += [
    "RecoveryGridConfig",
    "default_recovery_grid",
    "plot_recovery_grid",
    "run_recovery_grid",
    "write_recovery_grid",
]


from .linked_calibration import (
    LinkedCalibrationScenario,
    default_linked_calibration_scenarios,
    plot_linked_bootstrap_calibration,
    run_linked_bootstrap_calibration,
    simulate_correlated_block_counts,
    summarize_linked_bootstrap_calibration,
    write_linked_bootstrap_calibration,
)

__all__ += [
    "LinkedCalibrationScenario",
    "default_linked_calibration_scenarios",
    "plot_linked_bootstrap_calibration",
    "run_linked_bootstrap_calibration",
    "simulate_correlated_block_counts",
    "summarize_linked_bootstrap_calibration",
    "write_linked_bootstrap_calibration",
]


from .mechanistic_linkage import (
    MechanisticLinkageScenario,
    MechanisticGridConfig,
    default_mechanistic_linkage_scenarios,
    plot_mechanistic_linkage,
    run_mechanistic_linkage_calibration,
    simulate_msprime_counts,
    summarize_mechanistic_linkage,
    write_mechanistic_linkage_calibration,
    default_mechanistic_grid,
    run_mechanistic_linkage_grid,
    plot_mechanistic_linkage_grid,
    write_mechanistic_linkage_grid,
)

__all__ += [
    "MechanisticLinkageScenario",
    "MechanisticGridConfig",
    "default_mechanistic_linkage_scenarios",
    "plot_mechanistic_linkage",
    "run_mechanistic_linkage_calibration",
    "simulate_msprime_counts",
    "summarize_mechanistic_linkage",
    "write_mechanistic_linkage_calibration",
    "default_mechanistic_grid",
    "run_mechanistic_linkage_grid",
    "plot_mechanistic_linkage_grid",
    "write_mechanistic_linkage_grid",
]


from .thresholds import (
    CalibratedThresholds,
    load_calibrated_thresholds,
    calibrate_thresholds,
    evaluate_thresholds,
    scan_thresholds,
    select_thresholds,
    write_threshold_calibration,
)
from .decision_evidence import (
    DecisionEvidenceScenario,
    default_decision_evidence_scenarios,
    run_decision_evidence_benchmark,
    simulate_counts_from_expected_spectrum,
    write_decision_evidence_benchmark,
)

__all__ += [
    "CalibratedThresholds",
    "load_calibrated_thresholds",
    "calibrate_thresholds",
    "evaluate_thresholds",
    "scan_thresholds",
    "select_thresholds",
    "write_threshold_calibration",
    "DecisionEvidenceScenario",
    "default_decision_evidence_scenarios",
    "run_decision_evidence_benchmark",
    "simulate_counts_from_expected_spectrum",
    "write_decision_evidence_benchmark",
]


from .comparators import (
    compare_method_files,
    summarize_direction_only_comparison,
    standardize_diflow_benchmark,
    write_standardized_diflow_benchmark,
    validate_matched_comparison,
    summarize_method_comparison,
    validate_comparator_table,
)

__all__ += [
    "compare_method_files",
    "summarize_direction_only_comparison",
    "standardize_diflow_benchmark",
    "write_standardized_diflow_benchmark",
    "validate_matched_comparison",
    "summarize_method_comparison",
    "validate_comparator_table",
]


from .release_review import (
    ReleaseCriteria,
    review_external_comparison,
    review_linkage_coverage,
    review_recovery_grid,
    review_threshold_calibration,
    review_validation_campaign,
)

__all__ += [
    "ReleaseCriteria",
    "review_external_comparison",
    "review_linkage_coverage",
    "review_recovery_grid",
    "review_threshold_calibration",
    "review_validation_campaign",
]


from .data_requirements import (
    DataRequirementTargets,
    pareto_minimum_regimes,
    recommend_data_requirements,
    summarize_data_regimes,
    write_data_requirements,
)

__all__ += [
    "DataRequirementTargets",
    "pareto_minimum_regimes",
    "recommend_data_requirements",
    "summarize_data_regimes",
    "write_data_requirements",
]
