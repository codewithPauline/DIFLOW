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
