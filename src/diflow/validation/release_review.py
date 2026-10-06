"""Release-readiness review for completed DIFLOW validation campaigns."""

from __future__ import annotations

from dataclasses import dataclass, asdict
import json
from pathlib import Path

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class ReleaseCriteria:
    """Project policy targets used to review a validation campaign.

    These are release criteria, not universal population-genetic constants.
    They remain configurable so future releases can tighten them transparently.
    """

    min_direction_accuracy: float = 0.90
    max_false_direction_rate: float = 0.05
    min_directional_sensitivity: float = 0.80
    min_ci_coverage: float = 0.90
    max_ci_coverage: float = 0.99
    min_success_rate: float = 0.95


def _validate_probability(value: float, name: str) -> float:
    value = float(value)
    if not np.isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f"{name} must lie within [0, 1].")
    return value


def _criteria_checked(criteria: ReleaseCriteria) -> ReleaseCriteria:
    values = asdict(criteria)
    for name, value in values.items():
        _validate_probability(value, name)
    if criteria.min_ci_coverage > criteria.max_ci_coverage:
        raise ValueError("min_ci_coverage cannot exceed max_ci_coverage.")
    return criteria


def _metric_check(
    *,
    section: str,
    metric: str,
    observed: float | None,
    target: str,
    passed: bool | None,
    detail: str,
) -> dict:
    return {
        "section": section,
        "metric": metric,
        "observed": observed,
        "target": target,
        "passed": passed,
        "detail": detail,
    }


def review_recovery_grid(
    summary: pd.DataFrame,
    criteria: ReleaseCriteria,
) -> list[dict]:
    """Review recovery-grid accuracy, false-direction rate, and fit success."""
    required = {
        "scenario",
        "success_rate",
        "direction_accuracy",
        "false_directional_positive_rate",
    }
    if not required.issubset(summary.columns):
        missing = sorted(required - set(summary.columns))
        raise ValueError("recovery grid summary missing: " + ", ".join(missing))

    checks: list[dict] = []

    success = pd.to_numeric(summary["success_rate"], errors="coerce").dropna()
    observed_success = float(success.min()) if not success.empty else None
    checks.append(
        _metric_check(
            section="recovery_grid",
            metric="minimum_success_rate",
            observed=observed_success,
            target=f">= {criteria.min_success_rate:.3f}",
            passed=(
                observed_success >= criteria.min_success_rate
                if observed_success is not None
                else None
            ),
            detail="Worst recovery-grid optimization success rate across cells.",
        )
    )

    directional = summary[
        summary["scenario"].astype(str) != "symmetric"
    ].copy()
    accuracy = pd.to_numeric(
        directional["direction_accuracy"], errors="coerce"
    ).dropna()
    observed_accuracy = float(accuracy.min()) if not accuracy.empty else None
    checks.append(
        _metric_check(
            section="recovery_grid",
            metric="minimum_direction_accuracy",
            observed=observed_accuracy,
            target=f">= {criteria.min_direction_accuracy:.3f}",
            passed=(
                observed_accuracy >= criteria.min_direction_accuracy
                if observed_accuracy is not None
                else None
            ),
            detail="Worst directional recovery accuracy across non-symmetric grid cells.",
        )
    )

    symmetric = summary[
        summary["false_directional_positive_rate"].notna()
    ].copy()
    fpr = pd.to_numeric(
        symmetric["false_directional_positive_rate"], errors="coerce"
    ).dropna()
    observed_fpr = float(fpr.max()) if not fpr.empty else None
    checks.append(
        _metric_check(
            section="recovery_grid",
            metric="maximum_false_direction_rate",
            observed=observed_fpr,
            target=f"<= {criteria.max_false_direction_rate:.3f}",
            passed=(
                observed_fpr <= criteria.max_false_direction_rate
                if observed_fpr is not None
                else None
            ),
            detail="Worst false directional-positive rate under symmetric truth.",
        )
    )
    return checks


def review_threshold_calibration(
    selected: pd.DataFrame,
    criteria: ReleaseCriteria,
) -> list[dict]:
    """Review the selected empirical decision thresholds."""
    required = {
        "false_directional_positive_rate",
        "directional_sensitivity",
        "direction_accuracy_when_called",
    }
    if not required.issubset(selected.columns):
        missing = sorted(required - set(selected.columns))
        raise ValueError(
            "selected threshold table missing: " + ", ".join(missing)
        )
    if len(selected) != 1:
        raise ValueError("selected threshold table must contain exactly one row.")

    row = selected.iloc[0]
    fpr = float(row["false_directional_positive_rate"])
    sensitivity = float(row["directional_sensitivity"])
    accuracy = float(row["direction_accuracy_when_called"])

    return [
        _metric_check(
            section="threshold_calibration",
            metric="false_direction_rate",
            observed=fpr,
            target=f"<= {criteria.max_false_direction_rate:.3f}",
            passed=fpr <= criteria.max_false_direction_rate,
            detail="False-direction rate for the selected classifier thresholds.",
        ),
        _metric_check(
            section="threshold_calibration",
            metric="directional_sensitivity",
            observed=sensitivity,
            target=f">= {criteria.min_directional_sensitivity:.3f}",
            passed=sensitivity >= criteria.min_directional_sensitivity,
            detail="Fraction of directional truths correctly called.",
        ),
        _metric_check(
            section="threshold_calibration",
            metric="direction_accuracy_when_called",
            observed=accuracy,
            target=f">= {criteria.min_direction_accuracy:.3f}",
            passed=accuracy >= criteria.min_direction_accuracy,
            detail="Direction accuracy among directional calls.",
        ),
    ]


def review_linkage_coverage(
    summary: pd.DataFrame,
    criteria: ReleaseCriteria,
) -> list[dict]:
    """Review mechanistic linkage interval coverage and fit success."""
    required = {
        "method",
        "success_rate",
        "coverage_a_to_b",
        "coverage_b_to_a",
    }
    if not required.issubset(summary.columns):
        missing = sorted(required - set(summary.columns))
        raise ValueError(
            "mechanistic linkage summary missing: " + ", ".join(missing)
        )

    block = summary[
        summary["method"].astype(str).str.startswith("block_")
    ].copy()
    if block.empty:
        return [
            _metric_check(
                section="mechanistic_linkage",
                metric="block_results_present",
                observed=None,
                target="at least one block-bootstrap method",
                passed=False,
                detail="No block-bootstrap rows were found.",
            )
        ]

    success = pd.to_numeric(block["success_rate"], errors="coerce").dropna()
    min_success = float(success.min()) if not success.empty else None

    coverage_values = pd.concat(
        [
            pd.to_numeric(block["coverage_a_to_b"], errors="coerce"),
            pd.to_numeric(block["coverage_b_to_a"], errors="coerce"),
        ],
        ignore_index=True,
    ).dropna()
    min_coverage = (
        float(coverage_values.min()) if not coverage_values.empty else None
    )
    max_coverage = (
        float(coverage_values.max()) if not coverage_values.empty else None
    )

    return [
        _metric_check(
            section="mechanistic_linkage",
            metric="minimum_success_rate",
            observed=min_success,
            target=f">= {criteria.min_success_rate:.3f}",
            passed=(
                min_success >= criteria.min_success_rate
                if min_success is not None
                else None
            ),
            detail="Worst successful-fit rate across block-bootstrap linkage cells.",
        ),
        _metric_check(
            section="mechanistic_linkage",
            metric="minimum_interval_coverage",
            observed=min_coverage,
            target=f">= {criteria.min_ci_coverage:.3f}",
            passed=(
                min_coverage >= criteria.min_ci_coverage
                if min_coverage is not None
                else None
            ),
            detail="Worst empirical confidence-interval coverage across block methods.",
        ),
        _metric_check(
            section="mechanistic_linkage",
            metric="maximum_interval_coverage",
            observed=max_coverage,
            target=f"<= {criteria.max_ci_coverage:.3f}",
            passed=(
                max_coverage <= criteria.max_ci_coverage
                if max_coverage is not None
                else None
            ),
            detail=(
                "Upper coverage guardrail. Extreme overcoverage can indicate "
                "uninformatively wide uncertainty intervals."
            ),
        ),
    ]


def review_external_comparison(summary: pd.DataFrame) -> list[dict]:
    """Verify that an established-method benchmark is actually present."""
    required = {"method", "scenario", "replicates"}
    if not required.issubset(summary.columns):
        missing = sorted(required - set(summary.columns))
        raise ValueError(
            "method comparison summary missing: " + ", ".join(missing)
        )
    methods = sorted(set(summary["method"].astype(str)))
    external = [name for name in methods if name.lower() != "diflow"]
    return [
        _metric_check(
            section="external_comparison",
            metric="external_methods_present",
            observed=float(len(external)),
            target=">= 1 established external method",
            passed=len(external) >= 1,
            detail=(
                "External methods found: "
                + (", ".join(external) if external else "none")
            ),
        )
    ]




def review_data_requirements(path: str | Path) -> list[dict]:
    """Require empirical minimum-data guidance from the recovery campaign."""
    source = Path(path)
    data = json.loads(source.read_text(encoding="utf-8"))
    status = str(data.get("status", ""))
    passing = int(data.get("passing_regimes", 0))
    minima = data.get("pareto_minimum_regimes", [])
    passed = (
        status == "recommendations_available"
        and passing > 0
        and len(minima) > 0
    )
    return [
        _metric_check(
            section="data_requirements",
            metric="empirical_minimum_guidance",
            observed=float(passing),
            target=">= 1 tested passing regime with Pareto-minimum guidance",
            passed=passed,
            detail=(
                f"Status={status}; passing regimes={passing}; "
                f"Pareto minima={len(minima)}."
            ),
        )
    ]


def review_validation_campaign(
    results_dir: str | Path,
    *,
    output_dir: str | Path | None = None,
    criteria: ReleaseCriteria | None = None,
) -> tuple[pd.DataFrame, dict]:
    """Review completed validation outputs against explicit release criteria.

    Missing empirical studies are reported as blockers rather than silently
    ignored.
    """
    criteria = _criteria_checked(criteria or ReleaseCriteria())
    root = Path(results_dir)
    outdir = Path(output_dir) if output_dir is not None else root / "release_review"
    outdir.mkdir(parents=True, exist_ok=True)

    checks: list[dict] = []
    expected = {
        "recovery_grid": root / "recovery_grid" / "recovery_grid_summary.csv",
        "threshold_calibration": (
            root / "threshold_calibration" / "selected_thresholds.csv"
        ),
        "mechanistic_linkage": (
            root / "mechanistic_linkage" / "mechanistic_grid_summary.csv"
        ),
        "external_comparison": (
            root / "external_comparison" / "method_comparison_summary.csv"
        ),
        "data_requirements": (
            root / "data_requirements" / "data_requirements.json"
        ),
    }

    reviewers = {
        "recovery_grid": review_recovery_grid,
        "threshold_calibration": review_threshold_calibration,
        "mechanistic_linkage": review_linkage_coverage,
    }

    for section, path in expected.items():
        if section == "mechanistic_linkage" and not path.exists():
            legacy = root / "mechanistic_linkage" / "mechanistic_linkage_summary.csv"
            if legacy.exists():
                path = legacy
        if not path.exists():
            checks.append(
                _metric_check(
                    section=section,
                    metric="study_completed",
                    observed=None,
                    target=f"required file: {path.name}",
                    passed=False,
                    detail=f"Required empirical output not found: {path}",
                )
            )
            continue

        if section == "data_requirements":
            checks.extend(review_data_requirements(path))
            continue

        frame = pd.read_csv(path)
        if section == "external_comparison":
            checks.extend(review_external_comparison(frame))
        else:
            checks.extend(reviewers[section](frame, criteria))

    table = pd.DataFrame(checks)
    passed_values = table["passed"].fillna(False).astype(bool)
    overall_pass = bool(passed_values.all()) if len(table) else False
    blockers = table.loc[~passed_values, ["section", "metric", "detail"]]

    summary = {
        "release_ready": overall_pass,
        "criteria": asdict(criteria),
        "checks_total": int(len(table)),
        "checks_passed": int(passed_values.sum()),
        "checks_failed_or_missing": int((~passed_values).sum()),
        "blockers": blockers.to_dict(orient="records"),
    }

    table.to_csv(outdir / "release_review.csv", index=False)
    (outdir / "release_review.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    lines = [
        "# DIFLOW release review",
        "",
        f"**Release ready:** {'YES' if overall_pass else 'NO'}",
        "",
        f"Checks passed: {summary['checks_passed']} / {summary['checks_total']}",
        "",
        "## Acceptance criteria",
        "",
    ]
    for name, value in asdict(criteria).items():
        lines.append(f"- {name}: {value}")
    lines.extend(["", "## Checks", ""])
    for row in checks:
        state = "PASS" if row["passed"] is True else "FAIL"
        lines.append(
            f"- **{state}** — {row['section']} / {row['metric']}: "
            f"{row['detail']}"
        )
    (outdir / "release_review.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )
    return table, summary
