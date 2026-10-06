"""Generate reproducible HPC campaigns for DIFLOW validation studies."""

from __future__ import annotations

import json
from pathlib import Path

from .hpc import write_slurm_script


def write_validation_campaign(
    output_dir: str | Path,
    *,
    diflow_executable: str = "diflow",
    replicates: int = 50,
    starts: int = 10,
    bootstrap_replicates: int = 100,
    cpus: int = 8,
    mem_gb: int = 64,
    hours: int = 72,
    seed: int = 42,
    email: str | None = None,
) -> Path:
    """Write Slurm scripts and a manifest for release-grade validation work.

    The campaign intentionally separates the expensive studies into distinct
    jobs so failures can be rerun independently and resource use can be tuned
    per study.
    """
    if replicates < 1:
        raise ValueError("replicates must be at least 1.")
    if starts < 1:
        raise ValueError("starts must be at least 1.")
    if bootstrap_replicates < 2:
        raise ValueError("bootstrap_replicates must be at least 2.")

    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    results = root / "results"
    results.mkdir(exist_ok=True)

    jobs = {
        "recovery_grid": (
            f"{diflow_executable} benchmark --suite grid "
            f"--output {results / 'recovery_grid'} "
            f"--replicates {replicates} --starts {starts} --seed {seed}"
        ),
        "decision_evidence": (
            f"{diflow_executable} benchmark --suite decision "
            f"--output {results / 'decision_evidence'} "
            f"--replicates {replicates} --starts {starts} "
            f"--decision-bootstrap-replicates {bootstrap_replicates} "
            f"--seed {seed + 1000}"
        ),
        "linked_calibration": (
            f"{diflow_executable} benchmark --suite linked "
            f"--output {results / 'linked_calibration'} "
            f"--replicates {replicates} --starts {min(starts, 5)} "
            f"--linked-bootstrap-replicates {bootstrap_replicates} "
            f"--seed {seed + 2000}"
        ),
        "mechanistic_linkage": (
            f"{diflow_executable} benchmark --suite mechanistic "
            f"--output {results / 'mechanistic_linkage'} "
            f"--replicates {replicates} --starts {min(starts, 5)} "
            f"--linked-bootstrap-replicates {bootstrap_replicates} "
            f"--seed {seed + 3000}"
        ),
    }

    manifest = {
        "campaign_version": 1,
        "purpose": (
            "Release-grade DIFLOW validation campaign. Completion of scripts "
            "does not itself imply that acceptance criteria were met."
        ),
        "settings": {
            "replicates": replicates,
            "starts": starts,
            "bootstrap_replicates": bootstrap_replicates,
            "cpus": cpus,
            "mem_gb": mem_gb,
            "hours": hours,
            "seed": seed,
        },
        "jobs": {},
        "postprocessing": {
            "threshold_calibration": (
                f"{diflow_executable} calibrate "
                f"--evidence {results / 'decision_evidence' / 'decision_evidence.csv'} "
                f"--output {results / 'threshold_calibration'} --max-fpr 0.05"
            ),
            "release_review": [
                "Review recovery-grid direction accuracy and bias.",
                "Confirm symmetric false-direction rates meet the target.",
                "Review threshold sensitivity versus false-positive tradeoff.",
                "Review mechanistic linkage CI coverage across block sizes.",
                "Run and normalize external-method comparator studies.",
                "Document failure regimes before freezing release defaults.",
            ],
        },
    }

    for index, (name, command) in enumerate(jobs.items(), start=1):
        script = root / f"{index:02d}_{name}.slurm"
        write_slurm_script(
            script,
            command=command,
            job_name=f"DIFLOW_{name}",
            cpus=cpus,
            mem_gb=mem_gb,
            hours=hours,
            email=email,
        )
        manifest["jobs"][name] = {
            "script": str(script),
            "command": command,
            "results": str(results / name),
        }

    manifest_path = root / "campaign_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    submit_lines = ["#!/bin/bash", "set -euo pipefail", ""]
    for index, name in enumerate(jobs, start=1):
        submit_lines.append(f"sbatch {root / f'{index:02d}_{name}.slurm'}")
    submit_lines.append("")
    submit_path = root / "submit_all.sh"
    submit_path.write_text("\n".join(submit_lines), encoding="utf-8")

    return manifest_path
