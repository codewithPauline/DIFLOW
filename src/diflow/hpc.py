"""HPC/Slurm job-script generation for reproducible DIFLOW runs."""

from __future__ import annotations

from pathlib import Path
import re


def _safe_job_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    if not cleaned:
        raise ValueError("job_name must contain at least one valid character.")
    return cleaned[:64]


def write_slurm_script(
    output_path: str | Path,
    *,
    command: str,
    job_name: str = "DIFLOW",
    cpus: int = 8,
    mem_gb: int = 32,
    hours: int = 24,
    email: str | None = None,
) -> Path:
    """Write a conservative Slurm script for an existing DIFLOW CLI command."""
    if not command.strip():
        raise ValueError("command cannot be empty.")
    if cpus < 1:
        raise ValueError("cpus must be at least 1.")
    if mem_gb < 1:
        raise ValueError("mem_gb must be at least 1.")
    if hours < 1:
        raise ValueError("hours must be at least 1.")

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    name = _safe_job_name(job_name)

    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name={name}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        f"#SBATCH --cpus-per-task={cpus}",
        f"#SBATCH --mem={mem_gb}G",
        f"#SBATCH --time={hours}:00:00",
        f"#SBATCH --output={name}_%j.out",
        f"#SBATCH --error={name}_%j.err",
    ]
    if email:
        lines.extend(
            [
                f"#SBATCH --mail-user={email}",
                "#SBATCH --mail-type=END,FAIL",
            ]
        )

    lines.extend(
        [
            "",
            "set -euo pipefail",
            "",
            'echo "DIFLOW job started: $(date)"',
            'echo "Host: $(hostname)"',
            "",
            command.strip(),
            "",
            'echo "DIFLOW job finished: $(date)"',
            "",
        ]
    )

    path.write_text("\n".join(lines), encoding="utf-8")
    return path
