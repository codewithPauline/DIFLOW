"""Reproducibility metadata and input provenance for DIFLOW runs."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
from importlib import metadata
import json
from pathlib import Path
import platform
import sys
from typing import Any


def _version(distribution: str) -> str | None:
    try:
        return metadata.version(distribution)
    except metadata.PackageNotFoundError:
        return None


def sha256_file(path: str | Path, *, chunk_size: int = 1024 * 1024) -> str:
    """Return the SHA-256 digest of a file without loading it into memory."""
    source = Path(path)
    digest = hashlib.sha256()
    with source.open("rb") as handle:
        while True:
            chunk = handle.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def collect_provenance(
    *,
    inputs: dict[str, str | Path],
    settings: dict[str, Any],
) -> dict[str, Any]:
    """Collect software versions, input hashes, and resolved run settings."""
    input_records = {}
    for name, raw_path in inputs.items():
        path = Path(raw_path)
        input_records[name] = {
            "path": str(path),
            "sha256": sha256_file(path),
            "size_bytes": int(path.stat().st_size),
        }

    return {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "software": {
            "diflow": _version("diflow"),
            "python": sys.version.split()[0],
            "numpy": _version("numpy"),
            "pandas": _version("pandas"),
            "scipy": _version("scipy"),
            "matplotlib": _version("matplotlib"),
            "networkx": _version("networkx"),
            "pyproj": _version("pyproj"),
            "dadi": _version("dadi"),
            "msprime": _version("msprime"),
            "tskit": _version("tskit"),
        },
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
        },
        "inputs": input_records,
        "settings": settings,
    }


def write_provenance(
    output_path: str | Path,
    *,
    inputs: dict[str, str | Path],
    settings: dict[str, Any],
) -> Path:
    """Write a deterministic, human-readable JSON provenance record."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = collect_provenance(inputs=inputs, settings=settings)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path
