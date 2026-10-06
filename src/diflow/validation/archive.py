"""Archive manifests for DIFLOW validation and release evidence."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_validation_archive_manifest(
    results_dir: str | Path,
    *,
    exclude_manifest: bool = True,
) -> dict:
    """Build a deterministic checksum inventory of validation outputs."""
    root = Path(results_dir).resolve()
    if not root.exists() or not root.is_dir():
        raise ValueError(f"validation results directory not found: {root}")

    files = []
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        relative = path.relative_to(root).as_posix()
        if exclude_manifest and relative in {
            "validation_archive_manifest.json",
            "validation_archive_manifest.sha256",
        }:
            continue
        files.append(
            {
                "path": relative,
                "bytes": int(path.stat().st_size),
                "sha256": _sha256(path),
            }
        )

    if not files:
        raise ValueError("validation results directory contains no files.")

    return {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "root_name": root.name,
        "file_count": len(files),
        "total_bytes": int(sum(item["bytes"] for item in files)),
        "files": files,
    }


def write_validation_archive_manifest(
    results_dir: str | Path,
    *,
    output_path: str | Path | None = None,
) -> Path:
    """Write JSON and conventional SHA-256 checksum manifests."""
    root = Path(results_dir).resolve()
    manifest = build_validation_archive_manifest(root)

    output = (
        Path(output_path)
        if output_path is not None
        else root / "validation_archive_manifest.json"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    checksum_path = output.with_suffix(".sha256")
    checksum_lines = [
        f"{item['sha256']}  {item['path']}"
        for item in manifest["files"]
    ]
    checksum_path.write_text(
        "\n".join(checksum_lines) + "\n",
        encoding="utf-8",
    )
    return output
