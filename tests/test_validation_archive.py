import hashlib
import json

from diflow.validation.archive import (
    build_validation_archive_manifest,
    write_validation_archive_manifest,
)


def test_archive_manifest_hashes_files_deterministically(tmp_path):
    root = tmp_path / "results"
    root.mkdir()
    (root / "a.txt").write_text("alpha\n", encoding="utf-8")
    sub = root / "nested"
    sub.mkdir()
    (sub / "b.csv").write_text("x,y\n1,2\n", encoding="utf-8")

    manifest = build_validation_archive_manifest(root)

    assert manifest["file_count"] == 2
    paths = [item["path"] for item in manifest["files"]]
    assert paths == ["a.txt", "nested/b.csv"]

    expected = hashlib.sha256(b"alpha\n").hexdigest()
    assert manifest["files"][0]["sha256"] == expected


def test_write_archive_manifest_creates_json_and_sha256(tmp_path):
    root = tmp_path / "results"
    root.mkdir()
    (root / "result.csv").write_text("a\n1\n", encoding="utf-8")

    path = write_validation_archive_manifest(root)

    assert path.exists()
    checksum = path.with_suffix(".sha256")
    assert checksum.exists()

    data = json.loads(path.read_text(encoding="utf-8"))
    assert data["file_count"] == 1
    assert "result.csv" in checksum.read_text(encoding="utf-8")


def test_archive_manifest_excludes_existing_manifest_files(tmp_path):
    root = tmp_path / "results"
    root.mkdir()
    (root / "result.txt").write_text("result", encoding="utf-8")
    (root / "validation_archive_manifest.json").write_text("old", encoding="utf-8")
    (root / "validation_archive_manifest.sha256").write_text("old", encoding="utf-8")

    manifest = build_validation_archive_manifest(root)

    assert manifest["file_count"] == 1
    assert manifest["files"][0]["path"] == "result.txt"
