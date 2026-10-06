import json

from diflow.provenance import collect_provenance, sha256_file, write_provenance


def test_sha256_file_is_stable(tmp_path):
    path = tmp_path / "input.txt"
    path.write_text("DIFLOW\n", encoding="utf-8")
    assert sha256_file(path) == sha256_file(path)


def test_provenance_records_hashes_versions_and_settings(tmp_path):
    source = tmp_path / "input.txt"
    source.write_text("abc", encoding="utf-8")

    payload = collect_provenance(
        inputs={"example": source},
        settings={"seed": 42, "polarized": False},
    )
    assert payload["inputs"]["example"]["sha256"] == sha256_file(source)
    assert payload["settings"]["seed"] == 42
    assert "python" in payload["software"]

    output = tmp_path / "provenance.json"
    write_provenance(
        output,
        inputs={"example": source},
        settings={"seed": 42},
    )
    loaded = json.loads(output.read_text(encoding="utf-8"))
    assert loaded["settings"]["seed"] == 42
