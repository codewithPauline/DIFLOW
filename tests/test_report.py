import json

import pandas as pd

from diflow.report import write_html_report


def test_html_report_contains_summary_tables_and_provenance(tmp_path):
    results = tmp_path / "results"
    results.mkdir()

    pd.DataFrame(
        {
            "population_a": ["A", "A"],
            "population_b": ["B", "C"],
            "status": ["supported", "ambiguous"],
            "preferred_direction": ["A->B", "A->C"],
        }
    ).to_csv(results / "pairwise_results.csv", index=False)

    pd.DataFrame(
        {
            "population_a": ["A"],
            "population_b": ["B"],
            "model": ["asymmetric_migration"],
            "akaike_weight": [0.9],
        }
    ).to_csv(results / "model_rankings.csv", index=False)

    (results / "run_provenance.json").write_text(
        json.dumps({"settings": {"seed": 42}}),
        encoding="utf-8",
    )
    (results / "directional_map.png").write_bytes(b"fakepng")

    output = write_html_report(results)
    text = output.read_text(encoding="utf-8")

    assert output.name == "report.html"
    assert "DIFLOW results report" in text
    assert "Supported" in text
    assert "Pairwise inference" in text
    assert "Run provenance" in text
    assert "filterTable" in text
    assert "data:image/png;base64" in text
