"""Self-contained interactive HTML reporting for DIFLOW result directories."""

from __future__ import annotations

from html import escape
import base64
import json
from pathlib import Path

import pandas as pd


_TABLE_FILES = (
    ("Pairwise inference", "pairwise_results.csv"),
    ("Model rankings", "model_rankings.csv"),
    ("Directional flows", "directional_flows.csv"),
    ("Network summary", "network_summary.csv"),
    ("Run metadata", "run_metadata.csv"),
)


def _load_csv(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _status_counts(pairwise: pd.DataFrame | None) -> dict[str, int]:
    counts = {"supported": 0, "ambiguous": 0, "unsupported": 0}
    if pairwise is None or pairwise.empty or "status" not in pairwise.columns:
        return counts
    values = pairwise["status"].astype(str).value_counts()
    for key in counts:
        counts[key] = int(values.get(key, 0))
    return counts


def write_html_report(
    results_dir: str | Path,
    *,
    output_path: str | Path | None = None,
) -> Path:
    """Create a self-contained searchable HTML summary of a DIFLOW run."""
    root = Path(results_dir)
    if not root.exists():
        raise FileNotFoundError(root)
    if not root.is_dir():
        raise ValueError("results_dir must be a directory.")

    output = Path(output_path) if output_path else root / "report.html"
    output.parent.mkdir(parents=True, exist_ok=True)

    tables: list[tuple[str, str, pd.DataFrame]] = []
    pairwise = None
    for title, filename in _TABLE_FILES:
        frame = _load_csv(root / filename)
        if frame is not None:
            tables.append((title, filename, frame))
            if filename == "pairwise_results.csv":
                pairwise = frame

    counts = _status_counts(pairwise)
    pair_count = 0 if pairwise is None else len(pairwise)

    provenance_text = None
    provenance_path = root / "run_provenance.json"
    if provenance_path.exists():
        try:
            provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
            provenance_text = json.dumps(provenance, indent=2, sort_keys=True)
        except (json.JSONDecodeError, OSError):
            provenance_text = provenance_path.read_text(
                encoding="utf-8",
                errors="replace",
            )

    map_path = root / "directional_map.png"
    map_html = ""
    if map_path.exists():
        encoded = base64.b64encode(map_path.read_bytes()).decode("ascii")
        map_html = (
            '<section><h2>Directional map</h2>'
            f'<img class="map" src="data:image/png;base64,{encoded}" '
            'alt="DIFLOW directional migration map"></section>'
        )

    table_sections = []
    for index, (title, filename, frame) in enumerate(tables):
        table_id = f"table_{index}"
        search_id = f"search_{index}"
        html_table = frame.to_html(
            index=False,
            escape=True,
            table_id=table_id,
            classes=["dataframe", "diflow-table"],
            border=0,
        )
        table_sections.append(
            f"""
<section>
  <h2>{escape(title)}</h2>
  <p class="file">{escape(filename)} · {len(frame)} rows</p>
  <input id="{search_id}" class="search" type="search"
         placeholder="Filter rows..."
         oninput="filterTable('{table_id}', this.value)">
  <div class="table-wrap">{html_table}</div>
</section>
"""
        )

    provenance_html = ""
    if provenance_text is not None:
        provenance_html = (
            "<section><h2>Run provenance</h2>"
            "<details><summary>Show provenance JSON</summary>"
            f"<pre>{escape(provenance_text)}</pre></details></section>"
        )

    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>DIFLOW results report</title>
<style>
:root {{ color-scheme: light dark; }}
body {{ font-family: system-ui, sans-serif; margin: 0; background: #f5f5f5; color: #111; }}
main {{ max-width: 1400px; margin: 0 auto; padding: 28px; }}
h1 {{ margin-bottom: 4px; }}
.subtitle {{ color: #555; margin-top: 0; }}
.cards {{ display: grid; grid-template-columns: repeat(auto-fit,minmax(150px,1fr)); gap: 12px; margin: 22px 0; }}
.card, section {{ background: white; border: 1px solid #ddd; border-radius: 10px; padding: 16px; }}
.card strong {{ display: block; font-size: 1.8rem; }}
section {{ margin: 18px 0; }}
.file {{ color: #666; }}
.search {{ width: min(420px,100%); padding: 8px 10px; margin: 4px 0 12px; }}
.table-wrap {{ overflow-x: auto; }}
table {{ border-collapse: collapse; width: 100%; font-size: 0.88rem; }}
th, td {{ border-bottom: 1px solid #ddd; padding: 7px 9px; text-align: left; white-space: nowrap; }}
th {{ position: sticky; top: 0; background: #eee; }}
.map {{ max-width: 100%; height: auto; }}
pre {{ overflow-x: auto; padding: 12px; background: #f2f2f2; }}
@media (prefers-color-scheme: dark) {{
  body {{ background:#111; color:#eee; }}
  .card, section {{ background:#1b1b1b; border-color:#333; }}
  .subtitle,.file {{ color:#aaa; }}
  th {{ background:#292929; }}
  th,td {{ border-color:#333; }}
  pre {{ background:#151515; }}
}}
</style>
<script>
function filterTable(id, query) {{
  const q = query.toLowerCase();
  const rows = document.querySelectorAll('#' + id + ' tbody tr');
  rows.forEach(row => {{
    row.style.display = row.innerText.toLowerCase().includes(q) ? '' : 'none';
  }});
}}
</script>
</head>
<body>
<main>
<h1>DIFLOW results report</h1>
<p class="subtitle">{escape(str(root))}</p>
<div class="cards">
  <div class="card"><span>Pairwise rows</span><strong>{pair_count}</strong></div>
  <div class="card"><span>Supported</span><strong>{counts["supported"]}</strong></div>
  <div class="card"><span>Ambiguous</span><strong>{counts["ambiguous"]}</strong></div>
  <div class="card"><span>Unsupported</span><strong>{counts["unsupported"]}</strong></div>
</div>
{map_html}
{''.join(table_sections)}
{provenance_html}
</main>
</body>
</html>
"""
    output.write_text(html, encoding="utf-8")
    return output
