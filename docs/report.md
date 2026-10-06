# Interactive HTML report

DIFLOW can convert a completed results directory into a self-contained HTML
report:

    diflow report \
      --results results/ \
      --output results/report.html

If `--output` is omitted, the report is written to
`results/report.html`.

The report includes:

- summary cards for pairwise evidence status
- searchable pairwise inference results
- searchable model rankings
- directional-flow and network tables when available
- run metadata
- the directional map when available
- collapsible run provenance JSON

The report is deliberately self-contained and does not require a dashboard
server. It is intended for exploration, review, and sharing alongside the
underlying CSV/PDF outputs.

The HTML report does not replace the statistical tables. It is a presentation
layer over the same evidence.
