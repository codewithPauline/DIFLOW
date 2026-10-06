# External method comparison

DIFLOW provides a standardized comparison layer for simulation results from
other migration methods.

## Standard input schema

Each method CSV must contain:

- `scenario`
- `replicate`
- `truth_m_a_to_b`
- `truth_m_b_to_a`
- `estimated_m_a_to_b`
- `estimated_m_b_to_a`

Compare two or more methods with:

    diflow compare \
      --method DIFLOW=diflow_results.csv \
      --method OTHER=other_method_results.csv \
      --output comparison/

Outputs include combined replicate records, method/scenario summaries, and
direction-accuracy / false-direction plots when applicable.

## Scientific requirement

Only compare quantities with genuinely compatible interpretations.

A recent-migrant assignment estimate, a backward-time lineage migration rate,
a historical coalescent migration parameter, and a dadi-scaled demographic
migration parameter are not automatically interchangeable.

The comparison framework is implemented; the final established-method
benchmark study still requires actually running and normalizing the selected
external tools under matched simulated truth.
