import pandas as pd

from diflow.pipeline.outputs import pairwise_to_flow_table


def test_pairwise_to_flow_table_keeps_formal_edges_only():
    pairwise = pd.DataFrame(
        [
            {
                "population_a": "A",
                "population_b": "B",
                "status": "supported",
                "preferred_direction": "A->B",
                "m_a_to_b_scaled": 0.03,
                "m_b_to_a_scaled": 0.005,
                "directional_support": 0.99,
                "asymmetry_index": 0.71,
                "asymmetric_model_weight": 0.85,
            },
            {
                "population_a": "B",
                "population_b": "C",
                "status": "candidate",
                "preferred_direction": "B->C",
                "m_a_to_b_scaled": 0.02,
                "m_b_to_a_scaled": 0.01,
                "directional_support": 0.8,
                "asymmetry_index": 0.33,
                "asymmetric_model_weight": 0.75,
            },
            {
                "population_a": "C",
                "population_b": "D",
                "status": "ambiguous",
                "preferred_direction": "D->C",
                "m_a_to_b_scaled": 0.01,
                "m_b_to_a_scaled": 0.02,
                "directional_support": 0.90,
                "asymmetry_index": -0.33,
                "asymmetric_model_weight": 0.80,
            },
        ]
    )

    flows = pairwise_to_flow_table(pairwise)

    assert len(flows) == 2
    assert list(flows["source"]) == ["A", "D"]
    assert list(flows["destination"]) == ["B", "C"]
