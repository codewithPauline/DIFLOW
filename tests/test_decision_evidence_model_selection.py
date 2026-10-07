from types import SimpleNamespace

import pandas as pd

import diflow.validation.decision_evidence as decision_evidence
from diflow.demography import AsymmetricIMParams
from diflow.validation.decision_evidence import (
    DecisionEvidenceScenario,
    run_decision_evidence_benchmark,
)


def _scenario():
    return DecisionEvidenceScenario(
        "test",
        AsymmetricIMParams(
            nu_a=1.0,
            nu_b=1.0,
            split_time=0.5,
            m_a_to_b=1.0,
            m_b_to_a=0.25,
        ),
        "A->B",
    )


def _counts(*args, **kwargs):
    return pd.DataFrame(
        [
            ["1", 1, "A", "G", "A", 3, 1, 4],
            ["1", 1, "A", "G", "B", 2, 2, 4],
            ["1", 2, "A", "G", "A", 2, 2, 4],
            ["1", 2, "A", "G", "B", 3, 1, 4],
        ],
        columns=[
            "chrom",
            "pos",
            "ref",
            "alt",
            "population",
            "ref_count",
            "alt_count",
            "called_chromosomes",
        ],
    )


def _fit(model_name):
    params = {
        "isolation": {"nu_a": 1.0, "nu_b": 1.0, "split_time": 0.5},
        "symmetric_migration": {
            "nu_a": 1.0,
            "nu_b": 1.0,
            "split_time": 0.5,
            "migration": 0.4,
        },
        "asymmetric_migration": {
            "nu_a": 1.0,
            "nu_b": 1.0,
            "split_time": 0.5,
            "m_a_to_b": 0.8,
            "m_b_to_a": 0.2,
        },
        "secondary_contact_symmetric": {
            "nu_a": 1.0,
            "nu_b": 1.0,
            "isolation_time": 0.5,
            "contact_time": 0.2,
            "migration": 0.4,
        },
        "secondary_contact_asymmetric": {
            "nu_a": 1.0,
            "nu_b": 1.0,
            "isolation_time": 0.5,
            "contact_time": 0.2,
            "m_a_to_b": 0.2,
            "m_b_to_a": 1.1,
        },
    }
    return SimpleNamespace(
        best_parameters=params[model_name],
        best_log_likelihood=-10.0,
        stable=True,
        converged_fraction=1.0,
    )


def test_decision_evidence_records_no_direction_for_symmetric_winner(monkeypatch):
    monkeypatch.setattr(
        decision_evidence,
        "simulate_counts_from_expected_spectrum",
        _counts,
    )
    monkeypatch.setattr(
        decision_evidence,
        "fit_multistart",
        lambda spectrum, model_name, **kwargs: _fit(model_name),
    )
    monkeypatch.setattr(
        decision_evidence,
        "rank_models",
        lambda scores, use_aicc=False: pd.DataFrame(
            {
                "model": [
                    "secondary_contact_symmetric",
                    "secondary_contact_asymmetric",
                    "asymmetric_migration",
                    "symmetric_migration",
                    "isolation",
                ],
                "akaike_weight": [0.8, 0.1, 0.05, 0.03, 0.02],
            }
        ),
    )

    result = run_decision_evidence_benchmark(
        scenarios=(_scenario(),),
        replicates=1,
        chromosomes=4,
        segregating_sites=10,
        starts=1,
        bootstrap_replicates=2,
        bootstrap_starts=1,
    )

    row = result.iloc[0]
    assert row["best_model"] == "secondary_contact_symmetric"
    assert row["directional_model"] is None
    assert row["preferred_direction"] == "none"
    assert row["asymmetric_model_weight"] == 0.0
    assert row["directional_support"] == 0.0


def test_decision_evidence_bootstraps_selected_secondary_contact(monkeypatch):
    seen = []
    monkeypatch.setattr(
        decision_evidence,
        "simulate_counts_from_expected_spectrum",
        _counts,
    )
    monkeypatch.setattr(
        decision_evidence,
        "fit_multistart",
        lambda spectrum, model_name, **kwargs: _fit(model_name),
    )
    monkeypatch.setattr(
        decision_evidence,
        "rank_models",
        lambda scores, use_aicc=False: pd.DataFrame(
            {
                "model": [
                    "secondary_contact_asymmetric",
                    "asymmetric_migration",
                    "secondary_contact_symmetric",
                    "symmetric_migration",
                    "isolation",
                ],
                "akaike_weight": [0.82, 0.08, 0.05, 0.03, 0.02],
            }
        ),
    )

    def fake_bootstrap(*args, model_name, **kwargs):
        seen.append(model_name)
        return SimpleNamespace(
            m_a_to_b_lower=0.1,
            m_a_to_b_upper=0.3,
            m_b_to_a_lower=0.9,
            m_b_to_a_upper=1.3,
            probability_a_to_b_stronger=0.05,
            preferred_direction_support=0.95,
        )

    monkeypatch.setattr(
        decision_evidence,
        "bootstrap_asymmetric_jsfs",
        fake_bootstrap,
    )

    result = run_decision_evidence_benchmark(
        scenarios=(_scenario(),),
        replicates=1,
        chromosomes=4,
        segregating_sites=10,
        starts=1,
        bootstrap_replicates=2,
        bootstrap_starts=1,
    )

    row = result.iloc[0]
    assert seen == ["secondary_contact_asymmetric"]
    assert row["directional_model"] == "secondary_contact_asymmetric"
    assert row["estimated_m_a_to_b"] == 0.2
    assert row["estimated_m_b_to_a"] == 1.1
    assert row["preferred_direction"] == "B->A"
