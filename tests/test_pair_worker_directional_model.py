from types import SimpleNamespace

import numpy as np
import pandas as pd

import diflow.pipeline.pair_worker as pair_worker
from diflow.pipeline.pair_worker import PairInferenceTask, infer_pair_task


def _task():
    return PairInferenceTask(
        pair_index=0,
        population_a="A",
        population_b="B",
        distance_km=10.0,
        counts=pd.DataFrame(),
        projection_chromosomes=4,
        starts=2,
        maxiter=10,
        bootstrap_replicates=0,
        bootstrap_starts=2,
        bootstrap_block_bp=None,
        min_model_weight=0.7,
        min_directional_support=0.95,
        min_abs_asymmetry=0.25,
        polarized=False,
        prepare_only=False,
        seed=7,
    )


def _projected():
    spectrum = np.ones((5, 5), dtype=float)
    return SimpleNamespace(
        spectrum=spectrum,
        loci_used=10,
        loci_skipped=0,
    )


def _fake_fit_factory():
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
            "m_a_to_b": 0.9,
            "m_b_to_a": 0.1,
        },
        "secondary_contact_symmetric": {
            "nu_a": 1.0,
            "nu_b": 1.0,
            "isolation_time": 0.6,
            "contact_time": 0.2,
            "migration": 0.5,
        },
        "secondary_contact_asymmetric": {
            "nu_a": 1.0,
            "nu_b": 1.0,
            "isolation_time": 0.6,
            "contact_time": 0.2,
            "m_a_to_b": 0.2,
            "m_b_to_a": 1.1,
        },
    }

    def fake_fit(spectrum, model_name, **kwargs):
        return SimpleNamespace(
            best_parameters=params[model_name],
            best_log_likelihood=-10.0,
            stable=True,
            converged_fraction=1.0,
        )

    return fake_fit


def test_pair_worker_does_not_emit_direction_for_symmetric_winner(monkeypatch):
    monkeypatch.setattr(
        pair_worker,
        "pairwise_projected_jsfs",
        lambda *args, **kwargs: _projected(),
    )
    monkeypatch.setattr(pair_worker, "fit_multistart", _fake_fit_factory())
    monkeypatch.setattr(
        pair_worker,
        "rank_models",
        lambda scores, use_aicc=False: pd.DataFrame(
            {
                "model": [
                    "symmetric_migration",
                    "asymmetric_migration",
                    "secondary_contact_asymmetric",
                    "isolation",
                ],
                "akaike_weight": [0.75, 0.15, 0.08, 0.02],
            }
        ),
    )

    result = infer_pair_task(_task())

    assert result.row["best_model"] == "symmetric_migration"
    assert result.row["status"] == "unsupported"
    assert result.row["directional_model"] is None
    assert result.row["preferred_direction"] is None


def test_pair_worker_uses_secondary_contact_winner_parameters(monkeypatch):
    monkeypatch.setattr(
        pair_worker,
        "pairwise_projected_jsfs",
        lambda *args, **kwargs: _projected(),
    )
    monkeypatch.setattr(pair_worker, "fit_multistart", _fake_fit_factory())
    monkeypatch.setattr(
        pair_worker,
        "rank_models",
        lambda scores, use_aicc=False: pd.DataFrame(
            {
                "model": [
                    "secondary_contact_asymmetric",
                    "asymmetric_migration",
                    "symmetric_migration",
                    "isolation",
                ],
                "akaike_weight": [0.82, 0.10, 0.06, 0.02],
            }
        ),
    )

    result = infer_pair_task(_task())

    assert result.row["best_model"] == "secondary_contact_asymmetric"
    assert result.row["directional_model"] == "secondary_contact_asymmetric"
    assert result.row["m_a_to_b_scaled"] == 0.2
    assert result.row["m_b_to_a_scaled"] == 1.1
    assert result.row["preferred_direction"] == "B->A"



def test_pair_worker_does_not_emit_direction_for_symmetric_secondary_contact(
    monkeypatch,
):
    monkeypatch.setattr(
        pair_worker,
        "pairwise_projected_jsfs",
        lambda *args, **kwargs: _projected(),
    )
    monkeypatch.setattr(pair_worker, "fit_multistart", _fake_fit_factory())
    monkeypatch.setattr(
        pair_worker,
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
                "akaike_weight": [0.80, 0.10, 0.05, 0.03, 0.02],
            }
        ),
    )

    result = infer_pair_task(_task())

    assert result.row["best_model"] == "secondary_contact_symmetric"
    assert result.row["status"] == "unsupported"
    assert result.row["directional_model"] is None
    assert result.row["preferred_direction"] is None
