import math

import pytest

from diflow.demography import ModelScore, aic, aicc, rank_models


def test_aic_formula():
    assert aic(-100.0, 5) == pytest.approx(210.0)


def test_aicc_exceeds_aic_when_defined():
    assert aicc(-100.0, 5, 100) > aic(-100.0, 5)


def test_rank_models_and_weights_sum_to_one():
    scores = [
        ModelScore("isolation", -120.0, 3, 50),
        ModelScore("symmetric", -110.0, 4, 50),
        ModelScore("asymmetric", -103.0, 5, 50),
    ]
    ranked = rank_models(scores, use_aicc=True)

    assert ranked.iloc[0]["model"] == "asymmetric"
    assert ranked["akaike_weight"].sum() == pytest.approx(1.0)
    assert ranked.iloc[0]["delta"] == pytest.approx(0.0)


def test_aicc_returns_infinity_for_too_few_observations():
    assert math.isinf(aicc(-10.0, 5, 6))



def test_symmetric_secondary_contact_is_registered():
    from diflow.demography.fit_models import MODEL_SPECS

    spec = MODEL_SPECS["secondary_contact_symmetric"]
    assert spec["names"] == [
        "nu_a",
        "nu_b",
        "isolation_time",
        "contact_time",
        "migration",
    ]
    assert len(spec["initial"]) == 5
    assert len(spec["lower"]) == 5
    assert len(spec["upper"]) == 5
