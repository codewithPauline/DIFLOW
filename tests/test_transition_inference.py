import pytest

from diflow.inference import estimate_one_generation
from diflow.simulation import simulate_two_population


def test_recovers_asymmetric_migration_in_simple_case():
    truth_a_to_b = 0.03
    truth_b_to_a = 0.005

    sim = simulate_two_population(
        generations=1,
        loci=20000,
        ne_a=5000,
        ne_b=5000,
        m_a_to_b=truth_a_to_b,
        m_b_to_a=truth_b_to_a,
        seed=2026,
    )

    fit = estimate_one_generation(
        p_a_initial=sim["p_a_initial"],
        p_b_initial=sim["p_b_initial"],
        p_a_final=sim["p_a_final"],
        p_b_final=sim["p_b_final"],
        ne_a=5000,
        ne_b=5000,
    )

    assert fit.success
    assert fit.m_a_to_b == pytest.approx(truth_a_to_b, abs=0.002)
    assert fit.m_b_to_a == pytest.approx(truth_b_to_a, abs=0.002)
    assert fit.m_a_to_b > fit.m_b_to_a


def test_recovers_near_symmetric_migration():
    sim = simulate_two_population(
        generations=1,
        loci=15000,
        ne_a=4000,
        ne_b=6000,
        m_a_to_b=0.02,
        m_b_to_a=0.02,
        seed=11,
    )

    fit = estimate_one_generation(
        p_a_initial=sim["p_a_initial"],
        p_b_initial=sim["p_b_initial"],
        p_a_final=sim["p_a_final"],
        p_b_final=sim["p_b_final"],
        ne_a=4000,
        ne_b=6000,
    )

    assert fit.success
    assert fit.m_a_to_b == pytest.approx(0.02, abs=0.0025)
    assert fit.m_b_to_a == pytest.approx(0.02, abs=0.0025)


def test_rejects_mismatched_vectors():
    with pytest.raises(ValueError):
        estimate_one_generation(
            p_a_initial=[0.1, 0.2],
            p_b_initial=[0.2],
            p_a_final=[0.1, 0.2],
            p_b_final=[0.2, 0.3],
            ne_a=1000,
            ne_b=1000,
        )
