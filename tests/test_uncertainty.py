import pytest

from diflow.inference import bootstrap_one_generation
from diflow.simulation import simulate_two_population


def test_bootstrap_supports_known_stronger_direction():
    sim = simulate_two_population(
        generations=1,
        loci=3000,
        ne_a=3000,
        ne_b=3000,
        m_a_to_b=0.04,
        m_b_to_a=0.005,
        seed=99,
    )

    result = bootstrap_one_generation(
        p_a_initial=sim["p_a_initial"],
        p_b_initial=sim["p_b_initial"],
        p_a_final=sim["p_a_final"],
        p_b_final=sim["p_b_final"],
        ne_a=3000,
        ne_b=3000,
        replicates=40,
        seed=123,
    )

    assert result.m_a_to_b_lower < 0.04 < result.m_a_to_b_upper
    assert result.m_b_to_a_lower <= 0.005 <= result.m_b_to_a_upper
    assert result.probability_a_to_b_stronger > 0.95
    assert result.replicates == 40


def test_bootstrap_validates_arguments():
    with pytest.raises(ValueError):
        bootstrap_one_generation(
            p_a_initial=[0.1, 0.2],
            p_b_initial=[0.2, 0.3],
            p_a_final=[0.1, 0.2],
            p_b_final=[0.2, 0.3],
            ne_a=1000,
            ne_b=1000,
            replicates=1,
        )
