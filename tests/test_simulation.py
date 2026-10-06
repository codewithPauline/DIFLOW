import pandas as pd
import pytest

from diflow.simulation import simulate_two_population


def test_simulator_shape_and_columns():
    out = simulate_two_population(
        generations=5,
        loci=100,
        ne_a=1000,
        ne_b=1000,
        m_a_to_b=0.02,
        m_b_to_a=0.005,
        seed=42,
    )
    assert isinstance(out, pd.DataFrame)
    assert len(out) == 100
    assert {"p_a_final", "p_b_final"}.issubset(out.columns)


def test_simulator_is_reproducible():
    kwargs = dict(
        generations=3,
        loci=20,
        ne_a=500,
        ne_b=800,
        m_a_to_b=0.03,
        m_b_to_a=0.01,
        seed=7,
    )
    a = simulate_two_population(**kwargs)
    b = simulate_two_population(**kwargs)
    pd.testing.assert_frame_equal(a, b)


def test_invalid_migration_rate_rejected():
    with pytest.raises(ValueError):
        simulate_two_population(
            generations=1,
            loci=10,
            ne_a=100,
            ne_b=100,
            m_a_to_b=1.1,
            m_b_to_a=0.0,
        )
