import pytest

from diflow.demography import AsymmetricIMParams, to_dadi_split_asym_mig


def test_diflow_to_dadi_direction_mapping():
    params = AsymmetricIMParams(
        nu_a=1.2,
        nu_b=0.8,
        split_time=0.4,
        m_a_to_b=3.0,
        m_b_to_a=0.5,
    )

    dadi_params = to_dadi_split_asym_mig(params)

    assert dadi_params == pytest.approx((1.2, 0.8, 0.4, 0.5, 3.0))


def test_parameters_reject_zero_or_negative_values():
    with pytest.raises(ValueError):
        AsymmetricIMParams(
            nu_a=1.0,
            nu_b=1.0,
            split_time=0.5,
            m_a_to_b=0.0,
            m_b_to_a=0.1,
        )
