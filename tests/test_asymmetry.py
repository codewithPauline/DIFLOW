import math

import pytest

from diflow.core import asymmetry_index, directional_ratio


def test_asymmetry_is_zero_for_equal_rates():
    assert asymmetry_index(0.02, 0.02) == pytest.approx(0.0)


def test_asymmetry_direction():
    assert asymmetry_index(0.03, 0.01) == pytest.approx(0.5)
    assert asymmetry_index(0.01, 0.03) == pytest.approx(-0.5)


def test_asymmetry_bounds():
    assert asymmetry_index(0.01, 0.0) == pytest.approx(1.0)
    assert asymmetry_index(0.0, 0.01) == pytest.approx(-1.0)


def test_both_zero():
    assert asymmetry_index(0.0, 0.0) == pytest.approx(0.0)
    assert directional_ratio(0.0, 0.0) == pytest.approx(1.0)


def test_ratio_zero_denominator():
    assert math.isinf(directional_ratio(0.01, 0.0))


def test_negative_rate_rejected():
    with pytest.raises(ValueError):
        asymmetry_index(-0.01, 0.02)
