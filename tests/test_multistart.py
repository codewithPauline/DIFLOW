import numpy as np
import pytest

from diflow.demography.multistart import (
    _relative_spread,
    generate_starting_points,
)


def test_starting_points_reproducible():
    a = generate_starting_points("asymmetric_migration", starts=5, seed=42)
    b = generate_starting_points("asymmetric_migration", starts=5, seed=42)
    assert a == b


def test_first_start_is_canonical():
    points = generate_starting_points("isolation", starts=3, seed=1)
    assert points[0] == pytest.approx([1.0, 1.0, 0.5])


def test_generated_points_within_bounds():
    points = generate_starting_points("symmetric_migration", starts=20, seed=9)
    for point in points:
        nu_a, nu_b, split_time, migration = point
        assert 1e-3 <= nu_a <= 100
        assert 1e-3 <= nu_b <= 100
        assert 1e-4 <= split_time <= 20
        assert 1e-5 <= migration <= 50


def test_relative_spread_detects_stability():
    assert _relative_spread(np.array([1.0, 1.02, 0.98])) < 0.1
    assert _relative_spread(np.array([1.0, 4.0])) > 1.0
