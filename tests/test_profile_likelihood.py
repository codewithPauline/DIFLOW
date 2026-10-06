import pandas as pd
import pytest

from diflow.demography.profile_likelihood import (
    confidence_interval_from_profile,
    profile_grid,
)


def test_profile_grid_contains_mle_and_is_bounded():
    values = profile_grid(
        1.0,
        lower=1e-5,
        upper=50.0,
        points=9,
        fold_range=5.0,
    )
    assert 1.0 in values
    assert values.min() >= 0.2
    assert values.max() <= 5.0


def test_profile_interval_uses_likelihood_cutoff():
    table = pd.DataFrame(
        {
            "fixed_value": [0.25, 0.5, 1.0, 2.0, 4.0],
            "log_likelihood": [-12.5, -10.8, -10.0, -10.9, -13.0],
        }
    )
    lower, upper = confidence_interval_from_profile(table, confidence=0.95)
    assert lower == pytest.approx(0.5)
    assert upper == pytest.approx(2.0)
