import numpy as np
import pytest

dadi = pytest.importorskip("dadi")

from diflow.demography import AsymmetricIMParams, expected_spectrum, fit_multistart


def test_dadi_backend_generates_and_fits_folded_spectrum():
    truth = AsymmetricIMParams(
        nu_a=1.0,
        nu_b=1.0,
        split_time=0.5,
        m_a_to_b=1.0,
        m_b_to_a=0.25,
    )
    expected = expected_spectrum(truth, (6, 6))
    observed = np.asarray(expected, dtype=float)
    observed[0, 0] = 0.0
    observed[-1, -1] = 0.0

    fit = fit_multistart(
        observed,
        "asymmetric_migration",
        starts=2,
        maxiter=5,
        seed=7,
        polarized=False,
    )
    assert np.isfinite(fit.best_log_likelihood)
    assert fit.best_parameters["m_a_to_b"] > 0
    assert fit.best_parameters["m_b_to_a"] > 0
