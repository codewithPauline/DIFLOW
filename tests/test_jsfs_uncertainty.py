from types import SimpleNamespace

import pandas as pd
import pytest

from diflow.demography.jsfs_uncertainty import bootstrap_asymmetric_jsfs


def _counts():
    rows = []
    for pos, a_alt, b_alt in [
        (10, 1, 3),
        (20, 1, 2),
        (30, 0, 3),
        (40, 2, 3),
        (50, 1, 4),
        (60, 0, 2),
    ]:
        rows.append(["1", pos, "A", "G", "A", 4 - a_alt, a_alt, 4])
        rows.append(["1", pos, "A", "G", "B", 4 - b_alt, b_alt, 4])
    return pd.DataFrame(
        rows,
        columns=[
            "chrom",
            "pos",
            "ref",
            "alt",
            "population",
            "ref_count",
            "alt_count",
            "called_chromosomes",
        ],
    )


def fake_fit(spectrum, model_name, *, starts, seed, maxiter):
    # Deterministic test backend whose directional estimates depend on spectrum.
    axis_a = sum(i * spectrum[i, :].sum() for i in range(spectrum.shape[0]))
    axis_b = sum(j * spectrum[:, j].sum() for j in range(spectrum.shape[1]))
    scale = max(float(spectrum.sum()), 1.0)
    m_ab = 0.02 + 0.002 * axis_b / scale
    m_ba = 0.005 + 0.0005 * axis_a / scale
    return SimpleNamespace(
        best_parameters={"m_a_to_b": m_ab, "m_b_to_a": m_ba}
    )


def test_jsfs_bootstrap_returns_direction_support():
    result = bootstrap_asymmetric_jsfs(
        _counts(),
        "A",
        "B",
        chromosomes_a=4,
        chromosomes_b=4,
        replicates=20,
        seed=9,
        fit_function=fake_fit,
    )

    assert result.successful_replicates == 20
    assert result.attempted_replicates == 20
    assert result.m_a_to_b_mean > result.m_b_to_a_mean
    assert result.probability_a_to_b_stronger == pytest.approx(1.0)
    assert result.preferred_direction_support == pytest.approx(1.0)
    assert result.m_a_to_b_lower <= result.m_a_to_b_mean <= result.m_a_to_b_upper


def test_jsfs_bootstrap_requires_multiple_loci():
    tiny = _counts().query("pos == 10")
    with pytest.raises(ValueError):
        bootstrap_asymmetric_jsfs(
            tiny,
            "A",
            "B",
            chromosomes_a=4,
            chromosomes_b=4,
            replicates=10,
            fit_function=fake_fit,
        )



def test_block_bootstrap_records_resampling_units():
    result = bootstrap_asymmetric_jsfs(
        _counts(),
        "A",
        "B",
        chromosomes_a=4,
        chromosomes_b=4,
        replicates=12,
        seed=3,
        block_size_bp=25,
        fit_function=fake_fit,
    )

    assert result.resampling_unit == "25-bp genomic block"
    assert result.blocks_used == 3
    assert result.loci_used == 6
    assert result.successful_replicates == 12


def test_block_bootstrap_requires_multiple_blocks():
    with pytest.raises(ValueError, match="at least two resampling blocks"):
        bootstrap_asymmetric_jsfs(
            _counts(),
            "A",
            "B",
            chromosomes_a=4,
            chromosomes_b=4,
            replicates=10,
            block_size_bp=1000,
            fit_function=fake_fit,
        )


def test_block_bootstrap_rejects_nonpositive_window():
    with pytest.raises(ValueError, match="block_size_bp"):
        bootstrap_asymmetric_jsfs(
            _counts(),
            "A",
            "B",
            chromosomes_a=4,
            chromosomes_b=4,
            replicates=10,
            block_size_bp=0,
            fit_function=fake_fit,
        )



def test_directional_support_helper_is_public():
    from diflow.demography import directional_support_for_estimate

    assert directional_support_for_estimate(0.98, 1.0, 0.2) == pytest.approx(0.98)
    assert directional_support_for_estimate(0.98, 0.2, 1.0) == pytest.approx(0.02)
    assert directional_support_for_estimate(0.98, 0.5, 0.5) == pytest.approx(0.5)



def test_jsfs_bootstrap_uses_requested_directional_model():
    seen = []

    def recording_fit(spectrum, model_name, *, starts, seed, maxiter):
        seen.append(model_name)
        return SimpleNamespace(
            best_parameters={"m_a_to_b": 0.8, "m_b_to_a": 0.2}
        )

    result = bootstrap_asymmetric_jsfs(
        _counts(),
        "A",
        "B",
        chromosomes_a=4,
        chromosomes_b=4,
        replicates=6,
        seed=11,
        fit_function=recording_fit,
        model_name="secondary_contact_asymmetric",
    )

    assert result.successful_replicates == 6
    assert set(seen) == {"secondary_contact_asymmetric"}
