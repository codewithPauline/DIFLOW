import numpy as np
import pandas as pd
import pytest

from diflow.spectra import pairwise_projected_jsfs


def test_projection_retains_loci_with_extra_called_chromosomes():
    counts = pd.DataFrame(
        [
            ["1", 10, "A", "G", "A", 4, 2, 6],
            ["1", 10, "A", "G", "B", 1, 3, 4],
            ["1", 20, "C", "T", "A", 2, 2, 4],
            ["1", 20, "C", "T", "B", 4, 2, 6],
        ],
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

    result = pairwise_projected_jsfs(
        counts,
        "A",
        "B",
        chromosomes_a=4,
        chromosomes_b=4,
    )

    assert result.loci_used == 2
    assert result.loci_skipped == 0
    assert result.spectrum.shape == (5, 5)
    assert result.spectrum.sum() == pytest.approx(2.0)
    assert np.all(result.spectrum >= 0)


def test_projection_skips_locus_below_target():
    counts = pd.DataFrame(
        [
            ["1", 10, "A", "G", "A", 1, 1, 2],
            ["1", 10, "A", "G", "B", 2, 2, 4],
        ],
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

    result = pairwise_projected_jsfs(
        counts,
        "A",
        "B",
        chromosomes_a=4,
        chromosomes_b=4,
    )

    assert result.loci_used == 0
    assert result.loci_skipped == 1
    assert result.spectrum.sum() == pytest.approx(0.0)
