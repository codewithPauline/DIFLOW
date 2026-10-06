import pandas as pd

from diflow.spectra import pairwise_jsfs


def test_pairwise_jsfs_counts_cells():
    counts = pd.DataFrame(
        [
            ["1", 10, "A", "G", "A", 3, 1, 4],
            ["1", 10, "A", "G", "B", 1, 3, 4],
            ["1", 20, "C", "T", "A", 2, 2, 4],
            ["1", 20, "C", "T", "B", 4, 0, 4],
            ["1", 30, "G", "A", "A", 1, 1, 2],
            ["1", 30, "G", "A", "B", 2, 2, 4],
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

    result = pairwise_jsfs(
        counts,
        "A",
        "B",
        chromosomes_a=4,
        chromosomes_b=4,
    )

    assert result.spectrum.shape == (5, 5)
    assert result.spectrum[1, 3] == 1
    assert result.spectrum[2, 0] == 1
    assert result.loci_used == 2
    assert result.loci_skipped == 1
