"""Pairwise joint site-frequency spectrum utilities."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class PairwiseJSFS:
    """A two-population unfolded joint allele-count spectrum."""

    population_a: str
    population_b: str
    chromosomes_a: int
    chromosomes_b: int
    spectrum: np.ndarray
    loci_used: int
    loci_skipped: int


def pairwise_jsfs(
    counts: pd.DataFrame,
    population_a: str,
    population_b: str,
    *,
    chromosomes_a: int | None = None,
    chromosomes_b: int | None = None,
) -> PairwiseJSFS:
    """Build an unfolded 2D jSFS from long-form population allele counts.

    The current implementation uses only loci with exactly the requested number
    of called chromosomes in both populations. This avoids silently mixing
    different sample sizes. Statistical projection to smaller sample sizes will
    be implemented separately.

    The spectrum index [i, j] is the number of loci with i alternate alleles in
    population A and j alternate alleles in population B.
    """
    required = {
        "chrom",
        "pos",
        "ref",
        "alt",
        "population",
        "alt_count",
        "called_chromosomes",
    }
    if not required.issubset(counts.columns):
        raise ValueError("counts table is missing required allele-count columns.")
    if population_a == population_b:
        raise ValueError("population_a and population_b must be different.")

    sub = counts[counts["population"].isin([population_a, population_b])].copy()
    if sub.empty:
        raise ValueError("requested populations are absent from the counts table.")

    key = ["chrom", "pos", "ref", "alt"]
    alt = sub.pivot_table(index=key, columns="population", values="alt_count", aggfunc="first")
    called = sub.pivot_table(
        index=key,
        columns="population",
        values="called_chromosomes",
        aggfunc="first",
    )

    if population_a not in alt.columns or population_b not in alt.columns:
        raise ValueError("both requested populations must be represented.")

    complete = alt[[population_a, population_b]].notna().all(axis=1)
    complete &= called[[population_a, population_b]].notna().all(axis=1)

    if chromosomes_a is None:
        vals = called.loc[complete, population_a].astype(int)
        if vals.empty:
            raise ValueError("no loci have data in both populations.")
        chromosomes_a = int(vals.mode().iloc[0])
    if chromosomes_b is None:
        vals = called.loc[complete, population_b].astype(int)
        if vals.empty:
            raise ValueError("no loci have data in both populations.")
        chromosomes_b = int(vals.mode().iloc[0])

    if chromosomes_a < 1 or chromosomes_b < 1:
        raise ValueError("chromosome sample sizes must be positive.")

    usable = complete.copy()
    usable &= called[population_a].fillna(-1).astype(int) == chromosomes_a
    usable &= called[population_b].fillna(-1).astype(int) == chromosomes_b

    spectrum = np.zeros((chromosomes_a + 1, chromosomes_b + 1), dtype=int)
    for idx in alt.index[usable]:
        i = int(alt.at[idx, population_a])
        j = int(alt.at[idx, population_b])
        if 0 <= i <= chromosomes_a and 0 <= j <= chromosomes_b:
            spectrum[i, j] += 1

    total_candidates = int(len(alt.index))
    loci_used = int(spectrum.sum())
    return PairwiseJSFS(
        population_a=population_a,
        population_b=population_b,
        chromosomes_a=chromosomes_a,
        chromosomes_b=chromosomes_b,
        spectrum=spectrum,
        loci_used=loci_used,
        loci_skipped=total_candidates - loci_used,
    )
