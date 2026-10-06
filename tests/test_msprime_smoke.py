import pytest

msprime = pytest.importorskip("msprime")

from diflow.validation.mechanistic_linkage import (
    MechanisticLinkageScenario,
    simulate_msprime_counts,
)


def test_msprime_mechanistic_simulator_generates_two_population_counts():
    scenario = MechanisticLinkageScenario(
        name="smoke",
        m_a_to_b=5e-5,
        m_b_to_a=1.25e-5,
        expected_direction="A->B",
    )
    counts = simulate_msprime_counts(
        scenario,
        chromosomes_per_population=6,
        nref=1000,
        split_time_scaled=0.5,
        sequence_length=100000,
        recombination_rate=1e-8,
        mutation_rate=1e-7,
        seed=23,
    )
    assert not counts.empty
    assert set(counts["population"]) == {"A", "B"}
    assert counts["pos"].nunique() >= 2
    assert (counts["called_chromosomes"] == 6).all()
