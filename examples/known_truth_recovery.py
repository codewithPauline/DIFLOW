"""Run a simple known-truth directional migration recovery benchmark."""

from diflow.core import asymmetry_index
from diflow.inference import estimate_one_generation
from diflow.simulation import simulate_two_population


def main() -> None:
    truth_a_to_b = 0.03
    truth_b_to_a = 0.005

    sim = simulate_two_population(
        generations=1,
        loci=20000,
        ne_a=5000,
        ne_b=5000,
        m_a_to_b=truth_a_to_b,
        m_b_to_a=truth_b_to_a,
        seed=2026,
    )

    fit = estimate_one_generation(
        p_a_initial=sim["p_a_initial"],
        p_b_initial=sim["p_b_initial"],
        p_a_final=sim["p_a_final"],
        p_b_final=sim["p_b_final"],
        ne_a=5000,
        ne_b=5000,
    )

    print("Known truth")
    print(f"  A -> B: {truth_a_to_b:.6f}")
    print(f"  B -> A: {truth_b_to_a:.6f}")
    print()
    print("DIFLOW estimate")
    print(f"  A -> B: {fit.m_a_to_b:.6f}")
    print(f"  B -> A: {fit.m_b_to_a:.6f}")
    print(f"  asymmetry: {asymmetry_index(fit.m_a_to_b, fit.m_b_to_a):.6f}")
    print(f"  optimizer success: {fit.success}")


if __name__ == "__main__":
    main()
