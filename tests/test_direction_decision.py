import pytest

from diflow.decision import DirectionEvidence, classify_direction, decisions_to_flow_table


def test_supported_direction():
    evidence = DirectionEvidence(
        source="A",
        destination="B",
        migration_forward=0.03,
        migration_reverse=0.005,
        directional_support=0.99,
        asymmetric_model_weight=0.85,
        optimizer_stable=True,
    )
    result = classify_direction(evidence)
    assert result.status == "supported"
    assert result.preferred_direction == "A->B"
    assert result.asymmetry_index > 0


def test_low_model_weight_is_unsupported():
    evidence = DirectionEvidence(
        source="A",
        destination="B",
        migration_forward=0.03,
        migration_reverse=0.005,
        directional_support=0.99,
        asymmetric_model_weight=0.20,
        optimizer_stable=True,
    )
    result = classify_direction(evidence)
    assert result.status == "unsupported"


def test_unstable_fit_is_ambiguous():
    evidence = DirectionEvidence(
        source="A",
        destination="B",
        migration_forward=0.03,
        migration_reverse=0.005,
        directional_support=0.99,
        asymmetric_model_weight=0.90,
        optimizer_stable=False,
    )
    result = classify_direction(evidence)
    assert result.status == "ambiguous"
    assert "optimizer instability" in result.reason


def test_interval_separation_can_be_required():
    evidence = DirectionEvidence(
        source="A",
        destination="B",
        migration_forward=0.03,
        migration_reverse=0.01,
        directional_support=0.99,
        asymmetric_model_weight=0.90,
        optimizer_stable=True,
        forward_lower=0.02,
        forward_upper=0.04,
        reverse_lower=0.008,
        reverse_upper=0.015,
    )
    result = classify_direction(evidence, require_interval_separation=True)
    assert result.status == "supported"


def test_flow_table_keeps_supported_only_by_default():
    supported = classify_direction(
        DirectionEvidence(
            source="A",
            destination="B",
            migration_forward=0.03,
            migration_reverse=0.005,
            directional_support=0.99,
            asymmetric_model_weight=0.90,
            optimizer_stable=True,
        )
    )
    ambiguous = classify_direction(
        DirectionEvidence(
            source="B",
            destination="C",
            migration_forward=0.02,
            migration_reverse=0.015,
            directional_support=0.80,
            asymmetric_model_weight=0.90,
            optimizer_stable=True,
        )
    )
    table = decisions_to_flow_table([supported, ambiguous])
    assert len(table) == 1
    assert table.iloc[0]["source"] == "A"
    assert table.iloc[0]["destination"] == "B"
