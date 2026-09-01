"""Controls and separation metrics — the null machinery, with no model required."""

from __future__ import annotations

import pytest

from pepdesign.controls import (
    CONTROL_KINDS,
    composition_distance,
    make_controls,
    scramble,
)
from pepdesign.evaluate import cohens_d, roc_auc, separate, threshold_at_control_rate

PEPTIDES = [("p1", "ACDEFGHIK"), ("p2", "LMNPQRSTV"), ("p3", "WYACDEFGH")]


def test_scrambling_preserves_composition_exactly():
    original = "ACDEFGHIKLMNPQ"
    shuffled = scramble(original, seed=7)
    assert sorted(shuffled) == sorted(original)


def test_scrambled_controls_have_zero_composition_distance():
    """The check that 'composition-preserving' is true and not just claimed."""
    controls = make_controls(PEPTIDES, kinds=["scrambled"], seed=0)
    assert composition_distance(
        [sequence for _, sequence in PEPTIDES], [c.sequence for c in controls]
    ) == pytest.approx(0.0)


def test_length_matched_controls_match_length_and_not_composition():
    controls = make_controls(PEPTIDES, kinds=["length_matched"], seed=0)
    for control in controls:
        source = dict(PEPTIDES)[control.derived_from]
        assert len(control.sequence) == len(source)


def test_every_control_family_is_produced():
    controls = make_controls(PEPTIDES, kinds=CONTROL_KINDS, seed=0)
    assert {c.kind for c in controls} == set(CONTROL_KINDS)
    assert len(controls) == len(PEPTIDES) * len(CONTROL_KINDS)


def test_scrambles_of_different_peptides_use_different_seeds():
    """One shared shuffle would correlate the controls and shrink the effective null."""
    repeated = [("a", "ACDEFGHIK"), ("b", "ACDEFGHIK")]
    controls = make_controls(repeated, kinds=["scrambled"], seed=0)
    assert controls[0].sequence != controls[1].sequence


def test_unknown_control_kind_is_refused():
    with pytest.raises(ValueError, match="unknown control kind"):
        make_controls(PEPTIDES, kinds=["nonsense"])


def test_auc_is_one_for_perfect_separation():
    assert roc_auc([5.0, 6.0, 7.0], [1.0, 2.0, 3.0]) == pytest.approx(1.0)


def test_auc_is_half_for_identical_populations():
    """All ties. A trapezoid implementation rounds these; the rank one gives 0.5."""
    assert roc_auc([1.0] * 5, [1.0] * 5) == pytest.approx(0.5)


def test_auc_is_zero_when_the_order_is_reversed():
    assert roc_auc([1.0, 2.0], [8.0, 9.0]) == pytest.approx(0.0)


def test_auc_needs_both_populations():
    with pytest.raises(ValueError, match="both populations"):
        roc_auc([1.0], [])


def test_threshold_admits_the_requested_control_fraction():
    controls = [float(i) for i in range(100)]
    threshold = threshold_at_control_rate(controls, rate=0.05)
    admitted = sum(1 for value in controls if value >= threshold) / len(controls)
    assert admitted == pytest.approx(0.05, abs=0.02)


def test_cohens_d_is_zero_for_identical_distributions():
    assert cohens_d([1.0, 2.0, 3.0], [1.0, 2.0, 3.0]) == pytest.approx(0.0)


def test_separation_reports_a_confidence_interval_that_brackets_the_point():
    real = [float(i) for i in range(40)]
    control = [float(i) - 5 for i in range(40)]
    result = separate(real, control, "scrambled")
    assert result.auc_ci_low <= result.auc <= result.auc_ci_high
    assert result.control_kind == "scrambled"
