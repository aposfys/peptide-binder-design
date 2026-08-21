"""Guards against the two ways this pipeline could lie: circularity and missing nulls."""

from __future__ import annotations

import pytest

from pepdesign import GENERATOR_MODELS, SCORER_MODELS
from pepdesign.filters import Design, assert_no_shared_model, pass_rate, separation


def make_design(design_id: str, *, good: bool) -> Design:
    return Design(
        design_id=design_id,
        sequence="ACDEFGHIKL",
        self_consistency_rmsd=1.2 if good else 5.0,
        interface_ptm=0.85 if good else 0.30,
        predicted_aligned_error=6.0 if good else 20.0,
    )


def test_generator_and_scorer_pools_are_disjoint() -> None:
    """The constraint is structural, not a convention someone has to remember."""
    assert not GENERATOR_MODELS & SCORER_MODELS


def test_scoring_with_the_generating_model_is_refused() -> None:
    with pytest.raises(ValueError, match="self-agreement"):
        assert_no_shared_model("rfdiffusion", "rfdiffusion")


def test_a_valid_generator_scorer_pair_is_accepted() -> None:
    """A disjoint pair passes; the guard raises rather than returning a verdict."""
    assert_no_shared_model("rfdiffusion", "boltz2")


def test_unknown_models_are_rejected_rather_than_assumed_safe() -> None:
    with pytest.raises(ValueError):
        assert_no_shared_model("some-new-model", "boltz2")


def test_pass_rate_over_no_designs_is_an_error_not_zero() -> None:
    """A pass rate of 0.0 from an empty run would read as a real, terrible result."""
    with pytest.raises(ValueError):
        pass_rate([])


def test_separation_flags_an_uninformative_filter_stack() -> None:
    """If known binders do not beat scrambled controls, the filters measure nothing."""
    designs = [make_design("d1", good=True)]
    controls = [make_design("c1", good=True)]
    known = [make_design("k1", good=False)]
    result = separation(designs, controls, known)
    assert result["filters_are_informative"] == 0.0


def test_separation_confirms_an_informative_filter_stack() -> None:
    designs = [make_design("d1", good=True)]
    controls = [make_design("c1", good=False)]
    known = [make_design("k1", good=True)]
    result = separation(designs, controls, known)
    assert result["filters_are_informative"] == 1.0
    assert result["known_binders"] > result["controls"]
