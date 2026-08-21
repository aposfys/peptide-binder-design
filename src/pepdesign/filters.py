"""The filter stack, its thresholds, and the null distributions that give them meaning.

Every threshold here is conventional in the literature and arbitrary in isolation. A design
passing at ``rmsd < 2.0`` means nothing until scrambled controls have been pushed through
the identical stack, so :func:`separation` -- not :func:`passes` -- is the reportable output
of this module.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from pepdesign import GENERATOR_MODELS, SCORER_MODELS

#: Conventional thresholds. Present so they can be argued with, not because they are right.
DEFAULT_THRESHOLDS = {
    "self_consistency_rmsd": 2.0,  # angstrom, lower is better
    "interface_ptm": 0.75,  # higher is better
    "predicted_aligned_error": 10.0,  # angstrom, lower is better
}


@dataclass(frozen=True)
class Design:
    """One designed binder and the metrics computed for it."""

    design_id: str
    sequence: str
    self_consistency_rmsd: float
    interface_ptm: float
    predicted_aligned_error: float


def assert_no_shared_model(generator: str, scorer: str) -> None:
    """Refuse a configuration where the proposing and scoring models are the same.

    This is the central guard of the repository. It raises rather than warns, because a
    warning in a long pipeline run is a warning nobody reads.
    """
    # Equality is checked first, deliberately. The two pools are disjoint today, so an
    # identical pair would otherwise be reported as an unknown scorer -- a confusing
    # message for the one misconfiguration this function exists to catch.
    if generator == scorer:
        raise ValueError(
            f"generator and scorer are both {generator!r}: this measures self-agreement, "
            "not binding"
        )
    if generator not in GENERATOR_MODELS:
        raise ValueError(f"unknown generator {generator!r}")
    if scorer not in SCORER_MODELS:
        raise ValueError(f"unknown scorer {scorer!r}")


def passes(design: Design, thresholds: dict[str, float] | None = None) -> bool:
    """Whether one design clears every filter. Meaningless without a control distribution."""
    limits = dict(DEFAULT_THRESHOLDS if thresholds is None else thresholds)
    return (
        design.self_consistency_rmsd < limits["self_consistency_rmsd"]
        and design.interface_ptm > limits["interface_ptm"]
        and design.predicted_aligned_error < limits["predicted_aligned_error"]
    )


def pass_rate(designs: Sequence[Design], thresholds: dict[str, float] | None = None) -> float:
    """Fraction of designs clearing the stack."""
    if not designs:
        raise ValueError("cannot compute a pass rate over zero designs")
    return sum(1 for design in designs if passes(design, thresholds)) / len(designs)


def separation(
    designs: Sequence[Design],
    controls: Sequence[Design],
    known_binders: Sequence[Design],
    thresholds: dict[str, float] | None = None,
) -> dict[str, float]:
    """Pass rates for the three populations that make the filter stack interpretable.

    If ``known_binders`` does not clear the stack at a higher rate than ``controls``, the
    filters are not measuring binding and the design pass rate should not be reported as
    if they were. That comparison is returned, not left to the reader.
    """
    rates = {
        "designs": pass_rate(designs, thresholds),
        "controls": pass_rate(controls, thresholds),
        "known_binders": pass_rate(known_binders, thresholds),
    }
    rates["filters_are_informative"] = float(rates["known_binders"] > rates["controls"])
    return rates
