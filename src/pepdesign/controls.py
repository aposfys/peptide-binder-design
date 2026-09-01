"""Scrambled and length-matched decoy peptides.

Every threshold in this repository needs a null, and these are it. Three control families,
each removing something different:

``scrambled``
    The same residues in a different order. Composition is preserved exactly, so a filter
    that separates real peptides from these is responding to *order* -- which is the only
    thing that could carry binding information.
``length_matched``
    Random sequences at the background amino-acid frequency, matched only on length. The
    weakest control, and the one most filters pass trivially.
``composition_matched``
    Random sequences drawn from the pooled composition of the real set. Between the other
    two: shared bulk composition, no per-peptide correspondence.

A filter that cannot rank real peptides above *scrambled* ones is not measuring anything
about sequence, whatever it reports on novel designs.
"""

from __future__ import annotations

import random
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass

#: Background amino-acid frequencies (UniProt-wide, rounded). Used for length-matched
#: controls, which are deliberately the naive null.
BACKGROUND_FREQUENCIES = {
    "A": 0.0825,
    "R": 0.0553,
    "N": 0.0406,
    "D": 0.0545,
    "C": 0.0137,
    "Q": 0.0393,
    "E": 0.0675,
    "G": 0.0707,
    "H": 0.0227,
    "I": 0.0596,
    "L": 0.0966,
    "K": 0.0584,
    "M": 0.0242,
    "F": 0.0386,
    "P": 0.0470,
    "S": 0.0656,
    "T": 0.0534,
    "W": 0.0108,
    "Y": 0.0292,
    "V": 0.0687,
}

CONTROL_KINDS = ("scrambled", "length_matched", "composition_matched")


@dataclass(frozen=True)
class Control:
    """One control peptide and what it was derived from."""

    control_id: str
    sequence: str
    kind: str
    #: The real peptide this control corresponds to, where there is one.
    derived_from: str | None = None


def scramble(sequence: str, *, seed: int = 0) -> str:
    """Shuffle a sequence, preserving composition exactly."""
    residues = list(sequence)
    random.Random(seed).shuffle(residues)
    return "".join(residues)


def _draw(frequencies: dict[str, float], length: int, rng: random.Random) -> str:
    letters = list(frequencies)
    weights = [frequencies[letter] for letter in letters]
    return "".join(rng.choices(letters, weights=weights, k=length))


def make_controls(
    peptides: Sequence[tuple[str, str]],
    *,
    kinds: Sequence[str] = CONTROL_KINDS,
    per_peptide: int = 1,
    seed: int = 0,
) -> list[Control]:
    """Build every requested control family from a set of `(id, sequence)` peptides.

    Scrambling uses a per-peptide seed so a control set is reproducible without every
    peptide being shuffled the same way -- one shared shuffle would correlate the controls
    with each other and shrink the effective size of the null.
    """
    unknown = set(kinds) - set(CONTROL_KINDS)
    if unknown:
        raise ValueError(f"unknown control kind(s) {sorted(unknown)}")

    rng = random.Random(seed)
    pooled: Counter[str] = Counter()
    for _, sequence in peptides:
        pooled.update(sequence)
    total = sum(pooled.values())
    pooled_frequencies = (
        {letter: count / total for letter, count in pooled.items()} if total else {}
    )

    controls: list[Control] = []
    for index, (peptide_id, sequence) in enumerate(peptides):
        for replicate in range(per_peptide):
            offset = index * per_peptide + replicate
            if "scrambled" in kinds:
                controls.append(
                    Control(
                        control_id=f"scr_{peptide_id}_{replicate}",
                        sequence=scramble(sequence, seed=seed + offset),
                        kind="scrambled",
                        derived_from=peptide_id,
                    )
                )
            if "length_matched" in kinds:
                controls.append(
                    Control(
                        control_id=f"len_{peptide_id}_{replicate}",
                        sequence=_draw(BACKGROUND_FREQUENCIES, len(sequence), rng),
                        kind="length_matched",
                        derived_from=peptide_id,
                    )
                )
            if "composition_matched" in kinds and pooled_frequencies:
                controls.append(
                    Control(
                        control_id=f"cmp_{peptide_id}_{replicate}",
                        sequence=_draw(pooled_frequencies, len(sequence), rng),
                        kind="composition_matched",
                        derived_from=peptide_id,
                    )
                )
    return controls


def composition_distance(a: Sequence[str], b: Sequence[str]) -> float:
    """Total-variation distance between the pooled compositions of two peptide sets.

    Reported so "composition-matched" is a measured claim and not a label. Scrambled
    controls must score exactly 0.0 against their sources; anything else is a bug.
    """

    def profile(sequences: Sequence[str]) -> dict[str, float]:
        counts: Counter[str] = Counter()
        for sequence in sequences:
            counts.update(sequence)
        total = sum(counts.values())
        return {letter: count / total for letter, count in counts.items()} if total else {}

    left, right = profile(a), profile(b)
    letters = set(left) | set(right)
    return 0.5 * sum(abs(left.get(letter, 0.0) - right.get(letter, 0.0)) for letter in letters)
