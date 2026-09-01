"""Separation between real peptides and their controls, with the null that gives it meaning.

The reportable output is not a pass rate. It is how well a filter distinguishes real
peptides from each control family, and -- critically -- what threshold it would take to do
so, since a threshold quoted without its control distribution is not a threshold at all.
"""

from __future__ import annotations

import json
import random
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import mean, median, pstdev


@dataclass
class Separation:
    """How well a score separates real peptides from one control family."""

    control_kind: str
    n_real: int
    n_control: int
    mean_real: float
    mean_control: float
    auc: float
    auc_ci_low: float
    auc_ci_high: float
    #: Cohen's d. Reported alongside AUC because a large AUC on a tiny effect is possible
    #: and reads misleadingly.
    effect_size: float
    #: Score threshold that would keep 95% of controls out.
    threshold_at_5pct_control: float
    #: Fraction of real peptides that clear that threshold.
    real_recall_at_threshold: float


def roc_auc(positive: Sequence[float], negative: Sequence[float]) -> float:
    """AUC as the probability a random positive outscores a random negative.

    Computed by rank rather than by trapezoid, so ties are handled exactly -- half credit,
    which is what a tie is worth. Discretised scores produce many ties and a trapezoid
    implementation quietly rounds them in one direction.
    """
    if not positive or not negative:
        raise ValueError("need both populations to compute an AUC")
    combined = sorted([(value, 1) for value in positive] + [(value, 0) for value in negative])
    ranks: dict[int, float] = {}
    index = 0
    rank_sum_positive = 0.0
    while index < len(combined):
        stop = index
        while stop + 1 < len(combined) and combined[stop + 1][0] == combined[index][0]:
            stop += 1
        average_rank = (index + stop) / 2.0 + 1.0
        for position in range(index, stop + 1):
            if combined[position][1] == 1:
                rank_sum_positive += average_rank
        index = stop + 1
    del ranks
    n_pos, n_neg = len(positive), len(negative)
    return (rank_sum_positive - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


def bootstrap_auc(
    positive: Sequence[float],
    negative: Sequence[float],
    *,
    resamples: int = 1000,
    seed: int = 0,
) -> tuple[float, float]:
    """Percentile bootstrap interval for an AUC."""
    rng = random.Random(seed)
    values = []
    for _ in range(resamples):
        sampled_positive = [positive[rng.randrange(len(positive))] for _ in positive]
        sampled_negative = [negative[rng.randrange(len(negative))] for _ in negative]
        values.append(roc_auc(sampled_positive, sampled_negative))
    values.sort()
    return values[int(0.025 * len(values))], values[int(0.975 * len(values))]


def cohens_d(positive: Sequence[float], negative: Sequence[float]) -> float:
    """Standardised mean difference, pooled standard deviation."""
    if len(positive) < 2 or len(negative) < 2:
        return 0.0
    pooled = ((pstdev(positive) ** 2 + pstdev(negative) ** 2) / 2) ** 0.5
    if pooled == 0:
        return 0.0
    return (mean(positive) - mean(negative)) / pooled


def threshold_at_control_rate(negative: Sequence[float], rate: float = 0.05) -> float:
    """The score a filter would need to admit only ``rate`` of controls.

    This is the number a design paper should quote instead of a conventional cutoff: it is
    defined by the null rather than by tradition.
    """
    ordered = sorted(negative, reverse=True)
    index = max(0, min(len(ordered) - 1, int(rate * len(ordered))))
    return ordered[index]


def separate(
    real: Sequence[float], control: Sequence[float], kind: str, *, seed: int = 0
) -> Separation:
    """Measure one real-versus-control comparison."""
    threshold = threshold_at_control_rate(control)
    low, high = bootstrap_auc(real, control, seed=seed)
    return Separation(
        control_kind=kind,
        n_real=len(real),
        n_control=len(control),
        mean_real=round(mean(real), 4),
        mean_control=round(mean(control), 4),
        auc=round(roc_auc(real, control), 4),
        auc_ci_low=round(low, 4),
        auc_ci_high=round(high, 4),
        effect_size=round(cohens_d(real, control), 4),
        threshold_at_5pct_control=round(threshold, 4),
        real_recall_at_threshold=round(
            sum(1 for value in real if value >= threshold) / len(real), 4
        ),
    )


def build_findings(scored, *, notes: dict | None = None) -> dict:
    """Assemble findings from a list of :class:`~pepdesign.score.Scored` items."""
    by_kind: dict[str, list[float]] = {}
    for item in scored:
        by_kind.setdefault(item.kind, []).append(item.pseudo_log_likelihood)

    real = by_kind.get("real", [])
    if not real:
        raise ValueError("no real peptides were scored")

    separations = [
        asdict(separate(real, values, kind))
        for kind, values in sorted(by_kind.items())
        if kind != "real"
    ]
    return {
        "filter": "ESM-2 pseudo-log-likelihood, length-normalised",
        "not_run": (
            "self-consistency RMSD, interface pTM and predicted aligned error: all need a "
            "structure predictor on a GPU. No structure-based numbers are reported."
        ),
        "populations": {
            kind: {
                "n": len(values),
                "mean": round(mean(values), 4),
                "median": round(median(values), 4),
            }
            for kind, values in sorted(by_kind.items())
        },
        "separations": separations,
        "notes": notes or {},
    }


def write(findings: dict, out_path: Path) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(findings, indent=1))
    return out_path
