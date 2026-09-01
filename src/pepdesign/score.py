"""The one filter this repository can actually compute, and what it is worth.

**The structure-based stack in :mod:`pepdesign.filters` has not been run.** Self-consistency
RMSD, interface pTM and predicted aligned error all require a structure predictor on a GPU
that this project has never had. Those thresholds remain declared and unmeasured.

What can be computed on a laptop is a protein language model's pseudo-log-likelihood: how
unsurprising a sequence is to ESM-2. It is a real filter -- "designability" scores of this
shape are used to triage designs before anything expensive runs -- and it is a good
instrument for the repository's actual question, because it is *sequence-only*. If a
sequence-only filter separates known binders from their own scrambles, that separation
cannot be about the target, only about whether the sequence looks like a protein.

Which is the point. A filter that ranks real peptides above scrambles is measuring
proteinness. Reporting a pass rate from it as though it were evidence of binding is the
error, and the null distributions here are what make the error visible.
"""

from __future__ import annotations

import os
from collections.abc import Sequence
from dataclasses import dataclass

#: Small enough to run on CPU, large enough to be a real language model. Pinned, because
#: "an ESM score" is not a reproducible statement.
DEFAULT_CHECKPOINT = "facebook/esm2_t12_35M_UR50D"


@dataclass(frozen=True)
class Scored:
    """One sequence and its computed scores."""

    identifier: str
    sequence: str
    kind: str
    pseudo_log_likelihood: float
    length: int


class EsmScorer:
    """Pseudo-log-likelihood under a frozen ESM-2 checkpoint.

    Scores are length-normalised. Without that, the filter would rank short sequences above
    long ones for arithmetic reasons and every comparison between differently-sized
    populations would be meaningless -- and length is exactly what the length-matched
    controls hold constant, so the confound would be invisible in that arm and only in that
    arm.
    """

    def __init__(self, checkpoint: str = DEFAULT_CHECKPOINT, *, batch_size: int = 16):
        self.checkpoint = checkpoint
        self.batch_size = batch_size
        self._tokenizer = None
        self._model = None

    def _load(self) -> None:
        if self._model is not None:
            return
        os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")
        import warnings

        from transformers import AutoModelForMaskedLM, AutoTokenizer

        warnings.filterwarnings("ignore")
        self._tokenizer = AutoTokenizer.from_pretrained(self.checkpoint)
        self._model = AutoModelForMaskedLM.from_pretrained(self.checkpoint)
        self._model.eval()

    def score(self, sequence: str) -> float:
        """Length-normalised pseudo-log-likelihood of one sequence.

        Computed the exact way -- masking each position in turn -- rather than from a
        single unmasked forward pass. The cheap approximation lets the model see the
        residue it is predicting, which inflates the score for every sequence and, worse,
        inflates it unevenly.
        """
        import torch

        self._load()
        encoded = self._tokenizer(sequence, return_tensors="pt")
        input_ids = encoded["input_ids"]
        attention = encoded["attention_mask"]
        # Skip the leading <cls> and trailing <eos>.
        positions = list(range(1, input_ids.shape[1] - 1))
        if not positions:
            raise ValueError(f"nothing to score in {sequence!r}")

        total = 0.0
        with torch.no_grad():
            for start in range(0, len(positions), self.batch_size):
                chunk = positions[start : start + self.batch_size]
                batch = input_ids.repeat(len(chunk), 1)
                for row, position in enumerate(chunk):
                    batch[row, position] = self._tokenizer.mask_token_id
                logits = self._model(
                    input_ids=batch, attention_mask=attention.repeat(len(chunk), 1)
                ).logits
                log_probs = torch.log_softmax(logits, dim=-1)
                for row, position in enumerate(chunk):
                    true_token = input_ids[0, position]
                    total += float(log_probs[row, position, true_token])
        return total / len(positions)

    def score_many(self, items: Sequence[tuple[str, str, str]]) -> list[Scored]:
        """Score `(identifier, sequence, kind)` triples."""
        results: list[Scored] = []
        for identifier, sequence, kind in items:
            results.append(
                Scored(
                    identifier=identifier,
                    sequence=sequence,
                    kind=kind,
                    pseudo_log_likelihood=self.score(sequence),
                    length=len(sequence),
                )
            )
        return results
