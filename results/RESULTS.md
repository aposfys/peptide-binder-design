# Results

Run 2026-09-01.

## What was run, and what was not

**Not run: the structure-based filter stack.** Self-consistency RMSD, interface pTM and
predicted aligned error all need a structure predictor on a GPU this repository has never
had. Those thresholds are declared in `filters.py` and remain unmeasured. No structure-based
number appears below.

**Run instead:** ESM-2 pseudo-log-likelihood, length-normalised, computed exactly by masking
each position in turn. It is a real triage filter, and it is a *sequence-only* one — which
makes it the right instrument for this repository's question, because any separation it
achieves is necessarily about proteinness rather than about the target.

## Setup

| | |
| --- | --- |
| Real peptides | 120 chains, 8–30 residues, standard amino acids, from RCSB structures containing more than one protein entity |
| Controls | 120 each of scrambled, composition-matched and length-matched |
| Filter | ESM-2 `esm2_t12_35M_UR50D`, exact pseudo-log-likelihood, length-normalised |
| Composition distance, real vs scrambled | 0.000000 |

Scrambling preserves composition exactly — the measured distance of 0.000000 is the check
that it does — so the scrambled arm differs from the real one in residue *order* and nothing
else.

## The result

| Control family | Mean score | AUC | 95% CI | Cohen's *d* | Real recall at 5% control FPR |
| --- | ---: | ---: | --- | ---: | ---: |
| Real peptides | −2.804 | — | — | — | — |
| Composition-matched | −3.051 | 0.663 | [0.592, 0.730] | 0.60 | 28% |
| Length-matched | −3.026 | 0.632 | [0.554, 0.704] | 0.54 | 25% |
| **Scrambled** | **−2.876** | **0.547** | **[0.471, 0.614]** | **0.14** | **7%** |

**Against the naive controls the filter looks like it works. Against the proper one it does
not.** The confidence interval for the scrambled comparison spans 0.5, so at this sample
size the filter cannot be shown to distinguish a real peptide from the same residues in a
different order. Effect size falls from *d* = 0.60 to *d* = 0.14.

The recall column makes the practical version of the point. Set a threshold strict enough to
reject 95% of scrambled sequences and it also rejects 93% of the real peptides. That filter
is not a triage step; it is a coin flip with extra compute.

## Why this matters for the pass rates the field reports

Nothing above is a criticism of ESM. It is doing what a language model does — the
composition-matched and length-matched controls are drawn from background frequencies, real
peptides are not, and the model notices. That is a compositional signal, not a binding one,
and the scrambled arm removes it by construction.

Which is the whole argument of this repository. **The same filter reports AUC 0.66 or 0.55
depending only on which null you compare it against**, and the more flattering number is the
one produced by the easier control. A design pipeline that quotes a pass rate without saying
what its controls were has not reported a result; it has reported a choice of denominator.

## Limitations

- **The peptide set is a proxy.** Short protein chains in multi-protein structures are
  mostly peptides observed bound to a partner, but some are subunits of a complex and a few
  are crystallisation tags. It is noisy in a known direction and is not a curated binder set.
- **120 peptides is small.** It is enough to show the CI crosses 0.5 and not enough to
  bound the effect tightly.
- **One filter, not the stack.** This says nothing about whether self-consistency RMSD or
  interface pTM separate binders from scrambles. Those are the interesting ones and they
  remain unrun — but the null-distribution machinery they will need is built and tested.
