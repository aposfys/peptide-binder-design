# Results

Run 2026-09-29, over the peptide set retrieved 2026-09-01.

## What was run, and what was not

**Not run: the structure-based filter stack.** Self-consistency RMSD, interface pTM and
predicted aligned error all need a structure predictor on a GPU this repository has never
had. Those thresholds are declared in `filters.py` and remain unmeasured. RFdiffusion
generation is blocked for the same reason. No structure-based number appears below.

**Run instead:** ESM-2 pseudo-log-likelihood, length-normalised, computed exactly by masking
each position in turn. It is a real triage filter, and it is a *sequence-only* one, which
makes it the right instrument for this repository's question, because any separation it
achieves is necessarily about proteinness rather than about the target.

## Setup

| | |
| --- | --- |
| Real peptides | 190 chains, 8-30 residues, standard amino acids, from RCSB structures containing more than one protein entity, released on or before 2026-09-01 |
| Controls | 190 each of scrambled, composition-matched and length-matched |
| Filter | ESM-2 `esm2_t12_35M_UR50D`, exact pseudo-log-likelihood, length-normalised |
| Composition distance, real vs scrambled | 0.000000 |

Scrambling preserves composition exactly. The measured distance of 0.000000 is the check
that it does, so the scrambled arm differs from the real one in residue *order* and nothing
else.

## The result

The threshold column is set by the control distribution rather than by convention. 5% is
the target and 5.26% is what a 190-sequence sample can actually admit, so that is the
number quoted.

| Control family | Mean score | AUC | 95% CI | Cohen's *d* | Controls admitted | Real recall there |
| --- | ---: | ---: | --- | ---: | ---: | ---: |
| Real peptides | -2.823 | | | | | |
| Composition-matched | -3.056 | 0.684 | [0.632, 0.735] | 0.65 | 5.26% | 27.4% |
| Length-matched | -3.018 | 0.638 | [0.583, 0.696] | 0.55 | 5.26% | 26.8% |
| **Scrambled** | **-2.905** | **0.569** | **[0.512, 0.628]** | **0.18** | **5.26%** | **6.3%** |

**The same filter reports AUC 0.684 or 0.569 depending only on which null it is compared
against.** Effect size falls from *d* = 0.65 to *d* = 0.18, a factor of 3.6. At a matched
control admission rate of 5.26% the filter keeps 27.4% of real peptides against the
composition-matched null and 6.3% against the scrambled one, a factor of 4.3. The more
flattering number in every column is the one produced by the easier control.

The scrambled interval is `[0.512, 0.628]`. It excludes 0.5, so at 190 peptides there is a
small residual order signal and the honest statement is that it is small, not that it is
absent. What the comparison establishes is the size of the gap between nulls, not the
absence of a signal.

A threshold strict enough to admit 5.26% of scrambled sequences rejects 93.7% of the real
peptides. That is not a triage step.

## Sample size, and an earlier run that got this wrong

An earlier run scored only the first 120 of these 190 peptides and reported scrambled AUC
0.547 with a 95% CI of [0.471, 0.614]. That interval spans 0.5, and this repository
previously read it as the filter being indistinguishable from chance against a
composition-preserving null.

That reading was a sample-size artifact. The 70 unscored peptides were already in
`data/peptides.json` and needed no GPU. Scoring all 190 moves the scrambled interval off
0.5 without changing the argument about controls, which survives at both sample sizes and
is stronger at the larger one. Both runs reproduce from the committed code:
`pepdesign analysis --max-peptides 120 --force` returns the 120-peptide table exactly.

## Why this matters for the pass rates the field reports

Nothing above is a criticism of ESM. It is doing what a language model does. The
composition-matched and length-matched controls are drawn from background frequencies, real
peptides are not, and the model notices. That is a compositional signal, not a binding one,
and the scrambled arm removes it by construction.

Which is the whole argument of this repository. A design pipeline that quotes a pass rate
without saying what its controls were has not reported a result. It has reported a choice of
denominator.

## Reproducing this

- `results/scores.csv` holds every per-sequence score, 760 rows. Every AUC, bootstrap
  interval, effect size and recall figure above re-derives from it with
  `pepdesign.evaluate.separate`, with no model download and no network.
- `data/peptides.json` is the exact peptide set these numbers came from, and
  `pepdesign.targets.build` reads it in preference to the network. The RCSB query is now
  pinned to a release date and returns sorted identifiers, but the committed file predates
  that pin and a fresh retrieval returns a different set. The committed file, not the query,
  is what makes these numbers checkable.
- The bootstrap is seeded and the scorer runs in eval mode, so the whole pipeline is
  deterministic and runs on a laptop CPU.

## Limitations

- **The peptide set is a proxy.** Short protein chains in multi-protein structures are
  mostly peptides observed bound to a partner, but some are subunits of a complex and a few
  are crystallisation tags. It is noisy in a known direction and is not a curated binder set.
- **190 peptides is still small.** It is enough to separate the three nulls and to put the
  scrambled interval just off 0.5. It is not enough to bound the residual order signal
  tightly.
- **One filter, not the stack.** This says nothing about whether self-consistency RMSD or
  interface pTM separate binders from scrambles. Those are the interesting ones and they
  remain unrun. The null-distribution machinery they will need is built and tested.
