# Analysis

What was built, why it was built that way, and why a blocked experiment produced a better
result than the one that was blocked.

## What could not be run, and what replaced it

The structure-based filter stack — self-consistency RMSD, interface pTM, predicted aligned
error — needs a structure predictor on a GPU this project has never had. Those thresholds
are declared in `filters.py` and remain unmeasured. RFdiffusion generation is blocked for
the same reason.

Rather than report nothing, the repository measures the thing its filters were *for*:
whether a filter's apparent performance is a property of the filter or of the control set
it was compared against. That question needs no GPU, and it turns out to be the sharper one.

## Design decisions, and the reasoning

**A sequence-only filter is the right instrument, not a compromise.** ESM-2
pseudo-log-likelihood is a real triage score, and it has no idea what the target is. So any
separation it achieves is necessarily about *proteinness*, never about binding. That makes
it a clean probe: if it separates real peptides from controls, the separation cannot be
evidence of binding, and whatever the number is, it is a property of the controls.

**Three control families, removing different things.** Length-matched controls share only
length. Composition-matched controls share the pooled composition. **Scrambled controls
share composition exactly, per peptide, and differ only in residue order** — which is the
only place binding information could live. Measured composition distance between real and
scrambled is `0.000000`, and a test asserts it, so "composition-preserving" is a checked
claim rather than a label.

**Per-peptide scramble seeds.** One shared shuffle would correlate the controls with each
other and shrink the effective size of the null.

**Exact pseudo-log-likelihood, masking each position in turn.** The cheap single-pass
approximation lets the model see the residue it is predicting, which inflates every score
and — worse — inflates them unevenly.

**Length normalisation.** Without it the filter ranks short sequences above long ones for
arithmetic reasons. Length is exactly what the length-matched arm holds constant, so the
confound would have been invisible in that arm and only in that arm.

**Rank-based AUC, not trapezoid.** Discretised scores produce many ties, and a trapezoid
implementation rounds them in one direction. Ties get half credit, which is what a tie is
worth.

## What was measured

190 peptide chains observed bound in PDB complexes, against 190 of each control family:

| Control | AUC | 95% CI | Cohen's *d* | Controls admitted | Real recall there |
| --- | ---: | --- | ---: | ---: | ---: |
| Composition-matched | 0.684 | [0.632, 0.735] | 0.65 | 5.26% | 27.4% |
| Length-matched | 0.638 | [0.583, 0.696] | 0.55 | 5.26% | 26.8% |
| **Scrambled** | **0.569** | **[0.512, 0.628]** | **0.18** | **5.26%** | **6.3%** |

**The same filter reports 0.684 or 0.569 depending only on which null it is compared
against.** Effect size falls by a factor of 3.6 and recall at a matched control admission
rate by a factor of 4.3. A threshold strict enough to admit 5.26% of scrambles rejects
93.7% of the real peptides.

The scrambled interval is [0.512, 0.628], which excludes 0.5. So there is a small residual
order signal at this sample size, and the claim is about the size of the gap between nulls
rather than about the absence of a signal. An earlier run over the first 120 of these
peptides put the scrambled interval at [0.471, 0.614], spanning 0.5, and that was a
sample-size artifact rather than a null result. `results/RESULTS.md` has both runs.

None of this is a criticism of ESM. The naive controls are drawn from background
frequencies and real peptides are not, the model notices, and that is a compositional
signal. The scrambled arm removes it by construction, and what survives is a fraction of
what the naive controls suggested.

## What the peptide set actually is

Short protein chains (8-30 residues, standard amino acids) from RCSB structures containing
more than one protein entity, released on or before 2026-09-01. That is peptides observed bound to a protein partner — the
closest thing to a validated binder obtainable without a wet lab. It is **not** curated:
some short chains are subunits rather than ligands, and a few are crystallisation tags. The
population is noisy in a known direction, and that is stated here rather than after it has
been forgotten.

## What is not established

- Anything about self-consistency RMSD, interface pTM or interface confidence. Those are the
  interesting filters and they remain unrun.
- A tight bound on the residual order signal. 190 peptides puts the scrambled interval just
  off 0.5. It does not pin the size of what is left.
- Anything about designed sequences. Only observed ones were scored.

## What would change the conclusion

Two things, and only one of them needs hardware.

**More peptides, which is free.** The first version of this analysis scored 120 of the 190
peptides already committed to `data/peptides.json`, and the extra 70 moved the scrambled
interval off 0.5. Widening the length window or raising the search limit is CPU work and
would bound the residual signal better than anything else available here.

**A GPU, for the experiment this was designed for.** The same three control families through
the structure-based stack. The null machinery it needs is now built and tested. The
prediction implied by the result above is uncomfortable: a stack validated only against
length-matched decoys may be reporting the same kind of number.
