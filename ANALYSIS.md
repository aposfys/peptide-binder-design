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

120 peptide chains observed bound in PDB complexes, against 120 of each control family:

| Control | AUC | 95% CI | Cohen's *d* | Recall at 5% control FPR |
| --- | ---: | --- | ---: | ---: |
| Composition-matched | 0.663 | [0.592, 0.730] | 0.60 | 28% |
| Length-matched | 0.632 | [0.554, 0.704] | 0.54 | 25% |
| **Scrambled** | **0.547** | **[0.471, 0.614]** | **0.14** | **7%** |

**The same filter reports 0.663 or 0.547 depending only on which null it is compared
against, and the scrambled interval spans 0.5.** A threshold strict enough to reject 95% of
scrambles rejects 93% of the real peptides too.

None of this is a criticism of ESM. The naive controls are drawn from background
frequencies and real peptides are not; the model notices, and that is a compositional
signal. The scrambled arm removes it by construction, and what is left is not
distinguishable from chance at this sample size.

## What the peptide set actually is

Short protein chains (8–30 residues, standard amino acids) from RCSB structures containing
more than one protein entity. That is peptides observed bound to a protein partner — the
closest thing to a validated binder obtainable without a wet lab. It is **not** curated:
some short chains are subunits rather than ligands, and a few are crystallisation tags. The
population is noisy in a known direction, and that is stated here rather than after it has
been forgotten.

## What is not established

- Anything about self-consistency RMSD, interface pTM or interface confidence. Those are the
  interesting filters and they remain unrun.
- A bound on the effect. 120 peptides shows the CI crosses 0.5; it does not pin the size.
- Anything about designed sequences. Only observed ones were scored.

## What would change the conclusion

A GPU. The same three control families through the structure-based stack is the experiment
this repository was designed for, and the null machinery it needs is now built and tested.
The prediction implied by the result above is uncomfortable: a stack validated only against
length-matched decoys may be reporting the same kind of number.
