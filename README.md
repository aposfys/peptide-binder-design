# peptide-binder-design
The same filter scores 0.66 or 0.55 depending only on which null you compare it against.

[![CI](https://github.com/aposfys/peptide-binder-design/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/peptide-binder-design/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

A triage filter for peptide binders, evaluated against three control families to
show that the reported performance is a property of the control, not the filter.

```
make install
pepdesign analysis     # peptides, controls, scoring, separation
pepdesign controls     # write the control sets and their composition distance
pepdesign evaluate     # print the separation table from an existing run
make test              # 26 tests, no model and no network
```

### One filter, three nulls

120 peptide chains observed bound in PDB complexes, scored by ESM-2
pseudo-log-likelihood:

| Control family | AUC | 95% CI | Cohen's *d* | Real recall at 5% control FPR |
| --- | ---: | --- | ---: | ---: |
| Composition-matched | 0.663 | [0.592, 0.730] | 0.60 | 28% |
| Length-matched | 0.632 | [0.554, 0.704] | 0.54 | 25% |
| **Scrambled** | **0.547** | **[0.471, 0.614]** | **0.14** | **7%** |

Against naive controls the filter looks like it works. Against
composition-preserving scrambles **the confidence interval spans 0.5**. The
scrambled arm holds residue composition exactly fixed — measured distance
0.000000 — so it differs from the real peptides only in residue *order*, which is
the only place binding information could live. Set a threshold strict enough to
reject 95% of scrambles and it rejects 93% of the real peptides too.

This is not a criticism of ESM. Composition- and length-matched controls are drawn
from background frequencies and real peptides are not, so the model is reading a
compositional signal rather than a binding one. **A pipeline that quotes a pass
rate without saying what its controls were has reported a choice of denominator,
not a result.**

### Scope

The structure-based stack — self-consistency RMSD, interface pTM, predicted
aligned error — needs a structure predictor on a GPU this repo has never had.
Those thresholds are declared and unmeasured; `targets` and `generate` name the
GPU as the reason they are unimplemented. **No structure-based number appears
here**, and the result above does not depend on one.

### Where this sits in the literature

The general principle — that how you construct negatives drives apparent performance — is
not new, and this repository is not the first to say it. It is best established in TCR–pMHC
specificity prediction, where shuffling within a dataset is known to introduce leakage that
models exploit instead of learning recognition, and where tools such as STAPLER exist
specifically to mitigate it. Decoy selection has its own literature in virtual screening.

What I could not find published is this test applied to **protein-language-model scoring of
peptide binders**. The peptide-binder design literature — PepMLM and target-conditioned
masked language modelling, DiffPepBuilder, contrastive target-conditioned design — reports
performance against controls that differ from the positives in composition and length.
Composition-preserving scrambles are not standard practice there, and the result above is
what happens when you use them: a filter that reads as working at AUC 0.663 has a
confidence interval spanning 0.5 once the control keeps its residue census fixed.

There is also a known confound in the same direction worth naming: protein language models
transfer unevenly to peptide-length sequences, so a score calibrated on protein-length
input is already on uncertain ground before the control question arises.

### More

- [Analysis](ANALYSIS.md) — what was done and why, including what the peptide set is and isn't
- [Results](results/RESULTS.md) — full results and limitations
- [Design](docs/DESIGN.md) — the circularity problem in full, and the traps this avoids
