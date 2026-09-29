# peptide-binder-design
The same filter scores 0.68 or 0.57 depending only on which null you compare it against.

[![CI](https://github.com/aposfys/peptide-binder-design/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/peptide-binder-design/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

A triage filter for peptide binders, evaluated against three control families to
show that the reported performance is a property of the control, not the filter.

```
make install           # controls, evaluate and the tests
make install-score     # adds torch and transformers, which only the scorer needs
pepdesign evaluate     # print the separation table from the committed run
pepdesign controls     # write the control sets and their composition distance
pepdesign analysis     # rerun it all: peptides, controls, ESM-2 scoring, separation
make test              # 35 tests, no model and no network
```

`evaluate`, `controls`, `targets` and the tests run on a plain `make install`. Only
`analysis` needs the `score` extra, because only `analysis` loads ESM-2.

### One filter, three nulls

190 peptide chains observed bound in PDB complexes, scored by ESM-2
pseudo-log-likelihood. The threshold in the last two columns is set by the control
distribution, not by convention, and 5.26% is the closest a 190-sequence sample
gets to 5%:

| Control family | AUC | 95% CI | Cohen's *d* | Controls admitted | Real recall there |
| --- | ---: | --- | ---: | ---: | ---: |
| Composition-matched | 0.684 | [0.632, 0.735] | 0.65 | 5.26% | 27.4% |
| Length-matched | 0.638 | [0.583, 0.696] | 0.55 | 5.26% | 26.8% |
| **Scrambled** | **0.569** | **[0.512, 0.628]** | **0.18** | **5.26%** | **6.3%** |

The scrambled arm holds residue composition exactly fixed, measured distance
0.000000, so it differs from the real peptides only in residue *order*, which is
the only place binding information could live. Swap the naive control for that one
and the effect size falls by a factor of 3.6 and the recall by a factor of 4.3. A
threshold strict enough to admit 5.26% of scrambles rejects 93.7% of the real
peptides.

The scrambled interval excludes 0.5, so the residual order signal is small rather
than absent. The claim here is about the size of the gap between nulls.
`results/scores.csv` holds all 760 per-sequence scores, so every figure in the
table re-derives from committed files with no model and no network.

This is not a criticism of ESM. Composition- and length-matched controls are drawn
from background frequencies and real peptides are not, so the model is reading a
compositional signal rather than a binding one. **A pipeline that quotes a pass
rate without saying what its controls were has reported a choice of denominator,
not a result.**

### Scope

The structure-based stack (self-consistency RMSD, interface pTM, predicted
aligned error) needs a structure predictor on a GPU this repo has never had.
Those thresholds are declared and unmeasured, and `generate` names the GPU as the
reason it is unimplemented. **No structure-based number appears here**, and the
result above does not depend on one.

### More

- [Analysis](ANALYSIS.md), what was done and why, including what the peptide set is and is not
- [Results](results/RESULTS.md), the full table, how to recheck it, and the limitations
- [Design](docs/DESIGN.md), the circularity problem in full and the traps this avoids
- [Literature](docs/LITERATURE.md), what is and is not new here
- [References](docs/REFERENCES.md), every work the other four documents name
