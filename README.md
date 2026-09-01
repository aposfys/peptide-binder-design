# peptide-binder-design
Designing binders is easy; believing them is the hard part.

[![CI](https://github.com/aposfys/peptide-binder-design/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/peptide-binder-design/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

> **The structure-based filter stack has not been run.** Self-consistency RMSD, interface pTM and predicted aligned error all need a structure predictor on a GPU this repo has never had. Those thresholds are declared and unmeasured, and no structure-based number appears here.

What has been run is the part that does not need a GPU, and it turns out to demonstrate the repository's thesis on its own.

### The same filter, two nulls, two conclusions

120 peptide chains observed bound in PDB complexes, scored by ESM-2 pseudo-log-likelihood — a real triage filter — against three control families:

| Control family | AUC | 95% CI | Cohen's *d* | Real recall at 5% control FPR |
| --- | ---: | --- | ---: | ---: |
| Composition-matched | 0.663 | [0.592, 0.730] | 0.60 | 28% |
| Length-matched | 0.632 | [0.554, 0.704] | 0.54 | 25% |
| **Scrambled** | **0.547** | **[0.471, 0.614]** | **0.14** | **7%** |

**Against naive controls the filter looks like it works. Against composition-preserving scrambles the confidence interval spans 0.5.** The scrambled arm holds residue composition exactly fixed — measured distance 0.000000 — so it differs from the real peptides only in the order of the residues, which is the only place binding information could live.

Set a threshold strict enough to reject 95% of scrambles and it rejects 93% of the real peptides too.

None of this is a criticism of ESM. Composition-matched and length-matched controls are drawn from background frequencies, real peptides are not, and the model notices — that is a compositional signal, not a binding one. **The same filter reports 0.66 or 0.55 depending only on which null it is compared against, and the flattering number comes from the easier control.** A pipeline that quotes a pass rate without saying what its controls were has reported a choice of denominator, not a result.

### Running it

```
make install
python3 -m pepdesign.cli analysis      # peptides, controls, scoring, separation
make test
```

The generation and structure-scoring subcommands remain unimplemented and say so.

### What the peptide set is

Short protein chains (8–30 residues, standard amino acids only) from RCSB structures containing more than one protein entity — peptides observed bound to a protein partner, which is the closest thing to a validated binder obtainable without a wet lab. It is **not** a curated set: some short chains are subunits rather than ligands, and a few are crystallisation tags. Noisy in a known direction, and said here rather than after it has been forgotten.

### Layout

```
src/pepdesign/
  targets.py    RCSB peptide retrieval, filtering, deduplication
  controls.py   scrambled, composition-matched and length-matched nulls
  score.py      exact ESM-2 pseudo-log-likelihood, length-normalised
  filters.py    the structure-based stack and its no-shared-model guard (unrun)
  evaluate.py   AUC with bootstrap CI, effect size, threshold from the null
  cli.py        `python -m pepdesign.cli`
```

20 tests, none needing a model or a network.

### More

- [Full results and limitations](results/RESULTS.md)
- [The circularity problem stated in full, and the traps this pipeline avoids](docs/DESIGN.md)
