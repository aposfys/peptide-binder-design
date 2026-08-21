# peptide-binder-design — designing binders is easy, believing them is the hard part

[![CI](https://github.com/aposfys/peptide-binder-design/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/peptide-binder-design/actions/workflows/ci.yml)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Ruff](https://img.shields.io/badge/lint-ruff-261230)](https://docs.astral.sh/ruff/)

> **Status: skeleton.** No designs generated. The point of the repo is stated below before
> any design exists, so the evaluation cannot be adjusted to fit the results.

De novo binder design is the fastest-moving area in computational structural biology:
RFdiffusion now produces antibodies whose designed CDR loops match cryo-EM structures at
atomic accuracy, RFpeptides extends the same machinery to macrocyclic peptides, and a 2026
*Nature* review argues the binder-design problem is close to solved.

This repository designs peptide binders against a target with known binders — and then
spends most of its effort on the question the field finds least convenient:

**How much of a reported in-silico success rate survives when you remove the circularity?**

| | |
| --- | --- |
| **Backbones** | RFdiffusion / RFpeptides, target hotspot-conditioned |
| **Sequences** | ProteinMPNN |
| **Filters** | self-consistency RMSD, predicted alignment error, interface pTM |
| **Control** | scrambled and length-matched designs run through the identical filter stack |
| **Readout** | filter pass rate for designs vs controls vs *known experimental binders* |

## The circularity problem, stated plainly

The standard pipeline generates a backbone with a diffusion model trained on the PDB, then
validates it by checking that a structure predictor — trained on the same PDB, often sharing
components — predicts the same complex. Agreement between two models with shared training
data and shared inductive biases is not evidence of binding. It is evidence of agreement.

Published computational success rates and wet-lab hit rates are different numbers, and the
gap is the subject here. Without a lab, this repo cannot prove a design binds. What it *can*
do, and what is worth doing, is measure how well the standard filter stack separates:

1. designs generated for the target,
2. scrambled or length-matched decoy designs, and
3. **known experimental binders of the same target**, held out.

A filter stack that cannot rank known binders above scrambled sequences is not validating
anything, whatever pass rate it reports on novel designs. That result — if it appears — is
the finding.

## Traps this pipeline is built to avoid

- **Generator and scorer must not share weights.** If the model that proposes a backbone
  also scores it, the pipeline measures self-agreement. The scoring model is pinned before
  generation and the pin is asserted in the test suite.
- **Self-consistency RMSD has no null.** A threshold of 2 Å means nothing until you know
  what scrambled sequences score. Every threshold in this repo has a control distribution.
- **Hotspot leakage.** Conditioning on hotspot residues taken from a known complex, then
  evaluating against that same complex, is memorisation. Held-out targets never contribute
  hotspots.
- **Peptides are not small proteins.** Short, flexible binders are exactly where structure
  predictors are least reliable, and confidence metrics calibrated on globular domains do
  not transfer. Confidence is recalibrated on peptide-length controls.

## Layout

```
src/pepdesign/
  targets.py    target preparation, hotspot definition, held-out splits
  generate.py   backbone generation and sequence design
  filters.py    the filter stack, its thresholds, and their null distributions
  controls.py   scrambled and length-matched decoy designs
  evaluate.py   separation of designs / controls / known binders
  cli.py        `python -m pepdesign.cli`
```

```bash
make install && make data && make design && make analysis && make test
```

Generation backends are heavy optional extras and are GPU-bound; the filter and evaluation
layers run on CPU and are what CI exercises.
