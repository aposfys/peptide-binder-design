# peptide-binder-design
Designing peptide binders is easy; believing them is the hard part.

[![CI](https://github.com/aposfys/peptide-binder-design/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/peptide-binder-design/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

> **Status: skeleton.** No designs generated. The evaluation is stated in full before any design exists, so it cannot be adjusted to fit the results.

The standard pipeline generates a backbone with a diffusion model trained on the PDB, then validates it by checking that a structure predictor — trained on the same PDB, often sharing components — predicts the same complex. That is evidence of agreement, not of binding.

Without a lab this repo cannot prove a design binds. What it can do is measure how well the standard filter stack separates designs generated for the target, scrambled and length-matched decoys, and **known experimental binders of the same target**, held out. A filter stack that cannot rank known binders above scrambled sequences is not validating anything, whatever pass rate it reports on novel designs.

### Running it
```
make install && make data && make design && make analysis && make test
```
Generation backends are heavy optional extras and are GPU-bound; the filter and evaluation layers run on CPU and are what CI exercises.

### Layout
```
src/pepdesign/
  filters.py    the filter stack, its thresholds, and their null distributions
  cli.py        `python -m pepdesign.cli`
```
Planned: `targets.py` (target prep, hotspot definition, held-out splits), `generate.py`, `controls.py` (scrambled and length-matched decoys), `evaluate.py`.

### Design notes
[The circularity problem stated in full, and the traps the pipeline is built to avoid](docs/DESIGN.md)
