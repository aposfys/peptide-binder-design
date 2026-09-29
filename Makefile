.PHONY: install install-score data controls analysis evaluate design test clean clean-data all

PYTHON ?= python3
TARGET ?= 1KMV

all: analysis

## Install the package plus dev tooling. Enough for controls, evaluate and the tests.
install:
	$(PYTHON) -m pip install -e ".[dev]"

## Add torch and transformers, which only the ESM-2 scorer needs.
install-score:
	$(PYTHON) -m pip install -e ".[dev,score]"

## Peptides, controls, ESM-2 scoring and the separation table. Needs the score extra.
analysis:
	$(PYTHON) -m pepdesign.cli analysis

## The control sets and their composition distance. CPU only, no model.
controls:
	$(PYTHON) -m pepdesign.cli controls

## Print the separation table from results/findings.json.
evaluate:
	$(PYTHON) -m pepdesign.cli evaluate

## The peptide set, retrieved from RCSB and cached in data/peptides.json
data:
	$(PYTHON) -m pepdesign.cli targets

## Backbone generation and sequence design. GPU, unimplemented, exits with that message.
design:
	$(PYTHON) -m pepdesign.cli generate --target $(TARGET)

test:
	$(PYTHON) -m pytest -q

## Only the generated control dump. results/findings.json, results/scores.csv and
## results/RESULTS.md are committed artifacts and are left alone.
clean:
	rm -f results/controls.json
	find . -name __pycache__ -type d -exec rm -rf {} +

## Also delete cached structures and generated designs
clean-data: clean
	rm -f data/*.pdb data/*.cif data/*.fasta
