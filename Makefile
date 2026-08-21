.PHONY: install data design controls analysis test clean clean-data all

PYTHON ?= python3
TARGET ?= 1KMV

all: analysis

## Install the package plus dev tooling. Generation backends are GPU-bound extras:
##   pip install -e ".[design,folding]"
install:
	$(PYTHON) -m pip install -e ".[dev]"

## Target structures, hotspot definitions, and the held-out known-binder set
data:
	$(PYTHON) -m pepdesign.cli targets --target $(TARGET)

## Backbone generation and sequence design (GPU)
design: data
	$(PYTHON) -m pepdesign.cli generate --target $(TARGET)

## Scrambled and length-matched controls, through the identical filter stack
controls: design
	$(PYTHON) -m pepdesign.cli controls

## Separation of designs / controls / known binders -- the actual result
analysis: controls
	$(PYTHON) -m pepdesign.cli evaluate

test:
	$(PYTHON) -m pytest -q

clean:
	rm -rf results/*
	find . -name __pycache__ -type d -exec rm -rf {} +

## Also delete cached structures and generated designs
clean-data: clean
	rm -f data/*.pdb data/*.cif data/*.fasta
