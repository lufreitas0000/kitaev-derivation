SHELL := /bin/bash -o pipefail
REPO_ROOT := $(shell pwd)
PYTHON := $(if $(wildcard $(REPO_ROOT)/.venv/bin/python),$(REPO_ROOT)/.venv/bin/python,python3)

.PHONY: help verify test lint format typecheck benchmarks labbook manuscript blueprint clean

help:
	@echo "kitaev-derivation :: Symbolic Schrieffer-Wolff CAS & 66-State Fock Oracle"
	@echo "  make verify      : Run lint, typecheck, and full pytest suite"
	@echo "  make benchmarks  : Export ACID exchange_vs_fock_oracle.dat + .meta.json"
	@echo "  make labbook     : Compile doc/labbook/labbook_kitaev_derivation.pdf"
	@echo "  make manuscript  : Compile doc/manuscript/manuscript_kitaev_derivation.pdf"
	@echo "  make blueprint   : Generate build/blueprint/index.html"

lint:
	$(PYTHON) -m ruff check src tests scripts
	$(PYTHON) -m ruff format --check src tests scripts

format:
	$(PYTHON) -m ruff check --fix src tests scripts
	$(PYTHON) -m ruff format src tests scripts

typecheck:
	$(PYTHON) -m mypy --strict src tests scripts

test:
	$(PYTHON) -m pytest -v

benchmarks:
	$(PYTHON) scripts/generate_benchmarks.py

labbook: benchmarks
	latexmk -pdf -interaction=nonstopmode -halt-on-error -cd doc/labbook/labbook_kitaev_derivation.tex

manuscript: benchmarks
	latexmk -pdf -interaction=nonstopmode -halt-on-error -cd doc/manuscript/manuscript_kitaev_derivation.tex

blueprint: benchmarks
	$(PYTHON) scripts/build_blueprint.py

verify: lint typecheck test

clean:
	rm -rf .pytest_cache .mypy_cache .ruff_cache build
	latexmk -C -cd doc/labbook/labbook_kitaev_derivation.tex || true
	latexmk -C -cd doc/manuscript/manuscript_kitaev_derivation.tex || true
