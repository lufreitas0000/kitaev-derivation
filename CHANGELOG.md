# Changelog & Epistemic Unit-of-Work Ledger (`kitaev-derivation`)

## [Unreleased]

### [2026-11-19] [UOW-20261119-01] [Actor: antigravity] — Migrate KitaevDerivation.jl to uv pure-Python SymPy CAS and 66-state Fock oracle (`M2.1`)

#### Ideas & Contracts Evolved
- Replaced legacy Julia `KitaevDerivation.jl` (`archive/KitaevDerivation.jl/`) with unified `uv` pure-Python (`SymPy` + `scipy.sparse`) symbolic Schrieffer-Wolff reduction (`M2.1`) and decoupled formal Lean 4 CAR/projector proofs to `research/qmb/finite-QMB` (`doc/lean_handoff_finite_qmb.md`).
- Resolved legacy defects (`docs/critical_review_2026.md` in `0aed437`): exact `SymPy` rationals everywhere, right-to-left multiplet projection wired through `compute_bond_effective_hamiltonian`, fixed `Gamma/Gamma'` averaging on `y`-bonds, and added independent 66-state two-hole Fock cross-checking (`< 1e-13`).

#### Code & Artifacts Changed
- Implemented `src/kitaev_derivation/` (`conventions.py`, `algebra.py`, `kanamori_soc.py`, `kinetic_hopping.py`, `schrieffer_wolff.py`, `exchange_extraction.py`, `sparse_fock_oracle.py`), unit/oracle tests, `scripts/generate_benchmarks.py`, `scripts/build_blueprint.py`, `doc/labbook/labbook_kitaev_derivation.tex`, and `doc/manuscript/manuscript_kitaev_derivation.tex`.
