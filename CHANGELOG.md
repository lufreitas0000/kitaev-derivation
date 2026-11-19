# Changelog & Epistemic Unit-of-Work Ledger (`kitaev-derivation`)

## [Unreleased]

### [2026-11-19] [M2.2-WAVE3] [Actor: antigravity] — Fourth-order BDL/Takahashi exchange, 924-state hexagon ring exchange, and finite-QMB bridge PR (#52)

#### Ideas & Contracts Evolved
- Quantified fourth-order $O(t^4)$ bond exchange separating $\lambda_{\mathrm{SOC}} \to \infty$ charge renormalization from finite-$\lambda_{\mathrm{SOC}}$ neutral $j_{\mathrm{eff}}=3/2$ intermediate channels (`KD-003`), finding strong finite-SOC fourth-order shifts for idealized $\alpha\text{-RuCl}_3$ (`J4 = -0.001797`, `K4 = +0.000972` at $\lambda_{\mathrm{SOC}}=0.15$).
- Derived 6-site $\binom{12}{6}=924$-state Kramers hexagon 6th-order ring exchange cleanly via Peierls $\Phi=0 - \Phi=\pi$ subtraction (`KD-004`) matching 3-point Richardson BDL to $\le 7.31 \times 10^{-7}$ (`K_ring(W_p) = 1.277e-5` and `3.248e-5` on reference hexagons).
- Formalized `FiniteQMB.Models.KitaevKanamoriSW` in `research/qmb/finite-QMB` ([PR #52](https://github.com/lufreitas0000/finite-QMB/pull/52), `9ab7df5`, merged to `main` @ `7098b031`).

#### Code & Artifacts Changed
- Added `src/kitaev_derivation/services/fourth_order_exchange.py`, `tests/test_fourth_order_exchange.py`, extended `scripts/generate_benchmarks.py` and `scripts/build_blueprint.py`, exported `data/exports/extended_kitaev_material_couplings.{dat,meta.json}`, and updated `doc/lean_handoff_finite_qmb.md`, LaTeX deliverables, and `PROJECT_STATE.yaml`.

### [2026-11-19] [UOW-20261119-01] [Actor: antigravity] — Migrate KitaevDerivation.jl to uv pure-Python SymPy CAS and 66-state Fock oracle (`M2.1`)

#### Ideas & Contracts Evolved
- Replaced legacy Julia `KitaevDerivation.jl` (`archive/KitaevDerivation.jl/`) with unified `uv` pure-Python (`SymPy` + `scipy.sparse`) symbolic Schrieffer-Wolff reduction (`M2.1`) and decoupled formal Lean 4 CAR/projector proofs to `research/qmb/finite-QMB` (`doc/lean_handoff_finite_qmb.md`).
- Resolved legacy defects (`docs/critical_review_2026.md` in `0aed437`): exact `SymPy` rationals everywhere, right-to-left multiplet projection wired through `compute_bond_effective_hamiltonian`, fixed `Gamma/Gamma'` averaging on `y`-bonds, and added independent 66-state two-hole Fock cross-checking (`< 1e-13`).

#### Code & Artifacts Changed
- Implemented `src/kitaev_derivation/` (`conventions.py`, `algebra.py`, `kanamori_soc.py`, `kinetic_hopping.py`, `schrieffer_wolff.py`, `exchange_extraction.py`, `sparse_fock_oracle.py`), unit/oracle tests, `scripts/generate_benchmarks.py`, `scripts/build_blueprint.py`, `doc/labbook/labbook_kitaev_derivation.tex`, and `doc/manuscript/manuscript_kitaev_derivation.tex`.
