# kitaev-derivation

Exact symbolic Schrieffer-Wolff projector Computer Algebra System (SymPy) and independent 66-state two-hole Fock-space oracle (`scipy.sparse`) deriving microscopic $(J, K, \Gamma, \Gamma')$ effective spin Hamiltonians for $d^5$ edge-shared octahedral Kitaev materials.

Migrated cleanly from legacy Julia `KitaevDerivation.jl` (`archive/KitaevDerivation.jl/`, including checkpoint `0aed437`) into unified `uv`-managed pure Python + strict typing (`mypy --strict`, `ruff`), with formal CAR/projector theorems routed to Lean 4 in [`research/qmb/finite-QMB`](../../qmb/finite-QMB/README.md) (see [`doc/lean_handoff_finite_qmb.md`](doc/lean_handoff_finite_qmb.md)).

## Goal

Establish a bug-free, algebraically exact symbolic Schrieffer-Wolff projector CAS in Python (`SymPy`) verified to machine precision (`< 1e-13`) against an independent 66-state two-hole Fock-space numerical oracle (`scipy.sparse`) across $x, y, z$ bonds and $(U, J_H, t_1, t_2, t_3, t_4, \lambda_{\mathrm{SOC}})$ parameters, with ACID benchmark tables and self-contained Lab Book / Manuscript / Blueprint deliverables.

## Scope

- **Owns:** Symbolic CAR normal-ordering, exact $d^4$ Hubbard-Kanamori multiplet projectors ($P_{[L=0,1,2]}$), $t_{2g}$ $T\text{–}P$ SOC Kramers doublet projection, second-order (and future fourth-order ring-exchange) Schrieffer-Wolff reductions, independent 66-state Fock oracle, and ACID parameter exports for `kitaev-mf` / `ksl-sc-bilayer`.
- **Strictly Out of Scope:** Lattice-scale Bogoliubov-de Gennes mean-field loops (`kitaev-mf`), multi-site spin ED / 2-RDM SDP (`kitaev-numerics`), spectral walks (`kitaev-walk-geometry`), and formal Lean 4 proofs (`research/qmb/finite-QMB`).

## Definition of Done (DoD)

- `make verify` exits `0` (`ruff check`, `ruff format --check`, `mypy --strict src tests scripts`, and full `pytest` suite including randomized 66-state oracle parity `< 1e-13`).
- Benchmark table `data/processed/benchmarks/exchange_vs_fock_oracle.dat` + `.meta.json` generated deterministically with ACID `#` provenance headers.
- `doc/labbook/labbook_kitaev_derivation.pdf`, `doc/manuscript/manuscript_kitaev_derivation.pdf`, and `build/blueprint/index.html` compile with zero LaTeX errors or undefined references/citations.
- `scripts/uow.py check --repo kitaev-derivation`, `governance/pg10.py check`, and `scripts/fiscalize_state.py --repo kitaev-derivation` pass with zero warnings.

## Next Steps

Complete Wave 3 fourth-order ring-exchange extraction (`M2.2`) and $\alpha\text{-RuCl}_3$ / $d^5\text{–}d^7$ benchmark parameter tables for downstream `kitaev-mf` and `ksl-sc-bilayer` extended-$K\text{–}J\text{–}\Gamma\text{–}\Gamma'$ campaigns after Wave 1/Wave 2 QSL gates close.
