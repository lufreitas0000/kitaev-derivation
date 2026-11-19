# Project Horizon — kitaev-derivation

**Authority:** Intra-project ordered milestone horizon (`Phase -> Milestone`). Cross-project portfolio priorities live in [`/home/lu/Projects/HORIZON.md`](../../../HORIZON.md) and QSL state strategy in [`../GUIDING_QUESTIONS.md`](../GUIDING_QUESTIONS.md).

## Horizon 1 — Current Sprint (`M2.1` Complete at Foundation Cutover)
- `M2.1`: Migrate symbolic Schrieffer-Wolff projector CAS from legacy `KitaevDerivation.jl` (including `0aed437` Appendix A commutator oracles) to `uv` pure Python (`SymPy` + `mypy --strict`), independently verified against the 66-state `scipy.sparse` Fock-space oracle (`< 1e-13`).

## Horizon 2 — Next Ready Milestones (Wave 3)
- `M2.2`: Implement fourth-order two-plaquette ring-exchange & $K_2$ next-nearest-neighbor extraction and export material benchmark tables (`.dat` + `.meta.json`) for `kitaev-mf` and `ksl-sc-bilayer`.

## Horizon 3 — Cross-Repository Formal Bridge
- Formal Lean 4 CAR and multiplet-projector orthogonality/completeness proofs routed to `research/qmb/finite-QMB` per [`doc/lean_handoff_finite_qmb.md`](doc/lean_handoff_finite_qmb.md).

## Adaptation Ledger
| Date | Prior expectation | New evidence / obstacle | Changed plan and authority | Evidence link |
| :--- | :--- | :--- | :--- | :--- |
| 2026-11-19 | Julia `Symbolics.jl` CAS + embedded `lean/` stub | `Symbolics v5.5` precompile fragility on Julia 1.12, `ComplexF64` projection drift, and Rule 5 Lean/Julia coupling | Migrated cleanly to `uv` pure Python (`SymPy` + 66-state `scipy.sparse` oracle) and routed Lean proofs to `finite-QMB`; legacy `.jl` archived | [LOG_DISCOVERY_AND_PIVOTS.md](LOG_DISCOVERY_AND_PIVOTS.md) |
