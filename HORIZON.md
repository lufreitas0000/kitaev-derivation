# Project Horizon — kitaev-derivation

**Authority:** Intra-project ordered milestone horizon (`Phase -> Milestone`). Cross-project portfolio priorities live in [`/home/lu/Projects/HORIZON.md`](../../../HORIZON.md) and QSL state strategy in [`../GUIDING_QUESTIONS.md`](../GUIDING_QUESTIONS.md).

## Horizon 1 — Current Sprint (`M2.1` & `M2.2` Completed)
- `M2.1`: Migrated symbolic Schrieffer-Wolff projector CAS from legacy `KitaevDerivation.jl` to `uv` pure Python (`SymPy` + `mypy --strict`), independently verified against the 66-state `scipy.sparse` Fock-space oracle (`< 1e-13`).
- `M2.2`: Implemented fourth-order BDL/Takahashi bond exchange (`KD-003`) and 924-state 6-site Kramers hexagon ring exchange via Peierls $\pi$-flux subtraction (`KD-004`), exported idealized material tables (`data/exports/extended_kitaev_material_couplings.{dat,meta.json}`), and opened formal Lean 4 bridge PR `lufreitas0000/finite-QMB#52`.

## Horizon 2 — Downstream Ingestion & Wave 4 Refinement
- Ingest `data/exports/extended_kitaev_material_couplings.dat` into `kitaev-mf` / `ksl-sc-bilayer` Wave 3 extended-$J\text{-}K\text{-}\Gamma\text{-}\Gamma'$ sweeps and monitor `finite-QMB#52` review.

## Horizon 3 — Cross-Repository Formal Bridge
- Formal Lean 4 CAR and multiplet-projector orthogonality/completeness proofs routed to `research/qmb/finite-QMB` per [`doc/lean_handoff_finite_qmb.md`](doc/lean_handoff_finite_qmb.md).

## Adaptation Ledger
| Date | Prior expectation | New evidence / obstacle | Changed plan and authority | Evidence link |
| :--- | :--- | :--- | :--- | :--- |
| 2026-11-19 | Julia `Symbolics.jl` CAS + embedded `lean/` stub | `Symbolics v5.5` precompile fragility on Julia 1.12, `ComplexF64` projection drift, and Rule 5 Lean/Julia coupling | Migrated cleanly to `uv` pure Python (`SymPy` + 66-state `scipy.sparse` oracle) and routed Lean proofs to `finite-QMB`; legacy `.jl` archived | [LOG_DISCOVERY_AND_PIVOTS.md](LOG_DISCOVERY_AND_PIVOTS.md) |
