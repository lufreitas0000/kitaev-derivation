# Discovery & Pivot Log — `kitaev-derivation`

| Date | Module/Context | Obstruction | Pivot/Resolution | Status |
| :--- | :--- | :--- | :--- | :--- |
| 2026-11-19 | `KitaevDerivation.jl` -> `kitaev-derivation` | Legacy `Symbolics.jl v5.5` precompile breakage on Julia 1.12, `ComplexF64` ket evaluator drifting from symbolic AST, `Gamma/Gamma'` double-counting on `y`-bond, and embedded `lean/Axioms.lean` violating Rule 5 | Rebuilt in pure `uv` Python (`SymPy` exact arithmetic everywhere + fixed `Gamma/Gamma'` averages + right-to-left multiplet projection), added independent 66-state `scipy.sparse` Fock oracle (`< 1e-13`), archived `.jl` repo, and routed Lean proofs to `finite-QMB` | Resolved |
