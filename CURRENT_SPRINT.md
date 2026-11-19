# Current Sprint — kitaev-derivation / M2.1

Owning milestone: M2.1 (see [HORIZON.md](HORIZON.md))
Machine state: [PROJECT_STATE.yaml](PROJECT_STATE.yaml)

## Scientific / Engineering Objective

Complete the pure-Python `SymPy` migration of `KitaevDerivation.jl` (`M2.1`), verifying Appendix A Kanamori-hopping commutators (`[V_N, T_m]`, `[V_S, T_m]`, `[V_L, T_m]`) directly on 12-mode two-site Fock spaces (`N_tot=2,3`) and validating exact second-order $(J, K, \Gamma, \Gamma')$ extraction to `< 1e-13` against an independent 66-state two-hole Fock-space oracle (`scipy.sparse`), while decoupling formal Lean 4 proofs via [`doc/lean_handoff_finite_qmb.md`](doc/lean_handoff_finite_qmb.md).

## Scope and Non-Goals

Owns symbolic two-site $d^4/d^5/d^6$ Schrieffer-Wolff exchange extraction, unit-port conversion (`spin_12` $\leftrightarrow$ `pauli` per `../contracts/qsl_conventions.json`), ACID benchmark tables, and self-contained documentation. Fourth-order ring exchange (`M2.2`) remains staged for Wave 3.

## Exact Next Step

Wait for PI Wave 3 gate before executing fourth-order ring-exchange extraction (`M2.2`) and material $(J, K, \Gamma, \Gamma')$ sweeps for `kitaev-mf` / `ksl-sc-bilayer`.
