# Cross-Repository Formal Proof Handoff: `kitaev-derivation` $\to$ `research/qmb/finite-QMB`

Per **Workspace Rule 5** (Lean 4 vs. Python/Numericals separation) and **Decision D2**, `research/qsl/kitaev-derivation` strictly owns symbolic/numerical Python CAS and exact 66-state Fock cross-checks, while formal machine-checked algebraic proofs reside in [`research/qmb/finite-QMB`](../../../qmb/finite-QMB/README.md).

## Source Material Migrated from Legacy `archive/KitaevDerivation.jl/lean/Axioms.lean`

Legacy `KitaevDerivation.jl` contained a 30-line uncompiled Lean stub defining `FermionAlgebra` (extending `NonUnitalNonAssocRing`), idempotent `Projector`, and `OrthogonalProjectors`.

## Proposed Bounded Targets for `finite-QMB` (Not Executed by QSL Setup)

1. **Finite-Mode $t_{2g}$ CAR Algebra (`Cl(6)`):** Prove normal-ordering reduction, Pauli nilpotency $(h_{is})^2 = 0$, and `[h_{is}, h_j^\dagger B h_j] = \delta_{ij} (B h_i)_s` (Appendix A Eq. A1) by connecting finite two-site $t_{2g}$ Fock witnesses to `finite-QMB`'s abstract algebraic `CARSystem` (`FiniteQMB/Fermion/CAR/{Defs,Basic}.lean`, `noncomputable section`).
2. **Partition of Unity for $d^4$ Multiplets:** Prove idempotency $P_{[L]}^2 = P_{[L]}$, mutual orthogonality $P_{[L_a]} P_{[L_b]} = 0$ ($L_a \neq L_b$), and completeness $P_{[L=0]} + P_{[L=1]} + P_{[L=2]} = P^{(2)}$ on the 15-dimensional $\Lambda^2 \mathbb{C}^6$ two-hole space.
3. **Schrieffer–Wolff Denominator Resolvent Identity:** Prove $(E_0 - H_0)^{-1} P^{(2)} = -\sum_{L=0}^2 \Delta E_L^{-1} P_{[L]}$ under $[H_0, P_{[L]}] = 0$.

## Executed Bridge PR (`finite-QMB`)

- **Merged Pull Request:** [lufreitas0000/finite-QMB#52](https://github.com/lufreitas0000/finite-QMB/pull/52) (branch `feat/kitaev-kanamori-sw-bridge` @ `9ab7df5152f61f9c0769b5f1ba079eaa0ad3d932`, merged into `finite-QMB` `main` at `7098b03165662d626969bf3056dc31acded2df1a` with passing `Lean 4 CI` and remote/local branch deleted).
- **Module:** `FiniteQMB.Models.KitaevKanamoriSW` (wired into `FiniteQMB.lean`), compiled with `lake build FiniteQMB.Models.KitaevKanamoriSW FiniteQMB` (`9200` jobs, `-DwarningAsError=true`, zero `sorry`, standard axioms `[propext, Classical.choice, Quot.sound]`):
  - `car_hopping_commutator`: `[c_i, \sum_{j,k} B_{jk} c_j^\dagger c_k] = \sum_k B_{ik} c_k` for `CARSystem`.
  - `MultipletProjectorFamily`: orthogonal excited-multiplet projectors (`Q0_idem`, `P0_mul_Q0`) and reduced resolvent `reducedResolvent = \sum_L \Delta_L^{-1} P_L` (`multiplet_resolvent_mul`).
  - `heff2_multiplet_expansion`: channel decomposition of second-order Schrieffer–Wolff superexchange across Kanamori/SOC multiplets.
  - `witness_heff2_00_nonzero`: explicit `Fin 3` / `Fin 2` non-vacuity matrix witness (`\Delta = \{2, 5\}`, `H_{\mathrm{eff}, 00}^{(2)} = -7/10 \ne 0`).
