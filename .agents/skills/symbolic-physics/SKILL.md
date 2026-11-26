---
name: symbolic-physics
description: Pure-Python SymPy workflow for t2g orbital multiplet projectors, Kanamori d4 eigenvalues, right-to-left Schrieffer-Wolff exchange derivations, and sparse Fock oracle cross-checks in kitaev-derivation.
---
<!-- @LOCAL-OVERRIDE -->
# Symbolic Physics Workflow (`kitaev-derivation` Pure Python 3.12)

Exact analytical derivation skill for `research/qsl/kitaev-derivation` using pure Python 3.12 (`sympy`, `numpy`, `scipy`) with Lean 4 convention bridges—zero floating-point truncation and zero Julia legacy dependencies.

## Canonical Conventions (`src/kitaev_derivation/`)

1. **Exact SymPy Conventions (`domain/conventions.py`):**
   - Represent all amplitudes, angles, and coupling ratios via exact `Rational`, `sqrt`, and `I` on `ORBITALS = ("x", "y", "z")`, `SPINS = ("up", "dn")`, `BONDS = ("x", "y", "z")`.
   - Use `orbital_angular_momentum_matrices()` respecting $T$-$P$ equivalence $L_{t_{2g}} = -l_{\mathrm{eff}}$ ($[l_{\mathrm{eff}}^x, l_{\mathrm{eff}}^y] = i l_{\mathrm{eff}}^z$) and `pauli_matrices()` without silent basis rotations.
2. **Kanamori $d^4$ Projectors & Eigenvalues (`domain/kanamori_soc.py`):**
   - Build exact spectral projectors `d4_multiplet_projectors(site)` yielding `(P_[L=0], P_[L=1], P_[L=2])` with `kanamori_eigenvalues(U, J_H)` ($U + 2J_H$, $U - 3J_H$, $U - J_H$).
   - Canonicalize two-hole Slater determinants via `_ordered_pair(...)` so fermionic canonical anticommutation signs remain invariant.
3. **Schrieffer-Wolff Perturbation Pipeline (`domain/algebra.py`, `services/`):**
   - Assemble bond kinetic operators (`services/kinetic_hopping.py`: `build_t1_operator`, `build_tm1_operator`) and propagate intermediate states strictly right-to-left through `services/schrieffer_wolff.py` (`apply_resolvent`) and `services/fourth_order_exchange.py`.
   - Reduce normal-ordered strings solely with `apply_algebraic_rules` and `collect_terms`, then extract effective spin couplings via `services/exchange_extraction.py`.

## Mandatory Non-LLM Verification Gate

- Every derived analytical exchange prefactor ($J$, $K$, $\Gamma$, $\Gamma'$) MUST match `services/sparse_fock_oracle.py` across generic $(t_1..t_4, J_H/U, \lambda/U)$ parameter slices, hermiticity, and $SU(2)$ / cubic-$C_{3v}$ symmetry limits before promotion.
- Run focused pytest checks under `uv run pytest` with exact algebraic assertions; never coerce symbolic radicals into floats to pass tests.
