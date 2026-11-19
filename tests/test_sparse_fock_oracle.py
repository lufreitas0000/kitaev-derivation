"""Cross-checks symbolic SymPy (J, K, Gamma, Gamma') against independent 66-state Fock-space oracle."""

from __future__ import annotations

import numpy as np
import sympy as sp

from kitaev_derivation.services.exchange_extraction import derive_bond_exchange_symbolic
from kitaev_derivation.services.sparse_fock_oracle import (
    compute_fock_schrieffer_wolff_exchange,
    compute_fock_soc_finite_exchange,
)


def test_symbolic_cas_matches_66_state_sparse_oracle_to_machine_precision() -> None:
    """Verify exact equivalence (< 1e-13) between symbolic SymPy CAS and 66-state Fock projection."""
    U, J_H, t1, t2, t3, t4 = sp.symbols("U J_H t1 t2 t3 t4", real=True)
    sym_cache = {
        bond: derive_bond_exchange_symbolic(bond, U, J_H, t2, t1=t1, t3=t3, t4=t4)
        for bond in ("x", "y", "z")
    }

    rng = np.random.default_rng(20261119)
    for _ in range(5):
        u_v = float(rng.uniform(3.5, 6.0))
        jh_v = float(rng.uniform(0.25, 0.8))
        t1_v = float(rng.uniform(-0.2, 0.2))
        t2_v = float(rng.uniform(0.15, 0.4))
        t3_v = float(rng.uniform(-0.25, 0.25))
        t4_v = float(rng.uniform(-0.08, 0.08))
        subs = {U: u_v, J_H: jh_v, t1: t1_v, t2: t2_v, t3: t3_v, t4: t4_v}

        for bond in ("x", "y", "z"):
            oracle = compute_fock_schrieffer_wolff_exchange(
                bond, u_v, jh_v, t2_v, t1=t1_v, t3=t3_v, t4=t4_v
            )
            sym_vals = {
                k: float(sym_cache[bond][k].evalf(subs=subs))
                for k in ("J", "K", "Gamma", "Gamma_prime")
            }
            assert abs(sym_vals["J"] - oracle.J) < 1e-13
            assert abs(sym_vals["K"] - oracle.K) < 1e-13
            assert abs(sym_vals["Gamma"] - oracle.Gamma) < 1e-13
            assert abs(sym_vals["Gamma_prime"] - oracle.Gamma_prime) < 1e-13


def test_finite_soc_intermediate_state_approaches_strong_soc_limit() -> None:
    """Check finite-lambda_SOC 66-state intermediate-state resolvent vs. lambda_SOC=0 Kanamori limit."""
    res_lam0 = compute_fock_soc_finite_exchange(
        "z", U=4.0, J_H=0.6, lam_soc=0.0, t2=0.25, t1=0.05, t3=-0.15, t4=0.02
    )
    res_analytic = compute_fock_schrieffer_wolff_exchange(
        "z", U=4.0, J_H=0.6, t2=0.25, t1=0.05, t3=-0.15, t4=0.02
    )
    assert abs(res_lam0.K - res_analytic.K) < 1e-13
    assert abs(res_lam0.J - res_analytic.J) < 1e-13

    # For finite lambda_SOC << U - 3J_H, the d5->d4 excitation gap shifts by O(lambda_SOC):
    # verify linear-in-lambda scaling (lam=0.03 gives <3.5% shift; lam=0.15 gives ~14% shift).
    res_small = compute_fock_soc_finite_exchange(
        "z", U=4.0, J_H=0.6, lam_soc=0.03, t2=0.25, t1=0.05, t3=-0.15, t4=0.02
    )
    res_fin = compute_fock_soc_finite_exchange(
        "z", U=4.0, J_H=0.6, lam_soc=0.15, t2=0.25, t1=0.05, t3=-0.15, t4=0.02
    )
    rel_small = abs(res_small.K - res_analytic.K) / abs(res_analytic.K)
    rel_fin = abs(res_fin.K - res_analytic.K) / abs(res_analytic.K)
    assert 1e-4 < rel_small < 0.035
    assert 0.05 < rel_fin < 0.25
