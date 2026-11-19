"""Tests for symbolic (J, K, Gamma, Gamma') extraction, Jackeli-Khaliullin limit, and C2/C3 symmetries."""

from __future__ import annotations

import sympy as sp
from sympy import simplify

from kitaev_derivation.services.exchange_extraction import (
    convert_tensor_dict_to_pauli,
    derive_bond_exchange_symbolic,
)


def test_jackeli_khaliullin_exact_symbolic_formula_and_jh0_cancellation() -> None:
    """Verify exact closed-form Jackeli-Khaliullin Kitaev exchange and J_H=0 destructive interference."""
    U, J_H, t2 = sp.symbols("U J_H t2", positive=True)
    res_z = derive_bond_exchange_symbolic("z", U, J_H, t2, convention="spin_12")

    expected_K = -8 * J_H * t2**2 / (3 * (U - 3 * J_H) * (U - J_H))
    assert simplify(res_z["J"]) == 0
    assert simplify(res_z["Gamma"]) == 0
    assert simplify(res_z["Gamma_prime"]) == 0
    assert simplify(res_z["K"] - expected_K) == 0

    # Jackeli-Khaliullin destructive interference when J_H = 0
    assert simplify(res_z["K"].subs(J_H, 0)) == 0
    assert simplify(res_z["J"].subs(J_H, 0)) == 0

    # Pauli convention conversion (S = sigma / 2 -> factor 1/4) matching qsl_conventions.json
    pauli_z = convert_tensor_dict_to_pauli(res_z)
    assert simplify(pauli_z["K"] - expected_K / 4) == 0


def test_c2_and_c3_point_group_symmetry_general_hopping() -> None:
    """Verify C2 isotropy (J_xx = J_yy on z-bond) and C3 bond permutation invariance for (t1, t2, t3, t4)."""
    U, J_H, t1, t2, t3, t4 = sp.symbols("U J_H t1 t2 t3 t4", real=True)

    res_x = derive_bond_exchange_symbolic("x", U, J_H, t2, t1=t1, t3=t3, t4=t4)
    res_y = derive_bond_exchange_symbolic("y", U, J_H, t2, t1=t1, t3=t3, t4=t4)
    res_z = derive_bond_exchange_symbolic("z", U, J_H, t2, t1=t1, t3=t3, t4=t4)

    # C2 symmetry on z-bond: J_xx == J_yy and J_xz == J_yz symbolically
    assert sp.cancel(res_z["J_xx"] - res_z["J_yy"]) == 0
    assert sp.cancel(res_z["J_xz"] - res_z["J_yz"]) == 0

    # C3 rotational equivalence across x, y, z bonds: identical (J, K, Gamma, Gamma') invariants!
    for key in ("J", "K", "Gamma", "Gamma_prime"):
        assert sp.cancel(res_x[key] - res_z[key]) == 0
        assert sp.cancel(res_y[key] - res_z[key]) == 0

    # Verify textbook Rau-Lee-Kee PRL 112, 077204 (2014) Eq. (3) dependence when t4 = 0:
    # Gamma'(t4=0) == 0, while Gamma(t2 != 0, t1 != t3) = 16 J_H t2 (t1 - t3) / (9 (U - J_H)(U - 3 J_H))
    assert sp.cancel(res_z["Gamma_prime"].subs(t4, 0)) == 0
    gamma_t4_0 = sp.factor(sp.cancel(res_z["Gamma"].subs(t4, 0)))
    expected_gamma = 16 * J_H * t2 * (t1 - t3) / (9 * (U - J_H) * (U - 3 * J_H))
    assert sp.cancel(gamma_t4_0 - expected_gamma) == 0
    assert sp.cancel(gamma_t4_0.subs(t2, 0)) == 0
    assert sp.cancel(gamma_t4_0.subs(t1, t3)) == 0
    assert sp.cancel(gamma_t4_0.subs(J_H, 0)) == 0
