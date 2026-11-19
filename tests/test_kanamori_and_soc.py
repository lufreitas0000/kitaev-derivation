"""Tests for d4 Hubbard-Kanamori multiplet eigenvalues and t2g SOC Kramers doublets."""

from __future__ import annotations

import sympy as sp
from sympy import Matrix, Rational, simplify

from kitaev_derivation.domain.conventions import (
    kramers_doublet_vectors,
    orbital_angular_momentum_matrices,
    pauli_matrices,
)
from kitaev_derivation.domain.kanamori_soc import (
    d4_multiplet_basis,
    d4_multiplet_projectors,
    kanamori_eigenvalues,
    single_site_soc_matrix,
)
from kitaev_derivation.services.sparse_fock_oracle import compute_fock_kanamori_spectrum


def test_kanamori_d4_basis_dimensions_and_orthonormality() -> None:
    basis = d4_multiplet_basis()
    all_states = basis["P_[L=0]"] + basis["P_[L=1]"] + basis["P_[L=2]"]
    assert len(basis["P_[L=0]"]) == 1
    assert len(basis["P_[L=1]"]) == 9
    assert len(basis["P_[L=2]"]) == 5
    assert len(all_states) == 15

    # Exact SymPy orthonormality across all 15 two-hole states
    for i, bra in enumerate(all_states):
        for j, ket in enumerate(all_states):
            overlap = sp.Integer(0)
            for pair, val in ket.items():
                if pair in bra:
                    overlap += sp.conjugate(bra[pair]) * val
            expected = sp.Integer(1) if i == j else sp.Integer(0)
            assert simplify(overlap - expected) == 0


def test_kanamori_eigenvalues_match_exact_66_state_spectrum() -> None:
    U, J_H = sp.symbols("U J_H", real=True)
    p0, p1, p2 = d4_multiplet_projectors("i")
    eigs = kanamori_eigenvalues(U, J_H, "i")
    assert simplify(eigs[p0] - eigs[p1] - 5 * J_H) == 0
    assert simplify(eigs[p2] - eigs[p1] - 2 * J_H) == 0
    assert simplify(eigs[p0].subs(J_H, 0) - U) == 0

    u_num, jh_num = 4.2, 0.55
    e0, spectrum = compute_fock_kanamori_spectrum(u_num, jh_num)
    assert abs(e0 - (-5.0 * jh_num)) < 1e-12
    # Gaps sorted ascending: U-3JH (dim 9), U-JH (dim 5), U+2JH (dim 1)
    assert len(spectrum) == 3
    expected = [
        (float(eigs[p1].subs({U: u_num, J_H: jh_num})), 9),
        (float(eigs[p2].subs({U: u_num, J_H: jh_num})), 5),
        (float(eigs[p0].subs({U: u_num, J_H: jh_num})), 1),
    ]
    for (g_num, m_num), (g_sym, m_sym) in zip(spectrum, expected, strict=True):
        assert m_num == m_sym
        assert abs(g_num - g_sym) < 1e-12


def test_t2g_soc_and_kramers_effective_spin_algebra() -> None:
    lam = sp.Symbol("lambda", positive=True)
    h_so = single_site_soc_matrix(lam)
    plus, minus = kramers_doublet_vectors()

    # Exact eigenvectors at -lambda (j_eff = 1/2 Kramers doublet)
    assert (h_so * plus + lam * plus).applyfunc(simplify).is_zero_matrix
    assert (h_so * minus + lam * minus).applyfunc(simplify).is_zero_matrix
    assert simplify((plus.H * minus)[0, 0]) == 0
    assert simplify((plus.H * plus)[0, 0] - 1) == 0
    assert simplify((minus.H * minus)[0, 0] - 1) == 0

    # Spectrum is {-lambda (x2), +lambda/2 (x4)}, gap = 3*lambda/2
    evals = h_so.subs(lam, 2).eigenvals()
    assert evals == {sp.Integer(-2): 2, sp.Integer(1): 4}

    # Projection of j_eff = S + l_eff onto {|+>, |->} reproduces sigma / 2 identically
    l_mats = orbital_angular_momentum_matrices()
    paulis = pauli_matrices()
    basis = Matrix.hstack(plus, minus)
    for axis in ("x", "y", "z"):
        j_eff = sp.kronecker_product(
            Rational(1, 2) * paulis[axis], sp.eye(3)
        ) + sp.kronecker_product(sp.eye(2), l_mats[axis])
        proj = (basis.H * j_eff * basis).applyfunc(simplify)
        assert (proj - Rational(1, 2) * paulis[axis]).is_zero_matrix
