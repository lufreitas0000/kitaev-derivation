"""Exact symbolic Hubbard-Kanamori d4 multiplet basis and t2g spin-orbit operator."""

from __future__ import annotations

import sympy as sp
from sympy import Expr, Matrix, Rational, Symbol, sqrt

from kitaev_derivation.domain.algebra import MultipletProjector
from kitaev_derivation.domain.conventions import (
    orbital_angular_momentum_matrices,
    pauli_matrices,
)

TwoHoleAmplitudes = dict[tuple[tuple[str, str], tuple[str, str]], Expr]


def _ordered_pair(
    o1: str,
    s1: str,
    o2: str,
    s2: str,
    coeff: Expr,
) -> tuple[tuple[tuple[str, str], tuple[str, str]], Expr]:
    """Canonicalize two-hole Slater determinant order with fermionic exchange sign."""
    m1 = (o1, s1)
    m2 = (o2, s2)
    if m1 < m2:
        return (m1, m2), coeff
    return (m2, m1), -coeff


def d4_multiplet_projectors(
    site: str = "i",
) -> tuple[MultipletProjector, MultipletProjector, MultipletProjector]:
    """Return canonical multiplet projector triple (P_[L=0], P_[L=1], P_[L=2])."""
    return (
        MultipletProjector("P_[L=0]", site),
        MultipletProjector("P_[L=1]", site),
        MultipletProjector("P_[L=2]", site),
    )


def kanamori_eigenvalues(
    U: Expr | Symbol,
    J_H: Expr | Symbol,
    site: str = "i",
) -> dict[MultipletProjector, Expr]:
    """Map d4 multiplet projectors to exact Hubbard-Kanamori excitation energies Delta E_L.

    Uses the spherically symmetric Kanamori Hamiltonian V = (U - 3 J_H)/2 (N_i - 1)^2 - 2 J_H S_i^2 - (J_H/2) L_i^2
    relative to the two-site d5-d5 ground state (E_0 = -5 J_H):
        P_[L=0] (^1S, dim=1): U + 2 J_H
        P_[L=1] (^3P, dim=9): U - 3 J_H
        P_[L=2] (^1D, dim=5): U - J_H
    """
    p0, p1, p2 = d4_multiplet_projectors(site)
    return {
        p0: U + 2 * J_H,
        p1: U - 3 * J_H,
        p2: U - J_H,
    }


def d4_multiplet_basis() -> dict[str, list[TwoHoleAmplitudes]]:
    """Construct exact SymPy two-hole d4 multiplet eigenstates across L in {0, 1, 2}.

    Returns:
        Dict keyed by 'P_[L=0]' (1 state), 'P_[L=1]' (9 states), 'P_[L=2]' (5 states).
    """
    inv_sqrt2 = Rational(1, 1) / sqrt(2)
    inv_sqrt3 = Rational(1, 1) / sqrt(3)
    inv_sqrt6 = Rational(1, 1) / sqrt(6)

    l0_states: list[TwoHoleAmplitudes] = [
        dict(
            [
                _ordered_pair("x", "dn", "x", "up", inv_sqrt3),
                _ordered_pair("y", "dn", "y", "up", inv_sqrt3),
                _ordered_pair("z", "dn", "z", "up", inv_sqrt3),
            ]
        )
    ]

    l1_states: list[TwoHoleAmplitudes] = []
    l2_states: list[TwoHoleAmplitudes] = []

    for o1, o2 in [("x", "y"), ("x", "z"), ("y", "z")]:
        l1_states.append(dict([_ordered_pair(o1, "up", o2, "up", sp.Integer(1))]))
        l1_states.append(dict([_ordered_pair(o1, "dn", o2, "dn", sp.Integer(1))]))
        l1_states.append(
            dict(
                [
                    _ordered_pair(o1, "up", o2, "dn", inv_sqrt2),
                    _ordered_pair(o1, "dn", o2, "up", inv_sqrt2),
                ]
            )
        )
        l2_states.append(
            dict(
                [
                    _ordered_pair(o1, "up", o2, "dn", inv_sqrt2),
                    _ordered_pair(o1, "dn", o2, "up", -inv_sqrt2),
                ]
            )
        )

    l2_states.append(
        dict(
            [
                _ordered_pair("x", "dn", "x", "up", inv_sqrt2),
                _ordered_pair("y", "dn", "y", "up", -inv_sqrt2),
            ]
        )
    )
    l2_states.append(
        dict(
            [
                _ordered_pair("x", "dn", "x", "up", inv_sqrt6),
                _ordered_pair("y", "dn", "y", "up", inv_sqrt6),
                _ordered_pair("z", "dn", "z", "up", -2 * inv_sqrt6),
            ]
        )
    )

    return {
        "P_[L=0]": l0_states,
        "P_[L=1]": l1_states,
        "P_[L=2]": l2_states,
    }


def single_site_soc_matrix(lam: Expr | Symbol) -> Matrix:
    """Return exact 6x6 single-hole t2g SOC matrix H_SOC = +lambda * (S . l_eff).

    Eigenvalues are -lambda (2-fold j_eff=1/2 Kramers doublet) and +lambda/2 (4-fold j_eff=3/2).
    """
    l_mats = orbital_angular_momentum_matrices()
    paulis = pauli_matrices()
    h_so = sp.zeros(6, 6)
    for axis in ("x", "y", "z"):
        s_mat = Rational(1, 2) * paulis[axis]
        h_so += sp.kronecker_product(s_mat, l_mats[axis])
    return lam * h_so
