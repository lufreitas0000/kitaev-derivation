"""Directional t2g hopping matrices and Fock-sector hopping operators T_1, T_{-1}, T_0."""

from __future__ import annotations

import sympy as sp
from sympy import Expr, Matrix, Symbol

from kitaev_derivation.domain.algebra import (
    ZERO,
    AbstractQuantumOperator,
    FermionA,
    FermionC,
    OperatorString,
    OperatorSum,
    ParticleProjector,
    ScaledOperator,
)
from kitaev_derivation.domain.conventions import ORBITALS, SPINS


def directional_hopping_matrix(
    bond: str,
    t2: Expr | Symbol | float = 0,
    t1: Expr | Symbol | float = 0,
    t3: Expr | Symbol | float = 0,
    t4: Expr | Symbol | float = 0,
) -> Matrix:
    """Return the 3x3 real t2g orbital hopping matrix T_gamma on an edge-shared octahedra bond.

    Uses the standard Rau-Lee-Kee (PRL 112, 077204 (2014)) orbital hopping parameterization:
        z-bond (xy edge): [[t1, t2, t4], [t2, t1, t4], [t4, t4, t3]]
        x-bond (yz edge): [[t3, t4, t4], [t4, t1, t2], [t4, t2, t1]]
        y-bond (zx edge): [[t1, t4, t2], [t4, t3, t4], [t2, t4, t1]]
    Note that setting t1=t3=t4=0 recovers pure ligand-mediated Jackeli-Khaliullin hopping t2.
    """
    s1, s2, s3, s4 = map(sp.sympify, (t1, t2, t3, t4))
    if bond == "z":
        return Matrix([[s1, s2, s4], [s2, s1, s4], [s4, s4, s3]])
    if bond == "x":
        return Matrix([[s3, s4, s4], [s4, s1, s2], [s4, s2, s1]])
    if bond == "y":
        return Matrix([[s1, s4, s2], [s4, s3, s4], [s2, s4, s1]])
    raise ValueError(f"Unsupported bond '{bond}'; expected one of ('x', 'y', 'z').")


def build_t1_operator(
    bond: str,
    t2: Expr | Symbol | float,
    to_site: str = "i",
    from_site: str = "j",
    *,
    t1: Expr | Symbol | float = 0,
    t3: Expr | Symbol | float = 0,
    t4: Expr | Symbol | float = 0,
) -> OperatorSum:
    """Build directed doublon-creation hopping T_1(to <- from) = -sum P_to^(2) h^dagger T h P_from^(1)."""
    t_mat = directional_hopping_matrix(bond, t2=t2, t1=t1, t3=t3, t4=t4)
    p2_to = ParticleProjector(2, to_site)
    p1_from = ParticleProjector(1, from_site)
    terms: list[AbstractQuantumOperator] = []

    for ia, orb_a in enumerate(ORBITALS):
        for ib, orb_b in enumerate(ORBITALS):
            val = sp.simplify(t_mat[ia, ib])
            if val == 0:
                continue
            for spin in SPINS:
                op_str = OperatorString(
                    (
                        p2_to,
                        FermionC(to_site, orb_a, spin),
                        FermionA(from_site, orb_b, spin),
                        p1_from,
                    )
                )
                terms.append(ScaledOperator(-val, op_str))
    return OperatorSum(tuple(terms))


def build_tm1_operator(
    bond: str,
    t2: Expr | Symbol | float,
    to_site: str = "j",
    from_site: str = "i",
    *,
    t1: Expr | Symbol | float = 0,
    t3: Expr | Symbol | float = 0,
    t4: Expr | Symbol | float = 0,
) -> OperatorSum:
    """Build directed doublon-annihilation hopping T_{-1}(to <- from) = -sum P_to^(1) h^dagger T^dagger h P_from^(2)."""
    t_mat = directional_hopping_matrix(bond, t2=t2, t1=t1, t3=t3, t4=t4).H
    p1_to = ParticleProjector(1, to_site)
    p2_from = ParticleProjector(2, from_site)
    terms: list[AbstractQuantumOperator] = []

    for ia, orb_a in enumerate(ORBITALS):
        for ib, orb_b in enumerate(ORBITALS):
            val = sp.simplify(t_mat[ia, ib])
            if val == 0:
                continue
            for spin in SPINS:
                op_str = OperatorString(
                    (
                        p1_to,
                        FermionC(to_site, orb_a, spin),
                        FermionA(from_site, orb_b, spin),
                        p2_from,
                    )
                )
                terms.append(ScaledOperator(-val, op_str))
    return OperatorSum(tuple(terms))


def _drop_zero(op: AbstractQuantumOperator) -> AbstractQuantumOperator:
    return ZERO if op == ZERO else op
