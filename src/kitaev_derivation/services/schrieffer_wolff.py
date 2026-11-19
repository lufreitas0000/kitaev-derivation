"""Second-order Schrieffer-Wolff resolvent expansion across Kanamori d4 multiplet channels."""

from __future__ import annotations

from sympy import Expr, Symbol

from kitaev_derivation.domain.algebra import (
    AbstractQuantumOperator,
    MultipletProjector,
    OperatorString,
    OperatorSum,
    ParticleProjector,
    ScaledOperator,
    ZeroOp,
    apply_algebraic_rules,
    collect_terms,
)
from kitaev_derivation.domain.kanamori_soc import kanamori_eigenvalues
from kitaev_derivation.services.kinetic_hopping import build_t1_operator, build_tm1_operator


def apply_resolvent(
    h_t1: OperatorSum,
    U: Expr | Symbol,
    J_H: Expr | Symbol,
    site: str | None = None,
) -> OperatorSum:
    """Apply the exact d4 Kanamori resolvent R = -(H_0 - E_0)^{-1} to T_1.

    Replaces ParticleProjector(holes=2, site) with the spectral sum over L in {0, 1, 2}:
        P^(2)_site -> sum_L (-1 / Delta E_L) P_{[L=L]}_site
    """
    resolved_terms: list[AbstractQuantumOperator] = []

    for term in h_t1.terms:
        if not isinstance(term, ScaledOperator) or not isinstance(term.op, OperatorString):
            resolved_terms.append(term)
            continue

        factors = term.op.factors
        p2_idx = next(
            (i for i, f in enumerate(factors) if isinstance(f, ParticleProjector) and f.holes == 2),
            None,
        )
        if p2_idx is None:
            resolved_terms.append(term)
            continue

        p2 = factors[p2_idx]
        assert isinstance(p2, ParticleProjector)
        proj_site = site if site is not None else p2.site
        eig_map = kanamori_eigenvalues(U, J_H, proj_site)

        for label, delta_e in eig_map.items():
            new_factors = list(factors)
            new_factors[p2_idx] = MultipletProjector(label.L_label, proj_site)
            new_scalar = -term.scalar / delta_e
            resolved_terms.append(ScaledOperator(new_scalar, OperatorString(tuple(new_factors))))

    return OperatorSum(tuple(resolved_terms))


def _combine_directed_channel(tm1: OperatorSum, r_t1: OperatorSum) -> list[AbstractQuantumOperator]:
    """Multiply T_{-1} * (R T_1) and prune adjacent orthogonal projectors."""
    terms: list[AbstractQuantumOperator] = []
    for left in tm1.terms:
        if not isinstance(left, ScaledOperator):
            continue
        left_factors = left.op.factors if isinstance(left.op, OperatorString) else (left.op,)
        for right in r_t1.terms:
            if not isinstance(right, ScaledOperator):
                continue
            right_factors = (
                right.op.factors if isinstance(right.op, OperatorString) else (right.op,)
            )
            pruned = apply_algebraic_rules(left_factors + right_factors)
            if pruned and isinstance(pruned[0], ZeroOp):
                continue
            terms.append(ScaledOperator(left.scalar * right.scalar, OperatorString(pruned)))
    return terms


def compute_effective_hamiltonian(
    tm1: OperatorSum,
    t1: OperatorSum,
    U: Expr | Symbol,
    J_H: Expr | Symbol,
) -> OperatorSum:
    """Evaluate directed Schrieffer-Wolff contraction H_eff^(2) = T_{-1} R T_1."""
    r_t1 = apply_resolvent(t1, U, J_H)
    return collect_terms(OperatorSum(tuple(_combine_directed_channel(tm1, r_t1))))


def compute_bond_effective_hamiltonian(
    bond: str,
    U: Expr | Symbol,
    J_H: Expr | Symbol,
    t2: Expr | Symbol | float,
    site_i: str = "i",
    site_j: str = "j",
    *,
    t1: Expr | Symbol | float = 0,
    t3: Expr | Symbol | float = 0,
    t4: Expr | Symbol | float = 0,
) -> OperatorSum:
    """Compute the full undirected bond Schrieffer-Wolff operator H_eff^(2) = H_{i->d4} + H_{j->d4}.

    Includes both virtual double-hole channels:
        1) site_i goes d4, site_j goes d6 (T_{-1, j<-i} R_i T_{1, i<-j})
        2) site_j goes d4, site_i goes d6 (T_{-1, i<-j} R_j T_{1, j<-i})
    """
    t1_ij = build_t1_operator(bond, t2=t2, to_site=site_i, from_site=site_j, t1=t1, t3=t3, t4=t4)
    tm1_ji = build_tm1_operator(bond, t2=t2, to_site=site_j, from_site=site_i, t1=t1, t3=t3, t4=t4)

    t1_ji = build_t1_operator(bond, t2=t2, to_site=site_j, from_site=site_i, t1=t1, t3=t3, t4=t4)
    tm1_ij = build_tm1_operator(bond, t2=t2, to_site=site_i, from_site=site_j, t1=t1, t3=t3, t4=t4)

    terms = _combine_directed_channel(tm1_ji, apply_resolvent(t1_ij, U, J_H, site=site_i))
    terms.extend(_combine_directed_channel(tm1_ij, apply_resolvent(t1_ji, U, J_H, site=site_j)))
    return collect_terms(OperatorSum(tuple(terms)))
