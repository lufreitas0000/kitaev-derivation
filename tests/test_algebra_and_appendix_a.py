"""Tests for CAR normal ordering, projector algebra, and Appendix A commutator oracles."""

from __future__ import annotations

import sympy as sp

from kitaev_derivation.domain.algebra import (
    IDENTITY,
    ZERO,
    FermionA,
    FermionC,
    GenericOp,
    LowEnergyProjector,
    MultipletProjector,
    OperatorString,
    ParticleProjector,
    anticommutator,
    apply_algebraic_rules,
    apply_rewrite_rules,
    commutator,
    dagger,
    sort_normal_order,
)
from kitaev_derivation.domain.kanamori_soc import d4_multiplet_projectors


def test_car_anticommutator_and_normal_order() -> None:
    ci = FermionA("i", "x", "up")
    cid = FermionC("i", "x", "up")
    cj = FermionA("j", "y", "dn")
    cjd = FermionC("j", "y", "dn")

    assert anticommutator(ci, cid) == IDENTITY
    assert anticommutator(ci, cjd) == ZERO
    assert anticommutator(ci, cj) == ZERO

    # c_i c_i^\dagger -> I - c_i^\dagger c_i
    terms = sort_normal_order((ci, cid))
    assert terms == [((IDENTITY,), 1), ((cid, ci), -1)]

    # Pauli exclusion c_i c_i = 0
    assert sort_normal_order((ci, ci)) == []


def test_projector_idempotency_orthogonality_and_intertwining() -> None:
    p0, p1, p2 = d4_multiplet_projectors("i")
    p_half = LowEnergyProjector("1/2", "i")
    p2_hole = ParticleProjector(2, "i")
    p1_hole = ParticleProjector(1, "i")

    assert p0.L_label == "P_[L=0]"
    assert p1.L_label == "P_[L=1]"
    assert p2.L_label == "P_[L=2]"

    for idx_a, pa in enumerate((p0, p1, p2)):
        assert apply_algebraic_rules((pa, pa)) == (pa,)
        assert apply_algebraic_rules((pa, p2_hole)) == (pa,)
        assert apply_algebraic_rules((pa, p_half)) == (ZERO,)
        assert dagger(pa) == pa
        for idx_b, pb in enumerate((p0, p1, p2)):
            expected = (pa,) if idx_a == idx_b else (ZERO,)
            assert apply_algebraic_rules((pa, pb)) == expected

    # Intertwining: c^\dagger P^{(1)} = P^{(2)} c^\dagger ; c P^{(2)} = P^{(1)} c
    cid = FermionC("i", "z", "up")
    ci = FermionA("i", "z", "up")
    assert apply_algebraic_rules((cid, p1_hole)) == (p2_hole, cid)
    assert apply_algebraic_rules((ci, p2_hole)) == (p1_hole, ci)
    assert apply_algebraic_rules((cid, p2_hole)) == (ZERO,)


def test_appendix_a1_bilin_commutator_derived_via_car() -> None:
    """Derive Eq. (A1): [h_{i,s}, h_j^\\dagger B h_j] = delta_{ij} (B h_i)_s via normal ordering."""
    # Check [h_{i,x,up}, h_{i,x,up}^\dagger h_{i,y,up}] = h_{i,y,up} directly by expanding both products!
    h_ix = FermionA("i", "x", "up")
    hd_ix = FermionC("i", "x", "up")
    h_iy = FermionA("i", "y", "up")

    left = sort_normal_order((h_ix, hd_ix, h_iy))
    right = sort_normal_order((hd_ix, h_iy, h_ix))
    # Combine normal-ordered monomials with signs
    counts: dict[tuple[object, ...], int] = {}
    for seq, sgn in left:
        counts[seq] = counts.get(seq, 0) + sgn
    for seq, sgn in right:
        counts[seq] = counts.get(seq, 0) - sgn
    nonzero = {k: v for k, v in counts.items() if v != 0}
    assert nonzero == {(h_iy,): 1}


def test_appendix_a_macroscopic_commutator_oracles() -> None:
    """Ported Appendix A commutator rules from legacy checkpoint 0aed437 (Eqs. A8, A23, A25, A31, A54, A57)."""
    U, J_H = sp.symbols("U J_H", real=True)
    T_0, T_1, T_m1 = GenericOp("T_0"), GenericOp("T_1"), GenericOp("T_m1")
    T_1_core, T_m1_core, T_0_core = (
        GenericOp("T_1_core"),
        GenericOp("T_m1_core"),
        GenericOp("T_0_core"),
    )
    V_N, V_S, V_L = GenericOp("V_N"), GenericOp("V_S"), GenericOp("V_L")
    P_L0, P_L1, P_L2 = d4_multiplet_projectors("i")
    P_singlet = P_L0 + P_L2
    P_triplet = P_L1
    P_2 = ParticleProjector(2, "i")
    P_exc = GenericOp("P_exc")

    rules_vn = {
        commutator(V_N, T_1): (U - 3 * J_H) * T_1,
        commutator(V_N, T_0): ZERO,
        commutator(V_N, T_m1): -(U - 3 * J_H) * T_m1,
    }
    assert apply_rewrite_rules(commutator(V_N, T_1), rules_vn) == (U - 3 * J_H) * T_1
    assert apply_rewrite_rules(commutator(V_N, T_0), rules_vn) == ZERO
    assert apply_rewrite_rules(commutator(V_N, T_m1), rules_vn) == -(U - 3 * J_H) * T_m1

    vs_t1_rhs = -J_H * (3 * P_singlet - P_triplet) * T_1_core
    vs_tm1_rhs = J_H * T_m1_core * (3 * P_singlet - P_triplet)
    vs_t0_rhs = -2 * J_H * (P_2 * T_0_core * P_exc - P_exc * T_0_core * P_2)
    rules_vs = {
        commutator(V_S, T_1): vs_t1_rhs,
        commutator(V_S, T_m1): vs_tm1_rhs,
        commutator(V_S, T_0): vs_t0_rhs,
    }
    assert apply_rewrite_rules(commutator(V_S, T_1), rules_vs) == vs_t1_rhs
    assert apply_rewrite_rules(commutator(V_S, T_m1), rules_vs) == vs_tm1_rhs
    assert apply_rewrite_rules(commutator(V_S, T_0), rules_vs) == vs_t0_rhs

    vl_t1_rhs = -J_H * (2 * P_L0 + P_L1 - P_L2) * T_1_core
    vl_t0_rhs = -J_H * (
        3 * P_L0 * T_0_core * P_L2
        - 3 * P_L2 * T_0_core * P_L0
        - P_L1 * T_0_core * P_L0
        + P_L0 * T_0_core * P_L1
        + 2 * P_L1 * T_0_core * P_L2
        - 2 * P_L2 * T_0_core * P_L1
    )
    rules_vl = {commutator(V_L, T_1): vl_t1_rhs, commutator(V_L, T_0): vl_t0_rhs}
    assert apply_rewrite_rules(commutator(V_L, T_1), rules_vl) == vl_t1_rhs
    assert apply_rewrite_rules(commutator(V_L, T_0), rules_vl) == vl_t0_rhs
    assert isinstance(OperatorString((MultipletProjector("P_[L=0]", "i"),)), OperatorString)
