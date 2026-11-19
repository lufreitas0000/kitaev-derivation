"""Tests for CAR normal ordering, projector algebra, and Appendix A commutator oracles."""

from __future__ import annotations

import itertools
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
import sympy as sp
from numpy.typing import NDArray

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
from kitaev_derivation.services.sparse_fock_oracle import _single_site_operators


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


def test_ast_rewrite_rule_substitution() -> None:
    """Verify symbolic AST rewrite-table substitution on algebraic commutators."""
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

    vs_t1_rhs = J_H * (3 * P_singlet - P_triplet) * T_1_core
    vs_tm1_rhs = -J_H * T_m1_core * (3 * P_singlet - P_triplet)
    vs_t0_rhs = 2 * J_H * (P_2 * T_0_core * P_exc - P_exc * T_0_core * P_2)
    rules_vs = {
        commutator(V_S, T_1): vs_t1_rhs,
        commutator(V_S, T_m1): vs_tm1_rhs,
        commutator(V_S, T_0): vs_t0_rhs,
    }
    assert apply_rewrite_rules(commutator(V_S, T_1), rules_vs) == vs_t1_rhs
    assert apply_rewrite_rules(commutator(V_S, T_m1), rules_vs) == vs_tm1_rhs
    assert apply_rewrite_rules(commutator(V_S, T_0), rules_vs) == vs_t0_rhs

    vl_t1_rhs = J_H * (2 * P_L0 + P_L1 - P_L2) * T_1_core
    vl_t0_rhs = J_H * (
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


@dataclass(frozen=True)
class _FockSector:
    act: Callable[[int, int], NDArray[np.complex128]]
    v_n: NDArray[np.complex128]
    v_s: NDArray[np.complex128]
    v_l: NDArray[np.complex128]
    p1: tuple[NDArray[np.complex128], NDArray[np.complex128]]
    p2: tuple[NDArray[np.complex128], NDArray[np.complex128]]
    m0: tuple[NDArray[np.complex128], NDArray[np.complex128], NDArray[np.complex128]]
    m1: tuple[NDArray[np.complex128], NDArray[np.complex128], NDArray[np.complex128]]


def _build_fock_sector(ntot: int, u_val: float, jh_val: float) -> _FockSector:
    """Construct exact 12-mode two-site Fock operators on the N_tot-hole subspace."""
    basis = list(itertools.combinations(range(12), ntot))
    idx_map = {st: k for k, st in enumerate(basis)}
    dim = len(basis)

    def act(c_mode: int, a_mode: int) -> NDArray[np.complex128]:
        mat = np.zeros((dim, dim), dtype=np.complex128)
        for col, st in enumerate(basis):
            p = list(st)
            if a_mode not in p:
                continue
            pos = p.index(a_mode)
            sgn = -1 if pos % 2 else 1
            p.pop(pos)
            if c_mode in p:
                continue
            ins = 0
            while ins < len(p) and p[ins] < c_mode:
                ins += 1
            if ins % 2:
                sgn = -sgn
            p.insert(ins, c_mode)
            mat[idx_map[tuple(p)], col] += sgn
        return mat

    def bilinear(op6: NDArray[np.complex128], site: int) -> NDArray[np.complex128]:
        off = site * 6
        out = np.zeros((dim, dim), dtype=np.complex128)
        for r in range(6):
            for c in range(6):
                val = op6[r, c]
                if abs(val) > 1e-15:
                    out = out + val * act(off + r, off + c)
        return out

    n6, s6, l_hole6, _ = _single_site_operators()
    n0, n1 = bilinear(n6, 0), bilinear(n6, 1)
    s2_0 = sum(
        (bilinear(s, 0) @ bilinear(s, 0) for s in s6),
        start=np.zeros((dim, dim), dtype=np.complex128),
    )
    s2_1 = sum(
        (bilinear(s, 1) @ bilinear(s, 1) for s in s6),
        start=np.zeros((dim, dim), dtype=np.complex128),
    )
    l2_0 = sum(
        (bilinear(l_mat, 0) @ bilinear(l_mat, 0) for l_mat in l_hole6),
        start=np.zeros((dim, dim), dtype=np.complex128),
    )
    l2_1 = sum(
        (bilinear(l_mat, 1) @ bilinear(l_mat, 1) for l_mat in l_hole6),
        start=np.zeros((dim, dim), dtype=np.complex128),
    )

    eye = np.eye(dim, dtype=np.complex128)
    v_n = 0.5 * (u_val - 3.0 * jh_val) * ((n0 - eye) @ (n0 - eye) + (n1 - eye) @ (n1 - eye))
    v_s = -2.0 * jh_val * (s2_0 + s2_1)
    v_l = -0.5 * jh_val * (l2_0 + l2_1)

    n0d = np.real(np.diag(n0))
    n1d = np.real(np.diag(n1))
    p1 = (
        np.diag(np.isclose(n0d, 1.0).astype(np.complex128)),
        np.diag(np.isclose(n1d, 1.0).astype(np.complex128)),
    )
    p2 = (
        np.diag(np.isclose(n0d, 2.0).astype(np.complex128)),
        np.diag(np.isclose(n1d, 2.0).astype(np.complex128)),
    )

    def mult(
        s2: NDArray[np.complex128], l2: NDArray[np.complex128], p_two: NDArray[np.complex128]
    ) -> tuple[NDArray[np.complex128], NDArray[np.complex128], NDArray[np.complex128]]:
        p_l0 = p_two @ ((l2 - 2.0 * eye) @ (l2 - 6.0 * eye) / 12.0) @ p_two
        p_l1 = p_two @ (s2 / 2.0) @ p_two
        p_l2 = p_two @ (l2 @ (l2 - 2.0 * eye) / 24.0) @ p_two
        return p_l0, p_l1, p_l2

    return _FockSector(
        act=act,
        v_n=v_n,
        v_s=v_s,
        v_l=v_l,
        p1=p1,
        p2=p2,
        m0=mult(s2_0, l2_0, p2[0]),
        m1=mult(s2_1, l2_1, p2[1]),
    )


def test_appendix_a_kanamori_hopping_commutators_on_two_site_fock() -> None:
    """Verify Appendix A commutators [V_N, T_m], [V_S, T_m], [V_L, T_m] directly on Fock space.

    Uses spin-independent orbital hopping T_orb on two 6-mode sites in N_tot=2 (66 states)
    for T_{+1}, T_{-1} and N_tot=3 (220 states) for T_0. Verifies exact operator identities
    [V, T_1] = (E_final - E_initial) T_1, which also checks the consistent sign relative to
    [T_1, V] (note V_N + V_S + V_L yields the exact gaps U+2J_H, U-3J_H, U-J_H).
    """
    u_val, jh_val = 4.0, 0.6
    rng = np.random.default_rng(20261119)
    t_orb = rng.normal(size=(3, 3))
    t6 = np.kron(np.eye(2), t_orb)

    def _comm(op_a: NDArray[np.complex128], op_b: NDArray[np.complex128]) -> NDArray[np.complex128]:
        return np.asarray(op_a @ op_b - op_b @ op_a, dtype=np.complex128)

    # 1) N_tot = 2 (66 states): T_{+1} and T_{-1}
    s2 = _build_fock_sector(2, u_val, jh_val)
    hop01 = sum(
        (t6[a, b] * s2.act(a, 6 + b) for a in range(6) for b in range(6)),
        start=np.zeros((66, 66), dtype=np.complex128),
    )
    hop10 = hop01.conj().T
    p0_1, p1_1 = s2.p1
    p0_2, p1_2 = s2.p2
    p0_l0, p0_l1, p0_l2 = s2.m0
    p1_l0, p1_l1, p1_l2 = s2.m1

    t1 = p0_2 @ hop01 @ p1_1 + p1_2 @ hop10 @ p0_1
    tm1 = p1_1 @ hop10 @ p0_2 + p0_1 @ hop01 @ p1_2

    assert np.max(np.abs(_comm(s2.v_n, t1) - (u_val - 3.0 * jh_val) * t1)) < 1e-13
    assert np.max(np.abs(_comm(s2.v_n, tm1) + (u_val - 3.0 * jh_val) * tm1)) < 1e-13

    rhs_vs_t1 = jh_val * (
        (3.0 * (p0_l0 + p0_l2) - p0_l1) @ hop01 @ p1_1
        + (3.0 * (p1_l0 + p1_l2) - p1_l1) @ hop10 @ p0_1
    )
    rhs_vs_tm1 = -jh_val * (
        p1_1 @ hop10 @ (3.0 * (p0_l0 + p0_l2) - p0_l1)
        + p0_1 @ hop01 @ (3.0 * (p1_l0 + p1_l2) - p1_l1)
    )
    assert np.max(np.abs(_comm(s2.v_s, t1) - rhs_vs_t1)) < 1e-13
    assert np.max(np.abs(_comm(s2.v_s, tm1) - rhs_vs_tm1)) < 1e-13

    rhs_vl_t1 = jh_val * (
        (2.0 * p0_l0 + p0_l1 - p0_l2) @ hop01 @ p1_1 + (2.0 * p1_l0 + p1_l1 - p1_l2) @ hop10 @ p0_1
    )
    assert np.max(np.abs(_comm(s2.v_l, t1) - rhs_vl_t1)) < 1e-13

    # Sum check: [V_N + V_S + V_L, T_1] reproduces exact Delta E_{0,1,2} gaps!
    total_gap_rhs = (
        (u_val + 2.0 * jh_val) * p0_l0 + (u_val - 3.0 * jh_val) * p0_l1 + (u_val - jh_val) * p0_l2
    ) @ hop01 @ p1_1 + (
        (u_val + 2.0 * jh_val) * p1_l0 + (u_val - 3.0 * jh_val) * p1_l1 + (u_val - jh_val) * p1_l2
    ) @ hop10 @ p0_1
    assert np.max(np.abs(_comm(s2.v_n + s2.v_s + s2.v_l, t1) - total_gap_rhs)) < 1e-13

    # 2) N_tot = 3 (220 states): T_0 preserving doubly occupied count (2,1) <-> (1,2)
    s3 = _build_fock_sector(3, u_val, jh_val)
    h3_01 = sum(
        (t6[a, b] * s3.act(a, 6 + b) for a in range(6) for b in range(6)),
        start=np.zeros((220, 220), dtype=np.complex128),
    )
    h3_10 = h3_01.conj().T
    q0_2, q1_2 = s3.p2
    q0_l0, q0_l1, q0_l2 = s3.m0
    q1_l0, q1_l1, q1_l2 = s3.m1
    q0_exc = q0_l1 - (q0_l0 + q0_l2)
    q1_exc = q1_l1 - (q1_l0 + q1_l2)

    t0 = q0_2 @ h3_01 @ q1_2 + q1_2 @ h3_10 @ q0_2
    assert np.max(np.abs(_comm(s3.v_n, t0))) < 1e-13

    rhs_vs_t0 = (
        2.0
        * jh_val
        * (
            (q0_2 @ h3_01 @ q1_exc - q0_exc @ h3_01 @ q1_2)
            + (q1_2 @ h3_10 @ q0_exc - q1_exc @ h3_10 @ q0_2)
        )
    )
    assert np.max(np.abs(_comm(s3.v_s, t0) - rhs_vs_t0)) < 1e-13

    def _vl_cross(
        pa0: NDArray[np.complex128],
        pa1: NDArray[np.complex128],
        pa2: NDArray[np.complex128],
        h: NDArray[np.complex128],
        pb0: NDArray[np.complex128],
        pb1: NDArray[np.complex128],
        pb2: NDArray[np.complex128],
    ) -> NDArray[np.complex128]:
        return np.asarray(
            3.0 * pa0 @ h @ pb2
            - 3.0 * pa2 @ h @ pb0
            - pa1 @ h @ pb0
            + pa0 @ h @ pb1
            + 2.0 * pa1 @ h @ pb2
            - 2.0 * pa2 @ h @ pb1,
            dtype=np.complex128,
        )

    rhs_vl_t0 = jh_val * (
        _vl_cross(q0_l0, q0_l1, q0_l2, h3_01, q1_l0, q1_l1, q1_l2)
        + _vl_cross(q1_l0, q1_l1, q1_l2, h3_10, q0_l0, q0_l1, q0_l2)
    )
    assert np.max(np.abs(_comm(s3.v_l, t0) - rhs_vl_t0)) < 1e-13
