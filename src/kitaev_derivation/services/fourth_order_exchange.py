"""Wave 3 M2.2: Exact 66-state fourth-order BDL/Takahashi exchange and 6-site hexagon ring exchange."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np
from numpy.typing import NDArray

from kitaev_derivation.services.sparse_fock_oracle import (
    _build_hopping_66,
    _build_kanamori_66,
    _build_site_bilinear,
    _extract_j_matrix,
    _kramers_basis_66,
    _single_site_operators,
)

ComplexArray = NDArray[np.complex128]
FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class FourthOrderExchangeResult:
    """Exact 2nd + 4th order BDL/Takahashi bond exchange and channel decomposition."""

    J2: float
    K2: float
    Gamma2: float
    Gamma_prime2: float
    J4: float
    K4: float
    Gamma4: float
    Gamma_prime4: float
    J_total: float
    K_total: float
    Gamma_total: float
    Gamma_prime_total: float
    heff_2: ComplexArray
    heff_4: ComplexArray
    heff_4_chain_j32: ComplexArray
    heff_4_renorm: ComplexArray
    chain_norm: float
    renorm_norm: float
    richardson_bdl_max_diff: float


@dataclass(frozen=True)
class HexagonRingExchangeResult:
    """Exact 6-site half-filled Kramers hexagon 6th-order Peierls ring exchange."""

    u_eff: float
    t_iso: float
    t_kitaev: float
    k_ring_wp: float
    ring_frobenius_norm: float
    resolvent_vs_bdl_max_diff: float
    heff_ring: ComplexArray


def compute_nonperturbative_bdl_heff(
    h0: ComplexArray,
    v_hop: ComplexArray,
    p_basis: ComplexArray,
    eps: float,
    *,
    active_basis: ComplexArray | None = None,
) -> ComplexArray:
    """Compute exact non-perturbative BDL direct-rotation Heff on the 4-state low-energy manifold."""
    if active_basis is not None:
        h_eff_full = active_basis.conj().T @ (h0 + eps * v_hop) @ active_basis
        p_sub = active_basis.conj().T @ p_basis
    else:
        h_eff_full = h0 + eps * v_hop
        p_sub = p_basis

    h_eff_full = 0.5 * (h_eff_full + h_eff_full.conj().T)
    evals, evecs = np.linalg.eigh(h_eff_full)
    # Select the 4 eigenstates with maximum projection onto the P0 subspace (guaranteed low-energy at small eps)
    proj_weights = np.sum(np.abs(p_sub.conj().T @ evecs) ** 2, axis=0)
    idx4 = np.argsort(proj_weights)[-4:]
    # Sort the selected 4 eigenstates canonically for deterministic alignment
    sub_evals = evals[idx4]
    sub_vecs = evecs[:, idx4]

    overlap = p_sub.conj().T @ sub_vecs
    w_mat, _, vh_mat = np.linalg.svd(overlap)
    u_align = w_mat @ vh_mat
    return np.asarray(
        u_align @ np.diag(sub_evals) @ u_align.conj().T,
        dtype=np.complex128,
    )


def extract_richardson_fourth_order(
    h0: ComplexArray,
    v_hop: ComplexArray,
    p_basis: ComplexArray,
    e0: float,
    heff_2: ComplexArray,
    *,
    eps: float = 0.25,
    active_basis: ComplexArray | None = None,
) -> ComplexArray:
    """Eliminate O(eps^6, eps^8) contamination from non-perturbative BDL via 3-point Richardson."""

    def _d4(step: float) -> ComplexArray:
        h_bdl = compute_nonperturbative_bdl_heff(
            h0, v_hop, p_basis, step, active_basis=active_basis
        )
        return np.asarray(
            (h_bdl - e0 * np.eye(4, dtype=np.complex128) - (step**2) * heff_2) / (step**4),
            dtype=np.complex128,
        )

    s_vals = (eps, 0.78 * eps, 0.60 * eps)
    xs = [s * s for s in s_vals]
    fs = [_d4(s) for s in s_vals]
    rich = np.zeros((4, 4), dtype=np.complex128)
    for j in range(3):
        weight = 1.0
        for m in range(3):
            if m != j:
                weight *= -xs[m] / (xs[j] - xs[m])
        rich += weight * fs[j]
    return np.asarray(rich, dtype=np.complex128)


def _active_pq_partition(
    U: float,
    J_H: float,
    lam_soc: float,
) -> tuple[ComplexArray, ComplexArray, ComplexArray, ComplexArray, float]:
    """Return `(h0, p_basis, q_basis, active_basis, e0)` for finite-SOC (66) or infinite-SOC (34)."""
    v_kan, n0, n1 = _build_kanamori_66(U, J_H)
    p_basis = _kramers_basis_66()

    if lam_soc > 0.0:
        _, _, _, h_so_unit6 = _single_site_operators()
        h_so_66 = (
            _build_site_bilinear(h_so_unit6, 0) + _build_site_bilinear(h_so_unit6, 1)
        ).toarray()
        h0 = np.asarray(v_kan + lam_soc * h_so_66, dtype=np.complex128)
        e0 = float(np.real(p_basis[:, 0].conj().T @ h0 @ p_basis[:, 0]))
        comp_proj = np.eye(66, dtype=np.complex128) - p_basis @ p_basis.conj().T
        w, v = np.linalg.eigh(0.5 * (comp_proj + comp_proj.conj().T))
        q_basis = np.asarray(v[:, w > 0.5], dtype=np.complex128)
        active_basis = np.eye(66, dtype=np.complex128)
        return h0, p_basis, q_basis, active_basis, e0

    # Strict lambda_SOC -> infinity limit: j_eff=3/2 states are projected out (30 charge + 4 P0)
    n0_diag = np.real(np.diag(n0.toarray()))
    n1_diag = np.real(np.diag(n1.toarray()))
    q_idx = np.where(np.isclose(n0_diag, 2.0) | np.isclose(n1_diag, 2.0))[0]
    q_basis = np.eye(66, dtype=np.complex128)[:, q_idx]
    active_basis = np.hstack([p_basis, q_basis])
    h0 = np.asarray(v_kan, dtype=np.complex128)
    e0 = float(np.real(p_basis[:, 0].conj().T @ h0 @ p_basis[:, 0]))
    return h0, p_basis, q_basis, active_basis, e0


def compute_fourth_order_fock_exchange(
    bond: str,
    U: float,
    J_H: float,
    lam_soc: float,
    t2: float,
    *,
    t1: float = 0.0,
    t3: float = 0.0,
    t4: float = 0.0,
    convention: str = "spin_12",
    bdl_scale: float = 0.042,
) -> FourthOrderExchangeResult:
    """Compute exact 2nd and 4th order BDL/Takahashi exchange tensors and check Richardson BDL."""
    h0, p_basis, q_basis, active_basis, e0 = _active_pq_partition(U, J_H, lam_soc)
    v_hop = _build_hopping_66(bond, t1=t1, t2=t2, t3=t3, t4=t4)

    hq = q_basis.conj().T @ (h0 - e0 * np.eye(66, dtype=np.complex128)) @ q_basis
    r_q = np.linalg.inv(hq)
    r_q2 = r_q @ r_q

    v_pq = p_basis.conj().T @ v_hop @ q_basis
    v_qq = q_basis.conj().T @ v_hop @ q_basis
    v_qp = v_pq.conj().T

    heff_2 = np.asarray(-v_pq @ r_q @ v_qp, dtype=np.complex128)
    heff_4_chain = np.asarray(-v_pq @ r_q @ v_qq @ r_q @ v_qq @ r_q @ v_qp, dtype=np.complex128)
    # 4th-order wavefunction renormalization commutator +0.5 * {V_PQ R_Q^2 V_QP, V_PQ R_Q V_QP}
    w_op = v_pq @ r_q2 @ v_qp
    m_op = v_pq @ r_q @ v_qp
    heff_4_renorm = np.asarray(0.5 * (w_op @ m_op + m_op @ w_op), dtype=np.complex128)
    heff_4 = np.asarray(heff_4_chain + heff_4_renorm, dtype=np.complex128)

    # Scale dimensionless probe amplitude to keep |eps * V_PQ| / min_gap ~ bdl_scale, avoiding float64 eps^4 underflow
    hq_evals = np.linalg.eigvalsh(0.5 * (hq + hq.conj().T))
    # Geometric scale across V_PQ charge gap and full Q complement avoids underflow at small lam_SOC
    pq_weights = np.sum(np.abs(v_pq) ** 2, axis=0)
    eff_gap = float(np.sqrt(np.min(hq_evals) * np.max(hq_evals)))
    if float(np.sum(pq_weights)) > 0.0:
        eff_gap = float(np.sqrt(eff_gap * float(np.median(hq_evals))))
    v_scale = max(float(np.linalg.norm(v_pq, ord=2)), 1e-12)
    eps_probe = bdl_scale * eff_gap / v_scale
    d4_rich = extract_richardson_fourth_order(
        h0, v_hop, p_basis, e0, heff_2, eps=eps_probe, active_basis=active_basis
    )
    diff_bdl = float(np.max(np.abs(heff_4 - d4_rich)))

    j2, k2, g2, gp2, _ = _extract_j_matrix(heff_2, bond=bond, convention=convention)
    j4, k4, g4, gp4, _ = _extract_j_matrix(heff_4, bond=bond, convention=convention)
    jt, kt, gt, gpt, _ = _extract_j_matrix(heff_2 + heff_4, bond=bond, convention=convention)

    return FourthOrderExchangeResult(
        J2=j2,
        K2=k2,
        Gamma2=g2,
        Gamma_prime2=gp2,
        J4=j4,
        K4=k4,
        Gamma4=g4,
        Gamma_prime4=gp4,
        J_total=jt,
        K_total=kt,
        Gamma_total=gt,
        Gamma_prime_total=gpt,
        heff_2=heff_2,
        heff_4=heff_4,
        heff_4_chain_j32=heff_4_chain,
        heff_4_renorm=heff_4_renorm,
        chain_norm=float(np.linalg.norm(heff_4_chain)),
        renorm_norm=float(np.linalg.norm(heff_4_renorm)),
        richardson_bdl_max_diff=diff_bdl,
    )


def _hexagon_fock_basis() -> tuple[
    list[tuple[int, ...]],
    dict[tuple[int, ...], int],
    list[int],
    list[int],
    FloatArray,
]:
    """Enumerate all C(12, 6) = 924 6-site half-filled Kramers Fock states."""
    states = list(combinations(range(12), 6))
    state_to_idx = {st: k for k, st in enumerate(states)}
    p_indices: list[int] = []
    q_indices: list[int] = []
    diag_u = np.zeros(len(states), dtype=np.float64)

    for idx, st in enumerate(states):
        st_set = set(st)
        doubly = 0
        singly = True
        for site in range(6):
            occ = (1 if 2 * site in st_set else 0) + (1 if (2 * site + 1) in st_set else 0)
            if occ == 2:
                doubly += 1
            if occ != 1:
                singly = False
        diag_u[idx] = float(doubly)
        if singly:
            p_indices.append(idx)
        else:
            q_indices.append(idx)
    return states, state_to_idx, p_indices, q_indices, diag_u


def _act_single_hop_6site(
    st: tuple[int, ...], c_mode: int, a_mode: int, state_to_idx: dict[tuple[int, ...], int]
) -> tuple[int, int] | None:
    if a_mode not in st:
        return None
    p = list(st)
    pos = p.index(a_mode)
    sgn = -1 if (pos % 2) else 1
    p.pop(pos)
    if c_mode in p:
        return None
    ins = 0
    while ins < len(p) and p[ins] < c_mode:
        ins += 1
    if ins % 2:
        sgn = -sgn
    p.insert(ins, c_mode)
    return state_to_idx[tuple(p)], sgn


def _build_hexagon_hopping(
    states: list[tuple[int, ...]],
    state_to_idx: dict[tuple[int, ...], int],
    t_iso: float,
    t_kitaev: float,
    *,
    peierls_flux_pi: bool = False,
) -> ComplexArray:
    n = len(states)
    mat = np.zeros((n, n), dtype=np.complex128)
    paulis = (
        np.array([[0, 1], [1, 0]], dtype=np.complex128),
        np.array([[0, -1j], [1j, 0]], dtype=np.complex128),
        np.array([[1, 0], [0, -1]], dtype=np.complex128),
    )
    # 6 honeycomb plaquette edges in Kitaev bond cycle (x, y, z, x, y, z)
    bond_paulis = (0, 1, 2, 0, 1, 2)
    i2 = np.eye(2, dtype=np.complex128)

    for b_idx in range(6):
        s_i = b_idx
        s_j = (b_idx + 1) % 6
        phase = -1.0 if (peierls_flux_pi and b_idx == 5) else 1.0
        t_mat = phase * (t_iso * i2 + 1j * t_kitaev * paulis[bond_paulis[b_idx]])
        for k, st in enumerate(states):
            for r in range(2):
                for c in range(2):
                    val = t_mat[r, c]
                    if abs(val) < 1e-15:
                        continue
                    res1 = _act_single_hop_6site(st, 2 * s_i + r, 2 * s_j + c, state_to_idx)
                    if res1 is not None:
                        mat[res1[0], k] -= val * res1[1]
                    res2 = _act_single_hop_6site(st, 2 * s_j + c, 2 * s_i + r, state_to_idx)
                    if res2 is not None:
                        mat[res2[0], k] -= np.conj(val) * res2[1]
    return 0.5 * (mat + mat.conj().T)


def compute_hexagon_ring_exchange(
    u_eff: float,
    t_iso: float,
    *,
    t_kitaev: float = 0.0,
    bdl_scale: float = 0.08,
) -> HexagonRingExchangeResult:
    """Compute exact 6th-order 6-site hexagon ring exchange on 924-state Kramers Fock space."""
    states, state_to_idx, p_idx, q_idx, diag_u_unit = _hexagon_fock_basis()
    n_tot = len(states)
    h0 = np.diag(u_eff * diag_u_unit).astype(np.complex128)
    p_basis = np.eye(n_tot, dtype=np.complex128)[:, p_idx]
    r_q_vec = 1.0 / (u_eff * diag_u_unit[q_idx])

    v0 = _build_hexagon_hopping(states, state_to_idx, t_iso, t_kitaev, peierls_flux_pi=False)
    vpi = _build_hexagon_hopping(states, state_to_idx, t_iso, t_kitaev, peierls_flux_pi=True)

    # Exact 6th-order winding loop: -V_PQ (R_Q V_QQ)^4 R_Q V_QP (lower-order/renorm terms cancel in Phi=0 - Phi=pi)
    def _loop6(v_op: ComplexArray) -> ComplexArray:
        v_pq = v_op[np.ix_(p_idx, q_idx)]
        v_qq = v_op[np.ix_(q_idx, q_idx)]
        v_qp = v_op[np.ix_(q_idx, p_idx)]
        cur = (r_q_vec[:, None]) * v_qp
        for _ in range(4):
            cur = (r_q_vec[:, None]) * (v_qq @ cur)
        return np.asarray(-v_pq @ cur, dtype=np.complex128)

    heff_ring_resolvent = np.asarray(0.5 * (_loop6(v0) - _loop6(vpi)), dtype=np.complex128)

    # Rescale BDL probe dimensionless amplitude so (alpha * t_scale / u_eff) ~ bdl_scale,
    # keeping (alpha*t/U)^6 well above float64 roundoff while 3-point Richardson removes O(a^2, a^4).
    t_scale = max(float(np.hypot(t_iso, t_kitaev)), 1e-12)
    alpha_ref = bdl_scale * u_eff / t_scale

    def _bdl_flux_diff(step: float) -> ComplexArray:
        def _bdl(v_op: ComplexArray) -> ComplexArray:
            h_tot = h0 + step * v_op
            evals, evecs = np.linalg.eigh(h_tot)
            sub_e = evals[:64]
            sub_v = evecs[:, :64]
            ov = p_basis.conj().T @ sub_v
            w_m, _, vh_m = np.linalg.svd(ov)
            u_m = w_m @ vh_m
            return np.asarray(u_m @ np.diag(sub_e) @ u_m.conj().T, dtype=np.complex128)

        return np.asarray(0.5 * (_bdl(v0) - _bdl(vpi)) / (step**6), dtype=np.complex128)

    # Three even-power points in x = step^2 for 3-point Lagrange extrapolation to x = 0
    s_vals = (alpha_ref, 0.8 * alpha_ref, 0.65 * alpha_ref)
    xs = [s * s for s in s_vals]
    fs = [_bdl_flux_diff(s) for s in s_vals]
    d6_rich = np.zeros((64, 64), dtype=np.complex128)
    for j in range(3):
        weight = 1.0
        for m in range(3):
            if m != j:
                weight *= -xs[m] / (xs[j] - xs[m])
        d6_rich += weight * fs[j]
    diff_bdl = float(np.max(np.abs(heff_ring_resolvent - d6_rich)))

    # Project onto 6-spin Kitaev hexagon flux W_p = sigma_0^x sigma_1^y sigma_2^z sigma_3^x sigma_4^y sigma_5^z
    # Note: ordering inside p_idx matches single-hole mode occupation: (up=0, dn=1) at each site 0..5.
    # For each singly-occupied state, extract 6-bit spin configuration in site order 0..5.
    paulis = (
        np.array([[0, 1], [1, 0]], dtype=np.complex128),
        np.array([[0, -1j], [1j, 0]], dtype=np.complex128),
        np.array([[1, 0], [0, -1]], dtype=np.complex128),
    )
    # Construct permutation/diagonal isomorphism from p_idx ordering into tensor product of 6 qubits
    bits_map = np.zeros(64, dtype=np.int64)
    for k, fock_idx in enumerate(p_idx):
        st = states[fock_idx]
        code = 0
        for site in range(6):
            spin = st[site] - 2 * site
            code = (code << 1) | spin
        bits_map[k] = code

    wp_tensor = np.array([[1.0 + 0.0j]], dtype=np.complex128)
    for g in (0, 1, 2, 0, 1, 2):
        wp_tensor = np.asarray(np.kron(wp_tensor, paulis[g]), dtype=np.complex128)
    wp_in_p = wp_tensor[np.ix_(bits_map, bits_map)]

    k_ring_wp = float(np.real(np.trace(wp_in_p @ heff_ring_resolvent)) / 64.0)
    return HexagonRingExchangeResult(
        u_eff=u_eff,
        t_iso=t_iso,
        t_kitaev=t_kitaev,
        k_ring_wp=k_ring_wp,
        ring_frobenius_norm=float(np.linalg.norm(heff_ring_resolvent)),
        resolvent_vs_bdl_max_diff=diff_bdl,
        heff_ring=heff_ring_resolvent,
    )
