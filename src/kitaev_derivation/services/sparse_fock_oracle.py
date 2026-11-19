"""Independent scipy.sparse 66-state two-hole Fock-space oracle for cross-checking symbolic CAS."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy import sparse

from kitaev_derivation.domain.conventions import ORBITALS, SPINS

ComplexArray = NDArray[np.complex128]
FloatArray = NDArray[np.float64]

_BASIS_PAIRS: tuple[tuple[int, int], ...] = tuple(
    (a, b) for a in range(12) for b in range(a + 1, 12)
)
_PAIR_TO_IDX: dict[tuple[int, int], int] = {p: k for k, p in enumerate(_BASIS_PAIRS)}


@dataclass(frozen=True)
class FockOracleResult:
    """Independent numerical Fock-space exchange couplings on a single bond."""

    J: float
    K: float
    Gamma: float
    Gamma_prime: float
    J_matrix: FloatArray
    e0_ground: float
    min_d4_gap: float


def _act_hop(state_idx: int, c_mode: int, a_mode: int) -> tuple[int, int] | None:
    p = list(_BASIS_PAIRS[state_idx])
    if a_mode not in p:
        return None
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
    return _PAIR_TO_IDX[(p[0], p[1])], sgn


def _build_site_bilinear(op6: ComplexArray, site: int) -> sparse.csc_matrix:
    rows: list[int] = []
    cols: list[int] = []
    vals: list[complex] = []
    offset = site * 6
    for k in range(66):
        for r in range(6):
            for c in range(6):
                v = op6[r, c]
                if abs(v) < 1e-15:
                    continue
                res = _act_hop(k, offset + r, offset + c)
                if res is not None:
                    k2, sgn = res
                    rows.append(k2)
                    cols.append(k)
                    vals.append(v * sgn)
    return sparse.csc_matrix((vals, (rows, cols)), shape=(66, 66), dtype=np.complex128)


def _single_site_operators() -> tuple[
    ComplexArray, list[ComplexArray], list[ComplexArray], ComplexArray
]:
    lx = np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], dtype=np.complex128)
    ly = np.array([[0, 0, 1j], [0, 0, 0], [-1j, 0, 0]], dtype=np.complex128)
    lz = np.array([[0, -1j, 0], [1j, 0, 0], [0, 0, 0]], dtype=np.complex128)
    sx = 0.5 * np.array([[0, 1], [1, 0]], dtype=np.complex128)
    sy = 0.5 * np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
    sz = 0.5 * np.array([[1, 0], [0, -1]], dtype=np.complex128)

    n6: ComplexArray = np.eye(6, dtype=np.complex128)
    s6: list[ComplexArray] = [
        np.asarray(np.kron(s, np.eye(3)), dtype=np.complex128) for s in (sx, sy, sz)
    ]
    l_eff6: list[ComplexArray] = [
        np.asarray(np.kron(np.eye(2), l_mat), dtype=np.complex128) for l_mat in (lx, ly, lz)
    ]
    l_hole6: list[ComplexArray] = [-x for x in l_eff6]
    h_so_unit6: ComplexArray = np.zeros((6, 6), dtype=np.complex128)
    for s, l_mat in zip(s6, l_eff6, strict=True):
        h_so_unit6 = h_so_unit6 + (s @ l_mat)
    return n6, s6, l_hole6, h_so_unit6


def _build_kanamori_66(
    U: float, J_H: float
) -> tuple[ComplexArray, sparse.csc_matrix, sparse.csc_matrix]:
    n6, s6, l_hole6, _ = _single_site_operators()
    eye66 = sparse.identity(66, dtype=np.complex128, format="csc")
    n0 = _build_site_bilinear(n6, 0)
    n1 = _build_site_bilinear(n6, 1)

    v_tot = sparse.csc_matrix((66, 66), dtype=np.complex128)
    for site, ni in ((0, n0), (1, n1)):
        v_tot = v_tot + 0.5 * (U - 3.0 * J_H) * (ni - eye66) @ (ni - eye66)
        for smat in s6:
            si = _build_site_bilinear(smat, site)
            v_tot = v_tot - 2.0 * J_H * (si @ si)
        for lmat in l_hole6:
            li = _build_site_bilinear(lmat, site)
            v_tot = v_tot - 0.5 * J_H * (li @ li)
    return v_tot.toarray(), n0, n1


def _build_hopping_66(bond: str, t1: float, t2: float, t3: float, t4: float) -> ComplexArray:
    if bond == "z":
        t_orb = np.array([[t1, t2, t4], [t2, t1, t4], [t4, t4, t3]], dtype=np.complex128)
    elif bond == "x":
        t_orb = np.array([[t3, t4, t4], [t4, t1, t2], [t4, t2, t1]], dtype=np.complex128)
    elif bond == "y":
        t_orb = np.array([[t1, t4, t2], [t4, t3, t4], [t2, t4, t1]], dtype=np.complex128)
    else:
        raise ValueError(f"Unsupported bond '{bond}'.")

    t6 = np.kron(np.eye(2), t_orb)
    rows: list[int] = []
    cols: list[int] = []
    vals: list[complex] = []
    for k in range(66):
        for r in range(6):
            for c in range(6):
                v = t6[r, c]
                if abs(v) < 1e-15:
                    continue
                res1 = _act_hop(k, r, 6 + c)
                if res1 is not None:
                    rows.append(res1[0])
                    cols.append(k)
                    vals.append(-v * res1[1])
                res2 = _act_hop(k, 6 + c, r)
                if res2 is not None:
                    rows.append(res2[0])
                    cols.append(k)
                    vals.append(-np.conj(v) * res2[1])
    mat = sparse.csc_matrix((vals, (rows, cols)), shape=(66, 66), dtype=np.complex128).toarray()
    return np.asarray(mat, dtype=np.complex128)


def _kramers_basis_66() -> ComplexArray:
    # Single-hole mode idx: spin_idx*3 + orb_idx where SPINS=('up','dn'), ORBITALS=('x','y','z')
    def m(orb: str, spin: str) -> int:
        return SPINS.index(spin) * 3 + ORBITALS.index(orb)

    plus6 = np.zeros(6, dtype=np.complex128)
    plus6[m("z", "up")] = -1.0 / np.sqrt(3.0)
    plus6[m("y", "dn")] = -1j / np.sqrt(3.0)
    plus6[m("x", "dn")] = -1.0 / np.sqrt(3.0)

    minus6 = np.zeros(6, dtype=np.complex128)
    minus6[m("z", "dn")] = 1.0 / np.sqrt(3.0)
    minus6[m("y", "up")] = 1j / np.sqrt(3.0)
    minus6[m("x", "up")] = -1.0 / np.sqrt(3.0)

    kramers = (plus6, minus6)
    p_basis = np.zeros((66, 4), dtype=np.complex128)
    idx4 = 0
    for s0 in range(2):
        for s1 in range(2):
            for r in range(6):
                for c in range(6):
                    p_basis[_PAIR_TO_IDX[(r, 6 + c)], idx4] += kramers[s0][r] * kramers[s1][c]
            idx4 += 1
    return p_basis


def _extract_j_matrix(
    heff_4x4: ComplexArray, bond: str, convention: str
) -> tuple[float, float, float, float, FloatArray]:
    paulis = (
        np.array([[0, 1], [1, 0]], dtype=np.complex128),
        np.array([[0, -1j], [1j, 0]], dtype=np.complex128),
        np.array([[1, 0], [0, -1]], dtype=np.complex128),
    )
    scale = 4.0 if convention == "spin_12" else 1.0
    ops = [0.5 * p for p in paulis]
    j_mat = np.zeros((3, 3), dtype=np.float64)
    for a in range(3):
        for b in range(3):
            op_ab = np.kron(ops[a], ops[b])
            j_mat[a, b] = float(np.real(scale * np.trace(heff_4x4 @ op_ab)))

    if bond == "z":
        J = j_mat[0, 0]
        K = j_mat[2, 2] - J
        Gamma = 0.5 * (j_mat[0, 1] + j_mat[1, 0])
        Gamma_prime = 0.25 * (j_mat[0, 2] + j_mat[2, 0] + j_mat[1, 2] + j_mat[2, 1])
    elif bond == "x":
        J = j_mat[1, 1]
        K = j_mat[0, 0] - J
        Gamma = 0.5 * (j_mat[1, 2] + j_mat[2, 1])
        Gamma_prime = 0.25 * (j_mat[0, 1] + j_mat[1, 0] + j_mat[0, 2] + j_mat[2, 0])
    else:
        J = j_mat[0, 0]
        K = j_mat[1, 1] - J
        Gamma = 0.5 * (j_mat[0, 2] + j_mat[2, 0])
        Gamma_prime = 0.25 * (j_mat[0, 1] + j_mat[1, 0] + j_mat[1, 2] + j_mat[2, 1])
    return float(J), float(K), float(Gamma), float(Gamma_prime), j_mat


def compute_fock_kanamori_spectrum(U: float, J_H: float) -> tuple[float, list[tuple[float, int]]]:
    """Compute exact d5-d5 ground energy and unique d4-d6 excitation gaps + multiplicities."""
    v_tot, n0, n1 = _build_kanamori_66(U, J_H)
    d5_mask = np.isclose(np.diag(n0.toarray()), 1.0) & np.isclose(np.diag(n1.toarray()), 1.0)
    d4_0_mask = np.isclose(np.diag(n0.toarray()), 2.0)
    e0 = float(np.linalg.eigvalsh(v_tot[np.ix_(d5_mask, d5_mask)])[0])
    e_d4 = np.linalg.eigvalsh(v_tot[np.ix_(d4_0_mask, d4_0_mask)]) - e0
    vals, counts = np.unique(np.round(e_d4, 10), return_counts=True)
    return e0, [(float(v), int(c)) for v, c in zip(vals, counts, strict=True)]


def compute_fock_schrieffer_wolff_exchange(
    bond: str,
    U: float,
    J_H: float,
    t2: float,
    *,
    t1: float = 0.0,
    t3: float = 0.0,
    t4: float = 0.0,
    convention: str = "spin_12",
) -> FockOracleResult:
    """Evaluate exact 66-state second-order Schrieffer-Wolff projection into j_eff=1/2."""
    v_tot, n0, n1 = _build_kanamori_66(U, J_H)
    n0_diag = np.real(np.diag(n0.toarray()))
    n1_diag = np.real(np.diag(n1.toarray()))
    d5_mask = np.isclose(n0_diag, 1.0) & np.isclose(n1_diag, 1.0)
    q_mask = np.isclose(n0_diag, 2.0) | np.isclose(n1_diag, 2.0)

    e0 = float(np.linalg.eigvalsh(v_tot[np.ix_(d5_mask, d5_mask)])[0])
    vq = v_tot[np.ix_(q_mask, q_mask)] - e0 * np.eye(int(np.sum(q_mask)))
    min_gap = float(np.min(np.linalg.eigvalsh(vq)))
    vq_inv = np.linalg.inv(vq)

    ht = _build_hopping_66(bond, t1=t1, t2=t2, t3=t3, t4=t4)
    p_basis = _kramers_basis_66()
    heff_4x4 = -p_basis.conj().T @ ht[:, q_mask] @ vq_inv @ ht[q_mask, :] @ p_basis

    J, K, Gamma, Gamma_prime, j_mat = _extract_j_matrix(heff_4x4, bond=bond, convention=convention)
    return FockOracleResult(
        J=J,
        K=K,
        Gamma=Gamma,
        Gamma_prime=Gamma_prime,
        J_matrix=j_mat,
        e0_ground=e0,
        min_d4_gap=min_gap,
    )


def compute_fock_soc_finite_exchange(
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
) -> FockOracleResult:
    """Evaluate 66-state Schrieffer-Wolff WITH finite single-site SOC in intermediate states.

    Quantifies finite-lambda_SOC departures from the analytical strong-coupling limit
    lambda_SOC / (U - 3J_H) -> 0 in d^4 intermediate denominators (while retaining
    lambda_SOC >> t^2 / U so external states stay projected onto j_eff = 1/2).
    """
    v_kan, n0, n1 = _build_kanamori_66(U, J_H)
    _, _, _, h_so_unit6 = _single_site_operators()
    h_so_66 = (_build_site_bilinear(h_so_unit6, 0) + _build_site_bilinear(h_so_unit6, 1)).toarray()
    h0 = v_kan + lam_soc * h_so_66

    n0_diag = np.real(np.diag(n0.toarray()))
    n1_diag = np.real(np.diag(n1.toarray()))
    d5_mask = np.isclose(n0_diag, 1.0) & np.isclose(n1_diag, 1.0)
    q_mask = np.isclose(n0_diag, 2.0) | np.isclose(n1_diag, 2.0)

    p_basis = _kramers_basis_66()
    e0 = float(np.min(np.linalg.eigvalsh(h0[np.ix_(d5_mask, d5_mask)])))
    hq = h0[np.ix_(q_mask, q_mask)] - e0 * np.eye(int(np.sum(q_mask)))
    min_gap = float(np.min(np.linalg.eigvalsh(hq)))
    hq_inv = np.linalg.inv(hq)

    ht = _build_hopping_66(bond, t1=t1, t2=t2, t3=t3, t4=t4)
    heff_4x4 = -p_basis.conj().T @ ht[:, q_mask] @ hq_inv @ ht[q_mask, :] @ p_basis

    J, K, Gamma, Gamma_prime, j_mat = _extract_j_matrix(heff_4x4, bond=bond, convention=convention)
    return FockOracleResult(
        J=J,
        K=K,
        Gamma=Gamma,
        Gamma_prime=Gamma_prime,
        J_matrix=j_mat,
        e0_ground=e0,
        min_d4_gap=min_gap,
    )
