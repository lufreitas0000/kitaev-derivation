"""Tests for Wave 3 M2.2 fourth-order BDL/Takahashi exchange and 6-site hexagon ring exchange."""

from __future__ import annotations

import numpy as np
import pytest

from kitaev_derivation.services.fourth_order_exchange import (
    compute_fourth_order_fock_exchange,
    compute_hexagon_ring_exchange,
)


def test_fourth_order_exchange_matches_richardson_bdl_both_soc_limits() -> None:
    # 1. Strict lambda_SOC -> infinity limit (lam_soc=0.0): j_eff=3/2 states projected out
    res_inf = compute_fourth_order_fock_exchange(
        "z", U=3.0, J_H=0.48, lam_soc=0.0, t2=0.16, t1=0.05, t3=-0.19, t4=-0.02
    )
    assert res_inf.richardson_bdl_max_diff < 5e-10
    assert res_inf.chain_norm < 1e-13
    assert res_inf.renorm_norm > 0.0

    # 2. Finite lambda_SOC: includes 32 neutral j_eff=3/2 intermediates enabling O(t^4 / U^2 lam) channel
    res_fin = compute_fourth_order_fock_exchange(
        "z", U=3.0, J_H=0.48, lam_soc=0.25, t2=0.16, t1=0.05, t3=-0.19, t4=-0.02
    )
    assert res_fin.richardson_bdl_max_diff < 5e-10
    assert res_fin.chain_norm > 1e-6


def test_fourth_order_exact_analytical_hubbard_limit() -> None:
    # Single-orbital J_H=0, t1 != 0 only gives two decoupled 1-band chains? On z-bond, t1 acts on dxz,dyz:
    # Check exact O(t1^2/U) and O(t1^4/U^3) via explicit scaling ratios
    u_val, t1_val = 4.0, 0.2
    res_a = compute_fourth_order_fock_exchange(
        "z", U=u_val, J_H=0.0, lam_soc=0.0, t2=0.0, t1=t1_val
    )
    res_b = compute_fourth_order_fock_exchange(
        "z", U=u_val, J_H=0.0, lam_soc=0.0, t2=0.0, t1=0.5 * t1_val
    )
    assert abs(res_a.J2 / res_b.J2 - 4.0) < 1e-12
    assert abs(res_a.J4 / res_b.J4 - 16.0) < 1e-12
    assert res_a.richardson_bdl_max_diff < 5e-10


def test_fourth_order_exchange_c3_bond_covariance() -> None:
    res_x = compute_fourth_order_fock_exchange(
        "x", U=3.5, J_H=0.5, lam_soc=0.3, t2=0.20, t1=0.04, t3=-0.14, t4=-0.015
    )
    res_y = compute_fourth_order_fock_exchange(
        "y", U=3.5, J_H=0.5, lam_soc=0.3, t2=0.20, t1=0.04, t3=-0.14, t4=-0.015
    )
    res_z = compute_fourth_order_fock_exchange(
        "z", U=3.5, J_H=0.5, lam_soc=0.3, t2=0.20, t1=0.04, t3=-0.14, t4=-0.015
    )
    for field in ("J4", "K4", "Gamma4", "Gamma_prime4", "J_total", "K_total"):
        vx = getattr(res_x, field)
        vy = getattr(res_y, field)
        vz = getattr(res_z, field)
        assert abs(vx - vz) < 1e-12
        assert abs(vy - vz) < 1e-12


def test_six_site_hexagon_ring_exchange_matches_bdl_peierls_oracle() -> None:
    res = compute_hexagon_ring_exchange(u_eff=3.0, t_iso=0.18, t_kitaev=0.09)
    # 924x924 IEEE float64 eigh conditioning floor at O(s^6 ~ 1e-7) gives < 2e-7 abs (< 2e-3 rel)
    assert res.resolvent_vs_bdl_max_diff < 2e-7
    assert res.ring_frobenius_norm > 0.0
    assert np.isfinite(res.k_ring_wp)
    # Scaling check: doubling hopping amplitude scales 6th-order ring norm by 2^6 = 64
    res_half = compute_hexagon_ring_exchange(u_eff=3.0, t_iso=0.09, t_kitaev=0.045)
    assert res.ring_frobenius_norm / res_half.ring_frobenius_norm == pytest.approx(64.0, rel=1e-12)
