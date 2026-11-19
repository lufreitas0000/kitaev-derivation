"""Exact symbolic Kramers doublet projection and Pauli-trace exchange extraction."""

from __future__ import annotations

import sympy as sp
from sympy import Expr, I, Symbol

from kitaev_derivation.domain.algebra import (
    AbstractQuantumOperator,
    FermionA,
    FermionC,
    IdentityOp,
    MultipletProjector,
    OperatorString,
    OperatorSum,
    ParticleProjector,
    ScaledOperator,
    SpinComposite,
)
from kitaev_derivation.domain.conventions import (
    ORBITALS,
    SPINS,
    kramers_doublet_kets,
    pauli_matrix_element,
    spin12_to_pauli_coupling,
)
from kitaev_derivation.domain.kanamori_soc import d4_multiplet_basis
from kitaev_derivation.services.schrieffer_wolff import compute_bond_effective_hamiltonian

ModeTuple = tuple[str, str, str]  # (site, orbital, spin)
FockKet = dict[tuple[ModeTuple, ...], Expr]

_D4_BASIS = d4_multiplet_basis()
_KRAMERS = kramers_doublet_kets()


def _mode_order_key(mode: ModeTuple) -> tuple[str, int, int]:
    site, orb, spin = mode
    return (site, ORBITALS.index(orb), SPINS.index(spin))


def _canonicalize_kets(modes: list[ModeTuple]) -> tuple[tuple[ModeTuple, ...], int]:
    """Bubble-sort fermionic modes into canonical order and track parity sign."""
    arr = list(modes)
    sign = 1
    n = len(arr)
    for i in range(n):
        for j in range(n - 1 - i):
            if _mode_order_key(arr[j + 1]) < _mode_order_key(arr[j]):
                arr[j], arr[j + 1] = arr[j + 1], arr[j]
                sign *= -1
    return tuple(arr), sign


def _apply_fermion_c(state: FockKet, mode: ModeTuple) -> FockKet:
    out: FockKet = {}
    for kets, coeff in state.items():
        if mode in kets:
            continue
        num_less = sum(1 for m in kets if _mode_order_key(m) < _mode_order_key(mode))
        sign = -1 if (num_less % 2) else 1
        new_kets, sort_sgn = _canonicalize_kets([*kets, mode])
        # Inserting in sorted order avoids double-counting: num_less counts anticommutations past earlier kets
        del sort_sgn
        new_kets = tuple(sorted([*kets, mode], key=_mode_order_key))
        val = out.get(new_kets, sp.Integer(0)) + coeff * sign
        if val != 0:
            out[new_kets] = val
    return out


def _apply_fermion_a(state: FockKet, mode: ModeTuple) -> FockKet:
    out: FockKet = {}
    for kets, coeff in state.items():
        if mode not in kets:
            continue
        idx = kets.index(mode)
        sign = -1 if (idx % 2) else 1
        new_kets = kets[:idx] + kets[idx + 1 :]
        val = out.get(new_kets, sp.Integer(0)) + coeff * sign
        if val != 0:
            out[new_kets] = val
    return out


def _apply_particle_projector(state: FockKet, proj: ParticleProjector) -> FockKet:
    return {
        kets: coeff
        for kets, coeff in state.items()
        if sum(1 for m in kets if m[0] == proj.site) == proj.holes
    }


def _apply_multiplet_projector(state: FockKet, proj: MultipletProjector) -> FockKet:
    out: FockKet = {}
    basis_states = _D4_BASIS[proj.L_label]
    for kets, coeff in state.items():
        site_modes = [m for m in kets if m[0] == proj.site]
        if len(site_modes) != 2:
            continue
        # Ensure canonical (orbital, spin) pair order matching _ordered_pair
        m1 = (site_modes[0][1], site_modes[0][2])
        m2 = (site_modes[1][1], site_modes[1][2])
        if m1 < m2:
            pair, pair_sgn = (m1, m2), 1
        else:
            pair, pair_sgn = (m2, m1), -1

        other_modes = [m for m in kets if m[0] != proj.site]
        for bstate in basis_states:
            overlap = sp.conjugate(bstate.get(pair, sp.Integer(0))) * pair_sgn
            if overlap == 0:
                continue
            for (bp1, bp2), b_coeff in bstate.items():
                idx_s = [i for i, m in enumerate(kets) if m[0] == proj.site]
                rem_sgn = -1 if ((idx_s[0] + idx_s[1] - 1) % 2) else 1
                new_kets, ins_sgn = _canonicalize_kets(
                    [(proj.site, bp1[0], bp1[1]), (proj.site, bp2[0], bp2[1]), *other_modes]
                )
                total_coeff = coeff * rem_sgn * overlap * b_coeff * ins_sgn
                out[new_kets] = out.get(new_kets, sp.Integer(0)) + total_coeff

    return {k: sp.expand(v) for k, v in out.items() if sp.expand(v) != 0}


def evaluate_operator_string_on_ket(op_string: OperatorString, ket: FockKet) -> FockKet:
    """Apply an OperatorString right-to-left on an exact symbolic Fock state."""
    state = ket
    for factor in reversed(op_string.factors):
        if isinstance(factor, IdentityOp):
            continue
        if isinstance(factor, FermionC):
            state = _apply_fermion_c(state, (factor.site, factor.orbital, factor.spin))
        elif isinstance(factor, FermionA):
            state = _apply_fermion_a(state, (factor.site, factor.orbital, factor.spin))
        elif isinstance(factor, ParticleProjector):
            state = _apply_particle_projector(state, factor)
        elif isinstance(factor, MultipletProjector):
            state = _apply_multiplet_projector(state, factor)
        else:
            raise TypeError(f"Cannot evaluate factor {type(factor)} on Fock ket")
        if not state:
            return {}
    return state


def pseudospin_ket(site_i: str, spin_i: str, site_j: str, spin_j: str) -> FockKet:
    """Construct exact two-site Kramers j_eff=1/2 ket |spin_i, spin_j> in symbolic Fock space."""
    ket: FockKet = {}
    for (orb_i, rspin_i), c_i in _KRAMERS[spin_i].items():
        for (orb_j, rspin_j), c_j in _KRAMERS[spin_j].items():
            kets, sgn = _canonicalize_kets([(site_i, orb_i, rspin_i), (site_j, orb_j, rspin_j)])
            ket[kets] = sp.expand(ket.get(kets, sp.Integer(0)) + c_i * c_j * sgn)
    return {k: v for k, v in ket.items() if v != 0}


def _ket_inner_product(bra: FockKet, ket: FockKet) -> Expr:
    acc: Expr = sp.Integer(0)
    for kets, val in ket.items():
        if kets in bra:
            acc += sp.conjugate(bra[kets]) * val
    return sp.expand(acc)


def project_to_pseudospin(
    h_eff: OperatorSum,
    site_i: str = "i",
    site_j: str = "j",
) -> OperatorSum:
    """Project an effective two-site operator symbolically onto the 4-dim Kramers basis.

    Returns exact Pauli-basis OperatorSum sum_{mu,nu in {x,y,z}} C_{mu,nu} sigma^mu_i sigma^nu_j.
    """
    basis_kets = {(si, sj): pseudospin_ket(site_i, si, site_j, sj) for si in SPINS for sj in SPINS}

    matrix_elements: dict[tuple[str, str, str, str], Expr] = {
        (spi, spj, si, sj): sp.Integer(0)
        for spi in SPINS
        for spj in SPINS
        for si in SPINS
        for sj in SPINS
    }

    # Cache exact rational overlaps per OperatorString to avoid redundant ket evaluations.
    overlap_cache: dict[OperatorString, dict[tuple[str, str, str, str], Expr]] = {}

    for term in h_eff.terms:
        if isinstance(term, ScaledOperator):
            scalar = term.scalar
            op_str = term.op if isinstance(term.op, OperatorString) else OperatorString((term.op,))
        elif isinstance(term, OperatorString):
            scalar = sp.Integer(1)
            op_str = term
        else:
            continue

        cached = overlap_cache.get(op_str)
        if cached is None:
            cached = {}
            for si in SPINS:
                for sj in SPINS:
                    res = evaluate_operator_string_on_ket(op_str, basis_kets[(si, sj)])
                    if not res:
                        continue
                    for spi in SPINS:
                        for spj in SPINS:
                            ov = _ket_inner_product(basis_kets[(spi, spj)], res)
                            if ov != 0:
                                cached[(spi, spj, si, sj)] = ov
            overlap_cache[op_str] = cached

        for key, ov in cached.items():
            matrix_elements[key] = matrix_elements[key] + scalar * ov

    pauli_terms: list[AbstractQuantumOperator] = []
    for mu in ("x", "y", "z"):
        for nu in ("x", "y", "z"):
            coeff: Expr = sp.Integer(0)
            for si in SPINS:
                for sj in SPINS:
                    for spi in SPINS:
                        for spj in SPINS:
                            m_el = matrix_elements[(spi, spj, si, sj)]
                            if m_el == 0:
                                continue
                            # H = sum C_{mu,nu} sigma^mu sigma^nu => C_{mu,nu} = (1/4) Tr(H sigma^mu sigma^nu)
                            p_prod = pauli_matrix_element(mu, si, spi) * pauli_matrix_element(
                                nu, sj, spj
                            )
                            if p_prod != 0:
                                coeff += m_el * p_prod / 4
            coeff_simp = sp.cancel(sp.expand(coeff))
            if coeff_simp != 0:
                pauli_terms.append(
                    ScaledOperator(
                        coeff_simp,
                        OperatorString((SpinComposite(site_i, mu), SpinComposite(site_j, nu))),
                    )
                )
    return OperatorSum(tuple(pauli_terms))


def evaluate_pauli_trace(op: AbstractQuantumOperator) -> Expr:
    """Evaluate exact trace over tensor products of spin-1/2 Pauli operators."""
    if isinstance(op, ScaledOperator):
        return sp.expand(op.scalar * evaluate_pauli_trace(op.op))
    if isinstance(op, OperatorSum):
        return sp.expand(sum((evaluate_pauli_trace(t) for t in op.terms), sp.Integer(0)))
    if isinstance(op, SpinComposite):
        return sp.Integer(0)
    if not isinstance(op, OperatorString):
        return sp.Integer(0)
    if not all(isinstance(f, SpinComposite) for f in op.factors):
        return sp.Integer(0)

    by_site: dict[str, list[str]] = {}
    for f in op.factors:
        assert isinstance(f, SpinComposite)
        by_site.setdefault(f.site, []).append(f.component)

    total: Expr = sp.Integer(1)
    for comps in by_site.values():
        if len(comps) == 0:
            site_tr: Expr = sp.Integer(2)
        elif len(comps) == 2 and comps[0] == comps[1]:
            site_tr = sp.Integer(2)
        elif len(comps) == 3:
            tup = tuple(comps)
            if tup in {("x", "y", "z"), ("y", "z", "x"), ("z", "x", "y")}:
                site_tr = 2 * I
            elif tup in {("y", "x", "z"), ("x", "z", "y"), ("z", "y", "x")}:
                site_tr = -2 * I
            else:
                return sp.Integer(0)
        else:
            return sp.Integer(0)
        total *= site_tr
    return sp.expand(total)


def extract_exchange_tensors(
    h_pseudospin: OperatorSum,
    bond: str,
    site_i: str = "i",
    site_j: str = "j",
    convention: str = "spin_12",
) -> dict[str, Expr]:
    """Extract (J, K, Gamma, Gamma_prime, J_xx..J_zz) from a Pauli-expanded pseudospin Hamiltonian.

    Args:
        convention: 'spin_12' for H = sum J_{ab} S_i^a S_j^b (S = sigma/2, J_{ab} = 4 * C_{ab}),
                    or 'pauli' for H = sum K_{ab} sigma_i^a sigma_j^b (aligned with qsl_conventions.json).
    """

    def _coeff(a: str, b: str) -> Expr:
        probe = OperatorString((SpinComposite(site_i, a), SpinComposite(site_j, b)))
        tr_val = sp.Integer(0)
        for term in h_pseudospin.terms:
            prod = term * probe
            tr_val = tr_val + evaluate_pauli_trace(prod)
        c_pauli = sp.cancel(sp.re(tr_val / 4))
        if convention == "spin_12":
            return sp.cancel(4 * c_pauli)
        if convention == "pauli":
            return sp.cancel(c_pauli)
        raise ValueError(f"Unknown convention '{convention}'; expected 'spin_12' or 'pauli'.")

    j_mat: dict[tuple[str, str], Expr] = {
        (a, b): _coeff(a, b) for a in ("x", "y", "z") for b in ("x", "y", "z")
    }

    if bond == "z":
        J = j_mat[("x", "x")]
        K = sp.cancel(j_mat[("z", "z")] - J)
        Gamma = sp.cancel((j_mat[("x", "y")] + j_mat[("y", "x")]) / 2)
        Gamma_prime = sp.cancel(
            (j_mat[("x", "z")] + j_mat[("z", "x")] + j_mat[("y", "z")] + j_mat[("z", "y")]) / 4
        )
    elif bond == "x":
        J = j_mat[("y", "y")]
        K = sp.cancel(j_mat[("x", "x")] - J)
        Gamma = sp.cancel((j_mat[("y", "z")] + j_mat[("z", "y")]) / 2)
        Gamma_prime = sp.cancel(
            (j_mat[("x", "y")] + j_mat[("y", "x")] + j_mat[("x", "z")] + j_mat[("z", "x")]) / 4
        )
    elif bond == "y":
        J = j_mat[("x", "x")]
        K = sp.cancel(j_mat[("y", "y")] - J)
        Gamma = sp.cancel((j_mat[("x", "z")] + j_mat[("z", "x")]) / 2)
        Gamma_prime = sp.cancel(
            (j_mat[("x", "y")] + j_mat[("y", "x")] + j_mat[("y", "z")] + j_mat[("z", "y")]) / 4
        )
    else:
        raise ValueError(f"Unsupported bond '{bond}'.")

    return {
        "J": J,
        "K": K,
        "Gamma": Gamma,
        "Gamma_prime": Gamma_prime,
        "D": sp.Integer(0),
        "J_xx": j_mat[("x", "x")],
        "J_yy": j_mat[("y", "y")],
        "J_zz": j_mat[("z", "z")],
        "J_xy": j_mat[("x", "y")],
        "J_xz": j_mat[("x", "z")],
        "J_yz": j_mat[("y", "z")],
    }


def derive_bond_exchange_symbolic(
    bond: str,
    U: Expr | Symbol,
    J_H: Expr | Symbol,
    t2: Expr | Symbol | float,
    *,
    t1: Expr | Symbol | float = 0,
    t3: Expr | Symbol | float = 0,
    t4: Expr | Symbol | float = 0,
    convention: str = "spin_12",
) -> dict[str, Expr]:
    """End-to-end symbolic Schrieffer-Wolff -> Kramers projection -> (J, K, Gamma, Gamma') pipeline."""
    h_eff = compute_bond_effective_hamiltonian(
        bond=bond,
        U=U,
        J_H=J_H,
        t2=t2,
        t1=t1,
        t3=t3,
        t4=t4,
    )
    h_pseudo = project_to_pseudospin(h_eff)
    return extract_exchange_tensors(h_pseudo, bond=bond, convention=convention)


def convert_tensor_dict_to_pauli(spin12_tensors: dict[str, Expr]) -> dict[str, Expr]:
    """Convert a spin-1/2 exchange dict to Pauli sigma_i^a sigma_j^b units."""
    return {k: sp.cancel(spin12_to_pauli_coupling(v)) for k, v in spin12_tensors.items()}
