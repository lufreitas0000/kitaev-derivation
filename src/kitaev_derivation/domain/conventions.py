"""Physical conventions for t2g orbitals, T-P equivalence, Kramers doublets, and QSL unit ports.

All matrices and ket amplitudes use exact SymPy algebraic numbers (Rational, sqrt, I)
so symbolic derivations never introduce floating-point truncation errors.
"""

from __future__ import annotations

import sympy as sp
from sympy import Expr, I, Matrix, Rational, sqrt

ORBITALS: tuple[str, str, str] = ("x", "y", "z")
SPINS: tuple[str, str] = ("up", "dn")
BONDS: tuple[str, str, str] = ("x", "y", "z")


def orbital_angular_momentum_matrices() -> dict[str, Matrix]:
    """Return exact l_eff = 1 orbital matrices in the {|x>=|yz>, |y>=|zx>, |z>=|xy>} basis.

    Under T-P equivalence, L_{t2g} = -l_eff with [l_eff^x, l_eff^y] = i * l_eff^z.
    """
    lx = Matrix([[0, 0, 0], [0, 0, -I], [0, I, 0]])
    ly = Matrix([[0, 0, I], [0, 0, 0], [-I, 0, 0]])
    lz = Matrix([[0, -I, 0], [I, 0, 0], [0, 0, 0]])
    return {"x": lx, "y": ly, "z": lz}


def pauli_matrices() -> dict[str, Matrix]:
    """Return standard 2x2 Pauli matrices {I, x, y, z} with exact SymPy entries."""
    return {
        "I": Matrix([[1, 0], [0, 1]]),
        "x": Matrix([[0, 1], [1, 0]]),
        "y": Matrix([[0, -I], [I, 0]]),
        "z": Matrix([[1, 0], [0, -1]]),
    }


def pauli_matrix_element(axis: str, bra_spin: str, ket_spin: str) -> Expr:
    """Return <bra_spin | sigma^axis | ket_spin> for spin in {'up', 'dn'}."""
    row = 0 if bra_spin == "up" else 1
    col = 0 if ket_spin == "up" else 1
    val: Expr = pauli_matrices()[axis][row, col]
    return val


def kramers_doublet_kets() -> dict[str, dict[tuple[str, str], Expr]]:
    """Return exact single-hole j_eff = 1/2 Kramers doublet amplitudes in (orbital, spin) basis.

    Conventions (Jackeli & Khaliullin PRL 102, 017205 (2009); T-P equivalence L_{t2g} = -l_eff):
        |+> = (1/sqrt(3)) * (-|z, up> - i|y, dn> - |x, dn>)
        |-> = (1/sqrt(3)) * ( |z, dn> + i|y, up> - |x, up>)
    These satisfy (S . l_eff) |+/-> = -1 |+/-> and P_{1/2} (S + l_eff) P_{1/2} = sigma / 2.
    """
    inv_sqrt3 = Rational(1, 1) / sqrt(3)
    return {
        "up": {
            ("z", "up"): -inv_sqrt3,
            ("y", "dn"): -I * inv_sqrt3,
            ("x", "dn"): -inv_sqrt3,
        },
        "dn": {
            ("z", "dn"): inv_sqrt3,
            ("y", "up"): I * inv_sqrt3,
            ("x", "up"): -inv_sqrt3,
        },
    }


def single_hole_mode_index(orbital: str, spin: str) -> int:
    """Map (orbital, spin) to the canonical 6-mode single-site index.

    Ordering: [x up, y up, z up, x dn, y dn, z dn].
    """
    return SPINS.index(spin) * 3 + ORBITALS.index(orbital)


def kramers_doublet_vectors() -> tuple[Matrix, Matrix]:
    """Return 6-component SymPy column vectors |+> and |-> in canonical 6-mode order."""
    amps = kramers_doublet_kets()
    vec_p = sp.zeros(6, 1)
    vec_m = sp.zeros(6, 1)
    for (orb, spin), coeff in amps["up"].items():
        vec_p[single_hole_mode_index(orb, spin), 0] = coeff
    for (orb, spin), coeff in amps["dn"].items():
        vec_m[single_hole_mode_index(orb, spin), 0] = coeff
    return vec_p, vec_m


def spin12_to_pauli_coupling(spin12_value: Expr | float) -> Expr | float:
    """Convert S_i^a S_j^b exchange coefficient to sigma_i^a sigma_j^b Pauli convention.

    Since S_i^a = sigma_i^a / 2, we have J_{ab}^{(spin-1/2)} S_i^a S_j^b = (J_{ab}/4) sigma_i^a sigma_j^b.
    Aligns directly with research/qsl/contracts/qsl_conventions.json.
    """
    return spin12_value / 4
