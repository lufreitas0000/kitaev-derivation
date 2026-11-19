"""Exact symbolic Canonical Anti-commutation Relations (CAR) and projector algebra."""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp
from sympy import Expr

from kitaev_derivation.domain.conventions import ORBITALS, SPINS


class AbstractQuantumOperator:
    """Base class for non-commutative operator AST nodes."""

    def __mul__(self, other: object) -> AbstractQuantumOperator:
        if isinstance(other, AbstractQuantumOperator):
            return _mul_ops(self, other)
        if isinstance(other, (int, Expr)):
            return _scale_op(sp.sympify(other), self)
        return NotImplemented

    def __rmul__(self, other: object) -> AbstractQuantumOperator:
        if isinstance(other, (int, Expr)):
            return _scale_op(sp.sympify(other), self)
        return NotImplemented

    def __add__(self, other: object) -> AbstractQuantumOperator:
        if isinstance(other, AbstractQuantumOperator):
            return _add_ops(self, other)
        return NotImplemented

    def __sub__(self, other: object) -> AbstractQuantumOperator:
        if isinstance(other, AbstractQuantumOperator):
            return _add_ops(self, _scale_op(sp.Integer(-1), other))
        return NotImplemented

    def __neg__(self) -> AbstractQuantumOperator:
        return _scale_op(sp.Integer(-1), self)


@dataclass(frozen=True)
class ZeroOp(AbstractQuantumOperator):
    """Additive zero operator."""


@dataclass(frozen=True)
class IdentityOp(AbstractQuantumOperator):
    """Multiplicative identity operator."""


class FermionOp(AbstractQuantumOperator):
    """Base class for canonical fermionic creation/annihilation operators."""

    site: str
    orbital: str
    spin: str

    def sort_key(self) -> tuple[int, str, int, int]:
        raise NotImplementedError


@dataclass(frozen=True)
class FermionC(FermionOp):
    """Hole creation operator h^dagger_{site, orbital, spin}."""

    site: str
    orbital: str
    spin: str

    def sort_key(self) -> tuple[int, str, int, int]:
        return (0, self.site, ORBITALS.index(self.orbital), SPINS.index(self.spin))


@dataclass(frozen=True)
class FermionA(FermionOp):
    """Hole annihilation operator h_{site, orbital, spin}."""

    site: str
    orbital: str
    spin: str

    def sort_key(self) -> tuple[int, str, int, int]:
        return (1, self.site, ORBITALS.index(self.orbital), SPINS.index(self.spin))


class ProjectorOp(AbstractQuantumOperator):
    """Base class for Hermitian orthogonal subspace projectors."""

    site: str


@dataclass(frozen=True)
class ParticleProjector(ProjectorOp):
    """Particle/hole-number sector projector P^{(holes)}_{site}."""

    holes: int
    site: str


@dataclass(frozen=True)
class MultipletProjector(ProjectorOp):
    """Two-hole d4 orbital-angular-momentum multiplet projector P_{[L=0,1,2]}."""

    L_label: str
    site: str


@dataclass(frozen=True)
class LowEnergyProjector(ProjectorOp):
    """Single-hole spin-orbit j_eff projector P_{1/2} or P_{3/2}."""

    jeff: str
    site: str


@dataclass(frozen=True)
class SpinComposite(AbstractQuantumOperator):
    """Pseudospin Pauli operator sigma^{component}_{site}."""

    site: str
    component: str


@dataclass(frozen=True)
class GenericOp(AbstractQuantumOperator):
    """Named abstract macro operator (e.g., V_N, V_S, V_L, T_1) for commutator oracles."""

    name: str


@dataclass(frozen=True)
class ScaledOperator(AbstractQuantumOperator):
    """Exact SymPy scalar multiple of an operator."""

    scalar: Expr
    op: AbstractQuantumOperator


@dataclass(frozen=True)
class OperatorString(AbstractQuantumOperator):
    """Ordered non-commutative product of atomic operator factors."""

    factors: tuple[AbstractQuantumOperator, ...]


@dataclass(frozen=True)
class OperatorSum(AbstractQuantumOperator):
    """Linear sum of operator terms."""

    terms: tuple[AbstractQuantumOperator, ...]


ZERO = ZeroOp()
IDENTITY = IdentityOp()


def _is_zero_scalar(val: Expr) -> bool:
    return bool(sp.expand(val) == 0)


def _scale_op(scalar: Expr, op: AbstractQuantumOperator) -> AbstractQuantumOperator:
    if isinstance(op, ZeroOp) or _is_zero_scalar(scalar):
        return ZERO
    if scalar == 1:
        return op
    if isinstance(op, ScaledOperator):
        combined = sp.expand(scalar * op.scalar)
        if _is_zero_scalar(combined):
            return ZERO
        return op.op if combined == 1 else ScaledOperator(combined, op.op)
    if isinstance(op, OperatorSum):
        return OperatorSum(tuple(_scale_op(scalar, t) for t in op.terms))
    return ScaledOperator(scalar, op)


def _mul_ops(a: AbstractQuantumOperator, b: AbstractQuantumOperator) -> AbstractQuantumOperator:
    if isinstance(a, ZeroOp) or isinstance(b, ZeroOp):
        return ZERO
    if isinstance(a, IdentityOp):
        return b
    if isinstance(b, IdentityOp):
        return a
    if isinstance(a, OperatorSum):
        res: AbstractQuantumOperator = ZERO
        for t in a.terms:
            res = _add_ops(res, _mul_ops(t, b))
        return res
    if isinstance(b, OperatorSum):
        res = ZERO
        for t in b.terms:
            res = _add_ops(res, _mul_ops(a, t))
        return res
    if isinstance(a, ScaledOperator):
        return _scale_op(a.scalar, _mul_ops(a.op, b))
    if isinstance(b, ScaledOperator):
        return _scale_op(b.scalar, _mul_ops(a, b.op))
    af = a.factors if isinstance(a, OperatorString) else (a,)
    bf = b.factors if isinstance(b, OperatorString) else (b,)
    return OperatorString(af + bf)


def _add_ops(a: AbstractQuantumOperator, b: AbstractQuantumOperator) -> AbstractQuantumOperator:
    if isinstance(a, ZeroOp):
        return b
    if isinstance(b, ZeroOp):
        return a
    at = a.terms if isinstance(a, OperatorSum) else (a,)
    bt = b.terms if isinstance(b, OperatorSum) else (b,)
    return OperatorSum(at + bt)


def commutator(a: AbstractQuantumOperator, b: AbstractQuantumOperator) -> AbstractQuantumOperator:
    """Return algebraic commutator [a, b] = a*b - b*a."""
    return (a * b) - (b * a)


def anticommutator(
    a: AbstractQuantumOperator, b: AbstractQuantumOperator
) -> AbstractQuantumOperator:
    """Return exact CAR anticommutator {a, b} for fermionic operators or a*b + b*a otherwise."""
    if isinstance(a, ScaledOperator):
        return _scale_op(a.scalar, anticommutator(a.op, b))
    if isinstance(b, ScaledOperator):
        return _scale_op(b.scalar, anticommutator(a, b.op))
    if isinstance(a, FermionOp) and isinstance(b, FermionOp):
        if (
            type(a) is not type(b)
            and a.site == b.site
            and a.orbital == b.orbital
            and a.spin == b.spin
        ):
            return IDENTITY
        return ZERO
    return (a * b) + (b * a)


def dagger(op: AbstractQuantumOperator) -> AbstractQuantumOperator:
    """Return Hermitian conjugate of an operator AST node."""
    if isinstance(op, FermionC):
        return FermionA(op.site, op.orbital, op.spin)
    if isinstance(op, FermionA):
        return FermionC(op.site, op.orbital, op.spin)
    if isinstance(op, (ProjectorOp, IdentityOp, ZeroOp, GenericOp, SpinComposite)):
        return op
    if isinstance(op, ScaledOperator):
        return _scale_op(sp.conjugate(op.scalar), dagger(op.op))
    if isinstance(op, OperatorString):
        return OperatorString(tuple(dagger(f) for f in reversed(op.factors)))
    if isinstance(op, OperatorSum):
        return OperatorSum(tuple(dagger(t) for t in op.terms))
    raise TypeError(f"Unsupported operator type for dagger: {type(op)}")


def evaluate_projector_product(p1: ProjectorOp, p2: ProjectorOp) -> AbstractQuantumOperator:
    """Evaluate adjacent projector multiplication on the same or distinct sites."""
    if p1.site != p2.site:
        return OperatorString((p1, p2))
    # MultipletProjector(L=0,1,2) sits strictly inside ParticleProjector(holes=2).
    if isinstance(p1, MultipletProjector) and isinstance(p2, ParticleProjector):
        return p1 if p2.holes == 2 else ZERO
    if isinstance(p1, ParticleProjector) and isinstance(p2, MultipletProjector):
        return p2 if p1.holes == 2 else ZERO
    # LowEnergyProjector(1/2, 3/2) sits strictly inside ParticleProjector(holes=1).
    if isinstance(p1, LowEnergyProjector) and isinstance(p2, ParticleProjector):
        return p1 if p2.holes == 1 else ZERO
    if isinstance(p1, ParticleProjector) and isinstance(p2, LowEnergyProjector):
        return p2 if p1.holes == 1 else ZERO
    if type(p1) is not type(p2):
        return ZERO
    return p1 if p1 == p2 else ZERO


def apply_algebraic_rules(
    op_seq: list[AbstractQuantumOperator] | tuple[AbstractQuantumOperator, ...],
) -> tuple[AbstractQuantumOperator, ...]:
    """Prune projector orthogonality/idempotency, Pauli exclusion, and same-site intertwining."""
    if not op_seq:
        return (IDENTITY,)

    simplified: list[AbstractQuantumOperator] = []
    for op in op_seq:
        if isinstance(op, ZeroOp):
            return (ZERO,)
        if isinstance(op, IdentityOp):
            continue
        if not simplified:
            simplified.append(op)
            continue

        prev = simplified[-1]
        if isinstance(prev, ProjectorOp) and isinstance(op, ProjectorOp):
            res = evaluate_projector_product(prev, op)
            if isinstance(res, ZeroOp):
                return (ZERO,)
            if isinstance(res, ProjectorOp):
                simplified[-1] = res
                continue

        if isinstance(prev, FermionOp) and isinstance(op, FermionOp) and prev == op:
            return (ZERO,)

        # Same-site fermion-particle-projector intertwining:
        # h^\dagger P^{(n)} = P^{(n+1)} h^\dagger,   h P^{(n)} = P^{(n-1)} h
        if (
            isinstance(prev, FermionOp)
            and isinstance(op, ParticleProjector)
            and prev.site == op.site
        ):
            if isinstance(prev, FermionC):
                if op.holes >= 2:
                    return (ZERO,)
                simplified[-1] = ParticleProjector(op.holes + 1, op.site)
                simplified.append(prev)
                continue
            if op.holes <= 0:
                return (ZERO,)
            simplified[-1] = ParticleProjector(op.holes - 1, op.site)
            simplified.append(prev)
            continue

        simplified.append(op)

    return tuple(simplified) if simplified else (IDENTITY,)


def sort_normal_order(
    op_seq: tuple[AbstractQuantumOperator, ...] | list[AbstractQuantumOperator],
) -> list[tuple[tuple[AbstractQuantumOperator, ...], int]]:
    """Normal-order adjacent fermionic factors via exact CAR contractions.

    Projectors act as boundaries between fermion blocks and are preserved in order;
    fermions on distinct sites anticommute purely without contracting.
    """
    queue: list[tuple[list[AbstractQuantumOperator], int]] = [
        ([x for x in op_seq if not isinstance(x, IdentityOp)], 1)
    ]
    final_terms: list[tuple[tuple[AbstractQuantumOperator, ...], int]] = []

    while queue:
        seq, sign = queue.pop(0)
        if not seq:
            final_terms.append(((IDENTITY,), sign))
            continue
        if any(isinstance(x, ZeroOp) for x in seq):
            continue

        swapped = False
        for i in range(len(seq) - 1):
            op1, op2 = seq[i], seq[i + 1]
            if isinstance(op1, FermionOp) and isinstance(op2, FermionOp):
                if op1 == op2:
                    swapped = True
                    break  # Pauli nilpotency: c^2 = 0
                if op2.sort_key() < op1.sort_key():
                    same_mode = (
                        type(op1) is not type(op2)
                        and op1.site == op2.site
                        and op1.orbital == op2.orbital
                        and op1.spin == op2.spin
                    )
                    if same_mode:
                        # h_a h_a^\dagger = I - h_a^\dagger h_a
                        seq_id = seq[:i] + seq[i + 2 :]
                        queue.append((seq_id, sign))
                    seq_swap = list(seq)
                    seq_swap[i], seq_swap[i + 1] = op2, op1
                    queue.append((seq_swap, -sign))
                    swapped = True
                    break
        if not swapped:
            final_terms.append((tuple(seq), sign))

    return final_terms


def apply_rewrite_rules(
    op: AbstractQuantumOperator,
    rules: dict[AbstractQuantumOperator, AbstractQuantumOperator],
) -> AbstractQuantumOperator:
    """Recursively rewrite AST nodes matching exact structural rules (e.g. Appendix A)."""
    if op in rules:
        return rules[op]
    if isinstance(op, ScaledOperator):
        rewritten = apply_rewrite_rules(op.op, rules)
        return _scale_op(op.scalar, rewritten)
    if isinstance(op, OperatorString):
        new_factors = tuple(apply_rewrite_rules(f, rules) for f in op.factors)
        if any(isinstance(f, ZeroOp) for f in new_factors):
            return ZERO
        return OperatorString(new_factors)
    if isinstance(op, OperatorSum):
        acc: AbstractQuantumOperator = ZERO
        for t in op.terms:
            acc = _add_ops(acc, apply_rewrite_rules(t, rules))
        return acc
    return op


def collect_terms(op: AbstractQuantumOperator) -> OperatorSum:
    """Group terms sharing identical operator topologies and simplify exact SymPy scalars."""
    if isinstance(op, ZeroOp):
        return OperatorSum(())
    terms = op.terms if isinstance(op, OperatorSum) else (op,)
    buckets: dict[AbstractQuantumOperator, Expr] = {}
    for term in terms:
        if isinstance(term, ZeroOp):
            continue
        if isinstance(term, ScaledOperator):
            scalar = term.scalar
            base = term.op
        else:
            scalar = sp.Integer(1)
            base = term
        buckets[base] = buckets.get(base, sp.Integer(0)) + scalar

    cleaned: list[AbstractQuantumOperator] = []
    for base, scalar in buckets.items():
        s_simp = sp.cancel(sp.expand(scalar))
        if s_simp == 0:
            continue
        if s_simp == 1:
            cleaned.append(base)
        else:
            cleaned.append(ScaledOperator(s_simp, base))
    return OperatorSum(tuple(cleaned))
