"""Exact symbolic Schrieffer-Wolff CAS and two-site Fock oracle for Kitaev materials."""

from kitaev_derivation.domain.algebra import (
    FermionA,
    FermionC,
    GenericOp,
    IdentityOp,
    LowEnergyProjector,
    MultipletProjector,
    OperatorString,
    OperatorSum,
    ParticleProjector,
    ScaledOperator,
    SpinComposite,
    ZeroOp,
    anticommutator,
    collect_terms,
    commutator,
    dagger,
    sort_normal_order,
)
from kitaev_derivation.domain.conventions import (
    ORBITALS,
    SPINS,
    kramers_doublet_kets,
    orbital_angular_momentum_matrices,
    pauli_matrices,
    spin12_to_pauli_coupling,
)
from kitaev_derivation.domain.kanamori_soc import (
    d4_multiplet_basis,
    d4_multiplet_projectors,
    kanamori_eigenvalues,
    single_site_soc_matrix,
)
from kitaev_derivation.services.exchange_extraction import (
    derive_bond_exchange_symbolic,
    evaluate_pauli_trace,
    extract_exchange_tensors,
    project_to_pseudospin,
)
from kitaev_derivation.services.kinetic_hopping import (
    build_t1_operator,
    build_tm1_operator,
    directional_hopping_matrix,
)
from kitaev_derivation.services.schrieffer_wolff import (
    apply_resolvent,
    compute_effective_hamiltonian,
)
from kitaev_derivation.services.sparse_fock_oracle import (
    compute_fock_kanamori_spectrum,
    compute_fock_schrieffer_wolff_exchange,
    compute_fock_soc_finite_exchange,
)

__all__ = [
    "FermionA",
    "FermionC",
    "GenericOp",
    "IdentityOp",
    "LowEnergyProjector",
    "MultipletProjector",
    "ORBITALS",
    "OperatorString",
    "OperatorSum",
    "ParticleProjector",
    "SPINS",
    "ScaledOperator",
    "SpinComposite",
    "ZeroOp",
    "anticommutator",
    "apply_resolvent",
    "build_t1_operator",
    "build_tm1_operator",
    "collect_terms",
    "commutator",
    "compute_effective_hamiltonian",
    "compute_fock_kanamori_spectrum",
    "compute_fock_schrieffer_wolff_exchange",
    "compute_fock_soc_finite_exchange",
    "d4_multiplet_basis",
    "d4_multiplet_projectors",
    "dagger",
    "derive_bond_exchange_symbolic",
    "directional_hopping_matrix",
    "evaluate_pauli_trace",
    "extract_exchange_tensors",
    "kanamori_eigenvalues",
    "kramers_doublet_kets",
    "orbital_angular_momentum_matrices",
    "pauli_matrices",
    "project_to_pseudospin",
    "single_site_soc_matrix",
    "sort_normal_order",
    "spin12_to_pauli_coupling",
]
