# kitaev-derivation

Exact symbolic Schrieffer-Wolff projector Computer Algebra System (SymPy) and independent 66-state two-hole Fock-space oracle (`scipy.sparse`) deriving microscopic $(J, K, \Gamma, \Gamma')$ effective spin Hamiltonians for $d^5$ edge-shared octahedral Kitaev materials.

Migrated cleanly from legacy Julia `KitaevDerivation.jl` (`archive/KitaevDerivation.jl/`, including checkpoint `0aed437`) into unified `uv`-managed pure Python + strict typing (`mypy --strict`, `ruff`), with formal CAR/projector theorems routed to Lean 4 in [`research/qmb/finite-QMB`](../../qmb/finite-QMB/README.md) (see [`doc/lean_handoff_finite_qmb.md`](doc/lean_handoff_finite_qmb.md)).

## Goal

Develop an algebraically exact (within the stated model and approximation) symbolic Schrieffer-Wolff projector CAS in Python (`SymPy`) verified to machine precision (`< 1e-13`) against an independent 66-state two-hole Fock-space numerical oracle (`scipy.sparse`) across $x, y, z$ bonds and $(U,J_H,t_1,t_2,t_3,t_4)$ in the SOC-projector approximation, with numerical finite-SOC checks distinguished from symbolic results, with ACID benchmark tables and self-contained Lab Book / Manuscript / Blueprint deliverables.

## Scope

- **Owns:** Symbolic CAR normal-ordering, exact $d^4$ Hubbard-Kanamori multiplet projectors ($P_{[L=0,1,2]}$), $t_{2g}$ $T\text{–}P$ SOC Kramers doublet projection, second-order symbolic and fourth-order numerical isolated-bond Schrieffer-Wolff reductions, independent 66-state Fock oracle, separately modeled phenomenological Kramers hexagon winding exchange, and provenance-tagged parameter exports for `kitaev-mf` / `ksl-sc-bilayer`.
- **Strictly Out of Scope:** Lattice-scale Bogoliubov-de Gennes mean-field loops (`kitaev-mf`), multi-site spin ED / 2-RDM SDP (`kitaev-numerics`), spectral walks (`kitaev-walk-geometry`), and formal Lean 4 proofs (`research/qmb/finite-QMB`).

## Definition of Done (DoD)

- `make verify` exits `0` (`ruff check`, `ruff format --check`, `mypy --strict src tests scripts`, and full `pytest` suite including randomized 66-state oracle parity `< 1e-13`).
- Benchmark table `data/processed/benchmarks/exchange_vs_fock_oracle.dat` + `.meta.json` generated deterministically with ACID `#` provenance headers.
- `doc/labbook/labbook_kitaev_derivation.pdf`, `doc/manuscript/manuscript_kitaev_derivation.pdf`, and `build/blueprint/index.html` compile with zero LaTeX errors or undefined references/citations.
- `scripts/uow.py check --repo kitaev-derivation`, `governance/pg10.py check`, and `scripts/fiscalize_state.py --repo kitaev-derivation` pass with zero warnings.

## Scientific scope and next steps

The starting model is three-orbital low-spin $d^5$ ($t_{2g}$ only). Six spin-orbitals per site give $\binom{12}{2}=66$ states on an isolated two-hole bond, with 36 neutral states and four retained Kramers-product states. Hopping conserves isolated-cluster total hole number but changes local occupations; superexchange requires integrating out virtual excitations, not pre-projecting hopping. Full five-orbital $d$-shell/$e_g$ modeling is outside the active scope.

M2.1 and M2.2 are completed within their tested scope. The legacy "infSOC" / zero-lambda selector denotes a SOC-projector / charge-denominator approximation, not literal infinite SOC at fixed charge parameters; physical finite SOC is treated numerically in the bond oracle. The fourth-order coefficient cross-check does not certify accuracy at physical hopping. The 924-state hexagon calculation is an auxiliary two-mode Kramers Hubbard model, not yet a microscopic multiorbital sixth-order result; flux axes and operator interpretation are under audit (#8).

Near-term goals (planned, not implemented):
- Local electrostatic potentials: bias-dependent exchange and all DM components, including possible symmetry-protected zeros (#2).
- Direct orbital NNN hopping $t_2'$ and three-site connected processes (#5).
- Effective charge density including biased-bond and triangular terms, reproducing Pereira–Egger (2020) and Freitas et al. (2024) (#6).
- Convention correction (#4), scale-aware convergence diagnostics (#7), and hexagon audit (#8).

See [HORIZON.md](HORIZON.md) for the ordered todo list and issue links; [BACKLOG.md](BACKLOG.md) for explicit symbolic finite SOC, trigonal distortions, all-orders bounds and complete microscopic multicluster extensions. Current sprint and completed milestone evidence remain in [PROJECT_STATE.yaml](PROJECT_STATE.yaml).
