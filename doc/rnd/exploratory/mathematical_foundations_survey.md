---
id: RND-EXPL-MATHEMATICAL-FOUNDATIONS
title: "Mathematical foundations for controlled effective-model derivations"
status: draft
questions:
  - "Which mathematical structures genuinely compress, generalize or certify the current microscopic derivation?"
  - "Which useful exact results transfer from bonds to connected clusters and dressed observables, and under what assumptions?"
depends_on: []
contracts: []
claims: []
---

# Mathematical foundations — exploratory research charter and preliminary reading map

**Status:** preliminary orientation and research charter; [Issue #9](https://github.com/lufreitas0000/kitaev-derivation/issues/9). The comprehensive deep-research survey, guiding notes and worked verification examples are not completed. No candidate abstraction has been accepted for production. See the deferred implementation questions IDEA-009–017 in [BACKLOG.md](../../../BACKLOG.md).

## Scientific objective and Mental Image

Represent a microscopic model as a graph of local graded operator algebras, use symmetry to block-diagonalize atomic Hamiltonians and admissible effective operator spaces, generate only possible virtual paths, and reduce via resolvents/adjoint conjugation. The linked-cluster weights assemble the lattice result. The effective Hamiltonian and charge observables share a dressed embedding, not unrelated projection recipes.

Intelligence gain means discovering common structure that replaces special cases while preserving the content and rigor of the derivation. Read governance [Rule 03](https://github.com/lufreitas0000/project-governance/blob/main/core/03_editorial_rigor.md), [Rule 05](https://github.com/lufreitas0000/project-governance/blob/main/core/05_project_management.md) and [epistemology §2.6](https://github.com/lufreitas0000/project-governance/blob/main/docs/MODERN_GOVERNANCE_AND_EPISTEMOLOGY.md). Do not confuse compression with short notation, number of scripts, arbitrary merit scores or removing physical degrees of freedom without controlled approximations. Passive compression diagnostics are not optimization objectives.

## Preliminary knowledge map

These are established mathematical identities or candidate applications, not newly verified project invariants.

| Area | Structural principle | Possible concrete gain | First small test / limitation |
|---|---|---|---|
| Representation theory | Exterior powers and Schur multiplicity decomposition | Derive atomic multiplets and invert reduced blocks automatically; retain rational SOC dependence | 15-state two-hole multiplet decomposition, compared with full matrix; cubic rotational symmetry may be lowered |
| Noncommutative algebra | CAR quotient, normal forms and confluence; bilinear gl(M) relations | Canonical symbolic equality, reusable rewriting proof | Critical overlaps and independent full-Fock representation; a fixed N sector is not a faithful representation of the full CAR algebra |
| Symmetry/invariant theory | Equivariance of reduction, invariant subspaces of operator representations | Determine allowed terms and independent tensor coefficients before deriving them | Unbiased vs biased bond, TR-even triangle; allowed does not mean nonzero, antiunitarity requires real-linear treatment |
| Graph/gauge theory | Returning to fixed local occupation balances directed hops at every vertex | Enumerate admissible paths, prune exact zero processes, track loops and signed hopping monomials | Bond backtracking, NN–NNN triangle, NN hexagon; balanced flow is necessary, not sufficient for a nonzero ordered fermionic amplitude |
| Posets and incidence algebra | Möbius inversion of embedded subcluster weights | Reuse lower-cluster results without double counting | 3-site subtraction, disconnected additivity; requires consistent operator embeddings and SW gauge |
| Spectral/perturbation theory | Energy-dependent Schur complement, Riesz projectors and analytic continuation | Rational resolvents, independent reduction oracle, controlled error bounds | 66→36→4 nested elimination; Feshbach energy dependence differs from Hermitian SW, real-axis gaps alone do not certify Taylor convergence |
| Operator algebra | Isometric compression is unital completely positive, not multiplicative | Single interface for Hamiltonian, density, current and correlators | Dressed density and charge neutrality, product defect; derivative of gauge-dependent effective H needs a connection term |
| Lie algebra | BCH/adjoint grading; perturbation-order recursion, continuous unitary transformations | Unify higher-order H and observable generation | Reproduce coefficients through order 4 and density through order 3; isospectral flow is not Wilsonian RG |
| Category / algebraic QMB | Local CAR inclusions, graded tensor products, graph isomorphisms, compatible effective embeddings | Make locality and symmetry consistency explicit; quantify noncommuting reductions | Graph relabeling and disjoint union; enlarging an interacting cluster is not in general natural with respect to compression |

The first four areas and nested Schur reduction have particularly concrete applications to #2, #5, #6, #7 and IDEA-002. Categorical language should earn its place through a proved or testable compatibility statement, not a new abstraction layer for its own sake.

## Initial survey leads, checked on 2026-11-26

A general orientation search covered "Schrieffer Wolff transformation linked cluster", "Feshbach Schur map perturbation theory", "Diamond lemma ring theory", "Rota Möbius functions", "Stinespring positive functions C* algebras", "Haag Kastler local algebras", "Wegner flow equations", and "Kitaev materials trigonal spin orbit". The leads below are bibliographically checked; their detailed theorem hypotheses and relevance must be studied in the full survey.

- Bravyi, DiVincenzo & Loss (2011), *Schrieffer–Wolff transformation for quantum many-body systems*, Annals of Physics 326, 2793–2826, [arXiv:1105.0675](https://arxiv.org/abs/1105.0675). Start with finite-dimensional direct rotation, convergence radius and linked-cluster theorem.
- Dusson, Sigal & Stamm (2021), *The Feshbach–Schur map and perturbation theory*, [arXiv:2105.02058](https://arxiv.org/abs/2105.02058). Compare explicit eigenvalue/eigenvector error estimates to SW.
- Bergman (1978), *The diamond lemma for ring theory*, Advances in Mathematics 29, 178–218; [author's corrections and updates](https://math.berkeley.edu/~gbergman/papers/updates/diamond.html) and [publication list](https://math.berkeley.edu/~gbergman/papers/). Study terminating rewrite orders, critical overlaps and irreducible normal forms; consult original theorem, not merely a summary.
- Rota (1964), *On the foundations of combinatorial theory I. Theory of Möbius functions*, Z. Wahrscheinlichkeitstheorie 2, 340–368, [original paper](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/rota1.pdf). Adapt incidence-algebra inversion to operator-valued embedded cluster weights.
- Stinespring (1955), *Positive functions on C*-algebras*, [DOI:10.1090/S0002-9939-1955-0069403-4](https://doi.org/10.1090/S0002-9939-1955-0069403-4). Context for unital completely positive compression; finite-dimensional formula can be established directly.
- Wegner (1994), *Flow-equations for Hamiltonians*, [DOI:10.1002/andp.19945060203](https://doi.org/10.1002/andp.19945060203). Compare continuous unitary transformations, generators and truncation.
- Haag & Kastler (1964), *An algebraic approach to quantum field theory*, J. Math. Phys. 5, 848–861, [DOI:10.1063/1.1704187](https://doi.org/10.1063/1.1704187). Motivating local-algebra perspective, not a proof of categorical properties of our lattice SW reduction.
- Liu, Chaloupka & Khaliullin (2022), *Exchange interactions in d5 Kitaev materials: From Na2IrO3 to alpha-RuCl3*, [arXiv:2206.04326](https://arxiv.org/abs/2206.04326). Physical benchmark for spin–orbit and trigonal exchange; their extra e_g mechanisms are not in our active scope.

Additional survey directions: authoritative texts on exterior-power representation theory, angular momentum and Wigner–Eckart; Lie algebra filtrations and BCH; finite magnetic point-group representation and invariant operators; directed Eulerian multigraphs and rational invariant theory; block/resolvent analytic perturbation (Kato); physics linked-cluster and perturbative continuous unitary transformation literature; symbolic algebra and Lean proof feasibility. Search beyond these prompts, including failed generalizations and contrary examples. Reconcile with prior literature benchmarks in HORIZON.md.

## Concrete guiding-document deliverables (future)

1. A source-linked exploratory survey extending the above map: search log, bibliographic records, exact theorem statement/hypotheses, applicability to our model, competing approaches, open questions and negative findings.
2. A main analytical guiding note organized around algebraic foundations, general reductions, combinatorial/symmetry structure and physical applications (canonical Markdown; corresponding human-facing LaTeX/PDF as appropriate within existing document build). Keep the logical spine concise and full calculations in focused complements.
3. Complementary worked notes on multiplets and CAR; invariant operators and gauge/graph selection; linked cluster/incidence algebra; Feshbach/SW/Lie recursion and effective observables; formal/compositional perspective. Each should show a concrete minimal example and separately label exact facts, assumptions, finite-size evidence and conjectures.
4. A theory-to-implementation crosswalk: present module or repeated calculation, precise candidate abstraction, expected structure reduction, verification evidence, independently reviewable next step and rejection/rollback criterion. Promote only supported backlog ideas. No production code changes in the survey task.

## Evaluation before implementation

Use existing 66-state oracle, prospective 816-state cluster oracle, theoretical invariants and the published charge-operator limits. Test disjoint-cluster additivity, graph isomorphism covariance, all allowed/forbidden operator spaces, correct density dressing and convergence assumptions. Compare complexity/duplicate assumptions and reusable theorems against the current baseline; do not invent a scalar utility score. When a shared abstraction would add temporary complexity, name at least two actual consumer modules, provide an escrow and a bounded payoff/review horizon per governance. Record counterexamples and failures rather than hiding them. Keep formal theorem work in finite-QMB via the existing handoff and numerical checks in this repo.

This opening survey does not alter the active Wave 4 physics priorities or imply that any of the implementation ideas has been scheduled. 
