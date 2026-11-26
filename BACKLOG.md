# Backlog & Idea Incubation (`kitaev-derivation`)

Non-committal research queue, distinct from active milestones in PROJECT_STATE.yaml and planned near-future tasks in [HORIZON.md](HORIZON.md). IDs are stable; an entry promoted to the horizon points to its current authority rather than duplicating task checklists.

## IDEA-001 — Complete linked-cluster fourth-order lattice exchange
Status: incubating beyond the near-term three-site work ([#5](https://github.com/lufreitas0000/kitaev-derivation/issues/5)).  
Question: Extend the existing isolated-bond SW derivation to all connected lattice clusters at a given order, including the old multi-site $S_2$ generator proposal and four-site oracle. Enumerate connected hopping graphs, include reachable multiple-charge sectors (including $d^3$ when relevant), subtract all proper subcluster weights in compatible gauges, and verify graph embeddings/multiplicities. The isolated-bond fourth-order result is not the complete fourth-order lattice Hamiltonian.  
Review trigger: three-site symbolic/oracle parity and a demonstrably size-consistent linked-cluster expansion.

## IDEA-002 — Explicit symbolic finite-$\lambda_{\mathrm{SOC}}$ exchange and charge
Status: incubating; finite-SOC numerical bond oracle already exists in M2.2.  
Question: Retain $\lambda$ in the atomic resolvent and derive rational $J,K,\Gamma,\Gamma',\mathbf D$, dressed charge coefficients and higher-order connected processes as functions of $\lambda$, without calling the charge-denominator approximation literal $\lambda\to\infty$.  
Method: Within cubic $t_{2g}$, use fixed one-hole Kramers projector for $\lambda>0$ and block-diagonalize the 15-state two-hole SOC+Kanamori sector by total effective angular momentum $\mathcal J$: multiplicity blocks $2\times2$ for $\mathcal J=0,2$ and $1\times1$ for $\mathcal J=1$, with magnetic degeneracies. Invert exact blocks rather than expanding radical eigenvalues; include $+2\lambda$ reference shift where applicable and the directed onsite-potential shifts. Neutral denominators include $3\lambda/2$ and $3\lambda$. Preserve factored hopping monomials.  
First tests: reproduce the finite-$\lambda$ 66-state oracle, recover the charge-denominator approximation in its scale-separated second-order regime, distinguish the true fixed-$U,J_H$ large-$\lambda$ limit, check positivity and breakdown at physical $\lambda=0$. Compute fourth-order neutral channels rather than treating their absence as a controlled complete fourth-order result. Extend to the NNN and effective-density tasks after their initial cubic-projector versions.  
References: near-term convention audit #4 and validity benchmarks #7.

## IDEA-003 — Spin-$1$ / $j_{\mathrm{eff}}=3/2$ Quadrupolar Superexchange Comparator
Status: incubating (legacy idea, not part of the current $d^5$ program).  
Question: Derive microscopic $K$–$J$ vs. biquadratic exchange tensors for higher-spin Class II platforms intersecting `ksl-sc-bilayer` §5. Requires independent scope/model decision; no five-orbital implementation is implied.

## IDEA-004 — Electrostatic bias and DM superexchange
Status: promoted to near-term N1 in [HORIZON.md](HORIZON.md); authoritative task and tests in [#2](https://github.com/lufreitas0000/kitaev-derivation/issues/2).

## IDEA-005 — Trigonal distortions within $t_{2g}$
Status: incubating; separate from IDEA-002.  
Question: How does distortion of the octahedral environment modify the isolated doublet, SOC gap, exchange tensors and effective charge operator? Begin with a specified sign and axis convention for a local trigonal field, for example $H_{\mathrm{trig}}=\Delta_{\mathrm{trig}}[(\mathbf l_{\mathrm{eff}}\cdot\hat{\mathbf n})^2-2/3]$, $\hat{\mathbf n}\parallel[111]$. Diagonalize the local $6\times6$ SOC+crystal-field block and follow the lowest isolated Kramers pair as $\Delta_{\mathrm{trig}}/\lambda$ varies. Compute doublet admixture, gaps, anisotropic g tensor, full exchange/charge tensors and possible new symmetry-allowed terms. Local orbital splitting and distortion-induced changes in the hopping geometry are distinct models. Uniform inversion-preserving trigonal distortion is not by itself proof of a nonzero DM vector.  
Tests: cubic limit, Kramers degeneracy under time reversal, gap/domain checks, basis-gauge covariance and independent Fock oracle.  
Primary reference: Liu, Chaloupka & Khaliullin, PRB 105, 214411 (2022), https://doi.org/10.1103/PhysRevB.105.214411 (trigonal-field dependence; its extra e_g channels are outside our active scope).

## IDEA-006 — All-orders convergence and thermodynamic control
Status: incubating; quantitative diagnostics tracked in #7.  
Question: Establish a certified domain in $(U,J_H,\lambda,t_a,t_2',\delta)$ for microscopic downfolding, with geometric all-orders coefficient/remainder bounds where possible. Use analytic spectral projectors/direct rotation, complex-parameter gap and expansion radius, and the block connectivity $P_4\leftrightarrow Q_{\rm ch}\leftrightarrow Q_{\rm SOC}$ to exploit separate charge and SOC scales. A real-axis gap or agreement of low-order coefficients is not an all-orders proof. Investigate sharp two-stage resolvent/Schrieffer–Wolff bounds and locality-sensitive linked-cluster estimates for extensive systems. Distinguish microscopic convergence and finite-order approximation from RG relevance/infrared Kitaev phase stability; latter is a downstream theory problem, not established by convergence.  
Reference: Bravyi, DiVincenzo & Loss, Ann. Phys. 326, 2793–2826 (2011), https://arxiv.org/abs/1105.0675 (Lemma 3.4 and linked cluster theorem).

## IDEA-007 — Microscopic hexagon winding exchange
Status: incubating, after #8 audit and IDEA-001.  
Question: Derive the flux-sensitive sixth-order operator directly from the multiorbital $t_{2g}$ model, rather than assuming the phenomenological two-mode Kramers hopping. Use reachable sparse states/connected paths or exact symbolic multiplet blocks; the full six-hole space has $\binom{36}{6}=1947792$ states before pruning. Retain all relevant intermediate multiplets and terms, account for subtraction of lower-order subclusters, verify flux/axis and Pauli conventions and compare to the phenomenological model only with an explicitly derived matching prescription. Do not equate winding contribution to the entire sixth-order Hamiltonian or assume proportionality to $W_p$.

## IDEA-008 — Additional inversion-breaking mechanisms (conditional)
Status: incubating, triggered if scalar onsite bias in #2 does not generate DM or is inadequate.  
Question: Derive the most general symmetry-allowed electric-field-induced changes to orbital hopping matrices and local crystal fields, including antisymmetric orbital hopping in fixed oriented bonds. Distinguish scalar onsite bias from changes to hopping/ligand geometry and any external magnetic field. Determine which DM components can be nonzero and their field parity, with independently justified material parameters; no forced nonzero answer.

## Scope boundary
The active microscopic model is three-orbital low-spin $d^5$, not the full five-orbital $d$ shell. Explicit $e_g$ orbitals and high-spin extensions are outside the current program; this is a model-scope decision, not a universal claim of physical irrelevance. Keep formal CAR/projector/multiplet extensions in finite-QMB; no Lean proof or finite-cluster numerical test alone certifies the entire effective lattice model.
