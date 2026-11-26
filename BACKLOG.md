# Backlog & Idea Incubation (`kitaev-derivation`)

Non-committal idea queue separated from active milestones in `PROJECT_STATE.yaml`.

## IDEA-001 — Fourth-Order Ring-Exchange & Multi-Site $S_2$ Generator
Status: queued (Wave 3 / `M2.2`)  
Question: Extend `schrieffer_wolff.py` with `[S_1, T_0] + [S_2, V] = 0` third/fourth-order plaquette terms and validate against a $N_{\mathrm{sites}}=4$ sparse Fock cluster oracle.

## IDEA-002 — Finite-$\lambda_{\mathrm{SOC}}$ Perturbative Expansion & Trigonal Crystal Field $\Delta_{\mathrm{trig}}$
Status: incubating  
Question: Quantify departures from the infinite-$\lambda_{\mathrm{SOC}}$ Kramers manifold as $\lambda_{\mathrm{SOC}} / \Delta_{\mathrm{trig}}$ varies across $4d^5$ ($\alpha\text{-RuCl}_3$) vs. $5d^5$ ($\mathrm{Na}_2\mathrm{IrO}_3$) vs. $3d^7$ ($\mathrm{Co}^{2+}$).

## IDEA-003 — Spin-$1$ / $j_{\mathrm{eff}}=3/2$ Quadrupolar Superexchange Comparator
Status: incubating  
Question: Derive the microscopic $K\text{–}J$ vs. biquadratic exchange tensors for higher-spin Class II platforms intersecting `ksl-sc-bilayer` §5.

## IDEA-004 — Electrostatic Bias and Dzyaloshinskii–Moriya Superexchange
Status: inbox  
Question: How does a scalar onsite energy difference $\delta=\epsilon_i-\epsilon_j$ modify the second-order bond exchange tensor in the SOC-projector approximation, and is any DM component nonzero?  
Motivation: Determine electric-field tunability of $J,K,\Gamma,\Gamma'$ and antisymmetric exchange without conflating inversion breaking with a sufficient condition for DM.  
Provenance / scope: 2026-11-26 owner-requested investigation; [Issue #2](https://github.com/lufreitas0000/kitaev-derivation/issues/2). Begin with the existing real symmetric $T_\gamma$ and $H_E=\delta(n_i-n_j)/2$ in the hole representation; shift directed multiplet denominators by $\pm\delta$ and require positive charge gaps.  
First test / falsifier: Derive the full $3\times3$ tensor, recover $\delta=0$, check $\mathcal J(-\delta)=\mathcal J(\delta)^{\mathsf T}$ and pure-$t_2$ limits, test the residual-symmetry prediction $\mathbf D=0$, and independently cross-check the 66-state Fock oracle for admissible positive and negative biases. None of these finite-bias results is yet verified in this repository.  
Review trigger: Symbolic expressions and numerical parity evidence, including all three DM components and bias-dependent exchange tables. If DM vanishes in the scalar-bias model, investigate field-induced asymmetric hopping or local crystal fields separately. No change to the active sprint or completed milestones.
