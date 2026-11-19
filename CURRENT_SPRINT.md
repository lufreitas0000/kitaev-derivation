# Current Sprint — kitaev-derivation / M2.2

Owning milestone: M2.2 (see [HORIZON.md](HORIZON.md))
Machine state: [PROJECT_STATE.yaml](PROJECT_STATE.yaml)

## Scientific / Engineering Objective

Complete Wave 3 `M2.2`: exact 66-state fourth-order BDL/Takahashi bond exchange (separating $\lambda_{\mathrm{SOC}} \to \infty$ charge renormalization from finite-$\lambda_{\mathrm{SOC}}$ $j_{\mathrm{eff}}=3/2$ neutral intermediate channels), 924-state 6-site Kramers hexagon Peierls $\pi$-flux ring exchange ($K_{\mathrm{ring}}$), idealized $\alpha\text{-RuCl}_3$ / $\text{Na}_2\text{IrO}_3$ ACID exports (`data/exports/extended_kitaev_material_couplings.{dat,meta.json}`), and formal Lean 4 bridge PR (`lufreitas0000/finite-QMB#52`).

## Scope and Non-Goals

Exports represent idealized Hubbard–Kanamori–SOC reference models (not trigonal-field DFT fits). Follow-up formal proofs extending `FiniteQMB.Models.KitaevKanamoriSW` (multi-gap Sylvester `[H_0, S] = V_{od}` operator and concrete 66-state Fock isomorphism) are tracked on `research/qmb/finite-QMB` (`main` @ `7098b031`).

## Exact Next Step

Ingest `data/exports/extended_kitaev_material_couplings.dat` into downstream `kitaev-mf` / `ksl-sc-bilayer` Wave 2–3 sweeps (and track follow-up multi-gap Sylvester formalization in `finite-QMB`).
