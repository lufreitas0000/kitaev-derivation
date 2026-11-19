# Current Sprint — kitaev-derivation / M2.2

Owning milestone: M2.2 (see [HORIZON.md](HORIZON.md))
Machine state: [PROJECT_STATE.yaml](PROJECT_STATE.yaml)

## Scientific / Engineering Objective

Complete Wave 3 `M2.2`: exact 66-state fourth-order BDL/Takahashi bond exchange (separating $\lambda_{\mathrm{SOC}} \to \infty$ charge renormalization from finite-$\lambda_{\mathrm{SOC}}$ $j_{\mathrm{eff}}=3/2$ neutral intermediate channels), 924-state 6-site Kramers hexagon Peierls $\pi$-flux ring exchange ($K_{\mathrm{ring}}$), idealized $\alpha\text{-RuCl}_3$ / $\text{Na}_2\text{IrO}_3$ ACID exports (`data/exports/extended_kitaev_material_couplings.{dat,meta.json}`), and formal Lean 4 bridge PR (`lufreitas0000/finite-QMB#52`).

## Scope and Non-Goals

Exports represent idealized Hubbard–Kanamori–SOC reference models (not trigonal-field DFT fits). Formal Lean 4 theorem acceptance remains gated on maintainer review of `finite-QMB#52`.

## Exact Next Step

Await maintainer review/merge of [finite-QMB#52](https://github.com/lufreitas0000/finite-QMB/pull/52) and ingest `data/exports/extended_kitaev_material_couplings.dat` into downstream Wave 3 sweeps.
