"""Build static HTML Blueprint Portal for kitaev-derivation under build/blueprint/index.html."""

from __future__ import annotations

import html
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    out_dir = ROOT / "build" / "blueprint"
    out_dir.mkdir(parents=True, exist_ok=True)
    dat_file = ROOT / "data" / "processed" / "benchmarks" / "exchange_vs_fock_oracle.dat"
    dat_preview = (
        html.escape(dat_file.read_text(encoding="utf-8"))
        if dat_file.is_file()
        else "Run scripts/generate_benchmarks.py first."
    )
    exp_file = ROOT / "data" / "exports" / "extended_kitaev_material_couplings.dat"
    exp_preview = (
        html.escape(exp_file.read_text(encoding="utf-8"))
        if exp_file.is_file()
        else "Run scripts/generate_benchmarks.py first."
    )
    mermaid_code = "\n".join(
        [
            "flowchart LR",
            "  classDef ms fill:#1E3A5F,color:#FAF6EE,stroke:#1E3A5F,stroke-width:2px;",
            "  classDef inv fill:#E6F4EA,color:#137333,stroke:#137333,stroke-width:1.5px;",
            "  classDef handoff fill:#F5EFE2,color:#8C4A27,stroke:#8C4A27,stroke-width:1.5px;",
            '  M21["M2.1: Second-Order SymPy Projector CAS"]:::ms',
            '  M22["M2.2: 4th-Order BDL & 924-State Hexagon Ring"]:::ms',
            '  A1["12-Mode Fock / Appendix A Commutators (N=2,3)"]:::inv',
            '  CAS["Symbolic J-K-Gamma-Gamma\' (x,y,z C3 bonds)"]:::inv',
            '  F66["66-State scipy.sparse Two-Hole Fock Oracle (< 1e-13)"]:::inv',
            '  BDL4["4th-Order + j=3/2 Channel + 924-State Ring (< 7.4e-10 / 7.3e-7)"]:::inv',
            '  DAT["ACID Benchmark: exchange_vs_fock_oracle.dat"]:::inv',
            '  EXP["Export: extended_kitaev_material_couplings.dat"]:::inv',
            '  QMB["Merged PR #52: finite-QMB KitaevKanamoriSW.lean (7098b031)"]:::handoff',
            "  M21 --> A1",
            "  A1 --> CAS",
            "  M21 --> F66",
            "  CAS --> DAT",
            "  F66 --> DAT",
            "  M21 --> M22",
            "  M22 --> BDL4",
            "  BDL4 --> EXP",
            "  A1 --> QMB",
        ]
    )
    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>kitaev-derivation :: R&amp;D Blueprint Portal</title>
<script type="module">
  import mermaid from 'https://cdn.jsdelivr.net/npm/mermaid@10/dist/mermaid.esm.min.mjs';
  mermaid.initialize({{ startOnLoad: true, theme: 'neutral' }});
</script>
<style>
  body {{ background: #FAF6EE; color: #2B2B2B; font-family: Georgia, serif; max-width: 980px; margin: 2rem auto; padding: 0 1rem; line-height: 1.55; }}
  h1, h2 {{ color: #1E3A5F; }}
  code, pre {{ background: #EFE9DC; padding: 0.2rem 0.4rem; border-radius: 4px; font-size: 0.9rem; }}
  pre {{ padding: 0.8rem; overflow-x: auto; }}
  .card {{ border-left: 4px solid #8C4A27; padding: 0.6rem 1rem; background: #F5EFE2; margin: 1rem 0; }}
  .graph-box {{ background: #FFFFFF; border: 1px solid #D1C7B7; padding: 1rem; border-radius: 6px; overflow-x: auto; }}
</style>
</head>
<body>
<h1>kitaev-derivation — Symbolic Schrieffer–Wolff &amp; Fock Oracle Blueprint</h1>
<div class="card">
  <strong>Municipality Role (QSL Federation):</strong> Pure-Python (<code>SymPy</code>) CAR normal-ordering + Kanamori <code>P_[L=0,1,2]</code> multiplet projection + independent 66-state <code>scipy.sparse</code> Fock-space verification. Formal Lean 4 algebraic proofs route to <code>research/qmb/finite-QMB</code>.
</div>
<h2>1. Epistemic Dependency Graph (M2.1 &amp; Finite-QMB Handoff)</h2>
<div class="graph-box">
  <pre class="mermaid">
{html.escape(mermaid_code)}
  </pre>
</div>
<h2>2. Closed-Form Analytic Milestones (M2.1)</h2>
<ul>
  <li><strong>Jackeli–Khaliullin Limit:</strong> <code>J = 0, Gamma = 0, Gamma' = 0, K = -8 J_H t2^2 / (3 (U - 3 J_H)(U - J_H))</code> (spin-1/2), cancelling identically at <code>J_H = 0</code>.</li>
  <li><strong>Rau–Lee–Kee Symmetries:</strong> Exact <code>C2</code> isotropy (<code>J_xx = J_yy</code>), <code>C3</code> bond equivalence across <code>(x, y, z)</code>, and <code>Gamma = 16 J_H t2 (t1 - t3) / (9 (U - J_H)(U - 3 J_H))</code> at <code>t4 = 0</code>.</li>
  <li><strong>66-State Sparse Fock Oracle Parity:</strong> Symbolic <code>SymPy</code> and exact numerical 66-state Schrieffer–Wolff agree to <code>&lt; 1e-13</code> across randomized <code>(U, J_H, t1..t4)</code>.</li>
</ul>
<h2>3. Wave 3 Fourth-Order &amp; Ring Exchange (M2.2)</h2>
<p>Exact BDL/Takahashi 4th-order expansion separates infinite-SOC charge renormalization from finite-SOC <code>j_eff = 3/2</code> intermediate processes, while Peierls <code>Phi = 0 - Phi = pi</code> subtraction on the 924-state Kramers hexagon extracts the pure 6-site hexagon ring interaction.</p>
<h2>ACID Benchmark Tables</h2>
<pre>{dat_preview}</pre>
<pre>{exp_preview}</pre>
</body>
</html>
"""
    (out_dir / "index.html").write_text(page, encoding="utf-8")
    print(f"Built {out_dir / 'index.html'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
