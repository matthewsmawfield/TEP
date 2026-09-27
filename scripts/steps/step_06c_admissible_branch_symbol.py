#!/usr/bin/env python3
"""Principal symbol of the disformal perfect-fluid sector on B ≥ 0.

The matter metric is g̃_μν = A² g_μν + B ∂_μφ ∂_νφ. Its inverse exists
and is Lorentzian when A > 0 and 1 + (B/A²) X ≠ 0 with
X = g^{μν} ∂_μφ ∂_νφ, which is the paper's condition B X > -A².

For a perfect fluid of equation-of-state w coupled through g̃, the
scalar principal operator reduced on that background has coefficients

    Z_t = 1 + B (1+w) ρ / A²
    Z_s = 1 + B (1-w) ρ / (3 A²)

(these are the coefficients already identified in step_06). Strong
hyperbolicity of this sector requires Z_t > 0 and Z_s > 0 together with
a non-degenerate matter metric.

On the admissible branch B ≥ 0, for -1 ≤ w ≤ 1 and ρ ≥ 0, A² > 0:

    Z_t ≥ 1,    Z_s ≥ 1.

Both are strictly positive. The elliptic counterexamples in step_06 use
B < 0, which is the branch excluded by the null-cone condition. This
does not extend the claim to the noncanonical P(X) completion or to a
black-hole interior where |(B/A²) X| is not small; those remain separate
characteristic problems. It does close the disformal perfect-fluid
symbol on the branch the late-time theory uses.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np


def coefficients(B, w, rho, A2):
    Z_t = 1.0 + B * (1.0 + w) * rho / A2
    Z_s = 1.0 + B * (1.0 - w) * rho / (3.0 * A2)
    return float(Z_t), float(Z_s)


def main():
    ws = np.linspace(-1.0, 1.0, 9)
    rhos = np.array([0.0, 1e-30, 1e-6, 1.0, 1e6, 1e20])
    Bs = np.array([0.0, 1e-8, 1e-3, 1.0, 10.0])
    A2 = 1.0
    samples = []
    worst_Zt, worst_Zs = np.inf, np.inf
    for B in Bs:
        for w in ws:
            for rho in rhos:
                Zt, Zs = coefficients(float(B), float(w), float(rho), A2)
                worst_Zt = min(worst_Zt, Zt)
                worst_Zs = min(worst_Zs, Zs)
                if Zt <= 0 or Zs <= 0:
                    raise SystemExit(f"admissible branch not hyperbolic at B={B}, w={w}, rho={rho}")
    # Excluded branch really does go elliptic, so the sign is doing the work.
    Zt_bad, Zs_bad = coefficients(-5.0, 0.0, 1.0, 1.0)
    if not (Zt_bad < 0 and Zs_bad < 0):
        raise SystemExit("excluded-branch counterexample did not go elliptic; formula drifted")
    # Cone condition at the GW170817-sized deformation.
    cone_split = 1e-15
    out = {
        "step": "step_06c_admissible_branch_symbol",
        "branch": "B >= 0",
        "worst_Z_t": worst_Zt,
        "worst_Z_s": worst_Zs,
        "n_samples": int(len(Bs) * len(ws) * len(rhos)),
        "excluded_branch_example": {"B": -5.0, "w": 0.0, "rho": 1.0,
                                    "Z_t": Zt_bad, "Z_s": Zs_bad},
        "cone_nondegeneracy": (
            "1 + (B/A^2) X stays positive when |(B/A^2) X| is bounded by "
            f"the late-time cone split ~{cone_split:.0e}"
        ),
        "not_claimed": [
            "noncanonical P(X) characteristic speeds",
            "black-hole interior where |(B/A^2) X| is not small",
        ],
        "statement": (
            "On B>=0, -1<=w<=1, rho>=0 the disformal perfect-fluid "
            "principal coefficients satisfy Z_t>=1 and Z_s>=1. The "
            "elliptic examples sit on B<0."
        ),
    }
    dest = Path(__file__).resolve().parents[2] / "results" / "step_06c_admissible_branch_symbol.json"
    dest.write_text(json.dumps(out, indent=2))
    print(json.dumps({k: out[k] for k in ("worst_Z_t", "worst_Z_s", "n_samples",
                                          "excluded_branch_example")}, indent=2))
    print(f"wrote {dest}")


if __name__ == "__main__":
    main()
