#!/usr/bin/env python3
"""SN Ia luminosity-drift bound on the interior tracking fraction eps.

Step 71 showed that if interior fields track the ambient cosmological
drift by a fraction eps, then G_loc(z) ∝ (1+z)^{-2 eps} inside source
environments and characteristic mass scales drift as (1+z)^{3 eps}.
The sharpest cross-sector lever is the Chandrasekhar-mass channel: SN Ia
progenitors accumulate to the local M_Ch ∝ G^{-3/2}, so peak luminosity
drifts as

    L(z)/L(0) = (1+z)^{3 eps}   =>   dm(z) = -7.5 eps log10(1+z).

This signature is nearly linear in z at low z and cannot be absorbed by
the distance-ladder normalization (H0, M_B are global offsets). It is
therefore bounded directly by the SN Hubble-diagram residual floor:
Pantheon+ systematic/magnitude residual ~0.10-0.15 mag across z <= 1.

Result: eps <= ~0.05-0.10 (binding at the high-z end of the calibrated
baseline) — far below the critical eps_c = 2/3 at which the siren
discriminant vanishes and below full k-mouflage additive tracking
eps = 1 by an order of magnitude.

Consequences for the siren discriminant (step_70/71):
- Detector-mass exponent 3 eps - 1 stays in the range -0.70 to -0.85 at
  the SN bound — close to the naive -(1+z) map, not the (1+z)^2 flip.
- Separation vs GR is (1+z)^{1.7-1.85}: the discriminant survives,
  weakened modestly rather than reversed.

Interpretation for the architecture: interiors do NOT track the ambient
cosmological drift at anywhere near full strength over z <= 1-2. Two
consistent mechanisms in the corpus supply this decoupling:

- The canonical endpoint identity 1+z = A_o/A_e is a statement about
  ambient field VALUES at the endpoints of propagation; the deep-well
  interior field of a source galaxy is dominated by its own local
  depth, and the local well was SHALLOWER at earlier epochs (younger
  structures), pushing the effective interior shift below the ambient
  drift — eps < 1 structurally.
- The drift may sit preferentially in the ambient/void path sector
  (the path-integrated lapse channel already documented in section 9)
  rather than in the matter-hosted wells where sources live.

Either way, the SN ladder says interior G_loc is quasi-constant in
fractional terms over the SN baseline — which is also what the corpus's
own SN, galaxy-dynamics, and stellar-physics sectors tacitly require.
The remaining task is a derivation of eps from the actual interior
solution (the same compact-star/interior solver that supplies
sensitivities), not a free fit.

Outputs: results/step_73_sn_luminosity_bound.json
"""
import json
import os

import numpy as np


def run():
    zs = np.array([0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0])

    # Residual floors: Pantheon+ per-bin systematic ~0.1 mag; adopt
    # 0.10 and 0.15 as the bracket.
    bounds = {}
    for floor in (0.10, 0.15):
        eps_b = floor / (7.5 * np.log10(1 + zs))
        bounds[f"floor_{floor:.2f}_mag"] = {
            "z": zs.tolist(),
            "eps_bound": eps_b.tolist(),
            "tightest": float(eps_b.max() if False else eps_b[-1]),
        }
    # the operative bound is set by the deepest reliable SN baseline;
    # z=1 is the standard upper edge of the Pantheon+ calibrated range.
    i1 = int(np.where(zs == 1.0)[0][0])
    eps_010 = bounds["floor_0.10_mag"]["eps_bound"][i1]   # 0.044
    eps_015 = bounds["floor_0.15_mag"]["eps_bound"][i1]   # 0.066
    i05 = int(np.where(zs == 0.5)[0][0])
    eps_010_z05 = bounds["floor_0.10_mag"]["eps_bound"][i05]
    eps_015_z05 = bounds["floor_0.15_mag"]["eps_bound"][i05]

    eps_bracket = [0.04, 0.07, 0.10]
    discr = []
    for eps in eps_bracket:
        discr.append({
            "eps": eps,
            "M_det_exponent": 3 * eps - 1,
            "separation_vs_GR_exponent": 2 - 3 * eps,
            "dSN_mag_at_z1": -7.5 * eps * np.log10(2.0),
        })

    out = {
        "signature": "dm(z) = -7.5 eps log10(1+z); unabsorbable by "
                     "ladder normalization (z-shape, not offset)",
        "residual_floors_mag": [0.10, 0.15],
        "eps_bounds": {
            "z05": {"floor_0.10": eps_010_z05, "floor_0.15": eps_015_z05},
            "z1": {"floor_0.10": eps_010, "floor_0.15": eps_015},
        },
        "eps_allowed_range": ("eps <= ~0.05-0.10: the binding point is "
                              "the high-z end of the calibrated SN "
                              "baseline (eps <= 0.04-0.07 at z=1; "
                              "0.08-0.11 at z=0.5)"),
        "eps_full_tracking": 1.0,
        "eps_critical": 2.0 / 3.0,
        "eps_excluded_by_SN": ("eps = 1 predicts -2.26 mag at z=1, "
                               "~15-23x the residual floor: full "
                               "k-mouflage interior tracking is "
                               "falsified by the SN ladder"),
        "discriminant_at_bound": discr,
        "verdict": ("SN residuals bound the interior tracking fraction "
                    "to eps <= ~0.05-0.10 — far below the eps=2/3 "
                    "neutralization and the eps=1 flip. The siren "
                    "discriminant survives at near-naive strength: "
                    "M_det ~ (1+z)^{-0.70..-0.85} vs GR's (1+z), "
                    "separation ~(1+z)^{1.7-1.85}. Interiors are "
                    "effectively pinned against the ambient drift over "
                    "the SN baseline; deriving the small residual eps "
                    "from the interior solution (local well depth vs "
                    "ambient drift, plus shallower early wells) is the "
                    "remaining item, shared with the compact-star "
                    "solver."),
    }
    return out


if __name__ == "__main__":
    res = run()
    outdir = os.path.join(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))), "results")
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, "step_73_sn_luminosity_bound.json")
    with open(path, "w") as f:
        json.dump(res, f, indent=2)
    print(res["verdict"])
    for d in res["discriminant_at_bound"]:
        print(f"  eps={d['eps']}: M_det ~ (1+z)^{d['M_det_exponent']:+.2f}, "
              f"sep (1+z)^{d['separation_vs_GR_exponent']:.2f}, "
              f"dSN@z1={d['dSN_mag_at_z1']:+.2f} mag")
    print(f"Saved {path}")
