#!/usr/bin/env python3
"""Intrinsic-mass drift under the canonical endpoint map — is the siren
discriminant really M_det = M_tilde/(1+z) vs (1+z)M_tilde?

Step 70 showed the siren mapping carries G_loc ∝ A^2(phi)(1+alpha^2).
Inside source interiors the k-mouflage field profile is additive:
phi_int = phi_amb + (local well profile), so the ambient drift passes
through to the interior conformal factor multiplicatively,

    A_int(z)/A_int(0) = A_amb(z)/A_amb(0) = 1/(1+z)  (canonical map)

=> G_loc inside source galaxies drifts as (1+z)^{-2 eps}, where eps = 1
is the k-mouflage additive-tracking estimate and eps = 0 would be
chameleon-style pinning to a density-dependent minimum (NOT the
specified P(X) theory — but kept as the sanity bracket).

Characteristic gravitational mass scales in baryon units scale as
G^{-3/2} (Chandrasekhar, TOV maximum, Jeans, core-collapse and
pair-instability channels all share this scaling at fixed EOS/baryon
content). A binary of N baryons has M_tilde = N m_p fixed — what drifts
is the characteristic N at which collapse happens:

    M_tilde*(z) = M_tilde*_0 (1+z)^{3 eps}.

Detector-inferred mass (step_21, A_o = 1, screened):

    M_det(z) = M_tilde*(z)/(1+z) = M_0 (1+z)^{3 eps - 1}
    vs GR:  M_det^GR(z) = M_0 (1+z).

Consequences:
- eps = 1 (k-mouflage interior tracking): M_det ∝ (1+z)^2 — GROWING
  faster than GR's (1+z). The discriminant survives but FLIPS
  DIRECTION: TEP predicts detector-frame masses rising with redshift
  twice as fast as GR, separation factor (1+z) not (1+z)^2. High-z
  BHs appearing 'too massive' (GW190521-class) become a natural
  signature of the weaker-G epoch, not an anomaly.
- eps = 2/3: the discriminant vanishes exactly — TEP and GR detector
  masses coincide.
- eps = 0 (pinned interiors): the naive claim, M_det ∝ (1+z)^{-1}.

Binding cross-sector constraint — SN Ia. A white dwarf accumulates to
the local M_Ch ∝ G^{-3/2}; peak luminosity L ∝ M_Ni ∝ M_Ch gives

    L(z)/L(0) = (1+z)^{3 eps}  =>  Delta m = -2.5*3eps*log10(1+z)

    eps=1: -1.3 mag at z=0.5, -2.3 mag at z=1.

This is a LARGE intrinsic brightening that runs against the observed
relative dimness unless the TEP distance/propagation bookkeeping
overcompensates; the TEP-H0 magnitude-response channel (kappa_SN)
responds to environment spatially, not to a secular (1+z) drift. The
standardization relations themselves (stretch-luminosity) also shift.
The SN sector therefore sets an independent, tighter constraint on the
effective eps — a dedicated TEP-H0 re-analysis with intrinsic drift is
the required next step, alongside BBN-era G_loc consistency (which is
bounded separately and not extrapolated here).

Second-order caveat noted: galactic well depths at z were shallower
(younger structures), A_gal,e/A_gal,o > 1 partially counteracts the
ambient tracking — suppressing eps below 1, not enhancing it.

Output: results/step_71_intrinsic_mass_drift.json
"""
import json
import os

import numpy as np


def run():
    zs = np.array([0.1, 0.3, 0.5, 1.0, 2.0])
    eps_cases = {
        "eps_1_kmouflage_additive": 1.0,
        "eps_critical": 2.0 / 3.0,
        "eps_half": 0.5,
        "eps_0_pinned": 0.0,
    }
    table = {}
    for name, eps in eps_cases.items():
        g_ratio = (1 + zs) ** (-2 * eps)
        mstar = (1 + zs) ** (3 * eps)
        mdet_tep = (1 + zs) ** (3 * eps - 1)
        mdet_gr = 1 + zs
        sep = mdet_tep / mdet_gr
        dsn = -2.5 * 3 * eps * np.log10(1 + zs)
        table[name] = {
            "eps": eps,
            "z": zs.tolist(),
            "G_loc_ratio": g_ratio.tolist(),
            "Mstar_over_M0": mstar.tolist(),
            "M_det_TEP_over_M0": mdet_tep.tolist(),
            "M_det_GR_over_M0": mdet_gr.tolist(),
            "TEP_over_GR_separation": sep.tolist(),
            "delta_SN_mag_intrinsic": dsn.tolist(),
        }
    out = {
        "premise": ("G_loc ∝ A^2 inside source interiors tracks the "
                    "ambient drift (k-mouflage additive field profile); "
                    "characteristic mass scales ∝ G^{-3/2}"),
        "eps_definition": ("interior tracking fraction of the ambient "
                           "cosmological drift; 1 = full additive "
                           "tracking (k-mouflage), 0 = pinned "
                           "(chameleon-style, NOT the specified P(X))"),
        "eps_critical": 2.0 / 3.0,
        "eps_critical_meaning": ("at eps = 2/3 the TEP and GR "
                                 "detector-mass trends coincide exactly; "
                                 "the discriminant vanishes"),
        "eps_1_verdict": ("M_det ∝ (1+z)^2 vs GR (1+z): the "
                          "discriminant survives but flips direction — "
                          "TEP predicts faster-growing detector masses; "
                          "high-z 'too-massive' BH population is a "
                          "natural signature"),
        "binding_constraint": ("SN Ia intrinsic luminosity L ∝ "
                               "(1+z)^{3eps}; at eps=1, -1.3 mag at "
                               "z=0.5 and -2.3 mag at z=1 — runs "
                               "against observed dimness unless the "
                               "distance bookkeeping overcompensates; "
                               "requires a TEP-H0 re-analysis with "
                               "intrinsic drift before the eps=1 "
                               "endpoint is claimed"),
        "shallower_wells_caveat": ("galactic wells at z were shallower; "
                                   "A_gal,e/A_gal,o > 1 partially "
                                   "counteracts tracking (pushes eps "
                                   "below 1, never above)"),
        "table": table,
    }
    return out


if __name__ == "__main__":
    res = run()
    outdir = os.path.join(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))), "results")
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, "step_71_intrinsic_mass_drift.json")
    with open(path, "w") as f:
        json.dump(res, f, indent=2)
    print("eps_critical:", res["eps_critical"], "-", res["eps_critical_meaning"])
    for name, t in res["table"].items():
        print(f"\n{name} (eps={t['eps']:.3f})")
        print("  z     G/G0   M*/M0  Mdet_TEP  Mdet_GR  sep    dSN(mag)")
        for i, z in enumerate(t["z"]):
            print(f"  {z:<5.1f} {t['G_loc_ratio'][i]:.4f} "
                  f"{t['Mstar_over_M0'][i]:8.2f} "
                  f"{t['M_det_TEP_over_M0'][i]:9.2f} "
                  f"{t['M_det_GR_over_M0'][i]:8.2f} "
                  f"{t['TEP_over_GR_separation'][i]:.3f} "
                  f"{t['delta_SN_mag_intrinsic'][i]:+.2f}")
    print(f"\nSaved {path}")
