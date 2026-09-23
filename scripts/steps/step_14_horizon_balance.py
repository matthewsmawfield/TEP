#!/usr/bin/env python3
"""
step_14_horizon_balance.py
==========================
The residual geometry test.  The audit chain established:

  * archived step_27: A^-1-weighted identity 3<A^-1 Pi> + 4<A^-1 v.grad l>
    = 0 is the correct balance; conformal flow supplies ~1e-8 at
    corpus amplitudes; endpoint redshift survives (z = Pi Dt).
  * step_13: B<0 required; universal roll excluded (Q -> -1);
    partition structurally mandatory.
  * archived step_29: unforced tracking gives opposite-sign partition but
    ~1e-7 H0 -- the roll must be a V-driven homogeneous mode.
  * multimessenger bound: |Q| along photon paths <=~ 1e-15
    (GW170817); anisotropic light-speed shifts inside observed
    wells are excluded at comparable or tighter levels -- the
    disformal volume term cannot carry the balance anywhere
    exterior observers or endpoints live.

The only region remaining is the DEEP INTERIOR: behind the
temporal horizon (A -> 0), causally disconnected from exterior
measurement.  The corrected identity is A^-1-weighted, and A^-1
diverges at the horizon -- the balance may be saturated there.
Near the horizon all three amplifiers compound:
    A^-1 = e^{delta}     (temporal depth)
    |grad lnA| ~ delta/r  (steepening profile)
    v_in -> c            (relativistic infall)

This script integrates the weighted identity terms over a
composite well ensemble -- shallow Newtonian-tracking halo
(delta ~ 1e-6) matching to a power-law deep interior
(delta ~ delta_c (r_c/r)^p) truncated at the horizon r_h where
delta_H is set by the horizon asymptotics -- and reports the
depth-stratified cumulative of each term:

  drift term:  T_D(r)  = 3 Pi_bar e^{delta(r)}
  flow term:   T_F(r)  = 4 v_in(r) |grad lnA| e^{delta(r)}

Question: does the ratio T_F/T_D (the balance fraction carried)
accumulate to O(1) in the near-horizon region while staying
negligible at observable depths (delta <~ 1)?
"""
import json
import numpy as np
import os

C_SI = 2.998e8
H0 = 2.27e-18
MPC = 3.086e22
GYR = 3.156e16
G_SI = 6.674e-11
M_SUN = 1.989e30


def well_terms(r, delta_prof, v_prof, Pi_bar):
    """Per-unit-coordinate-volume integrands of the two weighted
    terms, as a function of radius inside one well."""
    Ainv = np.exp(delta_prof)                    # e^{-lnA}=e^{delta}
    g = np.abs(np.gradient(delta_prof, r))       # |grad lnA| [1/m]
    T_D = 3.0 * Pi_bar * Ainv
    T_F = 4.0 * v_prof * g * Ainv
    return T_D, T_F


def single_well(r_h, delta_H, p, r_c, delta_c,
                Pi_bar, v_ratio, n=4000):
    """Power-law interior delta(r) = delta_c (r_c/r)^p truncated
    at r_h where delta = delta_H (horizon asymptotics); v_infall
    ramping to v_ratio*c at the horizon (free-fall-like
    v ~ c sqrt(r_h/r)).  Returns cumulative integrals."""
    r_min = r_h * 1.0                            # integrate to r_h
    r_max = r_c * 10.0                           # outer matching zone
    r = np.logspace(np.log10(r_min), np.log10(r_max), n)
    delta = delta_c * (r_c / r) ** p
    delta = np.minimum(delta, delta_H)           # horizon floor
    v = v_ratio * C_SI * np.sqrt(np.clip(r_h / r, 0, 1))
    T_D, T_F = well_terms(r, delta, v, Pi_bar)
    dV = 4 * np.pi * r**2                        # coord volume elem
    cum_D = np.cumsum(T_D[::-1] * dV[::-1] *
                     np.abs(np.gradient(r)[::-1]))[::-1]
    cum_F = np.cumsum(T_F[::-1] * dV[::-1] *
                      np.abs(np.gradient(r)[::-1]))[::-1]
    return r, delta, T_D, T_F, cum_D, cum_F, dV


def main():
    res = {"units": "SI",
           "model": "composite well ensemble; depth-stratified "
                    "A^-1-weighted balance integral"}
    L_MPC = 100.0
    L3 = (L_MPC * MPC)**3

    out = {}
    # --- ensemble scan: how deep must interiors go? ---
    # galaxy-scale halo: delta_halo ~ 1e-6 shallow; deep core
    # delta_c at r_c; power-law index p; horizon at r_h with
    # delta_H ~ asymptotic depth (TEP-BH A->0).
    for tag, delta_H, p, rho_halo in [
            ("dH=20,p=2", 20.0, 2.0, 1e-6),
            ("dH=30,p=2", 30.0, 2.0, 1e-6),
            ("dH=30,p=1", 30.0, 1.0, 1e-6),
            ("dH=40,p=2", 40.0, 2.0, 1e-6)]:
        # a single representative well; weighting by well density
        # done below.  r_c: interior scale; r_h: horizon radius.
        r_c = 10.0 * MPC                          # 10 Mpc domain
        r_h = 1e-5 * MPC                          # ~0.3 pc horizon-ish
        delta_c = rho_halo                        # halo depth
        # solve the delta at r_c to be halo depth; profile
        # delta = delta_c*(r_c/r)^p reaches delta_H at r_h only if
        # delta_c*(r_c/r_h)^p = delta_H -> adjust normalization:
        delta_c_eff = delta_H * (r_h / r_c) ** p
        r, delta, T_D, T_F, cum_D, cum_F, dV = single_well(
            r_h, delta_H, p, r_c, delta_c_eff, H0, 0.1)
        # balance fraction interior to radius r: cum_F/cum_D at r_h
        # (cum arrays are cumulative from r_min inward)
        frac_in = cum_F / np.maximum(cum_D, 1e-300)
        # depth-stratified: fraction of total balance term
        # contributed at delta > threshold
        tot_F = cum_F[0]; tot_D = cum_D[0]
        def frac_below(depth):
            m = delta <= depth
            if not m.any():
                return 0.0
            return float(cum_F[m][-1] / tot_F) if tot_F else 0.0
        # observable-depth share of the flow term
        obs_share = frac_below(1e-3)
        out[tag] = {
            "balance_frac_integrated": float(tot_F / tot_D)
                if tot_D else None,
            "flow_term_share_at_delta<1e-3": obs_share,
            "depth_at_90pct_balance": float(
                delta[np.argmax(cum_F > 0.9 * tot_F)])
                if tot_F > 0 else None,
            "delta_H": delta_H, "p": p,
            "well_volume_frac": float(
                np.sum(dV * np.abs(np.gradient(r))) / L3),
        }
        print(tag, json.dumps(out[tag]), flush=True)

    res["scan"] = out
    res["verdict"] = {"statement":
        "closure by horizon saturation iff integrated "
        "flow/drift ratio -> O(1) with the overwhelming share at "
        "delta >> 1 while the observable-depth share stays "
        "consistent with endpoint redshift scatter"}
    outdir = os.path.join(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir,
                           "step_14_horizon_balance.json"), "w") as f:
        json.dump(res, f, indent=2)
    print(json.dumps(res["verdict"], indent=1))


if __name__ == "__main__":
    main()
