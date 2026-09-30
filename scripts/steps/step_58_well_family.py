#!/usr/bin/env python3
"""Temporal-well family solve (plan T-B1): settles Rule 20.

Does the central clock rate A(u_c) = e^{-u_c} of a compact source reach
a finite minimum N_min > 0 (floor-capped well) or fall asymptotically
toward zero (approach, never a completed halt)?

Family: fixed total scalar charge Q = 4 pi (same convention as
step_32/56), uniform-density sphere of radius R swept toward the
compact limit R -> 0.  Static spherical field equation

    (1/r^2) d/dr [ r^2 f(|u'|) u' ] = S(r;u)
    S = rho_shape(r) e^{-u} - C_V v'(u)

    rho_shape = 3/R^3 for r < R  (∫ rho d^3x = 4 pi), 0 outside.

e^{-u} is the source-starvation factor (the matter coupling shuts off
as clocks slow); v'(u) is the master-potential restoring term, which
dies at large u on both branches (quartic envelope e^{-(u/u_s)^4} kills
the matter term; the floor V_0 e^{-(u_s/u)^4} -> const), so in the
source-dominated deep well it is subleading and the asymptote question
is decided by the kinetic sector + starvation alone.

Method: Picard iteration with lagged flux law -- at fixed u-profile the
momentum equation is exactly integrable:

    f(q) q = F(r)/r^2 ,   F(r) = ∫_0^r s^2 S(s) ds ,  q = -u' >= 0
    u(r) = u_amb + ∫_r^{Rmax} q(s) ds

converges in a few iterations because e^{-u} damps the nonlinear
feedback.  Runs both sectors (baseline f = 1+q^2, two-branch
f = k q/sqrt(2) + q^2) and both potentials (none / master family).
"""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

K = 16.03
SQRT2 = np.sqrt(2.0)

# master potential (step_54 parametrisation), slope v'(u):
#   V = (l u^4/4) e^{-(u/u_s)^4} + V0 e^{-(u_s/u)^4}
# u_s = 10 per the step_54-declared family (R3); audited 2025-09-28:
# u_c ~ lnln(1/R) ~ 1-2 stays far below u_s over the whole family, so
# the result is insensitive to the choice (verified U_S 3 vs 10 ->
# identical u_c); keeping the declared value for corpus consistency.
U_S = 10.0
V0 = 1.0
LAM = 1.0


def vmatter(u):
    return (LAM * u ** 4 / 4.0) * np.exp(-(u / U_S) ** 4)


def vfloor(u):
    u = np.maximum(u, 1e-30)
    return V0 * np.exp(-(U_S / u) ** 4)


def vprime(u):
    """dV/du of the master family, analytic."""
    u = np.asarray(u, dtype=float)
    uc = np.maximum(u, 1e-30)
    e1 = np.exp(-(u / U_S) ** 4)
    e2 = np.exp(-(U_S / uc) ** 4)
    dm = LAM * u ** 3 * e1 * (1.0 - (u / U_S) ** 4)
    df = V0 * e2 * 4.0 * U_S ** 4 / uc ** 5
    return dm + df


def f_flux(q, sector):
    if sector == "baseline":
        return 1.0 + q * q
    return K * q / SQRT2 + q * q


def invert_flux(g, sector):
    """Solve q f(q) = g for q >= 0 (g >= 0)."""
    if g <= 0:
        return 0.0
    hi = max(1.0, g)
    while hi * f_flux(hi, sector) < g:
        hi *= 4.0
    return brentq(lambda q: q * f_flux(q, sector) - g, 0.0, hi,
                  xtol=1e-14, rtol=1e-13)


def solve_well(R, sector, use_V, C_V, nr=4000, rmax=40.0, itmax=200):
    """Fixed-charge source, radius R, r* units.  Returns central field.
    Log-spaced grid so the source interior is resolved at every R."""
    r = np.exp(np.linspace(np.log(R * 1e-4), np.log(rmax), nr))
    dr = np.diff(r)
    shape = np.where(r < R, 3.0 / R ** 3, 0.0)
    u = np.zeros(nr)
    for it in range(itmax):
        src = shape * np.exp(-u)
        if use_V:
            src = src - C_V * vprime(u)
        F = np.concatenate([[0.0],
                            np.cumsum(0.5 * (src[:-1] * r[:-1] ** 2
                                             + src[1:] * r[1:] ** 2)
                                      * dr)])
        g = np.maximum(F, 0.0) / np.maximum(r ** 2, 1e-30)
        q = np.array([invert_flux(gi, sector) for gi in g])
        # u(r) = ∫_r^rmax q ds  (suffix trapezoid; u(rmax) = 0)
        integ = 0.5 * (q[:-1] + q[1:]) * dr
        u_new = np.concatenate((np.cumsum(integ[::-1])[::-1], [0.0]))
        du = np.max(np.abs(u_new - u))
        u = u_new
        if du < 1e-10:
            break
    return float(u[0]), int(it), r, u


def family(sector, use_V, C_V):
    rows = []
    for R in [4.0, 2.0, 1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01,
              0.005, 0.002, 0.001]:
        uc, it, r, u = solve_well(R, sector, use_V, C_V)
        rows.append({"R_over_rstar": R, "u_central": uc,
                     "N_min": float(np.exp(-uc)), "iters": it})
    return rows


def growth_law(rows):
    """u_c vs ln(1/R) and vs ln ln(1/R): asymptotic-growth diagnostic.
    Double-log slope ~const with stable fit = ln ln divergence
    (approach, never a completed halt); saturating = finite floor."""
    x = np.log([1.0 / r["R_over_rstar"] for r in rows])
    y = np.array([r["u_central"] for r in rows])
    b, a = np.polyfit(x, y, 1)
    xp = x[-4:]
    yp = y[-4:]
    bp, ap = np.polyfit(xp, yp, 1)
    # double-log fit over the compact half
    xl = np.log(x[x > 0])
    yl = y[x > 0]
    bl, al = np.polyfit(xl, yl, 1)
    resid_ll = float(np.max(np.abs(yl - (al + bl * xl))))
    resid_lin = float(np.max(np.abs(yl - (np.polyval(
        np.polyfit(x[x > 0], yl, 1), x[x > 0])))))
    return {"log_slope_full": float(b), "log_slope_compact4": float(bp),
            "loglog_slope": float(bl),
            "resid_lnln_fit": resid_ll, "resid_ln_fit": resid_lin,
            "better": "lnln" if resid_ll < resid_lin else "ln"}


def main():
    out = {"k": K, "convention":
           "fixed charge Q=4pi, uniform sphere radius R -> 0; "
           "N_min = exp(-u_c); source starvation rho e^{-u}; "
           "units r_* = 1, g_t = 1"}
    for sector in ("two_branch", "baseline"):
        out[sector] = {
            "source_only": family(sector, use_V=False, C_V=0.0),
            "master_V_CV_1e-3": family(sector, use_V=True, C_V=1e-3)}
        out[sector]["growth_source_only"] = \
            growth_law(out[sector]["source_only"])
        out[sector]["growth_master_V"] = \
            growth_law(out[sector]["master_V_CV_1e-3"])
    dest = Path(__file__).resolve().parents[2] / "results" / \
        "step_58_well_family.json"
    dest.write_text(json.dumps(out, indent=2))
    for sector in ("two_branch", "baseline"):
        print(sector)
        for r in out[sector]["source_only"]:
            print("  R=%.2f  u_c=%.3f  N_min=%.3e" %
                  (r["R_over_rstar"], r["u_central"], r["N_min"]))
        print("  growth:", out[sector]["growth_source_only"])
    print("wrote", dest)


if __name__ == "__main__":
    main()
