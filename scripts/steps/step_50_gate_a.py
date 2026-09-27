#!/usr/bin/env python3
"""Gate A: one action, four closures.

1. One spherical profile phi(r). S_Sigma is the gradient projection,
   S_A is the clock projection. Both are read off that solution.
2. B(u) = B0 * u^2/(1+u^2) * exp(-u^4/2) is a family. B0 is not
   chosen. B0 >= 0 from the admissible branch. The upper end is the
   GW170817 cone split 1e-15 evaluated on this profile.
3. P(X) = X + X * sqrt(X^2 + eps) / Lambda^4. P_X > 0 at X=0.
   The sound speed stays positive down to eps -> 0.
4. kappa at the Cepheid radius is kappa_canonical * S_A(r_cep).
   That is the transfer map. No extra coefficient.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

C = 2.99792458e8
G = 6.67430e-11
MSUN = 1.98847e30
PC = 3.085677581e16
H0 = 70e3 / (1e6 * PC)
G_T = C * H0 / 2.0                 # beta_A^2 = 1
KAPPA_CANONICAL = 0.96e6
KAPPA_CEP_LO = 0.326e6
KAPPA_CEP_HI = 0.452e6
# Host used by step_05.
M_GAL = 1.0e11 * MSUN
R_GAL = 30.0e3 * PC
R_CEP = 8.0e3 * PC
# A Cepheid star, for the second projection of the same law.
M_STAR = 8.0 * MSUN
R_STAR = 50.0 * 6.957e8


def y_of_g(g):
    g = np.asarray(g, dtype=float)
    q = (g / G_T) ** 2
    y = np.empty_like(g)
    small = q < 1e-8
    y[small] = 1.0 - q[small]
    qq = q[~small]
    disc = (0.5 / qq) ** 2 + (1.0 / (3.0 * qq)) ** 3
    s = np.sqrt(np.maximum(disc, 0.0))
    y[~small] = np.cbrt(0.5 / qq + s) + np.cbrt(0.5 / qq - s)
    return y


def g_uniform(r, mass, radius):
    r = np.asarray(r, dtype=float)
    g = np.empty_like(r)
    inn = r < radius
    g[inn] = G * mass * r[inn] / radius ** 3
    g[~inn] = G * mass / r[~inn] ** 2
    return g


def phi_of_profile(r, mass, radius):
    """Screened and unscreened dimensionless fields, integrated inward
    from a flat infinity. du/dr = -y g/c^2 (u decreases outward)."""
    r = np.asarray(r, dtype=float)
    # integrate from the outside in so the boundary is u(inf)=0
    order = np.argsort(r)
    rr = r[order]
    # extend to a large outer radius
    r_out = max(rr[-1] * 20.0, radius * 50.0)
    grid = np.geomspace(rr[0], r_out, 20000)
    g = g_uniform(grid, mass, radius)
    y = y_of_g(g)
    # cumulative integral from the outside: u(r) = int_r^inf y g/c^2 dr
    # trap on reversed grid
    integrand = y * g / C ** 2
    integrand_u = g / C ** 2
    # reverse cumulative
    dr = np.diff(grid)
    du = 0.5 * (integrand[1:] + integrand[:-1]) * dr
    du_u = 0.5 * (integrand_u[1:] + integrand_u[:-1]) * dr
    u_rev = np.cumsum(du[::-1])[::-1]
    uu_rev = np.cumsum(du_u[::-1])[::-1]
    u = np.zeros_like(grid)
    uu = np.zeros_like(grid)
    u[:-1] = u_rev
    uu[:-1] = uu_rev
    # interpolate back
    u_r = np.interp(rr, grid, u)
    uu_r = np.interp(rr, grid, uu)
    y_r = y_of_g(g_uniform(rr, mass, radius))
    out_u = np.empty_like(r)
    out_uu = np.empty_like(r)
    out_y = np.empty_like(r)
    out_u[order] = u_r
    out_uu[order] = uu_r
    out_y[order] = y_r
    return out_u, out_uu, out_y


def projections(r, mass, radius):
    u, u_uns, y = phi_of_profile(np.atleast_1d(r), mass, radius)
    # S_Sigma: gradient projection. Screened gradient / Newtonian gradient.
    s_sigma = y
    # S_A: clock projection. (A-1)/(A_uns-1) -> u/u_uns for the small fields here.
    s_a = np.where(u_uns > 0, u / u_uns, 1.0)
    return u, u_uns, s_sigma, s_a


def regulate_px():
    """P = X + X*sqrt(X^2+eps)/Lam^4.
    P_X = 1 + (2 X^2 + eps) / (Lam^4 sqrt(X^2+eps)) >= 1.
    c_s^2 = P_X / (P_X + 2 X P_XX)."""
    lam4 = 1.0
    rows = []
    for eps in (1e-2, 1e-6, 1e-12, 0.0):
        xs = np.array([-2.0, -0.1, -1e-6, 0.0, 1e-6, 0.1, 2.0])
        px_min = np.inf
        cs_min = np.inf
        for X in xs:
            if eps == 0.0 and X == 0.0:
                px = 1.0
                # limiting P_XX at 0 is finite only for eps>0; the eps=0
                # theory is checked off zero and the limit of px is 1.
                cs = 1.0
            else:
                e = eps if eps > 0 else 0.0
                root = np.sqrt(X * X + e)
                px = 1.0 + (2.0 * X * X + e) / (lam4 * root)
                # P_XX = d/dX [(2X^2+e)/(Lam root)]
                # num = 2X^2+e, den = root
                # d(num/den) = (4X den - num * X/den) / den^2
                d = (4.0 * X * root - (2.0 * X * X + e) * (X / root)) / (root * root)
                pxx = d / lam4
                denom = px + 2.0 * X * pxx
                cs = px / denom if denom != 0 else np.nan
            px_min = min(px_min, px)
            cs_min = min(cs_min, cs)
        rows.append({
            "eps": eps,
            "P_X_min_on_sample": float(px_min),
            "c_s2_min_on_sample": float(cs_min),
            "ghost_free": bool(px_min > 0),
            "hyperbolic": bool(cs_min > 0),
        })
    return rows


def b_family(r, u, y, g):
    """Cone split sigma = (B/A^2) (du/dr * R_H)^2, A~1.
    B = B0 * shape(u). GW170817: |sigma| < 1e-15."""
    shape = u ** 2 / (1.0 + u ** 2) * np.exp(-0.5 * u ** 4)
    R_H = C / H0
    grad = y * g / C ** 2          # du/dr
    sigma_per_B0 = shape * (grad * R_H) ** 2
    peak = float(np.max(sigma_per_B0))
    b0_max = 1.0e-15 / peak if peak > 0 else None
    # Holonomy at B0=1 is the stored step_15 admissible-branch value.
    h_unit = 5.977e-13
    return {
        "form": "B(u)=B0*u^2/(1+u^2)*exp(-u^4/2), B0>=0",
        "sigma_per_B0_peak": peak,
        "B0_max_from_GW170817": b0_max,
        "holonomy_seconds_at_B0": {
            "0": 0.0,
            "B0_max": None if b0_max is None else h_unit * b0_max,
            "1_excluded_if_above_max": h_unit,
        },
    }


def main():
    r_cep = np.array([R_CEP])
    u, uu, s_sig, s_a = projections(r_cep, M_GAL, R_GAL)
    g_cep = float(g_uniform(r_cep, M_GAL, R_GAL)[0])
    # stellar surface, same law
    u_s, uu_s, ss_s, sa_s = projections(np.array([R_STAR]), M_STAR, R_STAR)
    kappa = KAPPA_CANONICAL * float(s_a[0])
    emp_mid = 0.5 * (KAPPA_CEP_LO + KAPPA_CEP_HI)
    # B family on a radial grid of the host
    r_grid = np.geomspace(0.05 * R_GAL, 5 * R_GAL, 400)
    ug, uug, yg = phi_of_profile(r_grid, M_GAL, R_GAL)[:3]
    gg = g_uniform(r_grid, M_GAL, R_GAL)
    fam = b_family(r_grid, ug, yg, gg)
    kin = regulate_px()
    if not all(row["ghost_free"] and row["hyperbolic"] for row in kin):
        raise SystemExit(f"kinetic completion failed: {kin}")
    out = {
        "step": "step_50_gate_a",
        "one_profile": {
            "configuration": "uniform sphere, 1e11 Msun, 30 kpc, Cepheid at 8 kpc",
            "g_m_s2": g_cep,
            "g_over_g_t": g_cep / G_T,
            "u": float(u[0]),
            "u_unscreened": float(uu[0]),
            "S_Sigma": float(s_sig[0]),
            "S_A": float(s_a[0]),
            "same_phi": True,
        },
        "stellar_surface_same_law": {
            "S_Sigma": float(ss_s[0]),
            "S_A": float(sa_s[0]),
            "u": float(u_s[0]),
        },
        "kappa": {
            "kappa_canonical": KAPPA_CANONICAL,
            "S_A_host": float(s_a[0]),
            "kappa_from_profile": kappa,
            "empirical_interval": [KAPPA_CEP_LO, KAPPA_CEP_HI],
            "empirical_midpoint": emp_mid,
            "canonical_over_empirical_mid": KAPPA_CANONICAL / emp_mid,
            "canonical_over_profile": KAPPA_CANONICAL / kappa if kappa else None,
        },
        "B_family": fam,
        "P_of_X": {
            "definition": "P=X+X*sqrt(X^2+eps)/Lambda^4",
            "samples": kin,
        },
    }
    dest = Path(__file__).resolve().parents[2] / "results" / "step_50_gate_a.json"
    dest.write_text(json.dumps(out, indent=2))
    print(json.dumps({
        "S_Sigma": out["one_profile"]["S_Sigma"],
        "S_A": out["one_profile"]["S_A"],
        "g_over_gt": out["one_profile"]["g_over_g_t"],
        "kappa_from_profile": kappa,
        "ratio_canon_over_emp": out["kappa"]["canonical_over_empirical_mid"],
        "B0_max": fam["B0_max_from_GW170817"],
        "kinetic_ok": True,
    }, indent=2))


if __name__ == "__main__":
    main()
