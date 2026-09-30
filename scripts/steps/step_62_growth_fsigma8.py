#!/usr/bin/env python3
"""AUD-3 / T3.1 (first pass) — quasi-static linear growth under the
derived operator.

The conformal-image bookkeeping keeps the ambient map
LambdaCDM-equivalent over the constrained range (AUD-4 branch a):
E^2(z) = Om (1+z)^3 + OL.  The operator enters through the
environmental suppression of the gravitational coupling,

    S_Sigma,cosmo(z) = 1 / P_X(xi_cosmo),   xi_cosmo = H^2/(2 H0^2)
    G_eff/G = 1 + 2 beta_A^2 S_Sigma

(plan AUD-3: phi_dot = -M_Pl H => X = M_Pl^2 H^2/2 => xi = H^2/(2H0^2)).

Growth ODE on the conformal-image background (the standard
modified-gravity projection; the full TEP constraint-slice growth
equation — static g, no physical Hubble friction — is AUD-6/T6.3,
flagged below):

    d^2 D / d ln a^2 + [2 + d ln H / d ln a] dD/d ln a
        = (3/2) Om (1+z)^3 / E^2 * (G_eff/G) D

Outputs D(z), f(z), f sigma_8(z), S_8 for baseline / two_branch /
exp_interp, confronted with representative RSD + weak-lensing values.

Outputs results/step_62_growth_fsigma8.json
"""
import json
import os

import numpy as np
from scipy.integrate import solve_ivp

OM, OL = 0.3, 0.7
S8_PLANCK = 0.834
BETA_A = -1.0
K = 16.03


def PX(xi, sector):
    xi = np.maximum(np.asarray(xi, dtype=float), 0.0)
    if sector == "lcdm":
        return np.inf * np.ones_like(xi)
    if sector == "baseline":
        return 1.0 + 2.0 * xi
    if sector == "two_branch":
        return K * np.sqrt(xi) + 2.0 * xi
    if sector == "exp_interp":
        return 1.0 - np.exp(-K * np.sqrt(xi)) + 2.0 * xi
    raise ValueError(sector)


def Ez(z):
    return np.sqrt(OM * (1.0 + z) ** 3 + OL)


def G_eff(z, sector):
    xi = Ez(z) ** 2 / 2.0
    return 1.0 + 2.0 * BETA_A ** 2 / PX(xi, sector)


def growth(sector):
    """Solve for D(a), normalized D(a_i) = a_i at z_i = 200.
    sector='lcdm' forces G_eff = 1 (LambdaCDM control)."""
    ai = 1.0 / 201.0

    def rhs(lna, y):
        a = np.exp(lna)
        z = 1.0 / a - 1.0
        dlnH = -(1.5 * OM * (1.0 + z) ** 3) / Ez(z) ** 2  # dlnH/dlna
        src = 1.5 * OM * (1.0 + z) ** 3 / Ez(z) ** 2 * G_eff(z, sector)
        return [y[1], src * y[0] - (2.0 + dlnH) * y[1]]

    sol = solve_ivp(rhs, [np.log(ai), 0.0], [ai, ai],
                    rtol=1e-10, atol=1e-14, dense_output=True)
    return sol


def gamma_fit(zgrid, f):
    """f = Om_m(z)^gamma: fit gamma on 0 < z < 2."""
    Om_z = OM * (1.0 + zgrid) ** 3 / Ez(zgrid) ** 2
    mask = (zgrid > 0.05) & (zgrid < 2.0)
    return float(np.polyfit(np.log(Om_z[mask]), np.log(f[mask]), 1)[0])


def main():
    out = {"framework": {
        "background": "LambdaCDM conformal image E^2 = Om(1+z)^3 + OL",
        "Om": OM, "OL": OL, "sigma8_0": S8_PLANCK, "beta_A": BETA_A,
        "caveat": ("standard modified-gravity projection of G_eff(z) "
                   "onto the conformal-image background; the full TEP "
                   "growth equation (static g, landscape X_env) is "
                   "AUD-6/T6.3 — this is the comparison channel against "
                   "the incumbent closure, not the final operator")}}

    zout = np.array([0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0])
    # LambdaCDM control (G_eff = 1): calibrates the reference-data
    # tension (the known S8 offset dominates the chi2 budget)
    for sector in ("lcdm", "baseline", "two_branch", "exp_interp"):
        sol = growth(sector)
        a = 1.0 / (1.0 + zout)
        D = sol.sol(np.log(a))[0]
        Dp = sol.sol(np.log(a))[1]
        f = Dp / D                       # d ln D / d ln a
        D0 = sol.sol([0.0])[0][0]
        s8 = S8_PLANCK * D / D0
        fs8 = f * s8
        S8 = s8 * np.sqrt(OM / 0.3)
        zg = np.linspace(0.0, 2.0, 200)
        Dg = sol.sol(np.log(1.0 / (1.0 + zg)))[0]
        fg = sol.sol(np.log(1.0 / (1.0 + zg)))[1] / Dg
        out[sector] = {
            "z": zout.tolist(),
            "D_over_D0": (D / D0).tolist(),
            "f": f.tolist(),
            "G_eff_over_G": [float(G_eff(z, sector)) for z in zout],
            "f_sigma8": fs8.tolist(),
            "sigma8": s8.tolist(),
            "S8": S8.tolist(),
            "growth_index_gamma": gamma_fit(zg, fg),
        }

    # --- reference data (representative, tracer-mixed; see note) ---
    out["reference_data"] = {
        "note": ("tracer-mixed f sigma8 compilations: BOSS DR12 "
                 "(Alam+17), eBOSS DR16, DES Y3 S8, KiDS-1000 S8"),
        "fs8_points": [
            {"z": 0.38, "fs8": 0.497, "err": 0.045, "src": "BOSS DR12"},
            {"z": 0.51, "fs8": 0.458, "err": 0.038, "src": "BOSS DR12"},
            {"z": 0.61, "fs8": 0.436, "err": 0.034, "src": "BOSS DR12"},
            {"z": 0.70, "fs8": 0.473, "err": 0.041, "src": "eBOSS DR16"},
            {"z": 1.48, "fs8": 0.462, "err": 0.045, "src": "eBOSS quasar"},
        ],
        "S8_points": [
            {"S8": 0.776, "err": 0.017, "src": "DES Y3 3x2pt"},
            {"S8": 0.759, "err": 0.024, "src": "KiDS-1000"},
        ],
    }

    # chi2 against the reference points per sector
    rd = out["reference_data"]
    for sector in ("lcdm", "baseline", "two_branch", "exp_interp"):
        zz = out[sector]["z"]
        fs8 = np.interp([p["z"] for p in rd["fs8_points"]], zz,
                        out[sector]["f_sigma8"])
        chi2_fs8 = float(np.sum(
            [(fi - p["fs8"]) ** 2 / p["err"] ** 2
             for fi, p in zip(fs8, rd["fs8_points"])]))
        chi2_s8 = float(np.sum(
            [(out[sector]["S8"][0] - p["S8"]) ** 2 / p["err"] ** 2
             for p in rd["S8_points"]]))
        out[sector]["chi2"] = {"fs8": chi2_fs8, "S8": chi2_s8,
                               "total": chi2_fs8 + chi2_s8}

    # --- AUD-6 / AUD-3 remainder: in-well X_env regime map ---
    # xi = (grad phi)^2/(2 (H0/c)^2), grad phi ~ a_phi/c^2 (beta_A = -1).
    # Extended regimes sit BELOW the ambient temporal xi = 0.5; compact
    # objects sit orders above it.  The growth vertex inside a well is
    # governed by the nonlinear pair map x_env = (1-y)/y (step_30),
    # not by a naive xi_spatial substitution: the sign question
    # X = X_t - X_s is resolved by the positive-operator response.
    H0_SI, C_L, G_SI, M_SUN = 2.27e-18, 3.0e8, 6.674e-11, 1.989e30
    G_T = 3.4e-10
    PC = 3.086e16

    def xi_spatial(a):
        return float((a / C_L**2 / (H0_SI / C_L))**2 / 2.0)

    def y_of_g(g):
        q = (g / G_T) ** 2
        disc = (0.5 / q) ** 2 + (1.0 / (3.0 * q)) ** 3
        s = np.sqrt(disc)
        return float(np.cbrt(0.5 / q + s) + np.cbrt(0.5 / q - s))

    regimes = {
        "cosmo_ambient_z0": {"kind": "temporal", "xi": 0.5},
        "disk_landscape_shear_u0.18": {"kind": "spatial",
                                     "xi": xi_spatial(0.18 * G_T)},
        "mw_halo_edge_30kpc": {"kind": "spatial",
                               "xi": xi_spatial((220e3) ** 2 / (30e3 * PC))},
        "solar_surface": {"kind": "spatial",
                          "xi": xi_spatial(G_SI * M_SUN / (6.96e8) ** 2)},
    }
    vertices = {
        "disk_pair_ambient": {"g_over_gt": 0.18},
        "halo_edge": {"g_over_gt":
                      float((220e3) ** 2 / (30e3 * PC) / G_T)},
        "solar_surface": {"g_over_gt":
                          float(G_SI * M_SUN / (6.96e8) ** 2 / G_T)},
    }
    for k, vv in vertices.items():
        y = y_of_g(vv["g_over_gt"] * G_T)
        vv["y"] = y
        vv["x_env_vertex"] = (1.0 - y) / y

    out["in_well_X_env_regime_map"] = {
        "spatial_xi_by_regime": regimes,
        "vertex_x_env": vertices,
        "reading": (
            "extended regions (disk landscape, halo bulk) carry spatial "
            "xi ~ 1e-3, below the ambient temporal xi = 0.5; compact "
            "objects carry xi ~ 1e20-1e23.  For the growth vertex the "
            "operative in-well variable is the positive-operator "
            "response x = (1-y)/y: y -> 0 inside wells gives x_env -> "
            "large (S_Sigma -> 0, G_eff -> 1), so collapsed regions "
            "grow Newtonianly as required.  The residual f sigma8 "
            "excess is therefore carried by the mode-weighted ambient "
            "fraction; the remaining quantitative step is convolving "
            "x_env over the perturbation support's collapsed fraction."),
    }

    out["verdict"] = {
        "lcdm_reference": ("LambdaCDM control: f sigma8(0.61) ~ 0.477, "
                           "S8 = 0.834, chi2 = 26.0 — the S8 tension "
                           "floor; two-branch sits at this level "
                           "(27.7), baseline does not (44.2)"),
        "note": ("two-branch G_eff(0) = 1.16 vs baseline 2.0 — the "
                 "small-X branch suppresses the late-time boost the "
                 "incumbent sector carried; compare chi2 columns."),
    }

    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, "..", "..", "results",
                        "step_62_growth_fsigma8.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=2)
    for sector in ("lcdm", "baseline", "two_branch", "exp_interp"):
        s = out[sector]
        print(f"{sector:12s} S8={s['S8'][0]:.3f} "
              f"f s8(0.61)~{np.interp(0.61, s['z'], s['f_sigma8']):.3f} "
              f"chi2={s['chi2']['total']:.1f} gamma={s['growth_index_gamma']:.3f}")
    print("wrote", os.path.abspath(path))


if __name__ == "__main__":
    main()
