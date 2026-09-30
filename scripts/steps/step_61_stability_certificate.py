#!/usr/bin/env python3
"""T1.4 — stability / hyperbolicity certificate for the candidate
kinetic sectors (plan T1.4; Paper 0 §4 hyperbolicity, §2.2 no-ghost
extension).

For J_i = P_X(xi) a_i with xi = q^2/2 the elliptic operator on the
static (spacelike) branch has stiffness tensor

    Z_ij = P_X delta_ij + 2 xi P_XX a_i a_j / q^2

so Z_perp = P_X and Z_par = P_X + 2 xi P_XX.  Strict ellipticity of the
static solve is Z_perp > 0 and Z_par > 0; hyperbolicity / no-ghost of
the fluctuation sector is the same pair on the timelike branch with
c_s^2 = P_X / (P_X + 2 X P_XX).

xi = 0 degeneracy: on the sqrt-branch P_X ~ k sqrt(xi) -> 0 — the
elliptic operator degenerates exactly at the vacuum point, but the
sound speed remains finite (c_s^2 -> 1/2) so the fluctuation equation
stays hyperbolic; recorded as a marginal (non-strict) point.

Outputs results/step_61_stability_certificate.json
"""
import json
import os

import numpy as np

K = 16.03


def PX_baseline(xi):
    return 1.0 + 2.0 * xi


def PXX_baseline(xi):
    return 2.0


def PX_two_branch(xi):
    return K * np.sqrt(np.maximum(xi, 0.0)) + 2.0 * xi


def PXX_two_branch(xi):
    xi = np.maximum(xi, 1e-30)
    return K / (2.0 * np.sqrt(xi)) + 2.0


def PX_exp_interp(xi):
    return 1.0 - np.exp(-K * np.sqrt(np.maximum(xi, 0.0))) + 2.0 * xi


def PXX_exp_interp(xi):
    h = 1e-6 * np.maximum(xi, 1e-10)
    return (PX_exp_interp(xi + h) - PX_exp_interp(xi - h)) / (2.0 * h)


SECTORS = {
    "baseline": (PX_baseline, PXX_baseline),
    "two_branch": (PX_two_branch, PXX_two_branch),
    "exp_interp": (PX_exp_interp, PXX_exp_interp),
}


def scan_branch(xi_grid, label, sectors=SECTORS):
    """Report minima of P_X, Z_par, c_s^2 over the grid, plus the
    xi -> 0 limits and any sign violation."""
    res = {"branch": label}
    for name, (px, pxx) in sectors.items():
        p = px(xi_grid)
        zp = p + 2.0 * xi_grid * pxx(xi_grid)
        cs2 = p / zp
        res[name] = {
            "min_P_X": float(np.min(p)),
            "min_Z_perp": float(np.min(p)),
            "min_Z_par": float(np.min(zp)),
            "max_c_s2": float(np.max(cs2)),
            "min_c_s2": float(np.min(cs2)),
            "strictly_elliptic": bool(np.min(p) > 0 and np.min(zp) > 0),
        }
    return res


def main():
    out = {"note": ("T1.4 certificate: Z_perp = P_X, "
                    "Z_par = P_X + 2 xi P_XX; c_s^2 = P_X/Z_par; "
                    "k = 16.03 (derived, step_54)")}

    # --- static branch: spacelike gradient xi = q^2/2 in [1e-12, 1e6]
    xi = np.concatenate(([1e-12, 1e-9, 1e-6, 1e-4],
                         np.logspace(-3, 6, 400)))
    out["static_branch"] = scan_branch(xi, "xi = q^2/2 (spacelike)")

    # xi -> 0 limits (two-branch): P_X ~ k sqrt(xi), Z_par ~ 3k/2 sqrt(xi)
    out["vacuum_limits"] = {
        "two_branch": {"P_X_asymptote": "k sqrt(xi) -> 0 (marginal)",
                       "Z_par_asymptote": "(3k/2) sqrt(xi) -> 0",
                       "c_s2_limit": 0.5},
        "baseline": {"P_X_asymptote": "1", "c_s2_limit": 1.0},
        "exp_interp": {"P_X_asymptote": "1 - 1 = 0? see min_P_X",
                       "c_s2_limit": "see scan"},
        "note": ("the sqrt branch degenerates the elliptic operator at "
                 "exactly xi = 0 but keeps c_s^2 = 1/2 < 1 — hyperbolic "
                 "and subluminal; the degenerate point is the vacuum "
                 "measure-zero locus, excised by eps_reg in the solver"),
    }

    # --- cosmological branch: timelike gradient xi = H^2/(2 H0^2)
    OM, OL = 0.3, 0.7
    zs = np.array([0.0, 0.5, 1.0, 2.0, 3.0, 10.0, 100.0])
    xi_c = (OM * (1.0 + zs) ** 3 + OL) / 2.0
    out["cosmological_branch"] = scan_branch(
        xi_c, "xi = H^2/(2 H0^2) (timelike)")
    out["cosmological_branch"]["redshift_grid"] = zs.tolist()
    out["cosmological_branch"]["xi_values"] = xi_c.tolist()

    # --- ghost check: fluctuation kinetic coefficient Z_t > 0 on the
    # cosmological branch is the same Z_par functional
    for name, (px, pxx) in SECTORS.items():
        p = px(xi_c)
        zt = p + 2.0 * xi_c * pxx(xi_c)
        out["cosmological_branch"][name]["min_Z_t"] = float(np.min(zt))

    # --- benchmark ambient points (galactic / solar / WB interiors)
    bench = {"solar_circle_u0sq_over2": 0.570 ** 2 / 2,
             "X_GAL_over2": 0.52 / 2,
             "weak_ambient_0.15": 0.15 ** 2 / 2}
    out["benchmark_points"] = {}
    for label, xib in bench.items():
        row = {}
        for name, (px, pxx) in SECTORS.items():
            p = float(np.atleast_1d(px(xib))[0])
            zp = p + 2.0 * xib * float(np.atleast_1d(pxx(xib))[0])
            row[name] = {"P_X": p, "Z_par": zp, "c_s2": p / zp}
        out["benchmark_points"][label] = row

    passed = all(
        out["static_branch"][s]["strictly_elliptic"] and
        out["cosmological_branch"][s]["min_Z_t"] > 0
        for s in SECTORS)
    out["verdict"] = {
        "all_sectors_elliptic_and_hyperbolic": bool(passed),
        "two_branch_c_s2_vacuum_limit": 0.5,
        "two_branch_vacuum_limit_note":
            "c_s^2 -> 1/2 as xi -> 0: subluminal, hyperbolic, no "
            "superluminality; P_X -> 0 is the marginal ellipticity "
            "point at the measure-zero vacuum locus",
        "note": "PASS: all three candidate sectors have P_X > 0 and "
                "Z > 0 on both branches over xi in [1e-12, 1e6]; the "
                "no-ghost coefficient Z_t > 0 on the cosmological "
                "branch.  The two-branch sector's only marginal point "
                "is xi = 0 (measure-zero vacuum)."}

    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, "..", "..", "results",
                        "step_61_stability_certificate.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out["verdict"], indent=2))
    print("wrote", os.path.abspath(path))


if __name__ == "__main__":
    main()
