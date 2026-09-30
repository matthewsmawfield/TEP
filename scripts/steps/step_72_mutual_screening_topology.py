#!/usr/bin/env python3
"""Dimensionless two-centre screening diagnostic, not an Earth–Moon prediction.

Reuse the step_56 finite-volume solver for unequal Gaussian sources and compare
its baseline and two-branch kinetic responses. The kinetic flux J=P_X grad(phi)
is not identified with either the conformal matter force or the photon
synchronization connection. No physical body-size, solar boundary condition,
clock-link model, or disformal coefficient is supplied here.
"""
import json
from pathlib import Path

import numpy as np

import step_56_embedded_two_branch as model

SEPARATION = 3.0
CHARGE_RATIO = 1.0 / 81.0
AMBIENT = 0.15


def source_density():
    norm = model.SRC_TOTAL / ((2.0 * np.pi) ** 1.5 * model.SIG_SRC ** 3)
    large = np.exp(-(model.SS ** 2 + (model.ZZ + SEPARATION / 2) ** 2)
                   / (2 * model.SIG_SRC ** 2))
    small = np.exp(-(model.SS ** 2 + (model.ZZ - SEPARATION / 2) ** 2)
                   / (2 * model.SIG_SRC ** 2))
    return norm * (large + CHARGE_RATIO * small)


def diagnostics(psi, sector, rho):
    a_s = np.gradient(psi, model.ds, axis=0)
    a_z = np.gradient(psi, model.dz, axis=1) + AMBIENT
    q = np.hypot(a_s, a_z)
    f = model.f_flux(q, sector)
    df_ds = np.gradient(f, model.ds, axis=0)
    df_dz = np.gradient(f, model.dz, axis=1)
    curl_j = df_ds * a_z - df_dz * a_s
    curl_j_direct = (np.gradient(f * a_z, model.ds, axis=0)
                     - np.gradient(f * a_s, model.dz, axis=1))
    curl_phi = (np.gradient(a_z, model.ds, axis=0)
                - np.gradient(a_s, model.dz, axis=1))
    bridge = (np.abs(model.z) < SEPARATION / 2 - 3 * model.SIG_SRC)
    bridge_idx = np.flatnonzero(bridge)
    j = bridge_idx[np.argmin(q[0, bridge])]
    mask = ((model.SS > 0.4) & (model.SS < 1.5)
            & (np.abs(model.ZZ) < 0.8)
            & (rho < 0.01 * rho.max()))
    residual = model.residual(psi, AMBIENT, sector)
    return {
        "residual_l2": float(np.linalg.norm(residual[2:-2, 2:-2])),
        "bridge_axis_min_gradient_over_gt": float(q[0, j]),
        "bridge_axis_min_z_over_rstar": float(model.z[j]),
        "bridge_axis_kinetic_coefficient": float(f[0, j]),
        "ambient_kinetic_coefficient": float(model.f_flux(np.array(AMBIENT), sector)),
        "max_abs_curl_dphi_source_free": float(np.max(np.abs(curl_phi[mask]))),
        "max_abs_curl_kinetic_flux_source_free": float(np.max(np.abs(curl_j[mask]))),
        "max_curl_identity_difference_source_free": float(
            np.max(np.abs((curl_j_direct - curl_j)[mask]))),
    }


def main():
    model.STARVE = False
    model.KPHI_C = 0.0
    model.KPHI_PERT = False
    model.PSI_FROZEN = None
    rho = source_density()
    model.rho = rho
    rows = {}
    for sector in ("baseline", "two_branch"):
        psi = model.solve_embedded(AMBIENT, sector)
        rows[sector] = diagnostics(psi, sector, rho)
        if rows[sector]["residual_l2"] > 1e-6:
            raise RuntimeError(f"{sector} finite-volume solve did not converge")
    result = {
        "scope": "dimensionless, unequal two-centre diagnostic; not an Earth-Moon or flyby prediction",
        "charge_ratio": CHARGE_RATIO,
        "separation_over_rstar": SEPARATION,
        "ambient_gradient_over_gt": AMBIENT,
        "source_width_over_rstar": model.SIG_SRC,
        "grid": [model.NS, model.NZ],
        "kinetic_sectors": rows,
        "interpretation": (
            "A curl of the kinetic flux is not a curl of the conformal force "
            "or of the disformal synchronization connection. Low local P_X "
            "does not establish an unscreened integrated force; the physical "
            "Earth-Moon-Sun boundary value problem, source charges, clock "
            "links, and LLR fit are not solved here. The finite-difference "
            "curl comparison is grid-limited and is qualitative only."
        ),
    }
    path = Path(__file__).resolve().parents[2] / "results" / "step_72_mutual_screening_topology.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
