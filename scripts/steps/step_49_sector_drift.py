#!/usr/bin/env python3
"""Sector-asymmetric drift of the R8 slice under the maximal-slicing lapse.

A homogeneous Pi^2 was a one-number device to put a zero in the lapse
operator. It is not the roll. The source is the matter trace, so the
normal momentum is supported only where rho_* exceeds the ambient floor.
Voids have pi = 0. Their field values therefore do not advance.

Time evolution on the frozen maximal slice is
    du/dt = N * pi
with N the lapse that satisfies int N f = 0 (mean preservation of K = 0).
Redshift and path length are read off that trajectory:
    A = exp(-u) on matter hosts,
    1 + z = A_obs / A_em,
    dL = <N> d tau
along the past-ward march. No H(z) template is inserted.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from step_20_landscape_existence import (  # noqa: E402
    M_PL, build_landscape, fd_laplacian, solve_landscape,
)

L = 1.0
N = 16
N_W = 4
U_WELL = 0.4
RHO0 = 0.5
RHO_AMB = 0.05
V0 = 60.0
RHOB_HALF = 30.0
U_S = 15.0
# Matter hosts sit above this fraction of the well weight.
HOST = 0.5


def rayleigh(ld):
    lam = -(float(ld["S0"].mean()) + 0.5 * float(ld["g"].mean()))
    f = (ld["rho_m"] - 2.0 * (ld["V"] + lam)) / (2.0 * M_PL ** 2)
    w = ld["wells"] / ld["wells"].max()
    phi = w ** 2
    dx = ld["dx"]
    g2 = sum(np.gradient(phi, dx, axis=a, edge_order=2) ** 2 for a in range(3))
    nrm = float((phi ** 2).sum())
    R = float((g2.sum()) + (f * phi ** 2).sum()) / nrm
    return R


def preserving_lapse(f):
    pos, neg = f > 0, f < 0
    P = float(f[pos].sum()) if pos.any() else 0.0
    Q = float((-f[neg]).sum()) if neg.any() else 0.0
    Nlap = np.ones_like(f)
    ok = P > 0 and Q > 0
    if ok:
        if P >= Q:
            Nlap = np.where(neg, P / Q, 1.0)
        else:
            Nlap = np.where(pos, Q / P, 1.0)
        Nlap = Nlap / Nlap.mean()
    return Nlap, ok, float((Nlap * f).mean())


def main():
    # Narrow wells so a void sector exists. The Rayleigh sign is checked
    # after the constraint solve, not used to reject the geometry.
    sigma, u_void = L / 8.0, -1.4
    ld = build_landscape(N, L, N_W, sigma, U_WELL, u_void, RHO0,
                         RHO_AMB, V0, RHOB_HALF, U_S)
    R = rayleigh(ld)
    # R8 slice: the constraint as solved, with no homogeneous kinetic term.
    rho = ld["rho_m"]
    sol = solve_landscape(ld["S0"], ld["g"], fd_laplacian(N, ld["dx"]),
                          tol=1e-8, max_newton=20, verbose=False)
    if not sol["converged"]:
        raise SystemExit("R8 constraint did not converge")
    lam = sol["lam"]
    # Momentum lives on the matter excess only. A homogeneous Pi^2 is not
    # added. The amplitude is kept small enough that 1/2 pi^2 is a
    # perturbation on |lambda*|, so the slice we evolve is still the
    # solved one.
    w = ld["wells"]
    host = w >= np.quantile(w, 0.80)
    void = w <= np.quantile(w, 0.05)
    excess = np.maximum(rho - rho[void].mean(), 0.0)
    pi_cap = np.sqrt(0.02 * abs(lam))
    shape = excess / np.maximum(excess[host].mean(), 1e-30)
    pi = np.where(host, pi_cap * shape / shape[host].max(), 0.0)
    f = (rho - 2.0 * (ld["V"] + lam)) / (2.0 * M_PL ** 2)
    Nlap, ok, drive = preserving_lapse(f)

    # Instantaneous rates. Voids have pi = 0, so they do not move.
    rate = Nlap * pi
    H_matter = float(rate[host].mean())
    H_void = float(rate[void].mean())
    u = ld["u"].copy()
    u_void0 = u[void].copy()
    void_band0 = float(np.max(np.abs(u_void0 - u_void0.mean())))

    rows = []
    tau = 0.0
    Lpath = 0.0
    u_obs = float(u[host].mean())
    dt = 0.002
    for step in range(0, 4001):
        if step % 40 == 0:
            u_m = float(u[host].mean())
            du = u_m - u_obs
            z = float(np.exp(du) - 1.0)  # A_obs/A_em - 1, A=exp(-u)
            dvoid = float(np.max(np.abs(u[void] - u_void0)))
            rows.append({
                "tau": tau,
                "L": Lpath,
                "z": z,
                "A_clock_emitter": float(np.exp(-u_m)),
                "A_clock_observer": float(np.exp(-u_obs)),
                "H_drift": H_matter,
                "delta_u_matter": du,
                "void_secular": dvoid,
                "void_in_band": bool(dvoid <= 1.0e-2),
            })
        u = u + dt * rate
        tau += dt
        Lpath += dt * float(Nlap.mean())
        if rows and rows[-1]["z"] > 1.5:
            break

    out = {
        "step": "step_49_sector_drift",
        "why_not_homogeneous": (
            "Homogeneous Pi^2 was only the lapse zero-mode closure. "
            "The source is rho_* A_phi, which is zero in the voids, "
            "so the normal momentum is zero there."
        ),
        "landscape": {
            "rayleigh": R, "sigma": sigma, "u_void": u_void,
            "lambda_star": float(lam),
            "lapse_preserves_mean_K": ok,
            "mean_Nf": drive,
            "constraint_converged": True,
            "ham_max": float(sol["history"][-1]),
            "pi_cap": float(pi_cap),
            "pi_energy_over_lambda": float((0.5 * pi ** 2).max() / abs(lam)),
        },
        "rates": {
            "H_drift_matter": H_matter,
            "H_drift_void": H_void,
            "void_initial_spatial_band": void_band0,
            "band_limit": 1.0e-2,
        },
        "z_of_L": rows,
    }
    dest = os.path.normpath(os.path.join(
        os.path.dirname(__file__), "..", "..", "results", "step_49_sector_drift.json"))
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    with open(dest, "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps({
        "H_matter": H_matter,
        "H_void": H_void,
        "void_band0": void_band0,
        "lapse_ok": ok,
        "zL": [(r["L"], r["z"], r["void_secular"], r["A_clock_emitter"]) for r in rows[::2]],
        "wrote": dest,
    }, indent=2))


if __name__ == "__main__":
    main()
