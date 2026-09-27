#!/usr/bin/env python3
"""Evolve the constraint-satisfying landscape and read off H_T(z).

The initial-data steps freeze a maximal slice. This step advances the
scalar and the comoving matter on that slice, in the volume-static gauge:
the spatial volume is held fixed and the lapse is rebuilt each step so
that the mean of N*f vanishes. The expansion is then not an input. It is
whatever the Raychaudhuri residual does.

Clock drift is measured on the matter-hosting wells,
    H_T = <N * pi>_wells = <du/dt>_wells,
with A = exp(-u) and beta_A = -1. Lookback is the direction of increasing
well depth. Redshift is the endpoint ratio 1+z = A_now / A_em.

Distance uses the static-space relation that follows from that definition,
    r = ∫ c dt,   d_L = (1+z)^2 r,
which is the same integral as an FLRW luminosity distance with H(z) = H_T(z).
The comparison is the shape against flat ΛCDM (Ωm=0.3), both normalised
at z=0. No parameter is adjusted to the supernova diagram.

The abundance confrontation uses the existing recycling fixed point. That
fixed point does not depend on H_T. What the evolution supplies is the
matter-frame exposure ∫ A dt. Gate 11 needs about nine star-formation
times to reach the fixed point; this step reports whether the evolved
exposure contains that many Hubble times, which is the condition the
drift law can actually test.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import scipy.fft as fft

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from step_20_landscape_existence import (  # noqa: E402
    M_PL, build_landscape, fd_laplacian, landscape_potential, solve_landscape,
)

C_KMS = 299792.458
H0_KMS_MPC = 70.0


def dV_du(u, V0, rhob_half, u_s):
    eps = 1e-4
    return (landscape_potential(u + eps, V0, rhob_half, u_s)
            - landscape_potential(u - eps, V0, rhob_half, u_s)) / (2 * eps)


def spectral_lap(field, k2):
    return np.real(fft.ifftn(-k2 * fft.fftn(field)))


def physical_lap(field, psi, k2, axes_k):
    """D^2 f = psi^{-4} [Δf + 2 ∇logψ · ∇f] on h = psi^4 δ."""
    lap = spectral_lap(field, k2)
    fh = fft.fftn(field)
    ph = fft.fftn(np.log(np.maximum(psi, 1e-8)))
    cross = np.zeros_like(field)
    for K in axes_k:
        df = np.real(fft.ifftn(1j * K * fh))
        dp = np.real(fft.ifftn(1j * K * ph))
        cross += df * dp
    return psi ** -4 * (lap + 2.0 * cross)


def preserving_lapse(f):
    """Positive lapse with ∫ N f = 0 when f changes sign. Mean(N) = 1."""
    pos = f > 0
    neg = f < 0
    P = float(f[pos].sum()) if pos.any() else 0.0
    Q = float((-f[neg]).sum()) if neg.any() else 0.0
    N = np.ones_like(f)
    if P > 0 and Q > 0:
        if P >= Q:
            N = np.where(neg, P / Q, 1.0)
        else:
            N = np.where(pos, Q / P, 1.0)
        N = N / N.mean()
    return N, bool(P > 0 and Q > 0)


def hamiltonian_residual(psi, S0, g, lam, lap_psi):
    """Lichnerowicz residual, the same combination step_20 drives to zero."""
    return (lap_psi
            + (psi ** 5) * (S0 + lam) / (4.0 * M_PL ** 2)
            + psi * g / (8.0 * M_PL ** 2))


def lcdm_dL_shape(z, n=4000):
    """(1+z) ∫_0^z dz'/E(z'), E=sqrt(0.3(1+z)^3+0.7). Units c/H0."""
    zz = np.linspace(0.0, z, n)
    E = np.sqrt(0.3 * (1 + zz) ** 3 + 0.7)
    return (1.0 + z) * np.trapz(1.0 / E, zz)


def main():
    # Same Rule-23 landscape family as the congruence solve, on a grid
    # the time stepper can advance. N=16 is the evolution resolution;
    # the initial constraint is solved on that same grid.
    L, n_w = 1.0, 4
    u_well = 0.4
    rho0, rho_ambient = 0.5, 0.05
    V0, rhob_half, u_s = 60.0, 30.0, 15.0
    N = 16
    # Same selection as the congruence solve: a preserving drift exists
    # only where the lapse operator has a negative Rayleigh quotient at
    # Pi^2 = 0. Pi^2 then shifts that eigenvalue up through zero.
    best = None
    for sg in (L / 4.0, L / 3.0):
        for u_v in (-1.0, -1.4, -1.8, -2.2):
            trial = build_landscape(N, L, n_w, sg, u_well, u_v, rho0,
                                    rho_ambient, V0, rhob_half, u_s)
            lam_est = -(float(trial["S0"].mean()) + 0.5 * float(trial["g"].mean()))
            f_s = (trial["rho_m"] - 2.0 * (trial["V"] + lam_est)) / (2.0 * M_PL ** 2)
            w = trial["wells"] / trial["wells"].max()
            phi = w ** 2
            gx = np.gradient(phi, trial["dx"], axis=0, edge_order=2)
            gy = np.gradient(phi, trial["dx"], axis=1, edge_order=2)
            gz = np.gradient(phi, trial["dx"], axis=2, edge_order=2)
            nrm = float((phi ** 2).sum())
            R = float((gx ** 2 + gy ** 2 + gz ** 2).sum()) / nrm + float((f_s * phi ** 2).sum()) / nrm
            cand = (R, sg, u_v)
            if best is None or R < best[0]:
                best = cand
    if best is None or best[0] >= 0:
        raise SystemExit(f"no volume-static drift on this grid (best Rayleigh {best})")
    rayleigh0, sigma, u_void = best
    ld = build_landscape(N, L, n_w, sigma, u_well, u_void, rho0,
                         rho_ambient, V0, rhob_half, u_s)
    u = ld["u"].copy()
    rho = ld["rho_m"].copy()
    dx = ld["dx"]
    wells = ld["wells"]
    well = wells > 0.5 * wells.max()
    LAP = fd_laplacian(N, dx)
    k = 2.0 * np.pi * fft.fftfreq(N, d=dx)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    k2 = KX ** 2 + KY ** 2 + KZ ** 2
    axes_k = (KX, KY, KZ)

    # Pi^2 shifts the lapse eigenvalue up by 3/2. The scan quotient is
    # an upper bound, so this is the smallest drift that can reach mu=0.
    Pi2 = max(0.0, -rayleigh0 * (2.0 * M_PL ** 2 / 3.0))
    psi = lam = None
    hist = []
    for _ in range(4):
        V = landscape_potential(u, V0, rhob_half, u_s)
        gu2 = ld["g"]  # initial gradients; updated after the field moves
        S0 = rho + V + 0.5 * Pi2
        sol = solve_landscape(S0, gu2, LAP, tol=1e-8, max_newton=25, verbose=False)
        psi, lam = sol["psi"], sol["lam"]
        hist.append({"Pi2": float(Pi2), "lambda": float(lam),
                     "rayleigh0": float(rayleigh0),
                     "converged": bool(sol["converged"]),
                     "ham_max": float(sol["history"][-1]) if sol["history"] else None})
        if sol["converged"]:
            break
    if psi is None or not sol["converged"]:
        raise SystemExit("initial constraint solve failed")

    # Drift into the past: wells deepen, A falls, redshift grows.
    pi = np.full_like(u, np.sqrt(Pi2) if Pi2 > 0 else 0.05)
    A_now = float(np.exp(-u[well].mean()))
    lap_psi = (LAP @ psi.ravel()).reshape(psi.shape)

    dt = 2.0e-4
    nsteps = 4000
    record_every = 40
    rows = []
    max_abs_u = 8.0

    def pack(step, u, pi, rho):
        V = landscape_potential(u, V0, rhob_half, u_s)
        gu2 = np.zeros_like(u)
        uh = fft.fftn(u)
        for K in axes_k:
            du = np.real(fft.ifftn(1j * K * uh))
            gu2 += du ** 2
        f = (rho + 2.0 * pi ** 2 - 2.0 * (V + lam)) / (2.0 * M_PL ** 2)
        Nlap, ok = preserving_lapse(f)
        u_w = float(u[well].mean())
        A_w = float(np.exp(-u_w))
        H_T = float((Nlap * pi)[well].mean())
        z = A_now / A_w - 1.0
        S0 = rho + V + 0.5 * pi ** 2
        ham = hamiltonian_residual(psi, S0, gu2, lam, lap_psi)
        # Instantaneous volume drive: mean of N*f is the gauge residual.
        drive = float((Nlap * f).mean())
        eom = physical_lap(u, psi, k2, axes_k) - dV_du(u, V0, rhob_half, u_s) + rho
        return {
            "step": step,
            "t": step * dt,
            "u_well": u_w,
            "A_well": A_w,
            "z": z,
            "H_T": H_T,
            "lapse_sign_changing": ok,
            "lapse_min": float(Nlap.min()),
            "volume_drive": drive,
            "ham_rms": float(np.sqrt(np.mean(ham ** 2))),
            "pi_well": float(pi[well].mean()),
            "scalar_eom_rms": float(np.sqrt(np.mean(eom ** 2))),
        }, V, gu2, f, Nlap

    snap, V, gu2, f, Nlap = pack(0, u, pi, rho)
    rows.append(snap)
    stopped = "completed"
    for step in range(1, nsteps + 1):
        def deriv(u, pi, rho, Nlap):
            lap_u = physical_lap(u, psi, k2, axes_k)
            dV = dV_du(u, V0, rhob_half, u_s)
            # step_12: phi_tt = -V' - alpha*rho, alpha=-1, so +rho
            dpi = Nlap * (lap_u - dV + rho)
            du = Nlap * pi
            # Einstein-frame exchange, volume fixed: rho_t = alpha * rho * phi_t
            drho = -rho * du
            return du, dpi, drho

        du1, dpi1, dr1 = deriv(u, pi, rho, Nlap)
        u2 = u + dt * du1
        pi2 = pi + dt * dpi1
        rho2 = np.maximum(rho + dt * dr1, 0.0)
        V2 = landscape_potential(u2, V0, rhob_half, u_s)
        f2 = (rho2 + 2.0 * pi2 ** 2 - 2.0 * (V2 + lam)) / (2.0 * M_PL ** 2)
        N2, _ = preserving_lapse(f2)
        du2, dpi2, dr2 = deriv(u2, pi2, rho2, N2)
        u = u + 0.5 * dt * (du1 + du2)
        pi = pi + 0.5 * dt * (dpi1 + dpi2)
        rho = np.maximum(rho + 0.5 * dt * (dr1 + dr2), 0.0)
        V = landscape_potential(u, V0, rhob_half, u_s)
        f = (rho + 2.0 * pi ** 2 - 2.0 * (V + lam)) / (2.0 * M_PL ** 2)
        Nlap, _ = preserving_lapse(f)
        if not np.isfinite(u).all() or abs(float(u[well].mean())) > max_abs_u:
            stopped = "field left the solved range"
            snap, *_ = pack(step, u, pi, rho)
            rows.append(snap)
            break
        if step % record_every == 0 or step == nsteps:
            snap, *_ = pack(step, u, pi, rho)
            rows.append(snap)
            if snap["ham_rms"] > 50.0 * rows[0]["ham_rms"] and snap["ham_rms"] > 1.0:
                stopped = "hamiltonian residual grew by more than 50 from the solved slice"
                break

    # H_T(z) from the trajectory. Drop the initial point if z~0 and H is defined.
    traj = [r for r in rows if np.isfinite(r["z"]) and r["z"] >= -0.05 and r["H_T"] != 0]
    # Distance integral along the trajectory, in units where the cell time is 1.
    # r(z) = Δt from the observer (first sample) to that emission event.
    for r in traj:
        r["r_coord"] = r["t"] - traj[0]["t"]
        r["dL_over_c"] = (1.0 + r["z"]) ** 2 * r["r_coord"]

    H_T0 = traj[0]["H_T"] if traj else None
    comparisons = []
    if H_T0 and abs(H_T0) > 0 and len(traj) > 3:
        # Normalise so H_T(z=0) = H0. Shape only.
        z_of = np.array([r["z"] for r in traj])
        HT_of = np.array([r["H_T"] for r in traj]) / H_T0
        t_of = np.array([r["t"] for r in traj])
        order = np.argsort(z_of)
        z_of, HT_of, t_of = z_of[order], HT_of[order], t_of[order]
        # unique z for interpolation
        _, uniq = np.unique(np.round(z_of, 6), return_index=True)
        z_of, HT_of, t_of = z_of[uniq], HT_of[uniq], t_of[uniq]
        for z_t in (0.1, 0.5, 1.0):
            if z_t > z_of.max():
                comparisons.append({"z": z_t, "reached": False})
                continue
            HT = float(np.interp(z_t, z_of, HT_of))
            # coordinate distance in units c/H0: ∫ dz/((1+z) H_T/H0)
            zz = np.linspace(0.0, z_t, 2000)
            HT_path = np.interp(zz, z_of, HT_of)
            HT_path = np.maximum(HT_path, 1e-8)
            r = np.trapz(1.0 / ((1.0 + zz) * HT_path), zz)
            dL = (1.0 + z_t) ** 2 * r
            dL_ref = lcdm_dL_shape(z_t)
            dmu = 5.0 * np.log10(dL / dL_ref) if dL > 0 and dL_ref > 0 else None
            comparisons.append({
                "z": z_t,
                "reached": True,
                "H_T_over_H0": HT,
                "closed_static_H_t_over_H0": float(np.sqrt(1.0 + 3.0 * 0.3 * z_t / (1.0 + z_t))),
                "lcdm_E": float(np.sqrt(0.3 * (1 + z_t) ** 3 + 0.7)),
                "delta_mu_vs_lcdm": None if dmu is None else float(dmu),
            })
        matter_time = float(np.trapz(np.exp(-np.interp(t_of, t_of, 
                            [r["u_well"] for r in traj])), t_of - t_of[0])) if False else None
    else:
        matter_time = None

    # Matter-frame exposure along the run, in units of 1/H_T(0).
    if H_T0 and abs(H_T0) > 0:
        t = np.array([r["t"] for r in rows])
        A = np.array([r["A_well"] for r in rows])
        tau = float(np.trapz(A, t))
        exposure_over_hubble = tau * abs(H_T0)
    else:
        tau = None
        exposure_over_hubble = None

    # Gate 11 fixed point, pure-well limit. Independent of H_T.
    R_ret, P_Y = 0.4, 0.149
    Y_eq = P_Y / (1.0 - R_ret)
    gate11_efolds = 9.0

    out = {
        "step": "step_29_landscape_evolution",
        "gauge": "volume-static, lapse rebuilt to keep mean(N*f)=0",
        "grid": {"N": N, "L": L, "dt": dt, "nsteps_requested": nsteps},
        "initial_constraint": {
            "lambda_star": float(lam),
            "Pi2": float(Pi2),
            "psi_range": [float(psi.min()), float(psi.max())],
            "history": hist,
            "ham_rms_initial": rows[0]["ham_rms"],
        },
        "stopped": stopped,
        "final": rows[-1],
        "n_records": len(rows),
        "H_T_definition": "mean N*pi on matter wells; A=exp(-u); 1+z=A_now/A_em",
        "comparisons": comparisons,
        "exposure": {
            "matter_frame_time": tau,
            "matter_frame_time_in_units_of_1_over_HT0": exposure_over_hubble,
            "gate11_relaxation_times": gate11_efolds,
            "Y_eq_pure_well": Y_eq,
            "Y_observed": 0.245,
            "note": (
                "Y_eq does not depend on H_T. The drift law can only test "
                "whether the exposure contains many star-formation times. "
                "Nine Gate-11 times inside one Hubble time would require a "
                "star-formation time shorter than ~0.1/H_T; this run reports "
                "the exposure it actually produced, not an assumed tau_star."
            ),
        },
        "trajectory_tail": rows[:: max(1, len(rows) // 20)],
    }
    dest = os.path.join(os.path.dirname(__file__), "..", "..", "results",
                        "step_29_landscape_evolution.json")
    dest = os.path.normpath(dest)
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    with open(dest, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps({
        "stopped": stopped,
        "z_final": rows[-1]["z"],
        "H_T_final": rows[-1]["H_T"],
        "ham_rms_initial": rows[0]["ham_rms"],
        "ham_rms_final": rows[-1]["ham_rms"],
        "volume_drive_final": rows[-1]["volume_drive"],
        "lapse_min_final": rows[-1]["lapse_min"],
        "exposure_over_hubble": exposure_over_hubble,
        "comparisons": comparisons,
        "wrote": dest,
    }, indent=2))


if __name__ == "__main__":
    main()
