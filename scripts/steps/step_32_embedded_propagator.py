#!/usr/bin/env python3
"""Embedded propagator: source on the Galactic nonlinear background.

Solves the quasi-static scalar equation

    div J(grad phi) = rho,     J_i = a_i (1 + |a|^2 / g_t^2)

for a point source embedded in a uniform background gradient
u_0 = g_0/g_t = sqrt(X_GAL), X_GAL = 0.52 (step_19 convention).
Units: g_t = 1, r_* = 1, so the Newtonian gradient is a_N = 1/r^2
and the isolated-source profile reduces to y_profile(s) of step_19.

Finite-volume discretization on an axisymmetric cylindrical (s, z)
grid: face fluxes J(a_face) with a_face the face gradient of
phi = u_0 z + psi; residual and Jacobian share the same stencil,
so Newton iteration converges. The perturbation psi is extracted
as y_emb(r) = |grad psi| / a_N along the axes parallel and
perpendicular to the Galactic shear for the Paper 13 estimator
check (TEP-WB step_017).
"""
import json
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import brentq

X_GAL = 0.52
U0_REF = np.sqrt(X_GAL)

# --- grid: s in [0, SMAX], z in [-ZMAX, ZMAX] ---
NS, NZ = 110, 220
SMAX, ZMAX = 10.0, 10.0
ds = SMAX / NS
dz = 2.0 * ZMAX / NZ
s = (np.arange(NS) + 0.5) * ds          # cell centres, s>0
z = (np.arange(NZ) + 0.5) * dz - ZMAX
SS, ZZ = np.meshgrid(s, z, indexing='ij')
RR = np.sqrt(SS ** 2 + ZZ ** 2)

# --- localized source: normalized Gaussian blob, total S = 4 pi ---
SIG_SRC = 0.15
SRC_TOTAL = 4.0 * np.pi                 # gives a_N = 1/r^2 unscreened
rho = (SRC_TOTAL / (2.0 * np.pi) ** 1.5 / SIG_SRC ** 3
       * np.exp(-RR ** 2 / (2.0 * SIG_SRC ** 2)))

# s-face radial coordinate: face between cell i-1 and i sits at s=i*ds
S_FACE = np.arange(NS + 1) * ds          # length NS+1
S_FACE[0] = 0.0


def face_fluxes(psi, u0):
    """Face gradients and fluxes of phi = u0 z + psi.

    s-face (i+1/2, j), i=0..NS-1 between cells i,i+1:
        a_s = (psi[i+1]-psi[i])/ds
        a_z = mean of cell d_z on both cells + u0
    z-face (i, j+1/2), j=0..NZ-1:
        a_z = (psi[:,j+1]-psi[:,j])/dz + u0
        a_s = mean of cell d_s on both cells
    Returns Js_sface (NS-1,NZ), Js_zface (NS,NZ-1) and face Z blocks.
    """
    d_ds = np.gradient(psi, ds, axis=0)
    d_dz = np.gradient(psi, dz, axis=1)

    # s-faces
    a_s_f = (psi[1:, :] - psi[:-1, :]) / ds                    # (NS-1,NZ)
    a_z_f = 0.5 * (d_dz[:-1, :] + d_dz[1:, :]) + u0
    am2 = a_s_f ** 2 + a_z_f ** 2
    Js_f = a_s_f * (1.0 + am2)
    Zss_f = 1.0 + am2 + 2.0 * a_s_f ** 2
    Zsz_f = 2.0 * a_s_f * a_z_f

    # z-faces
    a_z_g = (psi[:, 1:] - psi[:, :-1]) / dz + u0               # (NS,NZ-1)
    a_s_g = 0.5 * (d_ds[:, :-1] + d_ds[:, 1:])
    am2g = a_s_g ** 2 + a_z_g ** 2
    Jz_g = a_z_g * (1.0 + am2g)
    Zzz_g = 1.0 + am2g + 2.0 * a_z_g ** 2
    Zsz_g = 2.0 * a_s_g * a_z_g
    return (Js_f, Zss_f, Zsz_f), (Jz_g, Zzz_g, Zsz_g)


def residual(psi, u0):
    (Js_f, _, _), (Jz_g, _, _) = face_fluxes(psi, u0)
    R = np.zeros_like(psi)
    # s-divergence: (1/s) d_s(s J_s); faces at s=i*ds and (i+1)*ds
    sJp = S_FACE[1:-1][:, None] * Js_f          # flux at face i+1/2
    # cell i: outflow face i+1/2 (between i,i+1) -> Js_f[i-? ] careful:
    # Js_f[k] is the face between cell k and k+1, at s=(k+1)ds.
    # cell i has out-face Js_f[i] (s=(i+1)ds) and in-face Js_f[i-1].
    div_s = np.zeros_like(psi)
    s_out = S_FACE[1:NS]                        # face s for Js_f[k], k=0..NS-2 -> s=(k+1)ds
    div_s[:, :] = 0.0
    # accumulate
    div_s[:-1, :] += sJp / (SS[:-1, :] * ds)
    div_s[1:, :] -= sJp / (SS[1:, :] * ds)
    # s=0 axis limit handled by cell-0 out-flux only (in-flux = 0 at s=0):
    # for i=0: div_s = s_face(1) Js_f[0] / (s0 ds); JS_f[0] ~ a_s ~ 0 anyway
    # z-divergence
    div_z = np.zeros_like(psi)
    div_z[:, :-1] += Jz_g / dz
    div_z[:, 1:] -= Jz_g / dz
    return div_s + div_z - rho


def build_jacobian(psi, u0):
    (Js_f, Zss_f, Zsz_f), (Jz_g, Zzz_g, Zsz_g) = face_fluxes(psi, u0)
    N = NS * NZ
    rows, cols, vals = [], [], []

    def add(r, c, v):
        rows.append(np.asarray(r).ravel())
        cols.append(np.asarray(c).ravel())
        vals.append(np.asarray(v).ravel())

    JJ = np.arange(NZ)

    # ---- s-faces: face k between cells k,k+1, k=0..NS-2 ----
    for k in range(NS - 1):
        i, inb = k, k + 1
        sf = S_FACE[k + 1]
        # cell i: outflux +sf*Js/ds/s_i ; cell i+1: influx -sf*Js/ds/s_{i+1}
        for ci, sgn in ((i, +1.0), (inb, -1.0)):
            c_idx = ci * NZ + JJ
            denom = np.maximum(s[ci], 1e-12) * ds
            # d(Js)/d(psi[inb]-psi[i]) = Zss_f/ds
            coef = sgn * sf * Zss_f[k, :] / (denom * ds)
            add(c_idx, inb * NZ + JJ, coef)
            add(c_idx, i * NZ + JJ, -coef)
            # d(Js)/d d_z: Js depends on a_z on the face, which is
            # 0.5*(d_z[i]+d_z[inb]) with d_z centred -> j+/-1 taps
            coefz = sgn * sf * Zsz_f[k, :] / (denom * 4 * dz)
            for dj, w in ((1, 1.0), (-1, -1.0)):
                jj = np.clip(JJ + dj, 0, NZ - 1)
                mask = np.ones(NZ, bool)
                mask[-1 if dj > 0 else 0] = False
                add(c_idx[mask], (inb * NZ + jj[mask]),
                    coefz[mask] * w)
                add(c_idx[mask], (i * NZ + jj[mask]),
                    coefz[mask] * w)

    # ---- z-faces: face j between cells j,j+1, j=0..NZ-2 ----
    II = np.arange(NS)
    for j in range(NZ - 1):
        jnb = j + 1
        for cj, sgn in ((j, +1.0), (jnb, -1.0)):
            c_idx = II * NZ + cj
            coef = sgn * Zzz_g[:, j] / (dz * dz)
            add(c_idx, II * NZ + jnb, coef)
            add(c_idx, II * NZ + j, -coef)
            coefs = sgn * Zsz_g[:, j] / (4 * ds * dz)
            for di, w in ((1, 1.0), (-1, -1.0)):
                ii = np.clip(II + di, 0, NS - 1)
                mask = np.ones(NS, bool)
                mask[-1 if di > 0 else 0] = False
                add(c_idx[mask], (ii[mask] * NZ + jnb),
                    coefs[mask] * w)
                add(c_idx[mask], (ii[mask] * NZ + j),
                    coefs[mask] * w)

    A = sp.csr_matrix(
        (np.concatenate(vals),
         (np.concatenate(rows), np.concatenate(cols))), shape=(N, N))
    return A


def dirichlet_project(A):
    A = A.tolil()
    bd = ([(NS - 1) * NZ + j for j in range(NZ)]
          + [i * NZ + j for i in range(NS) for j in (0, NZ - 1)])
    for c in bd:
        A.rows[c] = [c]
        A.data[c] = [1.0]
    return A.tocsr()


def solve_embedded(u0):
    psi = np.zeros_like(SS)
    for it in range(100):
        R = residual(psi, u0)
        R = R.copy()
        # Dirichlet rows: correction must drive psi -> 0 on the boundary
        R[NS - 1, :] = psi[NS - 1, :]
        R[:, 0] = psi[:, 0]
        R[:, NZ - 1] = psi[:, NZ - 1]
        rn = np.linalg.norm(R[2:-2, 2:-2])
        if it % 10 == 0:
            print(f"  u0={u0:.3f} it={it} |R|={rn:.3e}", flush=True)
        if rn < 1e-8:
            break
        A = dirichlet_project(build_jacobian(psi, u0))
        dpsi = spla.spsolve(A, -R.reshape(-1)).reshape(NS, NZ)
        alpha = 1.0
        accepted = False
        for _ in range(40):
            cand = psi + alpha * dpsi
            rc = residual(cand, u0)
            if np.linalg.norm(rc[2:-2, 2:-2]) < rn:
                psi = cand
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            psi = psi + dpsi * 1e-6
    print(f"  u0={u0:.3f} done it={it} |R|="
          f"{np.linalg.norm(residual(psi, u0)[2:-2, 2:-2]):.3e}")
    return psi


def extract_profiles(psi):
    dpsi_ds = np.gradient(psi, ds, axis=0)
    dpsi_dz = np.gradient(psi, dz, axis=1)
    a_par = np.abs(dpsi_dz[0, :])
    j0 = np.argmin(np.abs(z))
    a_perp = np.abs(dpsi_ds[:, j0])
    return (np.abs(z), a_par), (s.copy(), a_perp)


def y_iso(x):
    if x >= 50.0:
        return 1.0
    return brentq(lambda yy: yy * (1.0 + yy * yy / x ** 4) - 1.0,
                  1e-30, 1.0, xtol=1e-14)


def main():
    out = {"X_GAL": X_GAL, "u0": U0_REF,
           "grid": {"NS": NS, "NZ": NZ, "SMAX": SMAX, "ZMAX": ZMAX}}

    print("--- isolated reference (u0=0) ---")
    psi0 = solve_embedded(0.0)
    (_, ap), (sp_, aq) = extract_profiles(psi0)
    sel = (sp_ > 0.3) & (sp_ < 4)
    y_num = aq[sel] * sp_[sel] ** 2
    y_ana = np.array([y_iso(x) for x in sp_[sel]])
    out["isolated_check"] = {
        "max_rel_dev_vs_yprofile": float(
            np.max(np.abs(y_num - y_ana) / np.maximum(y_ana, 1e-12))),
        "note": "solver vs step_19 y_profile in the u0=0 limit"}
    print(out["isolated_check"], flush=True)

    print("--- embedded (u0 = sqrt(0.52)) ---", flush=True)
    psiE = solve_embedded(U0_REF)
    (zp, ap), (sp_, aq) = extract_profiles(psiE)
    y_par = ap * zp ** 2
    y_perp = aq * sp_ ** 2

    Zpar = 1.0 + 3.0 * U0_REF ** 2
    Zperp = 1.0 + U0_REF ** 2
    mus = np.linspace(-1, 1, 4001)

    def y_lin(mu):
        r2 = mu ** 2 / Zpar + (1 - mu ** 2) / Zperp
        return (np.sqrt(mu ** 2 / Zpar ** 2 + (1 - mu ** 2) / Zperp ** 2)
                / (np.sqrt(Zpar) * Zperp * r2 ** 1.5))
    y_lin_mu = np.array([y_lin(m) for m in mus])
    ymean = float(np.mean(y_lin_mu))
    out["linear_propagator"] = {
        "Z_parallel": float(Zpar), "Z_perp": float(Zperp),
        "y_far_parallel": float(1.0 / Zperp),
        "y_far_perp": float(1.0 / np.sqrt(Zpar * Zperp)),
        "y_far_orient_mean": ymean,
        "implied_alpha_sat_pure_propagator":
            float(np.sqrt(1.0 + 2.0 * ymean) - 1.0),
        "implied_alpha_sat_propagator_x_vertex":
            float(np.sqrt(1.0 + 2.0 * ymean / (1 + X_GAL) ** 2) - 1.0),
        "observed_alpha_sat": 0.366}

    selr = (sp_ > 0.05) & (sp_ < 8)
    selz = (zp > 0.05) & (zp < 8)
    out["profiles"] = {
        "r_over_rstar": sp_[selr].tolist(),
        "y_emb_perp": y_perp[selr].tolist(),
        "z_over_rstar": zp[selz].tolist(),
        "y_emb_par": y_par[selz].tolist(),
        "y_iso_same_grid": [y_iso(x) for x in sp_[selr].tolist()]}

    # --- orientation-mean embedded profile y_emb(r) = <|grad psi|> r^2
    from scipy.ndimage import map_coordinates

    def orient_mean_profile(psi):
        g_s = np.gradient(psi, ds, axis=0)
        g_z = np.gradient(psi, dz, axis=1)
        gmag = np.sqrt(g_s ** 2 + g_z ** 2)
        radii = np.linspace(0.2, 6.0, 120)
        mus_p = np.linspace(-0.97, 0.97, 41)
        y_mu = np.zeros((len(mus_p), len(radii)))
        for k, mu in enumerate(mus_p):
            ss = radii * np.sqrt(1 - mu ** 2)
            zz = radii * mu
            ii = ss / ds - 0.5
            jj = (zz + ZMAX) / dz - 0.5
            y_mu[k] = map_coordinates(gmag, [ii, jj], order=1,
                                      mode='nearest') * radii ** 2
        return radii, y_mu.mean(axis=0)

    radii, y_emb_mean = orient_mean_profile(psiE)
    _, y_iso_mean = orient_mean_profile(psi0)
    yfar = float(y_emb_mean[-6:].mean())
    ih = int(np.argmin(np.abs(y_emb_mean - 0.5 * yfar)))
    yfar0 = float(y_iso_mean[-6:].mean())
    ih0 = int(np.argmin(np.abs(y_iso_mean - 0.5 * yfar0)))
    out["orientation_mean_profile"] = {
        "r_over_rstar": radii.tolist(),
        "y_emb": y_emb_mean.tolist(),
        "y_iso": y_iso_mean.tolist(),
        "y_far_embedded": yfar,
        "half_excess_r_over_rstar_embedded": float(radii[ih]),
        "half_excess_r_over_rstar_isolated": float(radii[ih0]),
        "scale_ratio_embedded_vs_isolated": float(radii[ih] / radii[ih0])}

    # --- direct solar-circle ambient (data-derived, not plateau-calibrated):
    # X_GAL = 0.52 was reverse-calibrated from the WB plateau under the
    # no-propagator decomposition (step_19 comment). The independent
    # dynamical estimate uses the measured centripetal acceleration
    # a_sun = v_c^2/R_0; under TEP the observed acceleration already
    # includes the scalar field's contribution, so u_sun = a_sun/g_t is
    # the field gradient at the solar circle.
    V_C = 220.0e3            # m/s, solar-circle circular speed (IAU-ish)
    R_0 = 8.1 * 3.085677581e19  # m, Galactocentric radius
    a_sun = V_C ** 2 / R_0    # ~1.94e-10 m/s^2
    u_sun = a_sun / 3.4e-10   # g_t = 3.4e-10 m/s^2 (constants.py)
    X_sun = u_sun ** 2
    out["direct_solar_circle_ambient"] = {
        "inputs": {"v_c_km_s": 220.0, "R_0_kpc": 8.1,
                   "g_t_m_s2": 3.4e-10},
        "a_sun_m_s2": float(a_sun), "u_sun": float(u_sun),
        "X_sun": float(X_sun),
        "S_sigma_direct": float(1.0 / (1.0 + X_sun)),
        "S_sigma_plateau_calibrated": float(1.0 / (1.0 + X_GAL)),
        "note": ("direct estimate S_Σ~0.75 (the 0.74 Paper 13 cites) vs "
                 "plateau-calibrated S_Σ~0.66 (X_GAL=0.52); the two "
                 "differ because X_GAL was tuned to reproduce the plateau "
                 "WITHOUT the propagator suppression this step derives")}

    # --- environmental sweep: same solve at weaker/stronger ambient ---
    out["environmental_sweep"] = {}
    for u0_sw in (0.3, float(u_sun), 1.1):
        print(f"--- sweep u0={u0_sw} ---", flush=True)
        psi_sw = solve_embedded(u0_sw)
        _, ym = orient_mean_profile(psi_sw)
        yfar_sw = float(ym[-6:].mean())
        ihs = int(np.argmin(np.abs(ym - 0.5 * yfar_sw)))
        out["environmental_sweep"][str(u0_sw)] = {
            "y_emb": ym.tolist(),
            "y_far": yfar_sw,
            "half_excess_r_over_rstar": float(radii[ihs])}
    out["environmental_sweep"]["0.721"] = {
        "y_emb": y_emb_mean.tolist(), "y_far": yfar,
        "half_excess_r_over_rstar": float(radii[ih])}
    dest = Path(__file__).resolve().parents[2] / "results" / \
        "step_32_embedded_propagator.json"
    dest.write_text(json.dumps(out, indent=2))
    print(json.dumps(out["linear_propagator"], indent=2))
    print("wrote", dest)


if __name__ == "__main__":
    main()
