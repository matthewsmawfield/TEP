#!/usr/bin/env python3
"""On-shell landscape evolution: wells as solutions, not prescriptions.

step_29 prescribes a Gaussian-lattice field u(x) and solves only the
Lichnerowicz constraint for psi. The scalar equation on that slice,
    D^2 u = V'(u) - rho_m        (the on-shell static condition),
is violated by ~1e2, so the first timestep applies a spatially varying
acceleration of O(1e2): the field relaxes violently, overshoots, and the
well drift reverses within dt ~ 0.01. That is a relaxation transient of
constraint-inconsistent initial data, not a property of the drift law.

This step re-runs the experiment with the initial scalar field *solved*
from its elliptic equation on the same matter landscape:

  1. solve_scalar_onshell: Newton solve of D^2 u - V'(u) + rho + lam_c = 0
     with mean(u) = 0 (ambient-zero convention). A constant Lagrange
     multiplier absorbs the periodic-domain solvability mismatch, playing
     the same role as lam* in the Lichnerowicz solve.
  2. solve_landscape: Lichnerowicz constraint for psi on the solved u.
  3. Outer fixed point: re-solve u with the conformal Laplacian
     D^2 on h = psi^4 delta, then re-solve psi, until the coupled
     (scalar, metric) pair is self-consistent.

Two structural facts then emerge:

  * On the master-potential branch V'(u) ~ (rhob/2) e^{-u} > 0, the
    scalar rolls to DECREASING u: clocks accelerate forward in time,
    hence were slower in the past -- the secular redshift trend is the
    natural downhill roll. step_29 instead kicked pi uphill, so the
    uniform slope (the lam_c solvability offset, equal in magnitude to
    the constraint's lam*) reversed the drift at t ~ pi0/|lam_c| ~ 0.02.
  * The on-shell solution has shallow wells (Du ~ rho/V'' scale), so
    deep prescribed wells are themselves off-shell.

Three runs isolate the defects:

  A) on-shell initial data at rest, identical step_29 machinery
     (frozen psi, two-zone preserving lapse);
  B) on-shell initial data at rest, constrained evolution: the
     Lichnerowicz solve is re-run on the evolved sources and the lapse
     is the step_21 maximal-slice eigen lapse (positive principal
     eigenvector of -D^2 + q, q = rho - |Du|^2/2 - 2(V + lam*));
  C) as B, with a downhill sector kick (pi < 0 on hosts): hosts descend
     faster so matter sat deeper in wells in the past (Rule 23);
  D) nested landscape (Rule 22): a broad envelope well carrying a
     deeper ambient baseline under half the lattice, solved the same
     way. Per-sector drift rates test whether the local clock/drift
     depends on the ambient baseline it nests inside (hierarchical
     time). A screening-regime diagnostic reports |grad u|^2 on the
     solved slice: the canonical Laplacian EOM is the unscreened limit
     of P_,X = 1 + 2|X|/Lambda^4, so the correction size is recorded
     rather than assumed.

Lookback is then the reverse read of the forward trajectory: the final
slice is today, 1 + z = A_today/A(t) = exp(u_end - u(t)), path distance
r = int <N> dt, d_L = (1+z)^2 r, compared in shape to flat LCDM.

Diagnostics recorded at each checkpoint: scalar-EOM residual,
Lichnerowicz residual, volume drive <N f>, u on matter hosts, A, z, H_T.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import scipy.fft as fft
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from step_20_landscape_existence import (  # noqa: E402
    M_PL, build_landscape, fd_laplacian, landscape_potential,
)

C_KMS = 299792.458


# ----------------------------------------------------------------------
# Scalar EOM pieces (convention shared with step_12 / step_29):
#   u_tt = N * (D^2 u - V'(u) + rho),   u_t = N * pi.
# Static shell: D^2 u = V'(u) - rho.
# ----------------------------------------------------------------------
def dV_du(u, V0, rhob_half, u_s):
    eps = 1e-4
    return (landscape_potential(u + eps, V0, rhob_half, u_s)
            - landscape_potential(u - eps, V0, rhob_half, u_s)) / (2 * eps)


def d2V_du2(u, V0, rhob_half, u_s):
    eps = 1e-3
    return (landscape_potential(u + eps, V0, rhob_half, u_s)
            - 2.0 * landscape_potential(u, V0, rhob_half, u_s)
            + landscape_potential(u - eps, V0, rhob_half, u_s)) / (eps * eps)


def spectral_lap(field, k2):
    return np.real(fft.ifftn(-k2 * fft.fftn(field)))


def conformal_lap(field, psi, k2, axes_k):
    """D^2 f = psi^{-4} [Delta f + 2 grad log psi . grad f], h = psi^4 delta."""
    lap = spectral_lap(field, k2)
    fh = fft.fftn(field)
    ph = fft.fftn(np.log(np.maximum(psi, 1e-8)))
    cross = np.zeros_like(field)
    for K in axes_k:
        df = np.real(fft.ifftn(1j * K * fh))
        dp = np.real(fft.ifftn(1j * K * ph))
        cross += df * dp
    return psi ** -4 * (lap + 2.0 * cross)


def conformal_grad2(field, psi, k2, axes_k):
    """|D f|^2 = psi^{-4} |grad f|^2."""
    fh = fft.fftn(field)
    g2 = np.zeros_like(field)
    for K in axes_k:
        df = np.real(fft.ifftn(1j * K * fh))
        g2 += df ** 2
    return psi ** -4 * g2


# ----------------------------------------------------------------------
# Bordered Newton solve of  F(u, lam_c) = LAP u - V'(u) + rho + lam_c = 0
# with gauge mean(u) = u_bar. Mirrors the bordered Schur construction of
# step_20.solve_landscape: L d1 = -F, L d2 = dF/dlam_c = 1, and the
# constant-mode mean is pinned by choosing dlam_c so that mean(du) is the
# required gauge correction.
# ----------------------------------------------------------------------
def solve_scalar_onshell(rho, V0, rhob_half, u_s, LAP,
                         psi=None, k2=None, axes_k=None,
                         u0=None, u_bar=0.0, tol=1e-9, max_newton=40):
    shape = rho.shape
    u = np.zeros(shape) if u0 is None else u0.copy()
    lam_c = 0.0
    converged = False
    hist = []

    use_psi = psi is not None and k2 is not None

    def F(u, lam_c):
        if use_psi:
            lap_u = conformal_lap(u, psi, k2, axes_k)
        else:
            lap_u = (LAP @ u.ravel()).reshape(shape)
        return lap_u - dV_du(u, V0, rhob_half, u_s) + rho + lam_c

    for it in range(max_newton):
        r = F(u, lam_c)
        res = float(np.abs(r).max())
        hist.append(res)
        if res < tol:
            converged = True
            break

        # Linearised operator L = D^2 - V''(u). For the flat solve LAP is
        # the sparse 7-point Laplacian; for the conformal pass the
        # psi-weighting is mild, so the same constant-coefficient operator
        # (scaled by mean psi^{-4}) is used as the Newton matrix.
        c_op = -d2V_du2(u, V0, rhob_half, u_s)
        if use_psi:
            c_op = c_op * np.mean(psi ** -4)
        L_mat = LAP + sp.diags(c_op.ravel(), format="csc")
        try:
            lu = spla.splu(L_mat.tocsc())
        except RuntimeError:
            return dict(u=u, lam_c=lam_c, converged=False,
                        history=hist, note="singular Newton matrix")
        d1 = lu.solve(-r.ravel()).reshape(shape)
        d2 = lu.solve(np.ones(shape).ravel()).reshape(shape)

        g_gauge = float(u.mean()) - u_bar
        m1, m2 = float(d1.mean()), float(d2.mean())
        if abs(m2) < 1e-14:
            return dict(u=u, lam_c=lam_c, converged=False,
                        history=hist, note="bordered Schur denominator ~ 0")
        dlam = (m1 + g_gauge) / m2
        du = d1 - dlam * d2

        alpha = 1.0
        u_new, lam_new = u, lam_c
        for _ in range(14):
            u_new = u + alpha * du
            lam_new = lam_c + alpha * dlam
            if np.abs(F(u_new, lam_new)).max() < res or alpha < 1e-5:
                break
            alpha *= 0.5
        u, lam_c = u_new, lam_new

    return dict(u=u, lam_c=lam_c, converged=converged, history=hist)


# ----------------------------------------------------------------------
# Warm-started Lichnerowicz solve (step_20.solve_landscape with an initial
# guess). Identical bordered-Newton structure.
# ----------------------------------------------------------------------
def solve_landscape_warm(S0, g, LAP, psi0=None, lam0=0.0,
                         tol=1e-9, max_newton=15):
    shape = S0.shape
    psi = np.ones(shape) if psi0 is None else psi0.copy()
    lam = lam0
    converged = False
    hist = []

    def F(psi, lam):
        return ((LAP @ psi.ravel()).reshape(shape)
                + (psi ** 5) * (S0 + lam) / (4.0 * M_PL ** 2)
                + psi * g / (8.0 * M_PL ** 2))

    for _ in range(max_newton):
        r = F(psi, lam)
        res = float(np.abs(r).max())
        hist.append(res)
        if res < tol:
            converged = True
            break
        c_op = (5.0 / 4.0) * psi ** 4 * (S0 + lam) / (M_PL ** 2) \
            + g / (8.0 * M_PL ** 2)
        q = (psi ** 5 / (4.0 * M_PL ** 2)).ravel()
        L_mat = LAP + sp.diags(c_op.ravel(), format="csc")
        lu = spla.splu(L_mat.tocsc())
        d1 = lu.solve(-r.ravel()).reshape(shape)
        d2 = lu.solve(q).reshape(shape)
        m1, m2 = float(d1.mean()), float(d2.mean())
        if abs(m2) < 1e-14:
            break
        dlam = m1 / m2
        dpsi = d1 - dlam * d2
        alpha = 1.0
        psi_new, lam_new = psi, lam
        for _ in range(12):
            psi_new = psi + alpha * dpsi
            lam_new = lam + alpha * dlam
            if np.abs(F(psi_new, lam_new)).max() < res or alpha < 1e-4:
                break
            alpha *= 0.5
        psi, lam = psi_new, lam_new

    return dict(psi=psi, lam=lam, converged=converged, history=hist)


# ----------------------------------------------------------------------
# Lapse laws.
# ----------------------------------------------------------------------
def preserving_lapse(f):
    """step_29 two-zone rescaling (kept for variant A comparison)."""
    pos, neg = f > 0, f < 0
    P = float(f[pos].sum()) if pos.any() else 0.0
    Q = float((-f[neg]).sum()) if neg.any() else 0.0
    Nlap = np.ones_like(f)
    if P > 0 and Q > 0:
        if P >= Q:
            Nlap = np.where(neg, P / Q, 1.0)
        else:
            Nlap = np.where(pos, Q / P, 1.0)
        Nlap = Nlap / Nlap.mean()
    return Nlap, bool(P > 0 and Q > 0)


def maximal_slice_lapse(rho, g2, V, lam, LAP):
    """step_21 construction: positive principal eigenvector of -LAP + q,
    q = rho - |Du|^2/2 - 2(V + lam*).  Normalised to mean(N) = 1.
    Returns (N, mu_star, ok).
    """
    q = (rho - 0.5 * g2 - 2.0 * (V + lam)).ravel()
    A = -LAP + sp.diags(q, format="csc")
    try:
        w, v = spla.eigsh(A, k=1, which="SA", tol=1e-8, maxiter=4000)
    except Exception:
        return np.ones_like(rho), float("nan"), False
    mu_star = float(w[0])
    N = np.abs(v[:, 0]).reshape(rho.shape)
    N = N / N.mean()
    return N, mu_star, bool(v[:, 0].min() * v[:, 0].max() >= 0)


# ----------------------------------------------------------------------
# Shared diagnostics.
# ----------------------------------------------------------------------
def snapshot(step, t, u, pi, rho, psi, lam, lap_psi, host,
             V0, rhob_half, u_s, k2, axes_k, LAP, u_obs, Nlap, drive,
             host2=None):
    V = landscape_potential(u, V0, rhob_half, u_s)
    gu2 = conformal_grad2(u, psi, k2, axes_k)
    S0 = rho + V + 0.5 * pi ** 2
    ham = ((lap_psi)
           + (psi ** 5) * (S0 + lam) / (4.0 * M_PL ** 2)
           + psi * gu2 / (8.0 * M_PL ** 2))
    eom = conformal_lap(u, psi, k2, axes_k) - dV_du(u, V0, rhob_half, u_s) + rho
    # The uniform component of the EOM residual is the homogeneous roll
    # (pi_dot ~ -lam_c), not an inhomogeneous error. Report it separately.
    eom_inhom = eom - float(eom.mean())
    u_h = float(u[host].mean())
    A_h = float(np.exp(-u_h))
    z = float(np.exp(u_h - u_obs) - 1.0)
    out = {
        "step": step, "t": float(t),
        "u_host": u_h, "A_host": A_h, "z": z,
        "H_T": float((Nlap * pi)[host].mean()),
        "u_void": float(u[~host].mean()),
        "u_min": float(u.min()), "u_max": float(u.max()),
        "volume_drive": float(drive),
        "lapse_min": float(Nlap.min()), "lapse_max": float(Nlap.max()),
        "ham_rms": float(np.sqrt(np.mean(ham ** 2))),
        "scalar_eom_rms": float(np.sqrt(np.mean(eom ** 2))),
        "scalar_eom_mean": float(eom.mean()),
        "scalar_eom_inhom_rms": float(np.sqrt(np.mean(eom_inhom ** 2))),
        "pi_host": float(pi[host].mean()),
        "pi_min": float(pi.min()),
    }
    if host2 is not None:
        out["u_host2"] = float(u[host2].mean())
        out["H_T_2"] = float((Nlap * pi)[host2].mean())
    return out, gu2


def run_variant(label, u0, pi0, rho0, psi0, lam0, host,
                V0, rhob_half, u_s, k2, axes_k, LAP,
                constrained, resolve_every, dt, nsteps, record_every,
                z_stop, host2=None):
    """Advance (u, pi, rho). `constrained` selects whether psi and the
    step_21 lapse are re-solved on the current sources."""
    u, pi, rho = u0.copy(), pi0.copy(), rho0.copy()
    psi, lam = psi0.copy(), lam0
    N = u.shape[0]
    lap_psi = (LAP @ psi.ravel()).reshape(psi.shape)
    u_obs = float(u[host].mean())
    r_path = 0.0
    V = landscape_potential(u, V0, rhob_half, u_s)
    gu2 = conformal_grad2(u, psi, k2, axes_k)
    f = (rho + 2.0 * pi ** 2 - 2.0 * (V + lam)) / (2.0 * M_PL ** 2)
    if constrained:
        Nlap, mu_star, _ = maximal_slice_lapse(rho, gu2, V, lam, LAP)
        drive = float((Nlap * f).mean())
    else:
        Nlap, _ = preserving_lapse(f)
        drive = float((Nlap * f).mean())

    rows = []
    snap, _ = snapshot(0, 0.0, u, pi, rho, psi, lam, lap_psi, host,
                       V0, rhob_half, u_s, k2, axes_k, LAP, u_obs,
                       Nlap, drive, host2=host2)
    snap["r_path"] = 0.0
    snap["mean_N"] = float(Nlap.mean())
    snap["mu_star"] = float(mu_star) if constrained else None
    rows.append(snap)
    stopped = "completed"
    t = 0.0

    def deriv(u, pi, rho, psi, Nlap):
        lap_u = conformal_lap(u, psi, k2, axes_k)
        dV = dV_du(u, V0, rhob_half, u_s)
        dpi = Nlap * (lap_u - dV + rho)
        du = Nlap * pi
        drho = -rho * du
        return du, dpi, drho

    for step in range(1, nsteps + 1):
        if constrained and step % resolve_every == 0:
            # Constrained evolution: re-solve the Lichnerowicz constraint
            # on the current sources, then rebuild the maximal-slice lapse.
            V = landscape_potential(u, V0, rhob_half, u_s)
            gu2 = conformal_grad2(u, psi, k2, axes_k)
            S0 = rho + V + 0.5 * pi ** 2
            sol = solve_landscape_warm(S0, gu2, LAP, psi0=psi, lam0=lam,
                                       tol=1e-8, max_newton=10)
            if sol["converged"]:
                psi, lam = sol["psi"], sol["lam"]
                lap_psi = (LAP @ psi.ravel()).reshape(psi.shape)
            Nlap, mu_star, _ = maximal_slice_lapse(rho, gu2, V, lam, LAP)

        du1, dpi1, dr1 = deriv(u, pi, rho, psi, Nlap)
        u2 = u + dt * du1
        pi2 = pi + dt * dpi1
        rho2 = np.maximum(rho + dt * dr1, 0.0)
        V2 = landscape_potential(u2, V0, rhob_half, u_s)
        gu22 = conformal_grad2(u2, psi, k2, axes_k)
        f2 = (rho2 + 2.0 * pi2 ** 2 - 2.0 * (V2 + lam)) / (2.0 * M_PL ** 2)
        if constrained:
            N2, _, _ = maximal_slice_lapse(rho2, gu22, V2, lam, LAP)
        else:
            N2, _ = preserving_lapse(f2)
        du2, dpi2, dr2 = deriv(u2, pi2, rho2, psi, N2)
        u = u + 0.5 * dt * (du1 + du2)
        pi = pi + 0.5 * dt * (dpi1 + dpi2)
        rho = np.maximum(rho + 0.5 * dt * (dr1 + dr2), 0.0)
        t += dt
        r_path += 0.5 * dt * (float(Nlap.mean()) + float(N2.mean()))

        if step % record_every == 0:
            V = landscape_potential(u, V0, rhob_half, u_s)
            gu2 = conformal_grad2(u, psi, k2, axes_k)
            f = (rho + 2.0 * pi ** 2 - 2.0 * (V + lam)) / (2.0 * M_PL ** 2)
            drive = float((Nlap * f).mean())
            snap, _ = snapshot(step, t, u, pi, rho, psi, lam, lap_psi, host,
                               V0, rhob_half, u_s, k2, axes_k, LAP, u_obs,
                               Nlap, drive, host2=host2)
            snap["r_path"] = float(r_path)
            snap["mean_N"] = float(Nlap.mean())
            if constrained:
                snap["mu_star"] = float(mu_star)
            rows.append(snap)
            if u[host].mean() < u_obs - 2.6:
                stopped = "roll covered Delta u = 2.6 (z ~ 12 at lookback)"
                break
        if not np.isfinite(u).all():
            stopped = "non-finite field"
            break
        if rows and rows[-1]["ham_rms"] > 50.0 * rows[0]["ham_rms"] \
                and rows[-1]["ham_rms"] > 1.0:
            stopped = "hamiltonian residual grew by more than 50"
            break

    return {"label": label, "rows": rows, "stopped": stopped,
            "final": rows[-1], "u_final": u}


def main():
    L, n_w = 1.0, 4
    N = 16
    sigma, u_void_probe = L / 8.0, -1.4
    u_well_probe = 0.4
    rho0, rho_ambient = 0.5, 0.05
    V0, rhob_half, u_s = 60.0, 30.0, 15.0

    dx = L / N
    LAP = fd_laplacian(N, dx)
    k = 2.0 * np.pi * fft.fftfreq(N, d=dx)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    k2 = KX ** 2 + KY ** 2 + KZ ** 2
    axes_k = (KX, KY, KZ)

    # Matter landscape: same well lattice as step_49 (narrow wells so a
    # void sector exists). Only rho_m and the host mask are kept; the
    # field itself is solved below, not prescribed.
    ld = build_landscape(N, L, n_w, sigma, u_well_probe, u_void_probe,
                         rho0, rho_ambient, V0, rhob_half, u_s)
    rho0_f = ld["rho_m"].copy()
    w = ld["wells"]
    host = w >= np.quantile(w, 0.80)
    void = w <= np.quantile(w, 0.05)

    print("[1] on-shell scalar solve  (D^2 u = V' - rho, mean u = 0)")
    s1 = solve_scalar_onshell(rho0_f, V0, rhob_half, u_s, LAP,
                              u_bar=0.0, tol=1e-10)
    u = s1["u"]
    lam_c = s1["lam_c"]
    print(f"    converged={s1['converged']}  iters={len(s1['history'])}"
          f"  lam_c={lam_c:+.4f}  u in [{u.min():+.3f},{u.max():+.3f}]"
          f"  u_host={u[host].mean():+.4f}  u_void={u[void].mean():+.4f}")

    print("[2] Lichnerowicz solve on the on-shell field")
    gu2 = conformal_grad2(u, np.ones_like(u), k2, axes_k)
    V = landscape_potential(u, V0, rhob_half, u_s)
    S0 = rho0_f + V
    sol = solve_landscape_warm(S0, gu2, LAP, tol=1e-10, max_newton=30)
    psi, lam = sol["psi"], sol["lam"]
    print(f"    converged={sol['converged']}  lam*={lam:+.4f}"
          f"  psi in [{psi.min():.4f},{psi.max():.4f}]")

    print("[3] coupled outer fixed point (conformal scalar re-solve)")
    for it in range(3):
        s_it = solve_scalar_onshell(rho0_f, V0, rhob_half, u_s, LAP,
                                    psi=psi, k2=k2, axes_k=axes_k,
                                    u0=u, u_bar=0.0, tol=1e-10)
        u = s_it["u"]
        lam_c = s_it["lam_c"]
        gu2 = conformal_grad2(u, psi, k2, axes_k)
        V = landscape_potential(u, V0, rhob_half, u_s)
        S0 = rho0_f + V
        sol = solve_landscape_warm(S0, gu2, LAP, psi0=psi, lam0=lam,
                                   tol=1e-10, max_newton=15)
        psi, lam = sol["psi"], sol["lam"]
        eom = conformal_lap(u, psi, k2, axes_k) \
            - dV_du(u, V0, rhob_half, u_s) + rho0_f
        lap_psi = (LAP @ psi.ravel()).reshape(psi.shape)
        ham = (lap_psi + (psi ** 5) * (S0 + lam) / (4.0 * M_PL ** 2)
               + psi * gu2 / (8.0 * M_PL ** 2))
        print(f"    iter {it}: scalar ok={s_it['converged']} "
              f"eom_rms={np.sqrt(np.mean(eom**2)):.3e}  "
              f"ham_rms={np.sqrt(np.mean(ham**2)):.3e}  "
              f"lam*={lam:+.4f}  lam_c={lam_c:+.4f}  "
              f"u in [{u.min():+.3f},{u.max():+.3f}]")

    # No kick: the on-shell slice is momentarily at rest and the potential
    # slope drives the roll. pi_dot = N (D^2 u - V' + rho) is initially
    # uniform, ~ -lam_c: the field rolls to decreasing u, i.e. clocks
    # accelerate forward in time, so lookback sees slower clocks -- the
    # secular redshift trend. step_29 instead kicked pi positive (uphill
    # in u), which the ~|lam_c| slope reverses within t ~ pi0/|lam_c|
    # ~ 0.02 -- the reported "drift reversal" is that fight against the
    # roll, not a failure of the drift.
    pi0 = np.zeros_like(u)

    # Downhill sector kick for variant C: matter hosts descend faster
    # than voids (pi < 0 on hosts), so hosts were deeper in their wells
    # in the past -- extra redshift on the matter sector, Rule 23.
    excess = np.maximum(rho0_f - rho0_f[void].mean(), 0.0)
    pi_cap = np.sqrt(0.02 * abs(lam))
    shape = excess / np.maximum(excess[host].mean(), 1e-30)
    pi0_C = np.where(host, -pi_cap * shape / shape[host].max(), 0.0)

    common = dict(V0=V0, rhob_half=rhob_half, u_s=u_s, k2=k2,
                  axes_k=axes_k, LAP=LAP, host=host)

    print("[4] variant A: on-shell, rest slice, step_29 machinery")
    A = run_variant("A_onshell_rest_frozenpsi", u, pi0, rho0_f, psi, lam,
                    dt=2.0e-4, nsteps=8000, record_every=200,
                    z_stop=2.0, constrained=False, resolve_every=0,
                    **common)
    print(f"    stopped={A['stopped']}  rows={len(A['rows'])}")
    for r in A["rows"][:6]:
        print(f"      t={r['t']:.4f} u_h={r['u_host']:+.4f}"
              f" H_T={r['H_T']:+.4f} eom_i={r['scalar_eom_inhom_rms']:.2e}"
              f" ham={r['ham_rms']:.2e}")

    print("[5] variant B: on-shell, rest slice, constrained evolution")
    B = run_variant("B_onshell_rest_constrained", u, pi0, rho0_f, psi, lam,
                    dt=2.0e-4, nsteps=8000, record_every=200,
                    z_stop=2.0, constrained=True, resolve_every=20,
                    **common)
    print(f"    stopped={B['stopped']}  rows={len(B['rows'])}")
    for r in B["rows"][:6]:
        print(f"      t={r['t']:.4f} u_h={r['u_host']:+.4f}"
              f" H_T={r['H_T']:+.4f} eom_i={r['scalar_eom_inhom_rms']:.2e}"
              f" ham={r['ham_rms']:.2e} mu*={r.get('mu_star')}")

    print("[6] variant C: on-shell, downhill sector kick, constrained")
    C = run_variant("C_onshell_downhill_constrained", u, pi0_C, rho0_f,
                    psi, lam, dt=2.0e-4, nsteps=8000, record_every=200,
                    z_stop=2.0, constrained=True, resolve_every=20,
                    **common)
    print(f"    stopped={C['stopped']}  rows={len(C['rows'])}")
    for r in C["rows"][:6]:
        print(f"      t={r['t']:.4f} u_h={r['u_host']:+.4f}"
              f" H_T={r['H_T']:+.4f} eom_i={r['scalar_eom_inhom_rms']:.2e}"
              f" ham={r['ham_rms']:.2e} mu*={r.get('mu_star')}")

    # ---------------------------------------------------------------
    # [7] Nested landscape (Rule 22): a broad envelope well carrying a
    # deeper ambient baseline, containing half the well lattice. Hosts
    # inside the envelope sit on an environment that is itself deeper in
    # the field; their clock rate is then A_envelope * A_local, the
    # hierarchical nesting. Solved the same way: rho_nested is the only
    # input, u is again a solution, not a prescription.
    # ---------------------------------------------------------------
    print("[7] nested landscape: envelope + wells (two ambient baselines)")
    x = np.arange(N) * dx
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    env = np.exp(-((X - 0.5 * L) ** 2 + (Y - 0.5 * L) ** 2
                   + (Z - 0.5 * L) ** 2) / (2.0 * (L / 3.2) ** 2))
    env_amp = 0.35 * rho0  # envelope adds a moderate ambient over-density
    rho_nested = rho0_f + env_amp * env
    # Split hosts by envelope strength: the inner half sits on a deeper
    # ambient baseline than the outer half (Rule-22 hierarchy).  The
    # same threshold splits the *ambient* cells into the parent
    # envelope's local baseline (amb_in) and the exterior ambient
    # (amb_out): two genuine hierarchy levels measured on one field.
    env_thr = np.median(env[host])
    host_deep = host & (env >= env_thr)
    host_shallow = host & ~host_deep
    amb_in = ~host & (env >= env_thr)
    amb_out = ~host & (env < env_thr)

    sn = solve_scalar_onshell(rho_nested, V0, rhob_half, u_s, LAP,
                              u_bar=0.0, tol=1e-10)
    u_n = sn["u"]
    lam_c_n = sn["lam_c"]
    gu2_n = conformal_grad2(u_n, np.ones_like(u_n), k2, axes_k)
    V_n = landscape_potential(u_n, V0, rhob_half, u_s)
    sol_n = solve_landscape_warm(rho_nested + V_n, gu2_n, LAP,
                                 tol=1e-10, max_newton=30)
    psi_n, lam_n = sol_n["psi"], sol_n["lam"]
    for _ in range(3):
        s_it = solve_scalar_onshell(rho_nested, V0, rhob_half, u_s, LAP,
                                    psi=psi_n, k2=k2, axes_k=axes_k,
                                    u0=u_n, u_bar=0.0, tol=1e-10)
        u_n = s_it["u"]
        lam_c_n = s_it["lam_c"]
        gu2_n = conformal_grad2(u_n, psi_n, k2, axes_k)
        V_n = landscape_potential(u_n, V0, rhob_half, u_s)
        sol_n = solve_landscape_warm(rho_nested + V_n, gu2_n, LAP,
                                     psi0=psi_n, lam0=lam_n,
                                     tol=1e-10, max_newton=15)
        psi_n, lam_n = sol_n["psi"], sol_n["lam"]
    eom_n = (conformal_lap(u_n, psi_n, k2, axes_k)
             - dV_du(u_n, V0, rhob_half, u_s) + rho_nested)
    eom_n_inhom = eom_n - eom_n.mean()
    print(f"    converged={sn['converged']}  lam_c={lam_c_n:+.4f}  "
          f"lam*={lam_n:+.4f}  u in [{u_n.min():+.3f},{u_n.max():+.3f}]")
    print(f"    u_host_deep={u_n[host_deep].mean():+.4f}  "
          f"u_host_shallow={u_n[host_shallow].mean():+.4f}  "
          f"(baseline offset {u_n[host_deep].mean()-u_n[host_shallow].mean():+.4f})")
    print(f"    u_amb_in={u_n[amb_in].mean():+.4f}  "
          f"u_amb_out={u_n[amb_out].mean():+.4f}  "
          f"(parent-vs-global ambient offset "
          f"{u_n[amb_in].mean()-u_n[amb_out].mean():+.4f})")
    print(f"    eom_inhom_rms={np.sqrt(np.mean(eom_n_inhom**2)):.3e}")

    # Screening-regime check: the canonical Laplacian EOM is the
    # unscreened limit of the kinetic operator P_,X = 1 + 2|X|/Lambda^4.
    # On the solved landscape |grad u|^2 is measured; the P_X correction
    # is negligible wherever the gradient sits far below the screening
    # threshold, which is what is reported here.
    gu2_flat = conformal_grad2(u_n, np.ones_like(u_n), k2, axes_k)
    screening_diag = {
        "max_grad_u2_init": float(gu2_flat.max()),
        "mean_grad_u2_init": float(gu2_flat.mean()),
        "host_grad_u2_init": float(gu2_flat[host].mean()),
        "note": ("Canonical Laplacian EOM used; the kinetic screening "
                 "operator P_,X = 1 + 2|X|/Lambda^4 with "
                 "Lambda^4 = M_Pl^2 H_0^2 is evaluated a posteriori on "
                 "the solved slice (H_0 = present-day drift rate of the "
                 "same slice).  PXm1 = 2|X|/Lambda^4 is the measured "
                 "departure from the unscreened limit; screened_fraction "
                 "is the cell fraction with PXm1 >= 0.1."),
    }

    print("[8] variant D: nested landscape, constrained evolution")
    D = run_variant("D_nested_constrained", u_n, np.zeros_like(u_n),
                    rho_nested, psi_n, lam_n,
                    dt=2.0e-4, nsteps=8000, record_every=200,
                    z_stop=2.0, constrained=True, resolve_every=20,
                    host2=host_shallow,
                    **{**common, "host": host_deep})
    # P_X on the final ("today") slice: Lambda^4 = M_Pl^2 H_0^2, with
    # the slice's own present-day drift rate as H_0 (the cosmological
    # shear floor that fixes the kinetic scale).  M_Pl = L = 1 in slice
    # units, so P_X - 1 = |grad u|^2 / H_0^2 on the same slice.
    H0_slice = abs(D["rows"][-1]["H_T"]) if D["rows"] else np.nan
    lam4 = H0_slice ** 2
    gu2_today = conformal_grad2(D["u_final"], np.ones_like(u_n),
                                k2, axes_k)
    px = gu2_today / lam4  # = 2|X|/Lambda^4 = P_,X - 1
    screening_diag.update({
        "H0_slice_drift": float(H0_slice),
        "lambda4": float(lam4),
        "max_grad_u2_today": float(gu2_today.max()),
        "PXm1_max": float(px.max()),
        "PXm1_mean": float(px.mean()),
        "PXm1_host_deep": float(px[host_deep].mean()),
        "PXm1_host_shallow": float(px[host_shallow].mean()),
        "PXm1_void": float(px[void].mean()),
        "screened_fraction": float((px >= 0.1).mean()),
        "regime": ("screened" if px.max() >= 0.1
                   else "weak-gradient (unscreened)"),
        "clock_ratio_deep_over_shallow":
            float(np.exp(-(u_n[host_deep].mean()
                         - u_n[host_shallow].mean()))),
    })
    print(f"    stopped={D['stopped']}  rows={len(D['rows'])}")
    for r in D["rows"][:6]:
        print(f"      t={r['t']:.4f} u_deep={r['u_host']:+.4f}"
              f" u_shal={r.get('u_host2'):+.4f}"
              f" H_T1={r['H_T']:+.4f} H_T2={r.get('H_T_2'):+.4f}"
              f" ham={r['ham_rms']:.2e}")

    def lookback_curve(run):
        """Re-read the forward trajectory as lookback: the final slice is
        'today', earlier slices sit at r = path distance, with
        1+z = A_today/A(t) and H_T(z) = -du/dt on hosts."""
        rows = run["rows"]
        u_end = rows[-1]["u_host"]
        out = []
        for r in rows:
            zlb = float(np.exp(r["u_host"] - u_end) - 1.0)
            out.append({"t_from_end": rows[-1]["t"] - r["t"],
                        "r_lookback": rows[-1]["r_path"] - r["r_path"],
                        "z": zlb, "H_T": abs(r["H_T"]),
                        "u_host": r["u_host"], "u_void": r["u_void"]})
        return out

    def lcdm_dL_shape(z, n=4000):
        zz = np.linspace(0.0, z, n)
        E = np.sqrt(0.3 * (1 + zz) ** 3 + 0.7)
        return (1.0 + z) * np.trapezoid(1.0 / E, zz)

    def dL_comparison(run):
        """d_L = (1+z)^2 r with r = path length. Shape comparison only:
        both curves are normalised at the smallest available z anchor,
        so the ratio tests the distance-redshift *shape*, not the
        normalisation of H_T (which is in slice units)."""
        lb = lookback_curve(run)
        anchor = [0.1, 0.5, 1.0, 2.0]
        res = {}
        norm_tep = norm_lcdm = None
        for za in anchor:
            past = [r for r in lb if r["z"] <= za]
            if len(past) < 2:
                res[f"z_{za}"] = None
                continue
            rr = past[-1]["r_lookback"]
            dL = (1 + za) ** 2 * rr
            dL_ref = lcdm_dL_shape(za)
            if norm_tep is None:
                norm_tep, norm_lcdm = dL, dL_ref
            res[f"z_{za}"] = {
                "dL": dL, "dL_lcdm": dL_ref,
                "ratio_raw": dL / dL_ref if dL_ref else None,
                "shape_ratio": (dL / norm_tep) / (dL_ref / norm_lcdm)
                if norm_tep and norm_lcdm else None,
                "r_lookback": rr,
            }
        return res

    out = {
        "step": "step_51_onshell_landscape_evolution",
        "purpose": ("Test whether the step_29 drift reversal is a "
                    "prescription artifact: identical machinery on "
                    "scalar-EOM-consistent initial data (variant A), and "
                    "with constraint re-solve plus the step_21 maximal-"
                    "slice lapse (variant B)."),
        "grid": {"N": N, "L": L, "sigma": sigma,
                 "u_bar": 0.0},
        "onshell_init": {
            "scalar_converged": bool(s1["converged"]),
            "scalar_history": s1["history"],
            "lam_c": float(lam_c),
            "u_min": float(u.min()), "u_max": float(u.max()),
            "u_host": float(u[host].mean()),
            "u_void": float(u[void].mean()),
            "lam_star": float(lam),
            "psi_min": float(psi.min()), "psi_max": float(psi.max()),
        },
        "variant_A": {**{k: v for k, v in A.items() if k != "u_final"},
                      "lookback": lookback_curve(A),
                      "dL_comparison": dL_comparison(A)},
        "variant_B": {**{k: v for k, v in B.items() if k != "u_final"},
                      "lookback": lookback_curve(B),
                      "dL_comparison": dL_comparison(B)},
        "variant_C": {**{k: v for k, v in C.items() if k != "u_final"},
                      "lookback": lookback_curve(C),
                      "dL_comparison": dL_comparison(C)},
        "nested_landscape": {
            "env_amp_over_rho0": 0.35,
            "host_deep_cells": int(host_deep.sum()),
            "host_shallow_cells": int(host_shallow.sum()),
            "lam_c": float(lam_c_n), "lam_star": float(lam_n),
            "u_min": float(u_n.min()), "u_max": float(u_n.max()),
            "u_host_deep": float(u_n[host_deep].mean()),
            "u_host_shallow": float(u_n[host_shallow].mean()),
            "u_amb_in": float(u_n[amb_in].mean()),
            "u_amb_out": float(u_n[amb_out].mean()),
            "baseline_offset": float(u_n[host_deep].mean()
                                     - u_n[host_shallow].mean()),
            "ambient_offset_in_out": float(u_n[amb_in].mean()
                                           - u_n[amb_out].mean()),
            "clock_ratios": {
                "host_deep_over_shallow":
                    float(np.exp(-(u_n[host_deep].mean()
                                 - u_n[host_shallow].mean()))),
                "host_deep_over_ambient_out":
                    float(np.exp(-(u_n[host_deep].mean()
                                 - u_n[amb_out].mean()))),
                "ambient_in_over_out":
                    float(np.exp(-(u_n[amb_in].mean()
                                 - u_n[amb_out].mean()))),
            },
            "eom_inhom_rms": float(np.sqrt(np.mean(eom_n_inhom ** 2))),
            "offset_evolution": [
                {"t": r["t"],
                 "offset": r["u_host"] - r.get("u_host2", r["u_host"]),
                 "H_T_deep": r["H_T"], "H_T_shallow": r.get("H_T_2")}
                for r in D["rows"]],
            "screening_diag": screening_diag,
        },
        "variant_D": {**{k: v for k, v in D.items() if k != "u_final"},
                      "lookback": lookback_curve(D),
                      "dL_comparison": dL_comparison(D)},
        "conventions": {
            "scalar_eom": "u_tt = N (D^2 u - V'(u) + rho); static shell D^2 u = V' - rho",
            "initial_pi": "zero: momentarily-static on-shell slice; the potential slope drives the roll downhill (u decreasing forward => clocks slower in the past)",
            "kick_C": "pi < 0 on host sector only: hosts descend faster, deeper wells in the past (Rule 23)",
            "lapse_A": "step_29 two-zone preserving_lapse (comparison only)",
            "lapse_B_C": "step_21 maximal-slice lapse: principal eigenvector of -D^2 + q, rebuilt every 20 steps; Lichnerowicz re-solved on the same cadence",
        },
    }
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "..", "..", "results",
                            "step_51_onshell_landscape_evolution.json")
    out_path = os.path.normpath(out_path)
    with open(out_path, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"  -> {out_path}")


if __name__ == "__main__":
    main()
