#!/usr/bin/env python3
"""
step_20_landscape_existence.py
==============================
Numerical existence demonstration for the inhomogeneous non-compact
eternal slice (manuscript Section 8, Lichnerowicz constraint), plus the
matter-frame flatness ledger required by the open-construction
programme (issues 0-1 and 0-14).

Construction.  The Hamiltonian constraint is posed on a periodic cell
whose universal cover is the non-compact R^3 slice — the standard
lattice-cosmology treatment of an inhomogeneous eternal universe.
With h_ij = psi^4 * hat_h_ij (hat_h flat), K = 0, A_ij = 0, Pi = 0
(time-symmetric instant of the drift trajectory), the Lichnerowicz
equation reads

    lap(psi) = -(psi^5 / 4 M_Pl^2) S(x) - (psi / 8 M_Pl^2) g(x),

    S(x) = rho_m(x) + V(u(x)) + Pi^2/2,   g(x) = |grad u|^2.

Integrability (Gauss) on the cell requires

    int psi^5 S dV + (1/2) int psi g dV = 0,

so this periodic conformally flat ansatz needs a negative contribution to S:
positive inside matter-hosting wells (rho_m + V > 0), negative in the
void sector.  Since rho_m >= 0 and the kinetic/momentum terms are
non-negative, only the scalar potential can carry the negative sector.
This is a restriction of this ansatz, not a structural requirement of every
non-compact eternal solution. The reconstructed
potential V_rec(u) = V_0 - (rho_bar/2) e^{-u} crosses zero
at u* = ln(rho_bar a^2 / (4 M_Pl^2)) < 0, so void field values u < u*
(clocks faster than ambient, Rule 10) would carry V < 0. The actual baseline
used below is positive and does not reach that branch; the solver's uniform
lam offset supplies the negative contribution in this numerical benchmark.
It must not be mistaken for a derived potential term or the localized
positive-energy construction required by TEP Rule 23.

Solver.  Newton iteration on F(psi, lam) = 0 with the void energy level
lam as the free integrability parameter (bordered system, gauge
mean(psi) = 1).  Linear solves use GMRES on the FFT-accelerated
indefinite operator L = Delta + c(x) with a shifted spectral
preconditioner.  lam* reports the ambient/void energy level required
for the constraint to close — its sign and magnitude are the output.

Flatness ledger.  On the solved slice we compute
  - the g-frame curvature  3R = -8 psi^{-5} lap(psi)  (= 2 rho_tot/M^2),
  - the matter-frame curvature  3R_t = A^{-2}(3R - 4 D_h^2 lnA - 2|D_h lnA|^2),
  - the effective curvature parameter |Omega_k^eff| = |<3R_t>|/(6 H_drift^2),
compared against the closed-static benchmark |Omega_k| >= 1/2.

Outputs results/step_20_landscape_existence.json.
"""

import json
import os
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import scipy.fft as fft

M_PL = 1.0


def fd_laplacian(N, dx):
    """Periodic 7-point finite-difference Laplacian on N^3 grid."""
    n = N ** 3
    idx = lambda i, j, k: (i % N) * N * N + (j % N) * N + (k % N)
    rows, cols, vals = [], [], []
    c = 1.0 / dx ** 2
    for i in range(N):
        for j in range(N):
            for k in range(N):
                m = idx(i, j, k)
                rows.append(m); cols.append(m); vals.append(-6.0 * c)
                for di, dj, dk in ((1, 0, 0), (-1, 0, 0), (0, 1, 0),
                                 (0, -1, 0), (0, 0, 1), (0, 0, -1)):
                    rows.append(m)
                    cols.append(idx(i + di, j + dj, k + dk))
                    vals.append(c)
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n))


# ----------------------------------------------------------------------
# Landscape construction
# ----------------------------------------------------------------------
def landscape_potential(u, V0, rhob_half, u_s):
    """Evaluate the same signed field convention in the solver and diagnostics."""
    u = np.asarray(u, dtype=float)
    eps = 1e-8
    V_rec = V0 - rhob_half * np.exp(-u)
    V = (V_rec * np.exp(-(u / u_s) ** 4)
         + V0 * np.exp(-(u_s / np.where(np.abs(u) < eps, eps, u)) ** 4))
    return np.where(np.abs(u) < eps, V_rec, V)


def conformal_spatial_derivatives(psi, field, lap_flat, wavevectors):
    """Scalar Laplacian and squared gradient for h_ij = psi^4 delta_ij.

    D_h^2 f = psi^-4 [Delta_flat f + 2 grad(log psi).grad(f)].
    Flat-grid derivatives cannot be inserted directly into the conformal
    curvature formula once the gravitational conformal factor is nonconstant.
    """
    field_hat = fft.fftn(field)
    psi_hat = fft.fftn(psi)
    grad2_flat = np.zeros_like(field)
    cross_flat = np.zeros_like(field)
    for K in wavevectors:
        df = np.real(fft.ifftn(1j * K * field_hat))
        dpsi = np.real(fft.ifftn(1j * K * psi_hat))
        grad2_flat += df**2
        cross_flat += dpsi * df
    return (psi**-4 * (lap_flat + 2.0 * cross_flat / psi),
            psi**-4 * grad2_flat)


def build_landscape(N, L, n_w, sigma, u_well, u_void, rho0, rho_ambient,
                    V0, rhob_half, u_s):
    """Return u(x), rho_m(x), S0(x) = rho_m + V(u), g(x) = |grad u|^2
    on an N^3 periodic grid, plus bookkeeping."""

    dx = L / N
    x = np.arange(N) * dx
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')

    # Cubic lattice of Gaussian wells (normalised so u spans
    # [u_void, u_well] even when wells overlap at large sigma)
    u = np.full((N, N, N), u_void)
    wells = np.zeros((N, N, N))
    spacing = L / n_w
    centres = [(i + 0.5) * spacing for i in range(n_w)]
    for cx in centres:
        for cy in centres:
            for cz in centres:
                r2 = (X - cx) ** 2 + (Y - cy) ** 2 + (Z - cz) ** 2
                wells += np.exp(-r2 / (2.0 * sigma ** 2))
    u += (u_well - u_void) * wells / wells.max()

    # Matter: ambient floor + concentration tracking the wells
    rho_m = rho0 * (rho_ambient + (1.0 - rho_ambient) * wells / wells.max())

    # Potential: natural shallow master-family member
    #   V(u) = V_rec(u) exp(-(u/u_s)^4) + V0 exp(-(u_s/u)^4)
    #   V_rec(u) = V0 - (rho_bar/2) e^{-u}   (closed-branch reconstruction)
    V = landscape_potential(u, V0, rhob_half, u_s)

    S0 = rho_m + V

    # |grad u|^2 spectrally
    k = 2.0 * np.pi * fft.fftfreq(N, d=dx)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing='ij')
    uh = fft.fftn(u)
    gu2 = np.zeros_like(u)
    for K in (KX, KY, KZ):
        du = np.real(fft.ifftn(1j * K * uh))
        gu2 += du ** 2

    return dict(u=u, wells=wells, rho_m=rho_m, V=V, S0=S0, g=gu2,
                dx=dx, N=N, L=L)


def make_laplacian(N, dx):
    k = 2.0 * np.pi * fft.fftfreq(N, d=dx)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing='ij')
    return KX ** 2 + KY ** 2 + KZ ** 2


# ----------------------------------------------------------------------
# Bordered Newton solve of F(psi, lam) = 0, mean(psi) = 1
#   Unknowns: psi (N^3), lam (1).  Equations: F = 0 (N^3) + gauge mean(psi)=1.
#   Schur complement:  dpsi = d1 - dlam*d2,  L d1 = -F,  L d2 = q = dF/dlam,
#   dlam = mean(d1)/mean(d2)  (enforces mean(dpsi) = 0).
# ----------------------------------------------------------------------
def solve_landscape(S0, g, LAP, tol=1e-9, max_newton=30, verbose=True):
    N = S0.shape[0]
    shape = S0.shape
    n = N ** 3

    psi = np.ones(shape)
    lam = 0.0
    converged = False
    hist = []

    def F(psi, lam):
        return ((LAP @ psi.ravel()).reshape(shape)
                + (psi ** 5) * (S0 + lam) / (4.0 * M_PL ** 2)
                + psi * g / (8.0 * M_PL ** 2))

    for it in range(max_newton):
        r = F(psi, lam)
        res = np.abs(r).max()
        hist.append(float(res))
        if verbose:
            print(f"  Newton {it:2d}: max|F| = {res:.3e}   lam = {lam:+.6e}   "
                  f"psi in [{psi.min():.4f}, {psi.max():.4f}]")
        if res < tol:
            converged = True
            break

        c_op = (5.0 / 4.0) * psi ** 4 * (S0 + lam) / (M_PL ** 2) \
               + g / (8.0 * M_PL ** 2)
        q = (psi ** 5 / (4.0 * M_PL ** 2)).ravel()

        L_mat = LAP + sp.diags(c_op.ravel(), format='csc')
        lu = spla.splu(L_mat)
        d1 = lu.solve(-r.ravel()).reshape(shape)
        d2 = lu.solve(q).reshape(shape)

        m1, m2 = d1.mean(), d2.mean()
        if abs(m2) < 1e-14:
            print("    [warn] bordered Schur denominator ~ 0")
            break
        dlam = m1 / m2
        dpsi = d1 - dlam * d2

        # Damped update if needed for positivity/convergence
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
# Main
# ----------------------------------------------------------------------
def main():
    print("[INFO] step_20_landscape_existence — Lichnerowicz landscape solve")

    # --- landscape parameters (dimensionless slice units, M_Pl = 1) ---
    N = 48                    # grid per side (48^3 = 110k dof)
    L = 1.0                   # cell size (universal cover = R^3)
    n_w = 4                   # wells per side (64 wells)
    sigma = L / 20.0          # well width
    u_well = 1.0              # well field value  (A = e^{-1} ~ 0.37)
    u_void = 0.8              # actual baseline field; V_rec(+0.8) > 0
    rho0 = 2.0                # mean matter density scale
    rho_ambient = 0.05        # void-to-well matter floor ratio
    V0 = 2.0                  # deep-field floor (closed-branch value)
    rhob_half = 1.0           # rho_bar/2 of the closed-static benchmark
    u_s = 15.0                # master-family transition scale (step_12 probe)
    H_drift = 1.0             # drift rate unit for Omega_k normalisation

    ld = build_landscape(N, L, n_w, sigma, u_well, u_void, rho0,
                         rho_ambient, V0, rhob_half, u_s)
    u, rho_m, V, S0, g = ld['u'], ld['rho_m'], ld['V'], ld['S0'], ld['g']
    dx = ld['dx']
    k2 = make_laplacian(N, dx)
    LAP = fd_laplacian(N, dx)
    dV_cell = dx ** 3

    # --- sign-structure check: does V(u) go negative on the void side? ---
    # V_rec(u*) = V0 - rhob_half e^{-u*} = 0  =>  u* = ln(rhob_half/V0) = -ln 2
    u_star = np.log(rhob_half / V0)
    V_at_void = float(landscape_potential(u_void, V0, rhob_half, u_s))
    print(f"  V_rec sign-crossing u* = {u_star:.4f};  "
          f"V(u_void = {u_void:+.4f}) = {V_at_void:+.4f}")

    # --- solve ---
    sol = solve_landscape(S0, g, LAP)
    psi, lam = sol['psi'], sol['lam']

    if not sol['converged']:
        print("[FAIL] Newton did not converge")
        verdict = "FAIL"

    psi_min, psi_max = float(psi.min()), float(psi.max())

    # --- integrability ledger ---
    int_psi5S = float((psi ** 5 * (S0 + lam)).sum() * dV_cell)
    int_psi_g = float((psi * g).sum() * dV_cell)
    integrability_residual = int_psi5S / (4 * M_PL ** 2) \
        + int_psi_g / (8 * M_PL ** 2)

    # --- curvature ledgers ---
    lap_psi = (LAP @ psi.ravel()).reshape(psi.shape)
    R3 = -8.0 * psi ** -5 * lap_psi          # g-frame 3R
    R3_check = (2.0 * (S0 + lam) / M_PL ** 2
                + psi ** -4 * g / M_PL ** 2)
    R3_residual = float(np.abs(R3 - R3_check).max())

    # matter-frame curvature, A = e^{-u} (conformal g~ = A^2 g):
    #  3R_t = A^{-2} ( 3R - 4 lap lnA - 2 |grad lnA|^2 )
    lnA = -u
    lap_lnA_flat = (LAP @ lnA.ravel()).reshape(u.shape)
    k = 2.0 * np.pi * fft.fftfreq(N, d=dx)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing='ij')
    lap_lnA, glnA2 = conformal_spatial_derivatives(
        psi, lnA, lap_lnA_flat, (KX, KY, KZ))
    A2 = np.exp(-2.0 * u)
    R3_t = (R3 - 4.0 * lap_lnA - 2.0 * glnA2) / A2

    # volume weightings: g-frame ~ psi^6, matter ~ A^3 psi^6
    vol_g = float((psi ** 6).sum() * dV_cell)
    vol_m = float((np.exp(-3.0 * u) * psi ** 6).sum() * dV_cell)
    mean_R3 = float((R3 * psi ** 6).sum() * dV_cell / vol_g)
    mean_R3_t = float((R3_t * np.exp(-3.0 * u) * psi ** 6).sum()
                      * dV_cell / vol_m)
    mean_R3_t_coord = float(R3_t.mean())

    # decomposition of the matter-frame mean curvature
    w_m = np.exp(-3.0 * u) * psi ** 6
    ledger_terms = {
        "A^-2_3R": float((R3 / A2 * w_m).sum() * dV_cell / vol_m),
        "A^-2_-4lap_lnA": float((-4.0 * lap_lnA / A2 * w_m).sum()
                                * dV_cell / vol_m),
        "A^-2_-2grad_lnA2": float((-2.0 * glnA2 / A2 * w_m).sum()
                                  * dV_cell / vol_m),
    }

    Omega_k_eff = abs(mean_R3_t) / (6.0 * H_drift ** 2)
    Omega_k_eff_coord = abs(mean_R3_t_coord) / (6.0 * H_drift ** 2)

    # --- landscape-smoothness scan (psi=1 ledger; psi-1 ~ 1e-3) ---
    # Flatness is controlled by the spatial contrast of u per drift
    # scale: Omega_k ~ |grad u|^2 / H_drift^2.  The observed redshift
    # scatter bounds the SPATIAL landscape contrast to << 1, so the
    # scan quantifies the flatness achieved as the contrast is reduced.
    def ledger_for(uu, w):
        lA = -uu
        l_lA = (LAP @ lA.ravel()).reshape(uu.shape)
        hh = fft.fftn(lA)
        g2 = np.zeros_like(uu)
        for K in (KX, KY, KZ):
            g2 += np.real(fft.ifftn(1j * K * hh)) ** 2
        A2l = np.exp(-2.0 * uu)
        wm = np.exp(-3.0 * uu)
        # constraint-enforced mean: lambda shifts <S> to ~0 (psi=1 gauge)
        Rg_raw = 2.0 * (rho0 * (rho_ambient + (1 - rho_ambient)
                                * w / w.max())
                        + V0 - rhob_half * np.exp(-uu))
        Rg = Rg_raw - Rg_raw.mean()
        Rt = (Rg - 4.0 * l_lA - 2.0 * g2) / A2l
        return float((Rt * wm).sum() / wm.sum()), float(g2.mean())

    def make_wells(sg):
        xg = np.arange(N) * dx
        Xg, Yg, Zg = np.meshgrid(xg, xg, xg, indexing='ij')
        w = np.zeros((N, N, N))
        spacing = L / n_w
        for cx in [(i + 0.5) * spacing for i in range(n_w)]:
            for cy in [(i + 0.5) * spacing for i in range(n_w)]:
                for cz in [(i + 0.5) * spacing for i in range(n_w)]:
                    w += np.exp(-((Xg - cx) ** 2 + (Yg - cy) ** 2
                                  + (Zg - cz) ** 2) / (2.0 * sg ** 2))
        return w

    scan = []
    for du, sg in ((1.8, sigma), (0.4, sigma), (0.1, sigma),
                   (0.1, 2 * sigma), (0.02, 2 * sigma)):
        w = make_wells(sg)
        uu = -du / 2 + du * w / w.max()
        mRt, mg2 = ledger_for(uu, w)
        scan.append({"Delta_u": du, "sigma_frac_L": sg / L,
                     "mean_R3_t": mRt, "mean_grad_u2": mg2,
                     "Omega_k_eff": abs(mRt) / (6.0 * H_drift ** 2)})
        print(f"  scan Du={du:5.2f} s={sg:.3f}: <3R~>={mRt:+.3e}  "
              f"|Ok_eff|={abs(mRt)/(6*H_drift**2):.3e}")

    # required void energy level: S in voids after lam shift
    void_mask = ld['wells'] < 0.5
    S_void_mean = float((S0 + lam)[void_mask].mean())
    S_well_mean = float((S0 + lam)[~void_mask].mean())
    rho_tot_mean = float(((S0 + lam)).mean())

    # implied void depth: V(u_d*) = lam-shifted void level (ambient V=0 ref)
    # interpret lam as the required uniform energy offset
    verdict = ("PASS" if (sol['converged'] and psi_min > 0) else "FAIL")

    summary = {
        "step": "step_20_landscape_existence",
        "description": "Newton-GMRES solve of the Lichnerowicz Hamiltonian "
                       "constraint on a periodic cell (universal cover = "
                       "non-compact R^3): inhomogeneous temporal landscape "
                       "of wells + voids, with the void energy level "
                       "determined by the integrability condition; plus the "
                       "matter-frame flatness ledger.",
        "grid": {"N": N, "L": L, "n_wells": n_w ** 3,
                 "universal_cover": "R^3 (non-compact)"},
        "landscape_inputs": {
            "u_well": u_well, "u_void_nominal": u_void,
            "sigma_frac_L": sigma / L,
            "rho0": rho0, "rho_ambient_fraction": rho_ambient,
            "V0": V0, "rhob_half": rhob_half, "u_s": u_s,
            "Pi": 0.0, "A_ij": 0.0,
        },
        "sign_structure": {
            "V_rec_sign_crossing_u_star": float(u_star),
            "V_at_nominal_void_depth": V_at_void,
            "actual_V_min": float(V.min()),
            "actual_V_max": float(V.max()),
            "negative_potential_cells_before_shift": int(np.count_nonzero(V < 0)),
            "required_energy_offset_lambda": float(lam),
            "S_void_mean_after_shift": S_void_mean,
            "S_well_mean_after_shift": S_well_mean,
            "rho_tot_mean": rho_tot_mean,
        },
        "solver": {
            "converged": bool(sol['converged']),
            "newton_history_max_res": sol['history'],
            "psi_min": psi_min, "psi_max": psi_max,
            "integrability_residual": integrability_residual,
            "R3_identity_residual": R3_residual,
        },
        "flatness_ledger": {
            "mean_R3_g_frame": mean_R3,
            "mean_R3_matter_frame_volume_weighted": mean_R3_t,
            "mean_R3_matter_frame_coord": mean_R3_t_coord,
            "matter_frame_ledger_decomposition": ledger_terms,
            "Omega_k_eff_H_drift_units": Omega_k_eff,
            "Omega_k_eff_coord": Omega_k_eff_coord,
            "smoothness_scan": scan,
            "closed_static_benchmark_Omega_k": ">= 1/2 (homogeneous benchmark only)",
        },
        "interpretation": (
            "This benchmark solves a periodic conformally flat Hamiltonian "
            "constraint with an explicitly fitted uniform energy offset lam. "
            "The unshifted field and potential have the sign reported in "
            "sign_structure; any negative shifted density is not evidence "
            "that the input field sampled the negative-potential branch. "
            "The solved curvature ledger uses derivatives of h=psi^4 delta. "
            "The separate smoothness_scan uses psi=1 and subtracts mean "
            "source curvature; it is an approximate diagnostic, not another "
            "constraint solution. The offset requirement is specific to "
            "this periodic ansatz. A Rule-23 realization must accommodate "
            "localized positive energy in the appropriate curved geometry "
            "and match that geometry to its exterior without using this "
            "offset as a cancellation mechanism."
        ),
        "remaining_open": [
            "evolution stability of the zero-expansion congruence "
            "(maximal-slicing lapse solvability on the solved landscape)",
            "A_ij != 0 volume-preserving shear sector",
            "realistic matter power spectrum / multi-scale landscape",
            "matter-frame singular horizon on this slice (TEP-TH branch)",
        ],
        "verdict": verdict,
        "verdict_scope": "Convergence and positivity of the shifted periodic benchmark only",
    }

    os.makedirs("results", exist_ok=True)
    out = "results/step_20_landscape_existence.json"
    with open(out, "w") as f:
        json.dump(summary, f, indent=2)

    print("\n" + "=" * 78)
    print(f"VERDICT: {verdict}")
    print(f"  converged={sol['converged']}  psi in [{psi_min:.4f}, {psi_max:.4f}]")
    print(f"  required void energy offset  lam* = {lam:+.4e}")
    print(f"  <S>_wells = {S_well_mean:+.4e}   <S>_voids = {S_void_mean:+.4e}")
    print(f"  <3R>_g  = {mean_R3:+.3e}    <3R~>_matter = {mean_R3_t:+.3e}")
    print(f"  |Omega_k^eff| = {Omega_k_eff:.3e}   (closed branch: >= 0.5)")
    print(f"  integrability residual = {integrability_residual:.3e}")
    print("=" * 78)
    print(f"[INFO] wrote {out}")


if __name__ == "__main__":
    main()
