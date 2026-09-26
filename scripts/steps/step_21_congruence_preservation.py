#!/usr/bin/env python3
"""
step_21_congruence_preservation.py
==================================
Evolution-side completion of the inhomogeneous-slice construction
(issues 0-1, continuation of step_20): whether the zero-expansion
congruence can be preserved on the constraint-satisfying landscape.

Preservation of K = 0 under evolution imposes the maximal-slicing
lapse equation

    D^2 N = N ( A_ij A^ij + (rho_tot + S_tot) / 2 M_Pl^2 ),

which on the periodic cell integrates to the mean condition
int N f = 0 with

    f(x) = A^2/M^2 + ( rho_m + 2 Pi^2 - 2 V ) / (2 M^2),

using rho + S = rho_m + 2 Pi^2 - 2 V for the canonical scalar sector
(the |grad phi|^2 terms cancel between rho and S).

Structural obstruction resolved here: the Hamiltonian integrability
condition forces <V> ~ -<rho_m> (negative mean potential, step_20),
while mean zero-acceleration wants <V> ~ +<rho_m>/2.  With uniform
weighting these are incompatible; the compatible class is a landscape
whose lapse-equation source f(x) CHANGES SIGN — f < 0 inside
potential-dominated wells (V > rho_m/2, the Rule-23 deep temporal
wells where the potential energy is extreme) and f > 0 in the
negative-V void sector.  This step therefore uses a deeper landscape
(V0 = 8) than step_20 so that the wells are potential-dominated.

Two preservation statements are computed on the SAME solved slice:

  (i)  Mean preservation: a positive lapse with int N f = 0 exists
       whenever f changes sign (explicit construction, verified).
  (ii) Exact preservation: 0 in spec(-Delta + f) is required.  Since
       Pi^2 shifts f uniformly, the eigenvalues shift uniformly and the
       required drift momentum is Pi*^2 = -M^2 mu_0(0).  The ground
       eigenfunction is then positive (Perron-Frobenius) and IS the
       preserving lapse.

The Hamiltonian constraint is re-solved with the tuned Pi^2 included
(the drift energy also enters rho_tot), iterating lambda* and Pi*
jointly.

Outputs results/step_21_congruence_preservation.json.
"""

import json
import os
import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import scipy.fft as fft

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from step_20_landscape_existence import (build_landscape, fd_laplacian,
                                       solve_landscape, M_PL)


def main():
    print("[INFO] step_21_congruence_preservation — lapse/congruence solve")

    # --- deeper landscape: wells must be potential-dominated (V > rho/2) ---
    # and the void sector must be deep enough to supply the integrability
    # slack natively (lam* ~ 0), so the uniform offset is not needed and
    # f = rho_m + 2 Pi^2 - 2 V changes sign across the landscape.
    L = 1.0
    n_w = 4
    sigma = L / 20.0
    u_well = 0.4
    rho0 = 0.5                # light matter loading keeps lapse window open
    rho_ambient = 0.05
    V0 = 60.0                 # deep floor (Rule-23 wells: V_w >> rho_w/2)
    rhob_half = 30.0          # keeps u* = ln(rhob_half/V0) = -ln 2
    u_s = 15.0

    # --- scan (void depth, well width) by variational bound ---
    # Preservation of K = 0 needs 0 in spec(-Delta + f) for a positive
    # lapse; Pi^2 >= 0 shifts the spectrum UP, so a physical drift
    # momentum exists iff mu_0(Pi^2 = 0) < 0, i.e. iff f supports a
    # bound state.  psi ~ 1 gives lam* ~ -<S0> - (1/2)<|grad u|^2 and
    # f ~ (rho_m - 2(V + lam*))/2M^2.  The Rayleigh quotient of the
    # well profile is an upper bound: R < 0 => mu_0 < 0 => Pi* exists.
    # Deep wells are exponentially cheap for the constraint (V_rec ~
    # -e^{-u}) but the lapse needs |f_w| sigma^2 >~ pi^2: the Rule-23
    # deep-well regime.
    def rayleigh(ld_s, sg):
        lam_est = -(float(ld_s['S0'].mean())
                    + 0.5 * float(ld_s['g'].mean()))
        f_s = (ld_s['rho_m'] - 2.0 * (ld_s['V'] + lam_est)) \
            / (2.0 * M_PL ** 2)
        w = ld_s['wells'] / ld_s['wells'].max()
        phi = w ** 2                     # smooth localised trial state
        gx = np.gradient(phi, ld_s['dx'], axis=0, edge_order=2)
        gy = np.gradient(phi, ld_s['dx'], axis=1, edge_order=2)
        gz = np.gradient(phi, ld_s['dx'], axis=2, edge_order=2)
        kin = float((gx ** 2 + gy ** 2 + gz ** 2).sum())
        pot = float((f_s * phi ** 2).sum())
        nrm = float((phi ** 2).sum())
        return kin / nrm + pot / nrm, lam_est, f_s

    scan = []
    u_v_star, sigma_star = None, None
    N_est = 40
    best_R = np.inf
    for sg in (L / 4.0, L / 3.0):
        for u_v in (-1.0, -1.4, -1.8, -2.2, -2.6):
            ld_s = build_landscape(N_est, L, n_w, sg, u_well, u_v, rho0,
                                   rho_ambient, V0, rhob_half, u_s)
            R, lam_est, f_s = rayleigh(ld_s, sg)
            scan.append({"u_void": u_v, "sigma_frac_L": sg,
                         "lambda_star_est": lam_est,
                         "rayleigh_bound": R,
                         "f_min": float(f_s.min()),
                         "f_max": float(f_s.max())})
            print(f"  scan u_v = {u_v:5.2f}  s = {sg:.3f}:  "
                  f"lam*~ {lam_est:+9.3f}  R = {R:+9.3f}  "
                  f"f in [{f_s.min():+.2f}, {f_s.max():+.2f}]")
            if R < best_R:
                best_R, u_v_star, sigma_star = R, u_v, sg

    # Mildest landscape that still binds: smallest |lam*| among R < 0
    # (the deep end is stiffer for Newton and not needed once mu_0 < 0).
    bound = [d for d in scan if d['rayleigh_bound'] < 0]
    if bound:
        best = min(bound, key=lambda d: abs(d['lambda_star_est']))
        u_v_star, sigma_star = best['u_void'], best['sigma_frac_L']
        best_R = best['rayleigh_bound']

    print(f"  selected u_v = {u_v_star}, sigma = {sigma_star:.3f} "
          f"(R = {best_R:+.3f})")

    # --- final solve on the self-consistent landscape ---
    N = 40
    sigma = sigma_star
    ld = build_landscape(N, L, n_w, sigma, u_well, u_v_star, rho0,
                         rho_ambient, V0, rhob_half, u_s)
    u, rho_m, V, g = ld['u'], ld['rho_m'], ld['V'], ld['g']
    wells = ld['wells']
    dx = ld['dx']
    LAP = fd_laplacian(N, dx)
    dVc = dx ** 3
    shape = u.shape
    u_void = u_v_star

    # --- joint (lambda*, Pi*^2) solve ---
    # mu_0(Pi^2) is affine to excellent accuracy: the constraint re-solve
    # gives lam*(Pi^2) ~ lam*_0 - Pi^2/2, so f shifts by 3 Pi^2/(2 M^2)
    # and d mu_0/d Pi^2 = 3/(2 M^2).  Newton on mu_0 = 0 converges in
    # one step up to the small psi^5-weighting residual.
    Pi2 = 0.0
    hist = []
    psi = lam = mu0 = Nvec = None
    f = np.zeros(shape)
    for outer in range(4):
        S0 = rho_m + V + 0.5 * Pi2
        sol = solve_landscape(S0, g, LAP, tol=1e-9, verbose=True)
        psi, lam = sol['psi'], sol['lam']
        if not sol['converged']:
            print("[FAIL] constraint solve did not converge")
            break

        # lapse source on the solved slice (A_ij = 0):
        f = (rho_m + 2.0 * Pi2 - 2.0 * (V + lam)) / (2.0 * M_PL ** 2)

        Aop = (-LAP + sp.diags(f.ravel(), format='csc')).tocsc()
        try:
            # shift-invert: sigma below the ground state (mu_0 >= min f
            # by Weyl) makes A - sigma*I SPD and the interior eigenvalue
            # dominant for the inverse iteration
            sig = float(f.min()) - 1.0
            ev = spla.eigsh(Aop, k=1, sigma=sig, which='LM', tol=1e-8,
                            return_eigenvectors=True)
            mu0, Nvec = float(ev[0][0]), ev[1][:, 0].reshape(shape)
        except Exception as e:
            print(f"[warn] eigsh failed ({e}); using variational bound")
            mu0, Nvec = float(f.min()), None

        # Newton update: d mu_0 / d Pi^2 = 3/(2 M^2)
        Pi2_new = Pi2 - mu0 * (2.0 * M_PL ** 2 / 3.0)
        hist.append({"Pi2": Pi2, "lambda": lam, "mu0": mu0,
                     "f_min": float(f.min()), "f_max": float(f.max())})
        print(f"  outer {outer}: lam={lam:+.4f}  Pi^2={Pi2:.4f}  "
              f"mu0={mu0:+.4e}  f in [{f.min():+.3f}, {f.max():+.3f}]")
        if abs(mu0) < 1e-6:
            Pi2 = Pi2_new
            break
        Pi2 = Pi2_new

    # --- final diagnostics on the tuned slice ---
    f_min, f_max = float(f.min()), float(f.max())
    sign_changing = (f_min < 0) and (f_max > 0)
    f_w = float(f[wells > 0.5].mean())
    f_v = float(f[wells < 0.5].mean())

    # (i) mean-preserving lapse: enhance N where f < 0 (or vice versa)
    P = float(f[f > 0].sum()) if (f > 0).any() else 0.0
    Q = float(-f[f < 0].sum()) if (f < 0).any() else 0.0
    if P >= Q and Q > 0:
        delta = (P - Q) / Q
        Nmean = np.where(f < 0, 1.0 + delta, 1.0)
    else:
        delta = (Q - P) / P if P > 0 else 0.0
        Nmean = np.where(f > 0, 1.0 + delta, 1.0)
    int_Nf = float((Nmean * f).sum() * dVc)
    mean_preserved = sign_changing and abs(int_Nf) < 1e-8

    # (ii) exact preservation: tuned Pi^2 => mu0 ~ 0; eigenfunction = lapse
    exact_lapse_positive = None
    if Nvec is not None:
        if Nvec.min() < 0:
            Nvec = -Nvec
        Nvec = Nvec / Nvec.mean()
        exact_lapse_positive = bool(Nvec.min() > 0)
        # residual of (Delta - f) N = 0
        res = (LAP @ Nvec.ravel()).reshape(shape) - f * Nvec
        eig_resid = float(np.abs(res).max())
    else:
        eig_resid = None

    # verification: lapse equation requires int N f = 0 for ANY zero mode
    int_eig_Nf = (float((Nvec * f).sum() * dVc) if Nvec is not None
                  else None)

    verdict = "PASS" if (sol['converged'] and sign_changing
                         and mean_preserved
                         and (mu0 is None or abs(mu0) < 1e-3)) else "FAIL"

    summary = {
        "step": "step_21_congruence_preservation",
        "description": "Zero-expansion congruence preservation on the "
                       "constraint-satisfying inhomogeneous slice: joint "
                       "solve of the Lichnerowicz integrability (lambda*) "
                       "and the maximal-slicing lapse zero-mode (mu_0 = 0 "
                       "via the drift momentum Pi*^2).",
        "landscape": {"V0": V0, "rhob_half": rhob_half, "u_well": u_well,
                      "u_void": u_void, "rho0": rho0, "N": N,
                      "variational_scan_N40": scan,
                      "note": "Rule-23 deep-well regime (V0 = 60): the "
                              "lapse source f changes sign only where "
                              "wells are potential-dominated, and a bound "
                              "state (mu_0 < 0 at Pi^2 = 0) requires "
                              "|f_w| sigma^2 >~ pi^2.  The scan minimises "
                              "the Rayleigh bound on mu_0."},
        "joint_solution": {
            "lambda_star": float(lam),
            "Pi_star2": float(Pi2),
            "mu_0_final": mu0,
            "psi_range": [float(psi.min()), float(psi.max())],
            "fixed_point_history": hist,
        },
        "lapse_source_structure": {
            "f_min": f_min, "f_max": f_max,
            "f_well_mean": f_w, "f_void_mean": f_v,
            "sign_changing": sign_changing,
            "interpretation": ("f < 0 in potential-dominated wells "
                               "(local expansion), f > 0 in negative-V "
                               "voids (local contraction): the eternal "
                               "slice is a steady-state circulation, not "
                               "a pointwise-static geometry") if sign_changing
                              else "f single-signed: congruence cannot be "
                                   "preserved on this landscape",
        },
        "mean_preservation": {
            "achieved": mean_preserved,
            "lapse_modulation_delta": float(delta),
            "int_Nf": int_Nf,
        },
        "exact_preservation": {
            "mu_0": mu0,
            "required_Pi2": float(Pi2),
            "eigenfunction_positive": exact_lapse_positive,
            "eigen_residual": eig_resid,
            "int_eigenfunction_f": int_eig_Nf,
        },
        "remaining_open": [
            "linearized inhomogeneous-mode growth (the landscape "
            "analogue of the Eddington mode) — requires perturbation "
            "evolution, not the initial-data level shown here",
            "A_ij != 0 sector",
            "long-time stationarity of the landscape statistics",
        ],
        "verdict": verdict,
    }

    os.makedirs("results", exist_ok=True)
    out = "results/step_21_congruence_preservation.json"
    with open(out, "w") as fo:
        json.dump(summary, fo, indent=2)

    print("\n" + "=" * 78)
    print(f"VERDICT: {verdict}")
    print(f"  lambda* = {lam:+.4f}   Pi*^2 = {Pi2:.4f}   mu0 = {mu0:+.3e}")
    print(f"  f: wells {f_w:+.3f} / voids {f_v:+.3f}  (sign-changing "
          f"= {sign_changing})")
    print(f"  mean-preserved = {mean_preserved}  (int Nf = {int_Nf:.2e})")
    if Nvec is not None:
        print(f"  eigenfunction lapse > 0: {exact_lapse_positive}")
    print("=" * 78)
    print(f"[INFO] wrote {out}")


if __name__ == "__main__":
    main()
