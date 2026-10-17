#!/usr/bin/env python3
"""
STEP 21 — Canonical landscape lapse: admissibility and evolving CMC branch
=======================================================================

On the step20 Hamiltonian slice (K_ij = 0, scalar momentum Pi = 0 and
pressureless matter at rest), the momentum constraint is identically zero.
The exact canonical source is q = (rho_m - 2 V_eff)/(2 M_Pl^2); spatial
scalar gradients cancel in rho_phi + S_phi. The covariant Laplacian is
D_h^2 N = psi^-6 div(psi^2 grad N).

A positive q on a periodic cell obstructs a positive lapse solving
(-D_h^2 + q) N = 0. Shifting the source to its principal eigenvalue is a
numerical diagnostic, not a physical completion: at fixed Hamiltonian
budget its required Pi^2 must be nonnegative.

Instead solve (-D_h^2 + q) N = C with physical-volume mean(N) = 1.
In the ADM convention d_t K = -D_h^2 N + N q on this K_ij = 0 slice,
C is a spatially constant initial Kdot. This is a positive lapse for an
evolving constant-mean-curvature foliation, without changing the source
or introducing negative kinetic energy. It establishes initial-slice
lapse consistency, not an eternal non-expanding cosmos, a rolling
noncanonical branch or long-time evolution. The prescribed scalar profile
is initial data; its subsequent acceleration is determined by evolution.
"""

import json
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from pathlib import Path

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "step_20", Path(__file__).parent / "step_20_landscape_existence.py")
_step20 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_step20)


def fast_fd_laplacian(N, dx):
    """Vectorized periodic 7-point Laplacian — same stencil as step_20's
    fd_laplacian but O(N^3) vectorized instead of Python loops."""
    n = N ** 3
    i0 = np.arange(n)
    ii = i0 // (N * N)
    jj = (i0 // N) % N
    kk = i0 % N
    c = 1.0 / dx ** 2
    rows = [i0]
    cols = [i0]
    vals = [-6.0 * c * np.ones(n)]
    for di, dj, dk in ((1, 0, 0), (-1, 0, 0), (0, 1, 0),
                       (0, -1, 0), (0, 0, 1), (0, 0, -1)):
        cols.append((((ii + di) % N) * N * N + ((jj + dj) % N) * N
                     + ((kk + dk) % N)))
        rows.append(i0)
        vals.append(c * np.ones(n))
    return sp.csr_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(n, n))


def build_covariant_neg_div_psi2_grad(N, dx, psi_3d):
    """Construct sparse symmetric positive semi-definite matrix for
    -div( psi^2(x) grad N ) on a periodic N^3 grid using face-averaged weights.
    Dividing by psi^6 gives the exact conformal 3-metric Laplacian -D_h^2 N."""
    n = N ** 3
    w = (psi_3d.reshape((N, N, N))) ** 2
    c = 1.0 / (dx ** 2)
    idx = np.arange(n).reshape((N, N, N))
    rows, cols, vals = [], [], []
    diag = np.zeros((N, N, N), dtype=float)
    for axis in (0, 1, 2):
        for shift in (1, -1):
            w_nb = np.roll(w, -shift, axis=axis)
            idx_nb = np.roll(idx, -shift, axis=axis)
            w_face = 0.5 * (w + w_nb) * c
            rows.append(idx.ravel())
            cols.append(idx_nb.ravel())
            vals.append(-w_face.ravel())
            diag += w_face
    rows.append(idx.ravel())
    cols.append(idx.ravel())
    vals.append(diag.ravel())
    return sp.csc_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(n, n))


def solve_unshifted_maximal_lapse(N=48, L=1.0, n_w=4):
    sigma, u_well, u_void = L / 20.0, 1.0, 0.8
    rho0, rho_ambient, V0, rhob_half, u_s = 2.0, 0.05, 2.0, 1.0, 15.0

    ld = _step20.build_landscape(N, L, n_w, sigma, u_well, u_void, rho0,
                                 rho_ambient, V0, rhob_half, u_s)
    S0, g, rho_m, V = ld["S0"], ld["g"], ld["rho_m"], ld["V"]
    LAP = fast_fd_laplacian(N, ld["dx"])
    sol = _step20.solve_landscape(S0, g, LAP, verbose=False)
    psi, lam_star = sol["psi"], sol["lam"]
    hamiltonian_residual = float(np.max(np.abs(
        (LAP @ psi.ravel()).reshape(psi.shape)
        + psi**5 * (S0 + lam_star) / (4 * _step20.M_PL**2)
        + psi * g / (8 * _step20.M_PL**2))))
    psi_flat = psi.ravel()
    psi6 = psi_flat ** 6

    # Covariant kinetic operator K_cov = -div(psi^2 grad)
    K_cov = build_covariant_neg_div_psi2_grad(N, ld["dx"], psi)

    # In the Hamiltonian constraint, S_ham = rho_m + V + Pi^2/2 + lam_H = S0 + lam_star,
    # so V_eff(x) = V(x) + lam_H = V(x) + lam_star - Pi^2/2.
    # In the exact canonical scalar stress tensor, rho_phi + S_phi = 2 Pi^2 - 2 V_eff(x)
    # (static spatial gradient terms +|D_h u|^2/2 and -|D_h u|^2/2 cancel identically!),
    # so the lapse source is:
    #   q_phys(x; Pi^2) = (1 / 2 M_Pl^2) [ rho_m(x) + 2 Pi^2 - 2 V_eff(x) ]
    #                   = (1 / 2 M_Pl^2) [ rho_m(x) - 2 (V(x) + lam_star) + 3 Pi^2 ].
    # Let delta_q = (3 / 2 M_Pl^2) Pi^2 - (1 / 2 M_Pl^2) * (shift) or solve directly for
    # the uniform scalar kinetic/potential combination delta_s = (3/2) Pi^2 - delta_V0
    # that sets the principal eigenvalue of (-D_h^2 + q_phys) N = 0 to mu = 0:
    q_base = 0.5 * (rho_m - 2.0 * (V + lam_star)).ravel()

    # For (-D_h^2 + q_base + delta_q) N = 0  <=>  [-psi^-6 div(psi^2 grad N) + q_base N] = (-delta_q) N:
    # Using the symmetric whitening similarity transform v = psi^3 N, the operator
    #   A_cov = diag(psi^-3) K_cov diag(psi^-3) + diag(q_base)
    # is standard symmetric (no M matrix needed!), so eigsh(A_cov, which="SA") converges in <1s!
    inv_psi3 = sp.diags(1.0 / (psi_flat ** 3), format="csc")
    A_cov_base = (inv_psi3 @ K_cov @ inv_psi3 + sp.diags(q_base, format="csc")).tocsc()
    rng = np.random.default_rng(20261017)
    v0_init = rng.normal(size=N ** 3)
    w_gen, v_gen = spla.eigsh(A_cov_base, k=1, which="SA", tol=1e-12, maxiter=20000, v0=v0_init)
    mu_unadjusted = float(w_gen[0])
    delta_q_req = -mu_unadjusted

    # Spectrally shifted diagnostic: physical admissibility is checked below.
    q_phys = q_base + delta_q_req
    H_phys = (K_cov + sp.diags(psi6 * q_phys, format="csc")).tocsc()
    alpha = v_gen[:, 0] / (psi_flat ** 3)
    alpha *= np.sign(alpha.mean())
    alpha /= float(alpha.mean())

    # Verify UNSHIFTED covariant lapse equation (-D_h^2 + q_phys) alpha = 0:
    lhs_unshifted_cov = (H_phys @ alpha) / psi6
    unshifted_rel_residual = float(np.linalg.norm(lhs_unshifted_cov) / np.linalg.norm(alpha))
    unshifted_max_residual = float(np.max(np.abs(lhs_unshifted_cov)))

    # Physical branch: retain the Hamiltonian source at Pi=0 and allow Kdot.
    # H_base = sqrt(h) (-D_h^2 + q_base). For q_base > 0 the
    # discrete maximum principle gives a positive inverse on this grid.
    H_base = (K_cov + sp.diags(psi6 * q_base, format="csc")).tocsc()
    preconditioner = sp.diags(1.0 / H_base.diagonal(), format="csc")
    lapse_raw, cg_info = spla.cg(H_base, psi6, M=preconditioner,
                                rtol=1e-11, atol=0.0, maxiter=20000)
    if cg_info != 0:
        raise RuntimeError(f"CMC lapse solve did not converge: {cg_info}")
    volume_mean = float(np.dot(psi6, lapse_raw) / psi6.sum())
    cmc_alpha = lapse_raw / volume_mean
    cmc_Kdot = 1.0 / volume_mean
    cmc_lhs = (H_base @ cmc_alpha) / psi6
    cmc_residual = float(np.linalg.norm(cmc_lhs - cmc_Kdot) /
                         max(np.linalg.norm(cmc_lhs), 1e-30))
    required_pi2 = 2.0 * _step20.M_PL ** 2 * delta_q_req / 3.0
    maximal_admissible = bool(sol["converged"] and required_pi2 >= 0
                              and alpha.min() > 0 and unshifted_rel_residual < 1e-9)
    cmc_valid = bool(sol["converged"] and q_base.min() > 0
                     and cmc_alpha.min() > 0 and cmc_residual < 1e-8)

    # Also record the legacy flat-Laplacian shifted eigenproblem for backward reference:
    q_legacy = (rho_m - 0.5 * g - 2.0 * (V + lam_star)).ravel()
    A_legacy = -LAP + sp.diags(q_legacy)
    w_leg, v_leg = spla.eigsh(A_legacy, k=1, which="SA", tol=1e-10, maxiter=20000, v0=v0_init)
    mu_legacy = float(w_leg[0])
    alpha_leg = v_leg[:, 0]
    alpha_leg *= np.sign(alpha_leg[np.argmax(np.abs(alpha_leg))])
    resid_legacy = float(np.abs(A_legacy @ alpha_leg - mu_legacy * alpha_leg).max())

    return {
        "N": N,
        "hamiltonian_max_residual": hamiltonian_residual,
        "required_Pi2_for_maximal_closure": float(required_pi2),
        "maximal_canonical_admissible": maximal_admissible,
        "q_base": q_base,
        "cmc_alpha": cmc_alpha,
        "cmc_Kdot": float(cmc_Kdot),
        "cmc_relative_residual": cmc_residual,
        "cmc_valid": cmc_valid,
        "cmc_volume_mean_lapse": float(np.dot(psi6, cmc_alpha) / psi6.sum()),
        "psi": psi,
        "lam_star": float(lam_star),
        "converged": bool(sol["converged"]),
        "q_phys": q_phys,
        "delta_q_req": float(delta_q_req),
        "mu_unadjusted": mu_unadjusted,
        "mu_phys_unshifted": 0.0,
        "alpha": alpha,
        "unshifted_rel_residual": unshifted_rel_residual,
        "unshifted_max_residual": unshifted_max_residual,
        "mu_legacy": mu_legacy,
        "resid_legacy": resid_legacy,
    }


def summarize(res):
    # Preferred-foliation khronon/cuscuton drift term P_drift = -mu_kh^2 sqrt(2X)
    # on the timelike cosmic roll (Pi > 0): its Hamiltonian energy density
    # 2 X P_{drift,X} - P_drift vanishes identically on the homogeneous roll
    # (leaving the Lichnerowicz solution psi and lambda* unchanged), while its
    # spatial stress trace S_drift = -3 mu_kh^2 Pi < 0 supplies delta_q = -mu_0 < 0
    # with strictly positive real drift mu_kh^2 Pi = +(2/3) M_Pl^2 mu_0 > 0,
    # solving (-D_h^2 + q_kh) N = 0 (d_t K = 0) with N(x) > 0 to machine precision.
    khronon_real_drift_density = float(2.0 * (_step20.M_PL ** 2) * res["mu_unadjusted"] / 3.0)
    khronon_static_valid = bool(
        res["converged"]
        and khronon_real_drift_density > 0.0
        and float(res["alpha"].min()) > 0.0
        and res["unshifted_rel_residual"] < 1e-9
    )
    return {
        "N": res["N"],
        "hamiltonian_converged": res["converged"],
        "hamiltonian_max_residual": res["hamiltonian_max_residual"],
        "psi_range": [float(res["psi"].min()), float(res["psi"].max())],
        "lambda_star": res["lam_star"],
        "maximal_slice_diagnostic": {
            "unadjusted_principal_eigenvalue": res["mu_unadjusted"],
            "source_shift_required": res["delta_q_req"],
            "required_Pi2": res["required_Pi2_for_maximal_closure"],
            "real_canonical_momentum": res["required_Pi2_for_maximal_closure"] >= 0,
            "adjusted_operator_relative_residual": res["unshifted_rel_residual"],
            "canonical_closure_established": res["maximal_canonical_admissible"],
        },
        "evolving_cmc_initial_slice": {
            "source_changed": False,
            "scalar_Pi": 0.0,
            "K_ij_initial": 0.0,
            "matter_momentum": 0.0,
            "momentum_constraint": "identically satisfied on this time-symmetric slice",
            "q_base_range": [float(res["q_base"].min()), float(res["q_base"].max())],
            "lapse_range": [float(res["cmc_alpha"].min()), float(res["cmc_alpha"].max())],
            "volume_mean_lapse": res["cmc_volume_mean_lapse"],
            "Kdot_initial_code_units": res["cmc_Kdot"],
            "lapse_equation_relative_residual": res["cmc_relative_residual"],
            "initial_slice_lapse_consistent": res["cmc_valid"],
        },
        "noncanonical_khronon_static_slice": {
            "operator_sector": "P_drift(X) = -mu_kh^2 sqrt(2X) (preferred-foliation khronon/cuscuton drift, 2X P_{,X} - P == 0)",
            "hamiltonian_source_unchanged": True,
            "Kdot_static_slice": 0.0,
            "required_positive_khronon_drift_mu2_Pi": khronon_real_drift_density,
            "real_timelike_drift_admissible": bool(khronon_real_drift_density > 0.0),
            "static_lapse_range": [float(res["alpha"].min()), float(res["alpha"].max())],
            "unshifted_static_lapse_relative_residual": res["unshifted_rel_residual"],
            "static_maximal_slice_established_on_khronon_branch": khronon_static_valid,
            "scope": "Lapse equation with imposed stress shift on one slice; scalar drift equation and time preservation not solved here.",
            "full_drift_equations_solved": False,
            "eternal_static_solution_established": False,
        },
    }


def main():
    grids = [summarize(solve_unshifted_maximal_lapse(N=n)) for n in (16, 48)]
    out = {
        "step": 21,
        "name": "canonical_lapse_admissibility_and_evolving_cmc_initial_slice",
        "scope": __doc__,
        "grids": grids,
        "verdict": "PASS_EVOLVING_CMC_INITIAL_SLICE" if all(
            x["evolving_cmc_initial_slice"]["initial_slice_lapse_consistent"] for x in grids
        ) else "INITIAL_SLICE_NOT_VERIFIED",
        "maximal_static_canonical_closure": all(
            x["maximal_slice_diagnostic"]["canonical_closure_established"] for x in grids),
        "interpretation": (
            "A positive lapse solves the unmodified canonical initial-slice CMC equation "
            "with a nonzero spatially uniform initial Kdot. The eigenvalue-shifted maximal "
            "diagnostic requires negative Pi squared on these backgrounds and is not a "
            "physical canonical completion. This result supports an evolving initial-data "
            "branch; it does not establish eternal zero expansion or the complete "
            "selected noncanonical action."
        ),
    }
    out_path = Path(__file__).resolve().parents[2] / "results/step_21_lapse_solvability.json"
    out_path.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
