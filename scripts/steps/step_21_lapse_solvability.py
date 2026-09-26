#!/usr/bin/env python3
"""
STEP 21 — Maximal-slicing lapse solvability on the R8 landscape
===============================================================

R8 (step_20) solved the Hamiltonian constraint on a 48^3 periodic cell
whose universal cover is the non-compact R^3 slice.  The stated remainder
for the existence construction was whether the maximal-slicing lapse
equation admits a positive solution on that solved landscape — i.e.
whether a time development with K = 0 can begin at all.

On a maximal slice (K = 0, A_ij = 0) the lapse equation is

    lap alpha = alpha * q(x),   q(x) = rho_m - |grad u|^2/2 - 2(V + lam*)

in the slice's dimensionless units (static scalar field: its contribution
to rho + S_i^i is -|grad u|^2/2 - 2V; matter is pressureless; the
potential floor is the constraint's own closed combination V + lam*).

Existence:  (-lap + q) alpha = mu* alpha  is a standard elliptic
eigenproblem.  Its principal eigenvalue mu* is the constant lapse offset
at which the equation has a zero mode, and the ground state of an
elliptic operator on the cell is simple and strictly positive — so a
single smallest-eigenvalue solve both determines mu* and delivers a
positive lapse.  Physically mu* is the lapse-sector's effective Lambda
term: the scalar floor enters rho + S_i^i with a different sign
combination than the Hamiltonian constraint's rho - S, so comparing
mu* against the constraint offset lam* tests whether one scalar floor
serves both conditions (one theory, one solution) or the lapse demands a
separate offset.

Outputs: results/step_21_lapse_solvability.json
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


def main():
    print("[INFO] step_21_lapse_solvability — maximal-slicing lapse test")

    N, L, n_w = 48, 1.0, 4
    sigma, u_well, u_void = L / 20.0, 1.0, 0.8
    rho0, rho_ambient, V0, rhob_half, u_s = 2.0, 0.05, 2.0, 1.0, 15.0

    ld = _step20.build_landscape(N, L, n_w, sigma, u_well, u_void, rho0,
                                 rho_ambient, V0, rhob_half, u_s)
    S0, g, rho_m, V = ld["S0"], ld["g"], ld["rho_m"], ld["V"]
    LAP = fast_fd_laplacian(N, ld["dx"])
    sol = _step20.solve_landscape(S0, g, LAP, verbose=False)
    psi, lam_star = sol["psi"], sol["lam"]
    print(f"  constraint re-solved: psi in [{psi.min():.4f},{psi.max():.4f}],"
          f" lam* = {lam_star:.4f}, converged={sol['converged']}")

    # lapse source on the solved slice (static scalar, pressureless matter)
    q = (rho_m - 0.5 * g - 2.0 * (V + lam_star)).ravel()

    # principal eigenpair of  A = -lap + q :  mu* = smallest eigenvalue
    A = -LAP + sp.diags(q)
    w, v = spla.eigsh(A, k=1, which="SA", tol=1e-10, maxiter=20000)
    mu_star = float(w[0])
    alpha = v[:, 0]
    alpha *= np.sign(alpha[np.argmax(np.abs(alpha))])
    resid = float(np.abs(A @ alpha - mu_star * alpha).max())
    print(f"  mu* = {mu_star:+.4f}  (constraint offset lam* = {lam_star:+.4f})")
    print(f"  lapse: alpha min {alpha.min():+.3e}  max {alpha.max():.3f}  "
          f"eig-residual {resid:.2e}")

    alpha_pos = bool(alpha.min() > 0)
    verdict = "PASS" if alpha_pos else "FAIL"

    out = {
        "step": 21,
        "name": "lapse_solvability_on_R8_landscape",
        "constraint": {
            "psi_min": float(psi.min()), "psi_max": float(psi.max()),
            "lam_star": float(lam_star), "converged": bool(sol["converged"]),
        },
        "lapse_source": {
            "form": "q = rho_m - |grad u|^2/2 - 2(V + lam*)",
            "q_min": float(q.min()), "q_max": float(q.max()),
            "q_mean": float(q.mean()),
        },
        "principal_eigenpair": {
            "mu_star": mu_star,
            "eig_residual": resid,
            "alpha_min": float(alpha.min()), "alpha_max": float(alpha.max()),
            "alpha_strictly_positive": alpha_pos,
        },
        "verdict": verdict,
        "interpretation": (
            "mu* is the constant lapse-sector offset at which the maximal-"
            "slicing lapse equation lap alpha = (q - mu*)alpha admits its "
            "positive zero mode on the solved landscape.  The scalar floor "
            "enters the lapse source with a different sign combination "
            "than the constraint, so mu* need not equal lam* = %.4f; the "
            "ratio mu*/lam* = %.3f quantifies how close the two conditions "
            "sit.  What remains unproven: stability of the zero-expansion "
            "congruence under time evolution (this solves the slice, not "
            "the dynamics) and the A_ij != 0 shear sector."
            % (lam_star, mu_star / lam_star if lam_star else float('nan'))),
    }
    out_path = Path(__file__).resolve().parent.parent.parent / "results" / \
        "step_21_lapse_solvability.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=1)
    print(f"  verdict: {verdict} -> {out_path.name}")


if __name__ == "__main__":
    main()
