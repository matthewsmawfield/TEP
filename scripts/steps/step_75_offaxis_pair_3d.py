#!/usr/bin/env python3
"""Off-axis two-centre pair solve on a 3D Cartesian grid.

Step_59 solved the pair axis PARALLEL to the ambient gradient (the
only axisymmetric configuration on the cylinder grid).  This step
removes that restriction: the same nonlinear flux law

    div J(a) = rho ,   J_i = a_i f(|a|),
    f(q) = eps_reg + k q/sqrt(2) + q^2   (two-branch)
    f(q) = 1 + q^2                        (baseline control)
    f(q) = 1                              (linear reference)

is solved on a Cartesian box for a pair whose separation vector sits
at angle theta to the ambient gradient u0 z_hat.  Dimensionless units
g_t = 1, d = 1, as in step_59; parameters matched (sigma_src = 0.15,
u0 = 0.15, charge = 4 pi per star) so the theta = 0 run is a direct
cross-check on step_59's axisymmetric y_pair(d=1, u0=0.15) = 0.1524.

Sources: two Gaussians of width SIG_SRC at r = +- (d/2) n_hat,
n_hat = (sin theta, 0, cos theta).  Ambient enters as the u0 z term in
a = grad(phi), phi = u0 z + psi, psi -> 0 at the box boundary.

Observables:
  * mutual force along n_hat from the momentum-conserving stress flux
    T_ij = f a_i a_j - delta_ij W,  W(q) = int_0^q f q' dq',
    through a sphere around star B, with the B-only solution
    subtracted (removes the ambient body force and self-terms, as in
    step_59); y_pair(theta) = F_pair / F_pair_linear;
  * bridge minimum |a| on the segment between the sources;
  * the parallel-orientation check: theta = 0 must reproduce the
    step_59 scaling to discretisation accuracy.

Coarse-grid first pass (N=64, sigma_src ~ 1.2 cells): the angular
dependence of the response is the deliverable, not the converged
asymptote; absolute values carry the discretisation error shared with
the sigma=0.15 source representation.
"""
import json
import os
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.ndimage import map_coordinates
from scipy.integrate import trapezoid

K = 16.03
EPS_REG = 1e-8
SQRT2 = np.sqrt(2.0)

# grid: cubic box, half-length LBOX, d = 1
N = 64
LBOX = 4.0
dx = 2.0 * LBOX / N
X = (np.arange(N) + 0.5) * dx - LBOX
XX, YY, ZZ = np.meshgrid(X, X, X, indexing='ij')

SIG_SRC = 0.15                  # matched to step_59
Q_STAR = 4.0 * np.pi            # unit-charge stars as in step_59
D_PAIR = 1.0
U0 = 0.15                       # step_59 lowest ambient: direct check
R_CAP = max(D_PAIR / 3.0, 4.0 * SIG_SRC)   # stress-sphere radius (step_59)


def f_flux(q, sector):
    if sector == "linear":
        return np.ones_like(q)
    if sector == "baseline":
        return 1.0 + q * q
    return EPS_REG + K * q / SQRT2 + q * q


def fp_flux(q, sector):
    if sector == "linear":
        return np.zeros_like(q)
    if sector == "baseline":
        return 2.0 * q
    return K / SQRT2 + 2.0 * q


def W_of_q(q, sector):
    """W(q) = int_0^q f(q') q' dq'; T_ij = f a_i a_j - delta_ij W."""
    q = np.asarray(q, dtype=float)
    out = np.empty_like(q)
    qcap = max(1.0, float(np.abs(q).max()) * 1.05)
    qs = np.linspace(0.0, qcap, 2001)
    ftab = f_flux(qs, sector)
    wtab = np.concatenate([[0.0],
                           np.cumsum(0.5 * (ftab[:-1] * qs[:-1]
                                            + ftab[1:] * qs[1:])
                                     * np.diff(qs))])
    out.ravel()[:] = np.interp(np.abs(q).ravel(), qs, wtab)
    return out


def pair_dir(theta):
    return np.array([np.sin(theta), 0.0, np.cos(theta)])


def make_rho(theta, stars=(+1, -1)):
    n = pair_dir(theta)
    rho = np.zeros_like(XX)
    for sgn in stars:
        c = np.stack([XX, YY, ZZ]) + sgn * 0.5 * D_PAIR * n[:, None,
                                                          None, None]
        r = np.sqrt(np.sum(c**2, axis=0))
        rho += np.exp(-r**2 / (2 * SIG_SRC**2))
    rho *= Q_STAR / (2.0 * np.pi) ** 1.5 / SIG_SRC**3
    return rho


def face_fluxes(psi, u0, sector):
    """Face fluxes on the Cartesian grid."""
    d = [np.gradient(psi, dx, axis=ax) for ax in range(3)]
    d[2] = d[2] + u0

    out = []
    for ax in range(3):
        a_ax_f = (np.roll(psi, -1, axis=ax) - psi) / dx
        if ax == 2:
            a_ax_f = a_ax_f + u0
        a_tr = []
        for bx in range(3):
            if bx == ax:
                continue
            a_tr.append(0.5 * (d[bx] + np.roll(d[bx], -1, axis=ax)))
        a1, a2 = a_tr
        am2 = a_ax_f**2 + a1**2 + a2**2
        q = np.sqrt(np.maximum(am2, 1e-30))
        out.append({
            "ax": ax,
            "J": a_ax_f * f_flux(q, sector),
            "f": f_flux(q, sector),
            "fp": fp_flux(q, sector),
            "q": q,
            "a_ax": a_ax_f,
            "a_tr": (a1, a2),
            "tr_axes": tuple(bx for bx in range(3) if bx != ax),
        })
    return out


def residual(psi, u0, sector, rho):
    faces = face_fluxes(psi, u0, sector)
    R = np.zeros_like(psi)
    for fc in faces:
        ax = fc["ax"]
        R += (fc["J"] - np.roll(fc["J"], 1, axis=ax)) / dx
    return R - rho


def build_jacobian(psi, u0, sector):
    """Newton Jacobian of div J.  Axial stiffness
    Z_aa = f + fp a_a^2/q multiplies the directional second difference;
    cross terms Z_ab = fp a_a a_b/q enter through the centred
    transverse derivative averaged over the two cells sharing each
    face.  Validated against finite differences of `residual`."""
    faces = face_fluxes(psi, u0, sector)
    Nv = psi.size
    idx = np.arange(Nv).reshape(psi.shape)

    rows, cols, vals = [], [], []

    def add(r, c, v):
        rows.append(np.asarray(r).ravel())
        cols.append(np.asarray(c).ravel())
        vals.append(np.asarray(v).ravel())

    for fc in faces:
        ax = fc["ax"]
        q, f_, fp = fc["q"], fc["f"], fc["fp"]
        a_ax, a_tr, tr_axes = fc["a_ax"], fc["a_tr"], fc["tr_axes"]

        # axial part: R_i = (J_i - J_{i-1})/dx
        Z_aa = f_ + fp * a_ax**2 / q
        coef_a = Z_aa / dx**2
        ci = idx
        ip = np.roll(idx, -1, axis=ax)
        add(ci, ip, +coef_a)
        add(ci, ci, -coef_a)
        Z_am = np.roll(Z_aa, 1, axis=ax)
        coef_am = Z_am / dx**2
        im = np.roll(idx, 1, axis=ax)
        add(ci, im, +coef_am)
        add(ci, ci, -coef_am)

        # cross terms: a_b(face i) = (psi_{i..b+1} - psi_{i..b-1}
        #   + psi_{i+1..b+1} - psi_{i+1..b-1}) / (4 dx)
        for bt, a_b in zip(tr_axes, a_tr):
            Z_ab = fp * a_ax * a_b / q
            w = Z_ab / (4.0 * dx)
            for shift_a, sgn_div in ((0, +1.0), (1, -1.0)):
                base = np.roll(idx, shift_a, axis=ax)
                wloc = np.roll(w, shift_a, axis=ax)
                for eps in (0, -1):
                    for sgn_b in (+1.0, -1.0):
                        # np.roll(arr, +1)[i] = arr[i-1], so a roll of
                        # -sgn_b places the column at the +sgn_b
                        # neighbour along bt
                        c = np.roll(np.roll(base, eps, axis=ax),
                                    -int(sgn_b), axis=bt)
                        coef = sgn_div * wloc * sgn_b / dx
                        add(ci, c, coef)

    A = sp.csr_matrix(
        (np.concatenate(vals),
         (np.concatenate(rows), np.concatenate(cols))),
        shape=(Nv, Nv))
    return A


def dirichlet_project(A):
    A = A.tolil()
    g = np.arange(N**3).reshape(N, N, N)
    bd = np.concatenate([
        g[0, :, :].ravel(), g[-1, :, :].ravel(),
        g[:, 0, :].ravel(), g[:, -1, :].ravel(),
        g[:, :, 0].ravel(), g[:, :, -1].ravel(),
    ])
    for c in bd:
        A.rows[c] = [int(c)]
        A.data[c] = [1.0]
    return A.tocsr()


def solve(u0, sector, rho, itmax=60, psi0=None, verbose=False):
    psi = np.zeros_like(XX) if psi0 is None else psi0.copy()
    bmask = np.zeros_like(psi, bool)
    bmask[0, :, :] = bmask[-1, :, :] = True
    bmask[:, 0, :] = bmask[:, -1, :] = True
    bmask[:, :, 0] = bmask[:, :, -1] = True
    rn = None
    for it in range(itmax):
        R = residual(psi, u0, sector, rho)
        R[bmask] = psi[bmask]
        rn = np.linalg.norm(R[~bmask])
        if verbose or it % 5 == 0:
            print(f"    [{sector}] it={it} |R|={rn:.3e}", flush=True)
        if rn < 1e-6 * np.linalg.norm(rho):
            break
        A = dirichlet_project(build_jacobian(psi, u0, sector))
        rhs = -R.reshape(-1)
        diag = A.diagonal()
        M = spla.LinearOperator(A.shape,
                                matvec=lambda v: v / np.maximum(
                                    np.abs(diag), 1e-12))
        dpsi, info = spla.cg(A, rhs, M=M, rtol=1e-7,
                             maxiter=4000)
        if info != 0:
            try:
                dpsi = spla.spsolve(A, rhs)
            except Exception:
                dpsi = spla.lsmr(A, rhs, atol=1e-8, btol=1e-8,
                                 maxiter=400)[0]
        dpsi = dpsi.reshape(psi.shape)
        dmax = np.max(np.abs(dpsi))
        if dmax > 1.0:
            dpsi *= 1.0 / dmax
        alpha = 1.0
        accepted = False
        for _ in range(30):
            cand = psi + alpha * dpsi
            rc = residual(cand, u0, sector, rho)
            rc[bmask] = 0.0
            if np.linalg.norm(rc[~bmask]) < rn:
                psi = cand
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            psi = psi + dpsi * 1e-7
    rn = np.linalg.norm(residual(psi, u0, sector, rho)[~bmask])
    if verbose:
        print(f"    [{sector}] done |R|={rn:.3e}", flush=True)
    return psi, rn / np.linalg.norm(rho)


def sphere_force(psi, u0, centre, sector, npts=61):
    """Momentum flux of T_ij = f a_i a_j - delta_ij W through a sphere
    of radius R_CAP centred at `centre`.  Returns the force vector."""
    grad = [np.gradient(psi, dx, axis=ax) for ax in range(3)]
    grad[2] = grad[2] + u0
    # uniform spherical sampling (equal-area via the golden-spiral
    # projection has uneven poles; use a theta-phi product grid)
    th = np.linspace(1e-6, np.pi - 1e-6, npts)
    ph = np.linspace(0.0, 2.0 * np.pi, 2 * npts, endpoint=False)
    TH, PH = np.meshgrid(th, ph, indexing='ij')
    nx = np.sin(TH) * np.cos(PH)
    ny = np.sin(TH) * np.sin(PH)
    nz = np.cos(TH)
    pts = np.stack([nx, ny, nz]) * R_CAP + centre[:, None, None]
    ii = (pts + LBOX) / dx - 0.5
    a = [map_coordinates(grad[ax], ii, order=1, mode='nearest')
         for ax in range(3)]
    q = np.sqrt(a[0]**2 + a[1]**2 + a[2]**2)
    fq = f_flux(q, sector)
    W = W_of_q(q, sector)
    # T.n where n = (nx, ny, nz):  (T.n)_i = f a_i (a.n) - n_i W
    an = a[0] * nx + a[1] * ny + a[2] * nz
    Tn = [fq * a[k] * an - [nx, ny, nz][k] * W for k in range(3)]
    dA = R_CAP**2 * np.sin(TH)          # (theta, phi) area element
    F = np.array([trapezoid(trapezoid(Tn[k] * dA, ph, axis=1),
                            th, axis=0)
                  for k in range(3)])
    return F


def bridge_min(psi, u0, theta, sector):
    g = [np.gradient(psi, dx, axis=ax) for ax in range(3)]
    g[2] = g[2] + u0
    amag = np.sqrt(g[0]**2 + g[1]**2 + g[2]**2)
    n = pair_dir(theta)
    tt = np.linspace(-0.48, 0.48, 97)
    pts = tt[None, :] * n[:, None]
    ii = (pts + LBOX) / dx - 0.5
    prof = map_coordinates(amag, ii, order=1, mode='nearest')
    imin = int(np.argmin(prof))
    return {"bridge_min_abs_a_over_gt": float(prof[imin]),
            "bridge_min_t_over_d": float(tt[imin]),
            "bridge_min_P_X": float(f_flux(np.asarray(prof[imin]),
                                           sector))}


def pair_force(theta, u0, sector, psi0=None):
    """Mutual pair force along n_hat: stress through a sphere around
    star B with the B-only solution subtracted (step_59 convention)."""
    n = pair_dir(theta)
    cB = 0.5 * D_PAIR * n
    rho_AB = make_rho(theta)
    rho_B = make_rho(theta, stars=(+1,))
    psi_AB, rAB = solve(u0, sector, rho_AB, psi0=psi0)
    psi_B, rB = solve(u0, sector, rho_B)
    F_AB = sphere_force(psi_AB, u0, cB, sector)
    F_B = sphere_force(psi_B, u0, cB, sector)
    F = F_AB - F_B
    return {"psi_AB": psi_AB, "F_vec": F, "F_par": float(F @ n),
            "res_AB": rAB, "res_B": rB}


def main():
    t0 = time.time()
    out = {"step": "step_75_offaxis_pair_3d",
           "note": ("3D Cartesian solve of div J = rho for a pair at "
                    "angle theta to the ambient gradient u0 z; "
                    "dimensionless g_t = d = 1; parameters matched to "
                    "step_59 (sigma=0.15, u0=0.15, Q=4pi)"),
           "grid": {"N": N, "LBOX_over_d": LBOX, "dx": dx},
           "u0_over_gt": U0, "sigma_src_over_d": SIG_SRC,
           "R_cap_over_d": R_CAP,
           "step59_reference_y_pair_theta0": 0.1524,
           "sectors": {}}

    thetas = (0, 30, 45, 60, 90)
    for sector in ("linear", "baseline", "two_branch"):
        res = {}
        psi_prev = None
        for theta_deg in thetas:
            theta = np.radians(theta_deg)
            print(f"  sector={sector} theta={theta_deg} deg",
                  flush=True)
            pf = pair_force(theta, U0, sector, psi0=psi_prev)
            psi_prev = pf["psi_AB"]
            row = {"theta_deg": float(theta_deg),
                   "F_pair_along_n": pf["F_par"],
                   "F_vec": [float(v) for v in pf["F_vec"]],
                   "resid_AB": pf["res_AB"], "resid_B": pf["res_B"]}
            row.update(bridge_min(pf["psi_AB"], U0, theta, sector))
            res[str(theta_deg)] = row
            print(f"    -> F_par={pf['F_par']:.4f} "
                  f"bridge_min={row['bridge_min_abs_a_over_gt']:.4f}",
                  flush=True)
        out["sectors"][sector] = res

    # y_pair normalised to the linear reference at the same theta
    yp = {}
    for sector in ("baseline", "two_branch"):
        yp[sector] = {
            th: (out["sectors"][sector][th]["F_pair_along_n"]
                 / out["sectors"]["linear"][th]["F_pair_along_n"])
            for th in out["sectors"]["linear"]}
    out["y_pair"] = yp
    out["interpretation"] = (
        "Orientation-dependence of the pair response: y_pair(theta) is "
        "the mutual force along the pair axis normalised to the linear "
        "solve at the same orientation.  theta=0 reproduces the "
        "step_59 axisymmetric geometry (reference y_pair = 0.1524 at "
        "u0=0.15, d=1); theta=90 is the first perpendicular solve.  "
        "Coarse grid (sigma_src ~ 1.2 cells): angular trend, not the "
        "converged asymptote.")

    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "..", "..", "results",
                        "step_75_offaxis_pair_3d.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"wrote {os.path.abspath(path)} "
          f"({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
