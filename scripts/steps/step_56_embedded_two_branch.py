#!/usr/bin/env python3
"""Embedded propagator under the two-branch kinetic sector
(plan T1.3 / T-W3) -- generalises step_32's flux law.

step_32 solves  div J(grad phi) = rho  for a point source embedded in a
uniform background gradient, with baseline J_i = a_i (1 + |a|^2)
(P_X = 1 + u^2 in u = |grad phi|/g_t-scale units).  The two-branch
sector of step_54 (P_X = k sqrt(xi) + 2 xi, xi = u^2/2) generalises to

    J_i = a_i f(|a|),   f(q) = k q / sqrt(2) + q^2  (+ eps_reg floor)

with the same finite-volume discretisation: face fluxes share the
Jacobian stencil and Newton iteration converges.  The baseline control
f(q) = 1 + q^2 is run first and must reproduce step_32's y_far ~ 0.56 /
half-excess ~0.6 r* to validate the port.

Observables extracted for the WB estimator (TEP-WB step_017):
  y_emb(r)  = |grad psi| * r^2        (scalar response vs canonical)
  y_far     = far-field suppression
  alpha_sat under three readings:
    (a) pure propagator:  v = sqrt(1 + 2 b^2 <y_emb>) - 1
    (b) consistent vertex: q^2 = P_X(u0^2/2)^{-2} factored in
    (c) baseline-fixed vertex q^2 = 0.433 for comparability
and the same for the fitted canonical R_s through the mass-convolved
profile (the estimator machinery lives in TEP-WB step_017; here the
profile-level outputs are emitted for it to consume).

Units: g_t = 1, r_* = 1 (a_N = 1/r^2), k = 16.03 (a_eff = a_0).
Degenerate-ellipticity regularisation: f(q) gains a floor eps_reg
(default 1e-8; reported) since P_X -> 0 as q -> 0 on the small-X branch.
"""
import json
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import brentq
from scipy.ndimage import map_coordinates

# ---------------------------------------------------------------------------
# sector definition:  f(q) multiplies a_i in the flux; f'(q) enters the
# anisotropic stiffness  Z_ij = f d_ij + f' a_i a_j / q
# ---------------------------------------------------------------------------
K = 16.03          # two-branch coefficient (a_eff = 4 sqrt(2) g_t / k)
EPS_REG = 1e-8
SQRT2 = np.sqrt(2.0)


def f_flux(q, sector):
    if sector == "baseline":
        return 1.0 + q * q
    if sector == "exp_interp":
        # P_X = 1 - exp(-k sqrt(xi)) + 2 xi,  sqrt(xi) = q/sqrt(2)
        return EPS_REG + (1.0 - np.exp(-K * q / SQRT2)) + q * q
    # two-branch: k |a|/sqrt(2) + |a|^2  ==  k sqrt(xi) + 2 xi
    return EPS_REG + K * q / SQRT2 + q * q


def fp_flux(q, sector):
    """df/dq."""
    if sector == "baseline":
        return 2.0 * q
    if sector == "exp_interp":
        return (K / SQRT2) * np.exp(-K * q / SQRT2) + 2.0 * q
    return K / SQRT2 + 2.0 * q


# --- grid: identical to step_32 ------------------------------------------
NS, NZ = 110, 220
SMAX, ZMAX = 10.0, 10.0
ds = SMAX / NS
dz = 2.0 * ZMAX / NZ
s = (np.arange(NS) + 0.5) * ds
z = (np.arange(NZ) + 0.5) * dz - ZMAX
SS, ZZ = np.meshgrid(s, z, indexing='ij')
RR = np.sqrt(SS ** 2 + ZZ ** 2)

SIG_SRC = 0.15
SRC_TOTAL = 4.0 * np.pi
rho = (SRC_TOTAL / (2.0 * np.pi) ** 1.5 / SIG_SRC ** 3
       * np.exp(-RR ** 2 / (2.0 * SIG_SRC ** 2)))

# source starvation: objects are wells -- the matter source couples
# through A(phi) = e^{-phi_tot}, phi_tot = PHI0 + u0 z + psi
# (ambient level + ambient gradient + perturbation).  Disabled by
# default so the baseline control is unchanged.
STARVE = False
PHI0 = 0.0
STARVE_LAM = 1.0          # continuation ramp on the exponent

# phi-dependent kinetic stiffness:  P(X,phi) = K(phi) * (X-part), so the
# flux law picks up a multiplicative K(phi_tot) per face.  K>1 in deep
# wells (KPHI_C > 0) stiffens the medium where the merged field is deep
# -- the two-centre bridge -- while leaving the far-field plateau
# (phi -> ambient) untouched.  Disabled by default (KPHI_C = 0).
KPHI_C = 0.0

# Picard freezing: when STARVE or KPHI_C is on, the phi that enters the
# source factor and K is taken from PSI_FROZEN (the previous iterate),
# so each inner Newton solve sees a standard spatially-varying-coefficient
# problem and converges cleanly.  solve_picard drives the outer loop.
PSI_FROZEN = None


def _phi_coeff_field(psi):
    return PSI_FROZEN if PSI_FROZEN is not None else psi


# KPHI_PERT: when True, K is built on the perturbation well depth psi
# only (the pair's own field), not the ambient ramp PHI0 + u0 z -- the
# "merged bridge" channel isolated from the box-scale ambient tilt.
KPHI_PERT = False


def kappa(phi):
    return np.exp(np.clip(KPHI_C * phi, -60.0, 60.0))


def source_term(psi, u0):
    if not STARVE:
        return rho
    phi_tot = PHI0 + u0 * ZZ + _phi_coeff_field(psi)
    return rho * np.exp(-np.clip(STARVE_LAM * phi_tot, -60.0, 60.0))

S_FACE = np.arange(NS + 1) * ds


def face_fluxes(psi, u0, sector):
    """Face gradients and fluxes of phi = u0 z + psi (same stencil as
    step_32; only the constitutive law f(q) changes)."""
    d_dz = np.gradient(psi, dz, axis=1)
    d_ds = np.gradient(psi, ds, axis=0)

    a_s_f = (psi[1:, :] - psi[:-1, :]) / ds
    a_z_f = 0.5 * (d_dz[:-1, :] + d_dz[1:, :]) + u0
    am2 = a_s_f ** 2 + a_z_f ** 2
    q = np.sqrt(np.maximum(am2, 1e-30))
    pf = _phi_coeff_field(psi)
    phi_f = (0.5 * (pf[:-1, :] + pf[1:, :]) if KPHI_PERT
             else PHI0 + u0 * z[None, :] + 0.5 * (pf[:-1, :] + pf[1:, :]))
    Ks = kappa(phi_f) if KPHI_C else 1.0
    Js_f = Ks * a_s_f * f_flux(q, sector)
    fp = Ks * fp_flux(q, sector)
    Zss_f = Ks * f_flux(q, sector) + fp * a_s_f ** 2 / q
    Zsz_f = fp * a_s_f * a_z_f / q

    a_z_g = (psi[:, 1:] - psi[:, :-1]) / dz + u0
    a_s_g = 0.5 * (d_ds[:, :-1] + d_ds[:, 1:])
    am2g = a_s_g ** 2 + a_z_g ** 2
    qg = np.sqrt(np.maximum(am2g, 1e-30))
    phi_g = (0.5 * (pf[:, :-1] + pf[:, 1:]) if KPHI_PERT
             else PHI0 + u0 * 0.5 * (z[:-1] + z[1:])[None, :]
             + 0.5 * (pf[:, :-1] + pf[:, 1:]))
    Kz = kappa(phi_g) if KPHI_C else 1.0
    Jz_g = Kz * a_z_g * f_flux(qg, sector)
    fpg = Kz * fp_flux(qg, sector)
    Zzz_g = Kz * f_flux(qg, sector) + fpg * a_z_g ** 2 / qg
    Zsz_g = fpg * a_s_g * a_z_g / qg
    return (Js_f, Zss_f, Zsz_f), (Jz_g, Zzz_g, Zsz_g)


def residual(psi, u0, sector):
    (Js_f, _, _), (Jz_g, _, _) = face_fluxes(psi, u0, sector)
    div_s = np.zeros_like(psi)
    sJp = S_FACE[1:-1][:, None] * Js_f
    div_s[:-1, :] += sJp / (SS[:-1, :] * ds)
    div_s[1:, :] -= sJp / (SS[1:, :] * ds)
    div_z = np.zeros_like(psi)
    div_z[:, :-1] += Jz_g / dz
    div_z[:, 1:] -= Jz_g / dz
    return div_s + div_z - source_term(psi, u0)


def build_jacobian(psi, u0, sector):
    (Js_f, Zss_f, Zsz_f), (Jz_g, Zzz_g, Zsz_g) = face_fluxes(psi, u0, sector)
    N = NS * NZ
    rows, cols, vals = [], [], []

    def add(r, c, v):
        rows.append(np.asarray(r).ravel())
        cols.append(np.asarray(c).ravel())
        vals.append(np.asarray(v).ravel())

    JJ = np.arange(NZ)
    for k in range(NS - 1):
        i, inb = k, k + 1
        sf = S_FACE[k + 1]
        for ci, sgn in ((i, +1.0), (inb, -1.0)):
            c_idx = ci * NZ + JJ
            denom = np.maximum(s[ci], 1e-12) * ds
            coef = sgn * sf * Zss_f[k, :] / (denom * ds)
            add(c_idx, inb * NZ + JJ, coef)
            add(c_idx, i * NZ + JJ, -coef)
            coefz = sgn * sf * Zsz_f[k, :] / (denom * 4 * dz)
            for dj, w in ((1, 1.0), (-1, -1.0)):
                jj = np.clip(JJ + dj, 0, NZ - 1)
                mask = np.ones(NZ, bool)
                mask[-1 if dj > 0 else 0] = False
                add(c_idx[mask], (inb * NZ + jj[mask]),
                    coefz[mask] * w)
                add(c_idx[mask], (i * NZ + jj[mask]),
                    coefz[mask] * w)

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
    if STARVE and PSI_FROZEN is None:
        # exact Newton: d(-rho e^{-lam phi_tot})/d psi = +lam * source
        A = A + sp.diags((STARVE_LAM * source_term(psi, u0)).reshape(-1),
                         offsets=0, format='csr')
    return A


def dirichlet_project(A):
    A = A.tolil()
    bd = ([(NS - 1) * NZ + j for j in range(NZ)]
          + [i * NZ + j for i in range(NS) for j in (0, NZ - 1)])
    for c in bd:
        A.rows[c] = [c]
        A.data[c] = [1.0]
    return A.tocsr()


def solve_embedded(u0, sector, itmax=120, psi0=None):
    psi = np.zeros_like(SS) if psi0 is None else psi0.copy()
    for it in range(itmax):
        R = residual(psi, u0, sector).copy()
        R[NS - 1, :] = psi[NS - 1, :]
        R[:, 0] = psi[:, 0]
        R[:, NZ - 1] = psi[:, NZ - 1]
        rn = np.linalg.norm(R[2:-2, 2:-2])
        if it % 10 == 0:
            print(f"  [{sector}] u0={u0:.3f} it={it} |R|={rn:.3e}",
                  flush=True)
        if rn < 1e-8:
            break
        A = dirichlet_project(build_jacobian(psi, u0, sector))
        try:
            dpsi = spla.spsolve(A, -R.reshape(-1)).reshape(NS, NZ)
        except Exception:
            dpsi = spla.lsmr(A, -R.reshape(-1))[0].reshape(NS, NZ)
        # globalisation: cap the step at the O(1) field scale so a
        # long Newton leg cannot overshoot the basin
        dmax = np.max(np.abs(dpsi))
        if dmax > 0.5:
            dpsi *= 0.5 / dmax
        alpha = 1.0
        accepted = False
        for _ in range(40):
            cand = psi + alpha * dpsi
            rc = residual(cand, u0, sector)
            if np.linalg.norm(rc[2:-2, 2:-2]) < rn:
                psi = cand
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            psi = psi + dpsi * 1e-6
    print(f"  [{sector}] u0={u0:.3f} done it={it} |R|="
          f"{np.linalg.norm(residual(psi, u0, sector)[2:-2, 2:-2]):.3e}",
          flush=True)
    return psi


def solve_picard(u0, sector, outer=40, omega=0.5, tol=1e-6):
    """Outer Picard loop for the phi-coupled terms (STARVE / KPHI_C):
    freeze phi in the source factor and K(phi) on the previous iterate,
    solve the standard Newton problem, damped-update the frozen field."""
    global PSI_FROZEN
    psi = np.zeros_like(SS)
    for it in range(outer):
        PSI_FROZEN = psi.copy()
        psi_new = solve_embedded(u0, sector, psi0=psi)
        delta = np.max(np.abs(psi_new - psi))
        psi = (1.0 - omega) * psi + omega * psi_new
        if delta < tol:
            break
    PSI_FROZEN = None
    print(f"  [{sector}] u0={u0:.3f} picard done outer={it} "
          f"dmax={delta:.3e}", flush=True)
    return psi


def solve_starved(u0, sector, lams=(0.1, 0.25, 0.45, 0.65, 0.8, 0.9,
                                    0.95, 1.0)):
    """Continuation in the starvation exponent: ramp STARVE_LAM with a
    warm-started exact-Newton solve at each stage.  The e^{-phi} source
    coupling is too stiff for Picard or cold Newton at lam = 1."""
    global STARVE_LAM
    psi = np.zeros_like(SS)
    for lam in lams:
        STARVE_LAM = lam
        psi = solve_embedded(u0, sector, psi0=psi)
    STARVE_LAM = 1.0
    return psi


def extract_profiles(psi):
    dpsi_ds = np.gradient(psi, ds, axis=0)
    dpsi_dz = np.gradient(psi, dz, axis=1)
    a_par = np.abs(dpsi_dz[0, :])
    j0 = np.argmin(np.abs(z))
    a_perp = np.abs(dpsi_ds[:, j0])
    return (np.abs(z), a_par), (s.copy(), a_perp)


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


def run_sector(sector, u0_list):
    out = {}
    for u0 in u0_list:
        psi = solve_embedded(u0, sector)
        (zp, ap), (sp_, aq) = extract_profiles(psi)
        radii, ym = orient_mean_profile(psi)
        yfar = float(ym[-6:].mean())
        ih = int(np.argmin(np.abs(ym - 0.5 * yfar)))
        # stiffness at the ambient (for the consistent-vertex reading)
        q0 = abs(u0)
        f0 = f_flux(q0, sector)
        fp0 = fp_flux(q0, sector)
        Zpar = f0 + fp0 * q0          # Z_ij a_i a_j / q along z
        Zperp = f0
        q2_vertex = (1.0 / f0) ** 2
        out[str(u0)] = {
            "u0": u0,
            "y_emb": ym.tolist(),
            "r_over_rstar": radii.tolist(),
            "y_far": yfar,
            "half_excess_r_over_rstar": float(radii[ih]),
            "ambient_PX": float(f0),
            "ambient_Z_parallel": float(Zpar),
            "ambient_Z_perp": float(Zperp),
            "q2_vertex": float(q2_vertex),
            "alpha_pure_propagator": float(np.sqrt(1 + 2 * yfar) - 1),
            "alpha_consistent_vertex":
                float(np.sqrt(1 + 2 * q2_vertex * yfar) - 1),
            "alpha_baseline_vertex_0433":
                float(np.sqrt(1 + 2 * 0.433 * yfar) - 1),
            # directional decomposition: y_par(z) along the ambient axis,
            # y_perp(s) transverse — the transverse profile bounds the
            # perpendicular two-centre pair response (its linearized limit)
            "y_emb_par": (ap * zp**2).tolist(),
            "z_over_rstar": zp.tolist(),
            "y_emb_perp": (aq * sp_**2).tolist(),
            "s_over_rstar": sp_.tolist(),
            "y_far_par": float((ap * zp**2)[zp > 7].mean()),
            "y_far_perp": float((aq * sp_**2)[sp_ > 7].mean()),
        }
    return out


def main():
    out = {"k": K, "eps_reg": EPS_REG, "grid":
           {"NS": NS, "NZ": NZ, "SMAX": SMAX, "ZMAX": ZMAX},
           "note": ("flux law J_i = a_i f(|a|); two-branch "
                    "f = k q/sqrt(2) + q^2; baseline f = 1 + q^2")}

    # sqrt(0.52) plateau-calibrated, 0.570 solar-circle direct, plus the
    # Gate-A discriminating sweep: the two-branch sector reproduces the
    # observed plateau alpha ~ 0.37 at ambient u0 ~ 0.14-0.15
    u0_list = [0.721, 0.570, 0.3, 0.15, 0.12, 1.1]

    print("=== baseline control (must reproduce step_32 ~0.56) ===")
    out["baseline_control"] = run_sector("baseline", [0.721])

    print("=== two-branch ===")
    out["two_branch"] = run_sector("two_branch", u0_list)

    print("=== exp_interp rival ===")
    out["exp_interp"] = run_sector("exp_interp", u0_list)

    # reference values from step_32 for validation
    try:
        s32 = json.loads((Path(__file__).resolve().parents[2] /
                          "results" /
                          "step_32_embedded_propagator.json").read_text())
        out["step32_reference"] = {
            "y_far": s32["orientation_mean_profile"]["y_far_embedded"],
            "half_excess": s32["orientation_mean_profile"]
                              ["half_excess_r_over_rstar_embedded"]}
    except Exception as e:
        out["step32_reference"] = f"unavailable: {e}"

    dest = Path(__file__).resolve().parents[2] / "results" / \
        "step_56_embedded_two_branch.json"
    dest.write_text(json.dumps(out, indent=2))
    print(json.dumps({k: {u: v["y_far"] for u, v in out[k].items()}
                      for k in ("baseline_control", "two_branch")},
                     indent=2))
    print("wrote", dest)


if __name__ == "__main__":
    main()
