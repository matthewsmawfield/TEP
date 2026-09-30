#!/usr/bin/env python3
"""Interacting two-centre solve (plan T1.3 follow-up): the pair is two
sources in the nonlinear medium, not one source in a fixed ambient.

Each star's scalar well sits inside the companion's field, so the
pair's effective scalar charge is modulated by the TOTAL local field
u_tot = |u0 zhat + grad psi_1 + grad psi_2| -- the monopole embedded
solve (step_32/56) misses this nonlinear superposition.

Geometry: Gaussian sources of total charge 4 pi each at z = +/-d/2 on
the axis -- pair axis PARALLEL to the ambient gradient, the only
axisymmetric configuration on the cylindrical grid (single-orientation
probe, flagged).

Mutual force by momentum conservation -- field stress through a sphere
enclosing only B:

    T_ij = f(q) a_i a_j - delta_ij P(q),   P(q) = -∫_0^q f(q') q' dq'

Pair-interaction part isolated by subtracting the B-only solution's
stress through the identical surface (removes the ambient body force
Q_B u0 and B's self-terms):

    F_pair(d) = F_stress[A+B+amb] - F_stress[B+amb]

Normalised response y_pair(d) = F_pair(d) / F_pair^linear(d), with the
linear reference from the f = 1 sector in identical geometry, so
y -> 1 is the unscreened limit.  Sweeps d over the transition region at
the ambients that mattered in step_019.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.integrate import trapezoid

sys.path.insert(0, str(Path(__file__).resolve().parent))
import step_56_embedded_two_branch as m

SIG = m.SIG_SRC
CHARGE = m.SRC_TOTAL          # 4 pi per star


def src_rho(centres):
    norm = CHARGE / (2.0 * np.pi) ** 1.5 / SIG ** 3
    rho = np.zeros_like(m.SS)
    for zc in centres:
        rho += norm * np.exp(-(m.SS ** 2 + (m.ZZ - zc) ** 2)
                             / (2 * SIG ** 2))
    return rho


def W_of_q(q, sector):
    """Flux potential W(q) = ∫_0^q f(q') q' dq' -- the Euclidean
    functional the solve extremises (dW/da_i = f a_i).  The
    momentum-flux tensor is T_ij = f a_i a_j - delta_ij W."""
    q = np.asarray(q, dtype=float)
    out = np.empty_like(q)
    qcap = max(1.0, float(np.abs(q).max()) * 1.05)
    qs = np.linspace(0.0, qcap, 2001)         # tabulate once
    ftab = m.f_flux(qs, sector)
    wtab = np.concatenate([[0.0],
                           np.cumsum(0.5 * (ftab[:-1] * qs[:-1]
                                            + ftab[1:] * qs[1:])
                                     * np.diff(qs))])
    flat = np.abs(q).ravel()
    out.ravel()[:] = np.interp(flat, qs, wtab)
    return out


def stress_z(psi, u0, zB, r_c, sector, nth=361):
    """Z-momentum flux through sphere centred on-axis at zB."""
    th = np.linspace(1e-6, np.pi - 1e-6, nth)
    sp_ = r_c * np.sin(th)
    zp_ = zB + r_c * np.cos(th)
    ii = sp_ / m.ds - 0.5
    jj = (zp_ + m.ZMAX) / m.dz - 0.5
    a_s = map_coordinates(np.gradient(psi, m.ds, axis=0), [ii, jj],
                          order=1, mode='nearest')
    a_z = map_coordinates(np.gradient(psi, m.dz, axis=1), [ii, jj],
                          order=1, mode='nearest') + u0
    q = np.sqrt(a_s ** 2 + a_z ** 2)
    fq = m.f_flux(q, sector)
    W = W_of_q(q, sector)
    if m.KPHI_C:
        # J_i = K(phi) f a_i  =>  T_ij = K(phi) (f a_i a_j - d_ij W)
        phipt = m.PHI0 + u0 * zp_ + map_coordinates(
            psi, [ii, jj], order=1, mode='nearest')
        Kp = m.kappa(phipt)
        fq = Kp * fq
        W = Kp * W
    T_zz = fq * a_z * a_z - W
    T_zs = fq * a_z * a_s
    dA = 2.0 * np.pi * r_c ** 2 * np.sin(th)
    return float(trapezoid((T_zz * np.cos(th) + T_zs * np.sin(th))
                           * dA, th))


def direct_force(psi, rho_B, u0):
    """F = |int rho_B w d_z psi dV| on a GIVEN field (used for the
    linear reference).  For the mutual pair force under K != 1 the
    B-only field must be subtracted first -- see pair_force_direct;
    B's own modulated well refracts the ambient and exerts a net
    self-force that does not cancel.  w = e^{-phi_tot} under STARVE,
    else 1 -- alpha = -1 folded into the magnitude."""
    dz_psi = np.gradient(psi, m.dz, axis=1)
    if m.STARVE:
        phi_tot = m.PHI0 + u0 * m.ZZ + psi
        w = np.exp(-np.clip(phi_tot, -60.0, 60.0))
    else:
        w = 1.0
    return float(np.abs((rho_B * w * dz_psi * m.SS
                         * (2.0 * np.pi)).sum()) * m.ds * m.dz)


def solve_with(rho_arr, u0, sector):
    m.rho = rho_arr
    if m.STARVE and not m.KPHI_C:
        return m.solve_starved(u0, sector)
    if m.KPHI_C:
        return m.solve_picard(u0, sector)
    return m.solve_embedded(u0, sector)


def solve_linear(rho_arr, u0):
    """Canonical Poisson reference (f = 1) via the same machinery.
    Saves and restores ALL coupling state -- KPHI_C / KPHI_PERT /
    STARVE -- so a reference computed after a coupled run is not
    silently contaminated by leftover modulation flags."""
    orig_f, orig_fp = m.f_flux, m.fp_flux
    orig_kc, orig_kp = m.KPHI_C, m.KPHI_PERT
    orig_starve, orig_phi0 = m.STARVE, m.PHI0
    m.f_flux = lambda q, s: np.ones_like(np.asarray(q, float))
    m.fp_flux = lambda q, s: np.zeros_like(np.asarray(q, float))
    m.KPHI_C, m.KPHI_PERT = 0.0, False
    m.STARVE = False
    m.rho = rho_arr
    psi = m.solve_embedded(u0, "baseline")
    m.f_flux, m.fp_flux = orig_f, orig_fp
    m.KPHI_C, m.KPHI_PERT = orig_kc, orig_kp
    m.STARVE, m.PHI0 = orig_starve, orig_phi0
    return psi


def pair_force(d, u0, sector):
    rho_AB = src_rho((-d / 2.0, +d / 2.0))
    rho_B = src_rho((+d / 2.0,))
    r_c = max(d / 3.0, 4.0 * SIG)
    psi_AB = solve_with(rho_AB, u0, sector)
    F_AB = stress_z(psi_AB, u0, d / 2.0, r_c, sector)
    psi_B = solve_with(rho_B, u0, sector)
    F_B = stress_z(psi_B, u0, d / 2.0, r_c, sector)
    return F_AB - F_B


def pair_force_direct(d, u0, sector):
    """Pair force by direct matter integral.  With K(phi) != 1, B's own
    modulated well refracts the ambient flux and exerts a net force on
    B itself (dielectric self-term), so the mutual part requires the
    B-only subtraction: F = |int rho_B w d_z(psi_AB - psi_B) dV|."""
    rho_AB = src_rho((-d / 2.0, +d / 2.0))
    rho_B = src_rho((+d / 2.0,))
    psi_AB = solve_with(rho_AB, u0, sector)
    psi_B = solve_with(rho_B, u0, sector)
    dz_diff = np.gradient(psi_AB - psi_B, m.dz, axis=1)
    if m.STARVE:
        phi_tot = m.PHI0 + u0 * m.ZZ + psi_AB
        w = np.exp(-np.clip(phi_tot, -60.0, 60.0))
    else:
        w = 1.0
    return float(np.abs((rho_B * w * dz_diff * m.SS
                         * (2.0 * np.pi)).sum()) * m.ds * m.dz)


def main():
    out = {"geometry": "pair axis || ambient gradient (axisymmetric)",
           "charge_per_star": CHARGE, "sigma_src": SIG,
           "note": ("y_pair(d) = F_pair/F_pair^linear; interaction "
                    "isolated by B-only stress subtraction; r_c = "
                    "max(d/3, 4 sigma_src)")}

    u0_list = [0.15, 0.57, 0.72]
    d_list = [1.0, 1.5, 2.0, 3.0, 4.0, 5.0]

    # linear reference forces (sector-independent)
    print("=== linear reference ===")
    ref = {}
    for d in d_list:
        rho_AB = src_rho((-d / 2.0, +d / 2.0))
        rho_B = src_rho((+d / 2.0,))
        r_c = max(d / 3.0, 4.0 * SIG)
        psi = solve_linear(rho_AB, 0.0)      # u0 irrelevant: uniform
        psiB = solve_linear(rho_B, 0.0)      # field integrates out
        F = (stress_z_linear(psi, 0.0, d / 2.0, r_c)
             - stress_z_linear(psiB, 0.0, d / 2.0, r_c))
        ref[d] = F
        print(f"  d={d}: F_lin={F:.5f}  (1/d^2={1/d**2:.5f})",
              flush=True)

    # linear reference via the direct matter integral (same convention
    # as the direct pair force: rho_B weighted total gradient)
    ref_d = {}
    for d in d_list:
        rho_AB = src_rho((-d / 2.0, +d / 2.0))
        rho_B = src_rho((+d / 2.0,))
        psi = solve_linear(rho_AB, 0.0)
        ref_d[d] = direct_force(psi, rho_B, 0.0)
        print(f"  d={d}: F_lin_direct={ref_d[d]:.5f}", flush=True)

    for sector in ("two_branch", "baseline"):
        out[sector] = {}
        for u0 in u0_list:
            rows = []
            for d in d_list:
                F = pair_force(d, u0, sector)
                Fd = pair_force_direct(d, u0, sector)
                rows.append({"d_over_rstar": d,
                             "F_pair": F, "F_linear": ref[d],
                             "F_pair_direct": Fd,
                             "F_linear_direct": ref_d[d],
                             "y_pair": F / ref[d] if ref[d] else None,
                             "y_pair_direct":
                                 Fd / ref_d[d] if ref_d[d] else None})
                print(f"  [{sector}] u0={u0} d={d}: "
                      f"y_pair={rows[-1]['y_pair']:.4f} "
                      f"y_dir={rows[-1]['y_pair_direct']:.4f}",
                      flush=True)
            out[sector][str(u0)] = rows

    # --- starved-source variant: rho -> rho e^{-phi_tot}, the objects
    # are wells, not bare charges.  phi_tot = PHI0 + u0 z + psi: ambient
    # level (phi0), ambient gradient (u0 z) and the pair's own field
    # (psi, which includes the companion's well at each star's core) all
    # suppress the effective scalar charge.  Sweep the (u0, phi0) plane
    # for the two-branch sector only. ---
    out["two_branch_starved"] = {}
    m.STARVE = True
    for u0 in (0.15, 0.57):
        for phi0 in (1.0,):
            m.PHI0 = phi0
            rows = []
            for d in d_list:
                F = pair_force_direct(d, u0, "two_branch")
                rows.append({"d_over_rstar": d, "F_pair": F,
                             "F_linear": ref_d[d],
                             "y_pair": F / ref_d[d] if ref_d[d] else None})
                print(f"  [tb starved] u0={u0} phi0={phi0} d={d}: "
                      f"y_pair={rows[-1]['y_pair']:.4f}", flush=True)
            out["two_branch_starved"][f"{u0}|{phi0}"] = {
                "u0": u0, "phi0": phi0, "rows": rows,
                "solver_note": ("continuation-Newton; residual floor "
                                "~1e2 where the e^{-lam phi} shell is "
                                "sub-cell -- partial convergence only")}
    m.STARVE = False
    m.PHI0 = 0.0

    # --- phi-dependent stiffness K(phi) = e^{c phi_tot}: the kinetic
    # sector's own field-level dependence.  c>0 stiffens the medium in
    # deep wells -- the pair's merged bridge is suppressed harder while
    # the far-field plateau (phi -> ambient) is untouched, so a sharp
    # transition is possible without touching the shallow-phi galactic
    # tail.  Quasi-Newton Jacobian neglects dK/dpsi (flagged). ---
    out["two_branch_kphi"] = {}
    for c in (0.5, 1.0, 2.0, -0.5):
        m.KPHI_C = c
        for u0 in (0.15, 0.57):
            rows = []
            for d in d_list:
                F = pair_force_direct(d, u0, "two_branch")
                rows.append({"d_over_rstar": d, "F_pair": F,
                             "F_linear": ref_d[d],
                             "y_pair": F / ref_d[d] if ref_d[d] else None})
                print(f"  [tb Kphi c={c}] u0={u0} d={d}: "
                      f"y_pair={rows[-1]['y_pair']:.4f}", flush=True)
            out["two_branch_kphi"][f"{c}|{u0}"] = {
                "kphi_c": c, "u0": u0, "rows": rows}
    m.KPHI_C = 0.0

    # --- bridge-depth variant: K(psi) responds to the pair's own merged
    # well, not the ambient ramp u0 z -- isolates the "molecules mixing"
    # channel the phi_tot version cannot separate from box-scale tilt. ---
    out["two_branch_kphi_pert"] = {}
    m.KPHI_PERT = True
    for c in (0.5, 1.0, 2.0, -0.5, -1.0):
        m.KPHI_C = c
        for u0 in (0.15, 0.3, 0.57):
            rows = []
            for d in d_list:
                F = pair_force_direct(d, u0, "two_branch")
                rows.append({"d_over_rstar": d, "F_pair": F,
                             "F_linear": ref_d[d],
                             "y_pair": F / ref_d[d] if ref_d[d] else None})
                print(f"  [tb Kpert c={c}] u0={u0} d={d}: "
                      f"y_pair={rows[-1]['y_pair']:.4f}", flush=True)
            out["two_branch_kphi_pert"][f"{c}|{u0}"] = {
                "kphi_c": c, "u0": u0, "rows": rows}
    m.KPHI_C = 0.0
    m.KPHI_PERT = False

    dest = Path(__file__).resolve().parents[2] / "results" / \
        "step_59_two_center.json"
    dest.write_text(json.dumps(out, indent=2))
    print("wrote", dest)


def stress_z_linear(psi, u0, zB, r_c, nth=361):
    """Momentum flux with f = 1 (canonical reference)."""
    th = np.linspace(1e-6, np.pi - 1e-6, nth)
    sp_ = r_c * np.sin(th)
    zp_ = zB + r_c * np.cos(th)
    ii = sp_ / m.ds - 0.5
    jj = (zp_ + m.ZMAX) / m.dz - 0.5
    a_s = map_coordinates(np.gradient(psi, m.ds, axis=0), [ii, jj],
                          order=1, mode='nearest')
    a_z = map_coordinates(np.gradient(psi, m.dz, axis=1), [ii, jj],
                          order=1, mode='nearest') + u0
    T_zz = a_z * a_z - 0.5 * (a_s ** 2 + a_z ** 2)
    T_zs = a_z * a_s
    dA = 2.0 * np.pi * r_c ** 2 * np.sin(th)
    return float(trapezoid((T_zz * np.cos(th) + T_zs * np.sin(th))
                           * dA, th))


if __name__ == "__main__":
    main()
