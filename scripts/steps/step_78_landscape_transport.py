#!/usr/bin/env python3
"""Sector-resolved drift structure of the evolving landscape: does the
noncanonical kinetic screen pin matter-hosting interiors while the
shallow ambient rolls?

Motivation
----------
step_76 measured the clock-channel transmission T = du_int/du_amb ~ 1 on
the leaking branch: interior field values ride ambient excursions
additively.  step_53 places real bodies on that branch.  If the ambient
drift itself were coherent across the landscape -- every cell rolling at
the same rate -- then interior clocks would track the cosmological roll,
G_loc inside sources would drift as (1+z)^{-2}, and the SN Ia luminosity
bound of step_73 (eps <= 0.05-0.10; eps=1 predicts -2.26 mag at z=1,
excluded ~15-23x) would falsify the tracking picture outright.  The only
architecture consistent with (i) a secular redshift trend, (ii) LLR
|Gdot/G| < 4e-13/yr, and (iii) the SN pinning bound is differential
drift: matter-hosting wells pinned (or weakly drifting) while the
shallow landscape sector rolls, with the redshift carried by transport
through the landscape rather than by interior clock drift.

The unscreened evolution of step_51 fails that precondition: the solved
canonical EOM rolls coherently (du_host/dt ~ du_void/dt to <0.2%), the
u-contrast flattens, and any interior drift would transmit at T ~ 1.
What was never tested is whether the corpus's own kinetic completion --
the same flux law J_i = a_i f(|a|) used in steps 56/59/75/76, with
f(q) = eps_reg + k q / sqrt(2) + q^2 (k = 16.03) on the two-branch sector
and f(q) = 1 + q^2 on the baseline -- breaks the coherence: in screened
(steep-gradient) regions the effective elliptic stiffness is f >> 1,
suppressing the spatial coupling that drags interiors along with the
ambient, while the local potential drive V'(u) is itself weaker at
deeper u on the reconstructed branch (V' ~ e^{-u}).

What this step computes
-----------------------
1. A nested landscape (well lattice + broad envelope, same construction
   as step_51) is solved on-shell under the SCREENED elliptic operator
       div( f(|grad u|) grad u ) = V'(u) - rho + lam_c
   via Picard iteration (linear screened Poisson solve per iterate,
   bordered gauge mean(u)=0).
2. Constrained evolution identical in form to step_51 variant B/D:
   u_t = N pi,  pi_t = N (div(f grad u) - V' + rho),  rho_t = -rho u_t,
   with psi and the step_21 maximal-slice lapse re-solved on the evolved
   sources.  Three kinetic branches are compared on identical initial
   data:
       canonical  f = 1                (reproduces step_51 coherence)
       baseline   f = 1 + q^2
       two-branch f = eps + k q/sqrt2 + q^2
   at a bracket of kinetic scales lam4 = Lambda_X^4 in slice units,
   including a self-consistent iteration lam4 = (slice drift rate)^2.
3. Sector decomposition: core (deepest u decile on hosts), flank,
   envelope-interior ambient, exterior ambient, void (lowest u decile).
   Reported per checkpoint: drift rate du/dt = <N pi> per sector, the
   u-contrast between sectors, and the pinning ratio
       R_pin = |du_core/dt| / |du_void/dt|.
   R_pin << 1 sustained = the transport architecture's precondition;
   R_pin ~ 1 = coherent roll (interior tracking, eps -> 1, SN-tension).
4. Transport bookkeeping: with the lapse re-normalised at a fixed
   observer cell (a well core), the relative-lapse drift rate
   <d ln N/dt> per sector and the path-averaged rate
   d ln omega_path/dt = <d ln N/dt - 2 d ln psi/dt>
   give the in-flight transport coefficient alongside the endpoint
   factor d u_host/dt.  The decomposition answers which channel carries
   the secular trend.

This is a diagnostic of drift *structure*, not a calibrated redshift
prediction: the slice is a periodic cell at toy parameters, drift rates
are in slice time units, and the self-consistent lam4 rescales with the
measured roll rate.
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
from step_51_onshell_landscape_evolution import (  # noqa: E402
    conformal_grad2, conformal_lap, dV_du, d2V_du2,
    maximal_slice_lapse, solve_landscape_warm,
)

OUT_PATH = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "..", "results", "step_78_landscape_transport.json"))

K_TWOB = 16.03
EPS_REG = 1e-8


# ----------------------------------------------------------------------
# Screened flux law f(q) and its derivative (for the Picard stiffness the
# exact Jacobian needs f + f'/q a a; Picard iterates the isotropic
# leading part, which converges for these mildly nonlinear profiles).
# ----------------------------------------------------------------------
def flux_f(q, branch):
    if branch == "canonical":
        return np.ones_like(q)
    if branch == "baseline":
        return 1.0 + q ** 2
    return EPS_REG + K_TWOB * q / np.sqrt(2.0) + q ** 2


def flux_fprime(q, branch):
    """f'(q) with the same regularisation floor as the flux solve."""
    qs = np.maximum(q, 1e-6)
    if branch == "canonical":
        return np.zeros_like(q)
    if branch == "baseline":
        return 2.0 * qs
    return K_TWOB / np.sqrt(2.0) + 2.0 * qs


def screened_div(fld, fvals, dx):
    """div( f grad fld ) by centred faces, periodic."""
    out = np.zeros_like(fld)
    for ax in range(3):
        d_fwd = (np.roll(fld, -1, axis=ax) - fld) / dx
        d_bwd = (fld - np.roll(fld, 1, axis=ax)) / dx
        f_fwd = 0.5 * (fvals + np.roll(fvals, -1, axis=ax))
        f_bwd = 0.5 * (fvals + np.roll(fvals, 1, axis=ax))
        out += (f_fwd * d_fwd - f_bwd * d_bwd) / dx
    return out


def cell_grad(u, dx):
    """Cell-centred centred-difference gradient components."""
    return [(np.roll(u, -1, axis=ax) - np.roll(u, 1, axis=ax))
            / (2.0 * dx) for ax in range(3)]


def div_vec(vec, dx):
    """Divergence of a cell-centred vector by centred differences."""
    return sum((np.roll(vec[ax], -1, axis=ax)
                - np.roll(vec[ax], 1, axis=ax)) / (2.0 * dx)
               for ax in range(3))


def jac_matvec_factory(u, dx, branch, lam4, V0, rhob_half, u_s):
    """Jacobian-vector action of F(u) = div(f(q) grad u) - V'(u).

    J v = div( f grad v ) + div( (f'/q)(a.grad v) a ) - V'' v
    evaluated at the cell-centred gradient a = grad u.  The product is
    assembled from stencil operations -- the sparse matrix is never
    formed; GMRES/LGMRES consumes the action directly.
    """
    a = cell_grad(u, dx)
    q = np.sqrt(sum(x ** 2 for x in a)) / np.sqrt(lam4)
    qsafe = np.maximum(q, 1e-8)
    fv = flux_f(qsafe, branch)
    fp = flux_fprime(qsafe, branch)
    # delta f = f'(q) (a.delta a)/(q lam4):  ratio = f'/(q lam4)
    ratio = fp / (qsafe * lam4)
    Vpp = d2V_du2(u, V0, rhob_half, u_s)

    def mv(v_flat):
        v = v_flat.reshape(u.shape)
        t1 = div_vec([fv * g for g in cell_grad(v, dx)], dx)
        av = sum(a[i] * g for i, g in enumerate(cell_grad(v, dx)))
        t2 = div_vec([ratio * av * a[i] for i in range(3)], dx)
        return (t1 + t2 - Vpp * v).ravel()
    return mv


def solve_scalar_screened(rho, V0, rhob_half, u_s, dx, branch, lam4,
                          u0=None, u_bar=0.0, tol=1e-9, max_newton=30):
    """Newton-Krylov solve of div(f grad u) - V' + rho + lam_c = 0 with
    bordered gauge mean(u) = u_bar.  The constant null direction of the
    periodic divergence is removed by the gauge constraint."""
    shape = rho.shape
    u = np.zeros(shape) if u0 is None else u0.copy()
    lam_c = 0.0
    history = []

    def F(u, lam_c):
        a = cell_grad(u, dx)
        q = np.sqrt(sum(x ** 2 for x in a)) / np.sqrt(lam4)
        fv = flux_f(np.maximum(q, 1e-8), branch)
        return (div_vec([fv * g for g in a], dx)
                - dV_du(u, V0, rhob_half, u_s) + rho + lam_c)

    mv = jac_matvec_factory(u, dx, branch, lam4, V0, rhob_half, u_s)

    for _ in range(max_newton):
        r = F(u, lam_c)
        res = float(np.abs(r).max())
        history.append(res)
        if res < tol:
            break
        mv = jac_matvec_factory(u, dx, branch, lam4,
                                V0, rhob_half, u_s)
        J = spla.LinearOperator((rho.size, rho.size), matvec=mv)
        du, info = spla.lgmres(J, -r.ravel(), rtol=1e-8,
                               maxiter=200, atol=0.0)
        du = du.reshape(shape)
        # bordered gauge: pick dlam so the update fixes the mean
        # J v = const -> solve J d2 = 1 for the gauge direction
        d2, _ = spla.lgmres(J, np.ones(rho.size), rtol=1e-8,
                            maxiter=200, atol=0.0)
        d2 = d2.reshape(shape)
        m2 = float(d2.mean())
        if abs(m2) < 1e-12:
            break
        dlam = (float((u + du).mean()) - u_bar) / m2
        du = du - dlam * d2
        alpha = 1.0
        for _ in range(16):
            un, ln = u + alpha * du, lam_c + alpha * dlam
            if np.abs(F(un, ln)).max() < res or alpha < 1e-5:
                break
            alpha *= 0.5
        u, lam_c = un, ln
    return dict(u=u, lam_c=lam_c,
                converged=bool(history and history[-1] < tol),
                history=history)


def flat_grad2(u, dx):
    fh = fft.fftn(u)
    g2 = np.zeros_like(u)
    k = 2.0 * np.pi * fft.fftfreq(u.shape[0], d=dx)
    for ax in range(3):
        kk = np.reshape(k, [1] * ax + [-1] + [1] * (2 - ax))
        g2 += np.real(fft.ifftn(1j * kk * fh)) ** 2
    return g2


def screened_eom_residual(u, rho, V0, rhob_half, u_s, dx, branch, lam4):
    """Same cell-centred convention as the Newton solve."""
    a = cell_grad(u, dx)
    q = np.sqrt(sum(x ** 2 for x in a)) / np.sqrt(lam4)
    fv = flux_f(np.maximum(q, 1e-8), branch)
    return (div_vec([fv * g for g in a], dx)
            - dV_du(u, V0, rhob_half, u_s) + rho)


def evolve(u0, rho0, psi, lam, V0, rhob_half, u_s, LAP, k2, axes_k,
           dx, branch, lam4, sectors, obs_mask,
           dt=2.0e-4, nsteps=6000, record_every=150, resolve_every=15,
           dt_snapshot=None):
    """Constrained evolution with the screened EOM; records full fields."""
    u, rho = u0.copy(), rho0.copy()
    pi = np.zeros_like(u)
    t = 0.0
    rows = []
    fields = []          # (t, u, N_norm_obs, psi) for transport integrals
    Nlap = np.ones_like(u)

    def N_norm(N):
        v = N / max(float(N[obs_mask].mean()), 1e-30)
        return v

    for step in range(1, nsteps + 1):
        if step % resolve_every == 0 or step == 1:
            V = landscape_potential(u, V0, rhob_half, u_s)
            gu2 = conformal_grad2(u, psi, k2, axes_k)
            S0 = rho + V + 0.5 * pi ** 2
            sol = solve_landscape_warm(S0, gu2, LAP, psi0=psi, lam0=lam,
                                       tol=1e-8, max_newton=8)
            if sol["converged"]:
                psi, lam = sol["psi"], sol["lam"]
            Nlap, mu_star, _ = maximal_slice_lapse(rho, gu2, V, lam, LAP)

        eom = screened_eom_residual(u, rho, V0, rhob_half, u_s, dx,
                                    branch, lam4)
        dpi = Nlap * eom
        du = Nlap * pi
        u2 = u + dt * du
        pi2 = pi + dt * dpi
        rho2 = np.maximum(rho + dt * (-rho * du), 0.0)
        eom2 = screened_eom_residual(u2, rho2, V0, rhob_half, u_s, dx,
                                     branch, lam4)
        u = u + 0.5 * dt * (du + Nlap * pi2)
        pi = pi + 0.5 * dt * (dpi + Nlap * eom2)
        rho = np.maximum(rho + 0.5 * dt * (-rho * du - rho2 * Nlap * pi2),
                         0.0)
        t += dt

        if step % record_every == 0 or step == 1:
            Nn = N_norm(Nlap)
            dudt = Nlap * pi
            row = {"t": float(t), "mu_star": float(mu_star)}
            for name, m in sectors.items():
                row[f"du_{name}"] = float(dudt[m].mean())
                row[f"u_{name}"] = float(u[m].mean())
                row[f"dlnN_{name}"] = float(Nn[m].mean())
                row[f"psi_{name}"] = float(psi[m].mean())
            row["u_min"], row["u_max"] = float(u.min()), float(u.max())
            row["pi_max_abs"] = float(np.abs(pi).max())
            eom_r = screened_eom_residual(u, rho, V0, rhob_half, u_s, dx,
                                          branch, lam4)
            row["eom_inhom_rms"] = float(np.sqrt(np.mean(
                (eom_r - eom_r.mean()) ** 2)))
            rows.append(row)
            fields.append((float(t), u.copy(), Nn.copy(), psi.copy()))
            if not np.isfinite(u).all():
                break
            if abs(row["du_void"]) > 0 and len(rows) > 3:
                pass
    return rows, fields


def transport_integral(fields, dx, n_rays=24):
    """Path integral of the metric-drift transport along null rays.

    d ln omega_g / dt = Ndot/N - 2 psidot/psi evaluated at the moving
    photon position, accumulated between recorded slices.  Rays sample
    the periodic box uniformly (light crossing time ~ L/mean N); the
    reported rate is per slice-time, i.e. the transport coefficient of
    the landscape as a medium.  Endpoint drift is reported separately.
    """
    if len(fields) < 3:
        return None
    rng = np.random.default_rng(0)
    L = fields[0][1].shape[0] * dx
    rays = rng.uniform(0, L, (n_rays, 3))
    dirs = rng.normal(size=(n_rays, 3))
    dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    path_rate = np.zeros(n_rays)
    endpoint = np.zeros(n_rays)
    for j in range(1, len(fields)):
        t0, u0, N0, p0 = fields[j - 1]
        t1, u1, N1, p1 = fields[j]
        dt_ep = t1 - t0
        if dt_ep <= 0:
            continue
        # advance rays one epoch at light speed (mean N ~ 1 normalised)
        rays = (rays + dirs * dt_ep) % L

        def sample(F, pos):
            ij = np.floor(pos / dx).astype(int) % F.shape[0]
            return F[ij[:, 0], ij[:, 1], ij[:, 2]]

        dlnN = np.log(np.maximum(sample(N1, rays), 1e-30)
                      / np.maximum(sample(N0, rays), 1e-30)) / dt_ep
        dlnp = np.log(np.maximum(sample(p1, rays), 1e-30)
                      / np.maximum(sample(p0, rays), 1e-30)) / dt_ep
        path_rate += (dlnN - 2.0 * dlnp)
        endpoint += (sample(u1, rays) - sample(u0, rays)) / dt_ep
    n_ep = len(fields) - 1
    return {"path_rate_per_ray": path_rate.tolist(),
            "path_rate_mean": float(path_rate.mean() / n_ep),
            "endpoint_u_drift_mean": float(endpoint.mean() / n_ep)}


def main():
    L, n_w = 1.0, 4
    N = 24
    sigma = L / 8.0
    u_well_probe, u_void_probe = 0.4, -1.4
    rho0, rho_ambient = 0.5, 0.05
    V0, rhob_half, u_s = 60.0, 30.0, 15.0

    dx = L / N
    LAP = fd_laplacian(N, dx)
    k = 2.0 * np.pi * fft.fftfreq(N, d=dx)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    k2 = KX ** 2 + KY ** 2 + KZ ** 2
    axes_k = (KX, KY, KZ)

    ld = build_landscape(N, L, n_w, sigma, u_well_probe, u_void_probe,
                         rho0, rho_ambient, V0, rhob_half, u_s)
    rho0_f = ld["rho_m"].copy()
    w = ld["wells"]
    host = w >= np.quantile(w, 0.80)
    void = w <= np.quantile(w, 0.05)

    # Nested envelope (same construction as step_51 variant D)
    x = np.arange(N) * dx
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    env = np.exp(-((X - 0.5 * L) ** 2 + (Y - 0.5 * L) ** 2
                   + (Z - 0.5 * L) ** 2) / (2.0 * (L / 3.2) ** 2))
    env_amp = 0.35 * rho0
    rho_nested = rho0_f + env_amp * env
    env_thr = np.median(env[host])
    sectors = {
        "core": host & (env >= env_thr),
        "host_shallow": host & (env < env_thr),
        "amb_in": ~host & (env >= env_thr),
        "amb_out": ~host & (env < env_thr),
        "void": void,
    }
    obs_mask = sectors["core"]  # observer convention: N = 1 on well cores

    results = {"step": "step_78_landscape_transport",
               "purpose": ("Sector-resolved drift structure under the "
                           "screened kinetic EOM: test whether matter-"
                           "hosting interiors pin while the ambient "
                           "rolls -- the precondition for carrying the "
                           "redshift trend in transport rather than in "
                           "interior clock drift (SN-bound consistency)."),
               "grid": {"N": N, "L": L, "sigma": sigma,
                        "nested_envelope_amp": env_amp},
               "kinetic_branches": {}}

    # lam4 bracket: the kinetic scale is self-consistent on the slice's
    # own drift rate (Lambda_X^4 = H_drift^2); scan the dependence.
    lam4_scan = [1e4, 1e2, 1.0, 1e-2]
    for lam4 in lam4_scan:
        for branch in ("canonical", "baseline", "two_branch"):
            tag = f"{branch}_lam4_{lam4:g}"
            print(f"[{tag}] on-shell screened solve")
            s = solve_scalar_screened(rho_nested, V0, rhob_half, u_s, dx,
                                      branch, lam4, u_bar=0.0,
                                      tol=1e-8, max_newton=30)
            u = s["u"]
            print(f"    conv={s['converged']} it={len(s['history'])}"
                  f" lam_c={s['lam_c']:+.4f}"
                  f" u=[{u.min():+.3f},{u.max():+.3f}]")

            gu2 = conformal_grad2(u, np.ones_like(u), k2, axes_k)
            V = landscape_potential(u, V0, rhob_half, u_s)
            sol = solve_landscape_warm(rho_nested + V, gu2, LAP,
                                       tol=1e-9, max_newton=20)
            psi0, lam0 = sol["psi"], sol["lam"]

            rows, fields = evolve(u, rho_nested, psi0, lam0,
                                  V0, rhob_half, u_s, LAP, k2, axes_k,
                                  dx, branch, lam4, sectors, obs_mask,
                                  nsteps=4500, record_every=150,
                                  resolve_every=15)
            if not rows:
                results["kinetic_branches"][tag] = {"failed": True}
                continue
            last = rows[-1]
            # drift coherence: ratio of sector drifts late-time
            core = abs(last["du_core"])
            vdr = abs(last["du_void"])
            amb = abs(last["du_amb_out"])
            pin_core = core / vdr if vdr else None
            pin_amb = core / amb if amb else None
            trans = transport_integral(fields, dx)
            results["kinetic_branches"][tag] = {
                "on_shell_converged": s["converged"],
                "lam_c": s["lam_c"],
                "lam_star": float(lam0),
                "u_range_init": [float(u.min()), float(u.max())],
                "n_checkpoints": len(rows),
                "final": last,
                "drift_history": [
                    {k: r[k] for k in r if k.startswith(("du_", "t"))}
                    for r in rows],
                "pinning_ratio_core_over_void": pin_core,
                "pinning_ratio_core_over_amb_out": pin_amb,
                "u_contrast_final": last["u_max"] - last["u_min"],
                "transport": trans,
            }
            print(f"    final drifts: core={last['du_core']:+.3f}"
                  f" amb_out={last['du_amb_out']:+.3f}"
                  f" void={last['du_void']:+.3f}"
                  f"  pin(core/void)={pin_core:.3f}")

    # verdict synthesis
    verdict = []
    for tag, r in results["kinetic_branches"].items():
        if r.get("failed"):
            continue
        pcv = r["pinning_ratio_core_over_void"]
        verdict.append({"branch": tag, "pin_core_over_void": pcv})
    results["verdict_table"] = verdict
    results["interpretation"] = (
        "pin_core_over_void ~ 1 => coherent roll: interiors track the "
        "ambient drift (with step_76 T ~ 1 this is the eps ~ 1 branch "
        "excluded by the SN luminosity floor).  pin << 1 sustained => "
        "differential drift: interiors pin, shallow landscape rolls, and "
        "the secular redshift can be carried by the transport integral "
        "between endpoints rather than by interior clock drift.")

    with open(OUT_PATH, "w") as fh:
        json.dump(results, fh, indent=1)
    print(f"  -> {OUT_PATH}")


if __name__ == "__main__":
    main()
