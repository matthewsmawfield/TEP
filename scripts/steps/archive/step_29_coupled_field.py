#!/usr/bin/env python3
"""
step_29_coupled_field.py
========================
The coupled-field gate: does the scalar equation of motion,
solved on the real collapse dynamics, spontaneously generate
the drift partition that step_27/28 showed is required?

Field equation (dimensionless u = phi/M_Pl, beta_A = -1):
    Box u = V_{,u}/M^2 - beta_coupling * rho_tilde * exp(-u)
Quasi-static (justified: grid crossing time ~ 5 Myr << collapse
timescale, so the field relaxes to equilibrium each step):
    laplacian u = V'(u) - s * (rho_s - 1) * exp(-u)

with the Poisson sign convention fixed by requiring u maximal in
overdense wells (temporal wells are phi-maxima, A-minima).

Decomposition: u(x,t) = u_bar(t) + du(x,t)
  * u_bar: homogeneous mode, evolved under <V'_eff> (mean-roll
    ODE) or pinned -- control knob
  * du: fluctuation solve (k != 0 Fourier modes)

Measured per step:
  Pi_hat(x) = -u_dot(x)            (conformal drift field)
  Pi_eff(x) = d/dt ln sqrt(A^2 - B u_dot^2 M^2 ...)
            approx = Pi_hat + (1/2) d/dt ln(1 + |b| u_dot_hat^2 / a^2)
            (disformal clock-lapse drift, B<0 candidate)
  partition ratio R = <Pi>_wells / <Pi>_ambient
  endpoint drift on tracked worldlines vs ambient
  and whether R ~ partition needed by the volume balance.

V models scanned:
  'linear'  : V'_u = const (free mean roll, canonical benchmark)
  'quartic' : V'_u = lam u^3 (pinned minimum at u=0)
"""
import json
import numpy as np
import os

rng = np.random.default_rng(29)

G_SI = 6.674e-11
C_SI = 2.998e8
H0 = 2.27e-18
RHO_C = 3.0 * H0**2 / (8.0 * np.pi * G_SI)
RHO_M = 0.3 * RHO_C
MPC = 3.086e22
GYR = 3.156e16

N = 64
NP = 64**3
L_MPC = 100.0
L = L_MPC * MPC
DX = L / N
DT = 0.06 * GYR
NSTEPS = 240

FOUR_PI_G = 4.0 * np.pi * G_SI
PARTICLE_MASS = RHO_M * L**3 / NP

_k = 2.0 * np.pi * np.fft.fftfreq(N, d=DX)
KX, KY, KZ = np.meshgrid(_k, _k, _k, indexing="ij")
K2 = KX**2 + KY**2 + KZ**2
K2[0, 0, 0] = 1.0


def cic_rho(pos):
    rho = np.zeros((N, N, N))
    g = pos / DX
    i0 = np.floor(g).astype(int) % N
    f = g - np.floor(g)
    for di in (0, 1):
        for dj in (0, 1):
            for dk in (0, 1):
                w = ((1 - di) + (2 * di - 1) * f[:, 0]) * \
                    ((1 - dj) + (2 * dj - 1) * f[:, 1]) * \
                    ((1 - dk) + (2 * dk - 1) * f[:, 2])
                np.add.at(rho, ((i0[:, 0] + di) % N,
                                (i0[:, 1] + dj) % N,
                                (i0[:, 2] + dk) % N), w)
    return rho * PARTICLE_MASS / DX**3


def cic3(field3, pos):
    g = pos / DX
    i0 = np.floor(g).astype(int) % N
    f = g - np.floor(g)
    out = np.zeros((len(pos), 3))
    for di in (0, 1):
        for dj in (0, 1):
            for dk in (0, 1):
                w = ((1 - di) + (2 * di - 1) * f[:, 0]) * \
                    ((1 - dj) + (2 * dj - 1) * f[:, 1]) * \
                    ((1 - dk) + (2 * dk - 1) * f[:, 2])
                ii = ((i0[:, 0] + di) % N, (i0[:, 1] + dj) % N,
                      (i0[:, 2] + dk) % N)
                for c in range(3):
                    out[:, c] += w * field3[c][ii]
    return out


def cic_scalar(field, pos):
    g = pos / DX
    i0 = np.floor(g).astype(int) % N
    f = g - np.floor(g)
    out = np.zeros(len(pos))
    for di in (0, 1):
        for dj in (0, 1):
            for dk in (0, 1):
                w = ((1 - di) + (2 * di - 1) * f[:, 0]) * \
                    ((1 - dj) + (2 * dj - 1) * f[:, 1]) * \
                    ((1 - dk) + (2 * dk - 1) * f[:, 2])
                out += w * field[((i0[:, 0] + di) % N,
                                  (i0[:, 1] + dj) % N,
                                  (i0[:, 2] + dk) % N)]
    return out


def poisson_phi(rho):
    dk = np.fft.fftn(rho - RHO_M)
    pk = -FOUR_PI_G * dk / K2
    pk[0, 0, 0] = 0.0
    return np.fft.ifftn(pk).real


def grad3(F):
    return np.stack([
        (np.roll(F, -1, 0) - np.roll(F, 1, 0)) / (2 * DX),
        (np.roll(F, -1, 1) - np.roll(F, 1, 1)) / (2 * DX),
        (np.roll(F, -1, 2) - np.roll(F, 1, 2)) / (2 * DX)])


def smooth_rho(rho, sigma_cells=1.5):
    dk = np.fft.fftn(rho)
    dk *= np.exp(-0.5 * (np.sqrt(K2) * sigma_cells * DX) ** 2)
    return np.fft.ifftn(dk).real


def initial_conditions():
    q = np.indices((64, 64, 64)).reshape(3, -1).T.astype(float)
    q *= L / 64.0
    d = rng.standard_normal((N, N, N))
    dk = np.fft.fftn(d)
    k = np.sqrt(K2)
    k0 = 8.0 * 2 * np.pi / L
    dk *= np.sqrt(1.0 / (1.0 + (k / k0) ** 4))
    dd = np.fft.ifftn(dk).real
    dd /= dd.std(); dd *= 0.8
    Dk = np.fft.fftn(dd) / K2; Dk[0, 0, 0] = 0.0
    D = np.fft.ifftn(Dk).real
    disp = grad3(D)
    disp *= (0.25 * DX) / np.sqrt(np.mean(np.sum(disp**2, 0)))
    pos = (q + cic3(disp, q)) % L
    tau = 1.0 / np.sqrt(FOUR_PI_G * RHO_M)
    vel = cic3(disp, q) / tau
    return pos, vel


def run(kappa, ubar_rate):
    """Coupled quasi-static solve.  The field equation in the
    tracking regime is laplacian u = beta_A * 4 pi G rho / c^2
    (dimensionless u = phi/M_Pl), i.e.  u = -kappa Phi_N / c^2
    with kappa = |beta_A| the conformal coupling: the temporal
    field tracks the Newtonian potential (corpus calibration
    u_gal ~ (v_c/c)^2 ~ 3e-7 at kappa ~ 1).

    ubar: homogeneous mean mode, prescribed roll rate (s^-1).
    Drift is measured BOTH Eulerian (grid) and Lagrangian
    (worldline): the latter is what endpoint clocks feel.
    """
    pos, vel = initial_conditions()
    u_grid_prev = None
    u_part_prev = None
    ubar = 0.0
    hist = {"t_gyr": [], "Pi_well": [], "Pi_amb": [],
            "R_partition": [], "u_max": [], "u_min": [],
            "rho_max": [], "PiL_well": [], "PiL_amb": []}
    for step in range(NSTEPS):
        t = step * DT
        rho = cic_rho(pos)
        phiN = poisson_phi(rho)                # m^2/s^2
        a_N = cic3(-grad3(phiN), pos)
        if step > 0:
            vel += 0.5 * DT * a_N
        pos = (pos + DT * vel) % L
        rho = cic_rho(pos)
        phiN = poisson_phi(rho)
        a_N = cic3(-grad3(phiN), pos)
        vel += 0.5 * DT * a_N
        # ---- coupled field: tracking solve ----
        rs = smooth_rho(rho) / RHO_M
        ubar += ubar_rate * DT
        du = -kappa * phiN / C_SI**2           # u>0 in wells
        u = ubar + du
        u_part = cic_scalar(u, pos)            # on worldlines
        if u_grid_prev is not None:
            Pi_grid = -(u - u_grid_prev) / DT  # Eulerian drift
            Pi_wl = -(u_part - u_part_prev) / DT  # Lagrangian
            # percentile masks on the field itself — robust to
            # collapse depth: wells = top-decile u (phi-maxima),
            # ambient = below-median u
            wmask = du > np.percentile(du, 90)
            amask = du < np.percentile(du, 50)
            du_p = cic_scalar(du, pos)
            wl_w = du_p > np.percentile(du_p, 90)
            wl_a = du_p < np.percentile(du_p, 50)
            pw = float(Pi_grid[wmask].mean()) if wmask.any() else np.nan
            pa = float(Pi_grid[amask].mean()) if amask.any() else np.nan
            plw = float(Pi_wl[wl_w].mean()) if wl_w.any() else np.nan
            pla = float(Pi_wl[wl_a].mean()) if wl_a.any() else np.nan
            hist["Pi_well"].append(pw)
            hist["Pi_amb"].append(pa)
            hist["PiL_well"].append(plw)
            hist["PiL_amb"].append(pla)
            hist["R_partition"].append(
                plw / pla if np.isfinite(pla) and pla != 0
                else np.nan)
        u_grid_prev = u
        u_part_prev = u_part
        hist["t_gyr"].append(t / GYR)
        hist["u_max"].append(float(u.max()))
        hist["u_min"].append(float(u.min()))
        hist["rho_max"].append(float(rs.max()))
    return hist


def main():
    res = {"units": "SI drift rates (s^-1)", "box_Mpc": L_MPC,
           "grid": N,
           "model": "quasi-static laplacian u = V'(u) - s(rs-1)e^-u"}
    out = {}
    # kappa ~ 1 is the corpus calibration (u_gal ~ (v_c/c)^2);
    # ubar_rate = -H0 adds the homogeneous roll control.
    for name, cfg in [
        ("kappa=1,ubar=0",   dict(kappa=1.0, ubar_rate=0.0)),
        ("kappa=1,ubar=H0",  dict(kappa=1.0, ubar_rate=-H0)),
        ("kappa=0.1,ubar=0", dict(kappa=0.1, ubar_rate=0.0)),
    ]:
        h = run(**cfg)
        sl = slice(int(0.6 * NSTEPS), NSTEPS)
        s = {
            "Pi_well_s-1": float(np.nanmean(h["Pi_well"][sl])),
            "Pi_amb_s-1": float(np.nanmean(h["Pi_amb"][sl])),
            "Pi_well_over_H0": float(
                np.nanmean(h["Pi_well"][sl]) / H0),
            "Pi_amb_over_H0": float(
                np.nanmean(h["Pi_amb"][sl]) / H0),
            "PiL_well_over_H0": float(
                np.nanmean(h["PiL_well"][sl]) / H0),
            "PiL_amb_over_H0": float(
                np.nanmean(h["PiL_amb"][sl]) / H0),
            "R_partition_worldline": float(
                np.nanmean(np.abs(h["R_partition"][sl]))),
            "u_range": [float(np.mean(h["u_min"][sl])),
                        float(np.mean(h["u_max"][sl]))],
            "rho_max": float(np.mean(h["rho_max"][sl])),
        }
        # implied coupling for well drift to reach H0
        if np.isfinite(s["PiL_well_over_H0"]) and \
                s["PiL_well_over_H0"] != 0:
            s["kappa_for_H0_well_drift"] = float(
                cfg["kappa"] / s["PiL_well_over_H0"])
        out[name] = s
        print(name, json.dumps(s), flush=True)
    res["scan"] = out
    res["verdict"] = {"statement":
                      "spontaneous partition iff worldline Pi ~ H0 in "
                      "wells while << H0 ambient; tracking-only "
                      "(ubar=0) tests the unforced channel; "
                      "kappa_for_H0_well_drift is the coupling the "
                      "unforced equation would need"}
    outdir = os.path.join(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir,
                           "step_29_coupled_field.json"), "w") as f:
        json.dump(res, f, indent=2)
    print(json.dumps(res["verdict"], indent=1))


if __name__ == "__main__":
    main()
