#!/usr/bin/env python3
"""
step_26_aligned_domains.py
==========================
Aligned bounded-domain landscape test — the hybrid that step_24/25
pointed to:

  * step_24: matter-sourced landscape needs beta_A ~ 2900 (rejected).
  * step_25: frozen random landscape gives ~2-30% correlation
    efficiency and Pi/H0 <= 0.23 in the physical window; large
    coherence blows up the lnA scatter bound.

The physical TEP picture is neither: V(phi) sets a BOUNDED landscape
contrast delta_lnA (domain structure of the scalar potential), while
density seeds the domain POSITIONS — matter sits in its own wells, so
the infall direction is systematically anti-aligned with the local
lnA gradient by construction.  Model:

    ell(x,t) = -delta_lnA * tanh((rho_s(x,t) - rho_t)/w_rho)

recomputed each step from the smoothed density: wells of bounded depth
delta_lnA live on structures; gradients concentrate on the isodensity
walls that accretion streams cross — the correlation is supported
exactly where inflow happens.

Force:  a = g_N - S c^2 grad ell     (screened conformal response)
Test:   Pi_inferred = -<v . grad ell>  vs  H0,
        with v_rms ~ 10^-3 c and delta_lnA <=~ 1e-2 bounded.

A self-consistent window (Pi ~ H0 at physical v with small delta_lnA)
is the concrete requirement that V(phi) must deliver.
"""
import json
import numpy as np
import os

rng = np.random.default_rng(23)

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
NSTEPS = 260

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


def landscape_grad(rho, delta_ell, rho_t, w_rho):
    """grad ell for ell = -delta_ell tanh((rho_s - rho_t)/w_rho)."""
    rs = smooth_rho(rho)
    grs = grad3(rs)
    Fp = np.cosh((rs - rho_t) / w_rho) ** -2 / w_rho
    return -delta_ell * Fp[None] * grs, rs


def run(pos0, vel0, delta_ell, rho_t, w_rho, S):
    pos = pos0.copy(); vel = vel0.copy()
    hist = {"t": [], "Pi_inf": [], "vrms": [], "rho_max": [],
            "gell_rms_Mpc": [], "wall_fraction": []}
    for step in range(NSTEPS):
        rho = cic_rho(pos)
        gell, rs = landscape_grad(rho, delta_ell, rho_t, w_rho)
        gvec = -grad3(poisson_phi(rho))
        a = cic3(gvec, pos) - S * C_SI**2 * cic3(gell, pos)
        if step > 0:
            vel += 0.5 * DT * a
        pos = (pos + DT * vel) % L
        rho = cic_rho(pos)
        gell, rs = landscape_grad(rho, delta_ell, rho_t, w_rho)
        gvec = -grad3(poisson_phi(rho))
        a = cic3(gvec, pos) - S * C_SI**2 * cic3(gell, pos)
        vel += 0.5 * DT * a
        # diagnostics
        gl = cic3(gell, pos)
        Pi_inf = -np.mean(np.sum(vel * gl, 1))
        hist["t"].append(step * DT / GYR)
        hist["Pi_inf"].append(float(Pi_inf))
        hist["vrms"].append(float(np.sqrt(np.mean(vel**2))))
        hist["rho_max"].append(float(rho.max() / RHO_M))
        hist["gell_rms_Mpc"].append(float(np.sqrt(
            np.mean(np.sum(gell**2, 0))) * MPC))
        # fraction of particles sitting on walls (gradient-bearing)
        rho_p = np.zeros(len(pos))
        gr = np.sqrt(np.sum(gell**2, 0))
        hist["wall_fraction"].append(float(
            np.mean(gr > 0.1 * gr.max())))
    return hist


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


def late(hist):
    sl = slice(int(0.6 * NSTEPS), NSTEPS)
    return {
        "Pi_inf_s-1": float(np.mean(hist["Pi_inf"][sl])),
        "Pi_over_H0": float(np.mean(hist["Pi_inf"][sl]) / H0),
        "vrms_over_c": float(np.mean(hist["vrms"][sl]) / C_SI),
        "rho_max_over_mean": float(np.mean(hist["rho_max"][sl])),
        "gell_rms_Mpc-1": float(np.mean(hist["gell_rms_Mpc"][sl])),
        "wall_fraction": float(np.mean(hist["wall_fraction"][sl])),
    }


def main():
    pos0, vel0 = initial_conditions()
    res = {"units": "SI", "box_Mpc": L_MPC, "grid": N, "NP": NP,
           "H0": H0,
           "model": "ell = -delta_lnA tanh((rho_s-rho_t)/w); "
                    "a = g_N - S c^2 grad ell; bounded domains "
                    "aligned with structure"}
    scan = {}
    for delta_ell in [1e-3, 3e-3, 1e-2]:
        for rho_t_f in [1.0, 3.0]:
            rho_t = rho_t_f * RHO_M
            w_rho = 0.5 * rho_t
            for S in [3e-4, 1e-3]:
                h = run(pos0, vel0, delta_ell, rho_t, w_rho, S)
                s = late(h)
                s["delta_lnA"] = delta_ell
                key = (f"dlA={delta_ell:g},rho_t={rho_t_f}rho,"
                       f"S={S:g}")
                scan[key] = s
                print(key, json.dumps(s), flush=True)
    res["scan"] = scan
    good = [k for k, v in scan.items()
            if (0.3 < v["Pi_over_H0"] < 3.0 and
                3e-4 < v["vrms_over_c"] < 3e-3 and
                v["delta_lnA"] <= 1e-2)]
    res["verdict"] = {
        "self_consistent_windows": good,
        "statement": (
            "aligned bounded domains sustain Pi ~ H0 at physical "
            "velocities with delta_lnA <= 1e-2 -> the V(phi)-domain "
            "channel closes the magnitude" if good else
            "aligned bounded domains still short of Pi ~ H0 in the "
            "scanned window -- report the residual gap"),
    }
    outdir = os.path.join(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "step_26_aligned_domains.json"),
              "w") as f:
        json.dump(res, f, indent=2)
    print(json.dumps(res["verdict"], indent=1))


if __name__ == "__main__":
    main()
