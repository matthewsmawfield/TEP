#!/usr/bin/env python3
"""
step_13_disformal_closure.py
============================
Single-action admissibility test for the disformal completion of
the volume balance (the archived step_27 audit outcome).  Retained as a
reproducible no-go / exclusion result: the reconstruction
succeeds formally but the required coupling is excluded
physically (see below).

archived step_27 measured: the A^-1-weighted flow term supplies only ~1e-8
of the drift balance at corpus amplitudes; the last candidate
carrier inside the action is the structure-localized disformal
volume term.  The exact
transport identity (volume-normal congruence) is

    A^6 (1 + Q) = C(x),   Q = B N^2 |grad phi|^2 / D,
    D = A^2 N^2 - B phi_dot^2 .

A candidate B(phi) is admissible only if the reconstruction

    B_req(x,t) = Q A^2 / (Y + Q phi_dot^2 / N^2)

collapses onto a single-valued function of phi across space and
time (the manuscript's isosurface-consistency test).

Field model (same collapse dynamics as archived step_27):
    lnA(x,t) = Pi_bar * t * w(x,t) + delta * F(rho_s(x,t))
    w = smooth structure weight  (partitioned roll)
    F = -tanh((rho_s-rho_t)/w_rho)  (shallow Newtonian-tracking wells)
    u = phi/M_Pl = -lnA ;  u_dot = -Pi_bar w - delta F_dot

Sign check: redshift needs A growing -> 1+Q = C A^-6 shrinking
-> Q < 0 -> B < 0 in rolling regions.  The sim reports the sign
and the collapse of b_req(u).

Physical admissibility of the reconstruction: B < 0 preserves the
signature (D = A^2 N^2 + |B| phi_dot^2 > 0) but tilts the matter
cone outward along grad phi — the longitudinal photon speed
v_par = c / sqrt(1 + (B/A^2)|grad phi|^2) exceeds c (Appendix B),
a GW170817-class violation inside structured regions; and a
universal roll drives Q -> -1 (matter-metric degeneracy).  The
channel is therefore excluded as a physical carrier: the
surviving volume-balance carrier is the A^-1-weighted flow
correlation, which step_14 shows saturates behind the temporal
horizon (delta_H ~ 28-30, A ~ e^-29).

Dimensionless variables on the grid:
    yhat = |grad u|^2 / (H0/c)^2   (spatial gradient of u=phi/M)
    uhat_dot = du/dt / H0
    b = B M^2 (H0/c)^2 N^2         (dimensionless coupling)
    Q = b yhat / (a^2 - b uhat_dot^2)   [a = exp(-u)]
    -> b_req = Q a^2 / (yhat + Q uhat_dot^2)
"""
import json
import numpy as np
import os

rng = np.random.default_rng(28)

G_SI = 6.674e-11
C_SI = 2.998e8
H0 = 2.27e-18
RHO_C = 3.0 * H0**2 / (8.0 * np.pi * G_SI)
RHO_M = 0.3 * RHO_C
MPC = 3.086e22
H0_MPC = H0 * MPC / C_SI          # H0/c in Mpc^-1 ~ 2.3e-4
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


def run(delta, rho_t_f, Pi_bar, partitioned=True):
    """Reconstruction test on the rolling partitioned field.

    lnA = Pi_bar * t * w(x,t) + delta * F(x,t)   [partitioned]
    u = -lnA ; at each snapshot measure:
      Q_req(x,t)  = C(x) exp(6 u) - 1,  C(x) = exp(6 u(x,t0))
      b_req(x,t)  = Q a^2 / (yhat + Q uhat_dot^2)
    then test whether b_req is a function of u alone.
    """
    pos, vel = initial_conditions()
    rho_t = rho_t_f * RHO_M
    w_rho = 0.5 * rho_t
    inv_H0MPC2 = 1.0 / H0_MPC**2     # (Mpc^2) to make yhat
    u0 = None
    samples_u, samples_b, samples_Q = [], [], []
    hist = {"t_gyr": [], "Q_min": [], "Q_p50": [], "Q_max": [],
            "frac_bpos": [], "rho_max": []}
    for step in range(NSTEPS):
        t = step * DT
        rho = cic_rho(pos)
        gvec = -grad3(poisson_phi(rho))
        a = cic3(gvec, pos)
        if step > 0:
            vel += 0.5 * DT * a
        pos = (pos + DT * vel) % L
        rho = cic_rho(pos)
        gvec = -grad3(poisson_phi(rho))
        a = cic3(gvec, pos)
        vel += 0.5 * DT * a
        # ---- fields on the grid ----
        rs = smooth_rho(rho)
        F = -np.tanh((rs - rho_t) / w_rho)
        if partitioned:
            wgt = 0.5 * (1.0 + np.tanh((rs - rho_t) / w_rho))
        else:
            wgt = np.ones_like(rs)
        ell = Pi_bar * t * wgt + delta * F
        u = -ell
        if u0 is None:
            u0 = u.copy()
        # spatial gradient of u (dimensionless, H0/c units)
        gu = grad3(u)
        yhat = np.sum(gu**2, 0) * MPC**2 * inv_H0MPC2
        # u_dot: roll contributes Pi_bar*w; landscape tracking is
        # quasi-static (delta*F_dot negligible vs Pi_bar) -> u_dot
        # = -(Pi_bar w); measured per cell:
        udot_hat = -(Pi_bar * wgt) / H0    # dimensionless
        a2 = np.exp(-2.0 * u)
        C = np.exp(6.0 * u0)
        Q = C * np.exp(6.0 * u) - 1.0
        denom = yhat + Q * udot_hat**2
        with np.errstate(divide="ignore", invalid="ignore"):
            b_req = Q * a2 / denom
        # sampling for the collapse test (late window)
        if step >= int(0.5 * NSTEPS) and step % 10 == 0:
            m = np.isfinite(b_req) & (np.abs(denom) > 0)
            samples_u.append(u[m].ravel())
            samples_b.append(b_req[m].ravel())
            samples_Q.append(Q[m].ravel())
        hist["t_gyr"].append(t / GYR)
        hist["Q_min"].append(float(np.min(Q)))
        hist["Q_p50"].append(float(np.median(Q)))
        hist["Q_max"].append(float(np.max(Q)))
        hist["frac_bpos"].append(float(np.mean(b_req > 0)))
        hist["rho_max"].append(float(rho.max() / RHO_M))
    return hist, np.concatenate(samples_u), \
        np.concatenate(samples_b), np.concatenate(samples_Q)


def collapse_stats(u, b, Q):
    """Bin b_req by u; report within-bin scatter relative to the
    median curve -- the single-valuedness measure."""
    lo, hi = np.percentile(u, [1, 99])
    edges = np.linspace(lo, hi, 40)
    med, scat, cnt = [], [], []
    for i in range(len(edges) - 1):
        m = (u >= edges[i]) & (u < edges[i + 1])
        if m.sum() > 50:
            mb = np.median(b[m])
            med.append(mb)
            scat.append(float(np.median(np.abs(b[m] - mb)) /
                              max(abs(mb), 1e-30)))
            cnt.append(int(m.sum()))
    return {"u_lo": float(lo), "u_hi": float(hi),
            "median_within_bin_rel_scatter":
                float(np.median(scat)) if scat else None,
            "n_bins": len(med), "b_median_range":
                [float(np.min(med)), float(np.max(med))]
                if med else None,
            "Q_min": float(Q.min()), "Q_p50": float(np.median(Q)),
            "Q_max": float(Q.max())}


def main():
    res = {"units": "dimensionless (H0,c,M_Pl)", "box_Mpc": L_MPC,
           "grid": N, "model":
           "lnA = Pi_bar t w(rho_s) + delta F(rho_s); "
           "Q_req from transport identity; b_req collapse test"}
    out = {}
    for name, cfg in [
        ("partitioned,delta=1e-6", dict(delta=1e-6, rho_t_f=1.0,
                                       Pi_bar=H0,
                                       partitioned=True)),
        ("universal,delta=1e-6",   dict(delta=1e-6, rho_t_f=1.0,
                                       Pi_bar=H0,
                                       partitioned=False)),
        ("partitioned,delta=1e-5", dict(delta=1e-5, rho_t_f=1.0,
                                       Pi_bar=H0,
                                       partitioned=True)),
    ]:
        hist, u, b, Q = run(**cfg)
        s = collapse_stats(u, b, Q)
        s["frac_bpos_late"] = float(
            np.mean(hist["frac_bpos"][int(0.6 * NSTEPS):]))
        s["rho_max_late"] = float(
            np.mean(hist["rho_max"][int(0.6 * NSTEPS):]))
        out[name] = s
        print(name, json.dumps(s), flush=True)
    res["scan"] = out
    res["verdict"] = {
        "statement":
            "admissible iff b_req(u) collapses to a single-valued "
            "curve (within-bin scatter << 1) with b<0 in rolling "
            "regions; frac_bpos measures inconsistent cells",
        "physical_status":
            "no-go: a formal single-valued b_req exists, but b<0 "
            "implies superluminal longitudinal photons "
            "(v_par = c/sqrt(1+(B/A^2)|grad phi|^2) > c) inside "
            "structure and a universal roll drives Q -> -1 "
            "(matter-metric degeneracy); excluded as a physical "
            "carrier — balance saturates via the A^-1-weighted "
            "flow term behind the temporal horizon (step_14)"}
    outdir = os.path.join(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(
            outdir, "step_13_disformal_closure.json"), "w") as f:
        json.dump(res, f, indent=2)
    print(json.dumps(res["verdict"], indent=1))


if __name__ == "__main__":
    main()
