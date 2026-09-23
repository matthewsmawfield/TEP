#!/usr/bin/env python3
"""
step_27_partition_audit.py
==========================
Corrected-measure closure test.  After the reviewer audit:

  * the average of the tilted-congruence expansion scalar is
    A^{-1}-weighted, NOT coordinate-volume weighted:
        <theta~> = <A^{-1}(3 Pi_hat + 4 v . grad lnA)> = 0
    since  <A^{-1} div v> = <A^{-1} v . grad lnA>  (parts).
  * redshift survival requires the drift to act on matter
    worldlines (well interiors must roll), tested by measuring
    D lnA / Dt along simulated particle trajectories and by
    pairwise endpoint redshifts between worldlines.
  * the landscape must pass a gradient-energy audit:
        rho_grad / rho_c = <|grad lnA|^2> / (6 H0^2)   (<=~ 1)

Model under test (corpus-consistent, drift-partitioned):

    lnA(x,t) = Pi_bar * t  +  delta * F(rho_s(x,t)),
    F = -tanh((rho_s - rho_t)/w_rho)     (wells on structure)

with Pi_bar = H0 the universal roll rate and delta ~ 1e-6 the
shallow Newtonian-tracking contrast.  The landscape is passive
(no scalar force on particles): consistent with screened wells.

Measured:
  balance:   B = -(4/3)<A^-1 v.grad l> / <A^-1 Pi_hat>
             (B -> 1 means flow supplies the balance)
  energy:    rho_grad/rho_c on the grid
  worldline: <D l/Dt> split into roll + d(delta F)/dt parts,
             for well-resident vs ambient particles
  redshift:  z_ij = exp(l_i(t_o) - l_j(t_e))-1 for particle
             pairs at separation r = c(t_o-t_e), binned vs r —
             does the Hubble slope Pi_bar survive the inflow?
"""
import json
import numpy as np
import os

rng = np.random.default_rng(27)

G_SI = 6.674e-11
C_SI = 2.998e8
H0 = 2.27e-18            # s^-1 ; in c=1 units  H0/c ~ Mpc^-1 scale
H0_MPC = H0 * C_SI / 3.086e22   # ~2.2e-4 Mpc^-1
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
NSTEPS = 300

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


def cic_scalar(field, pos):
    """interpolate a scalar grid field at particle positions"""
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


def landscape(rho, delta, rho_t, w_rho):
    """ell_spatial = -delta tanh((rho_s-rho_t)/w); returns F grid
    (in [-1,1]-ish, wells negative) and its gradient."""
    rs = smooth_rho(rho)
    F = -np.tanh((rs - rho_t) / w_rho)
    return F, grad3(F)


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


def run(delta, rho_t_f, Pi_bar):
    """delta: well depth in lnA; Pi_bar: universal roll rate (s^-1)."""
    pos, vel = initial_conditions()
    rho_t = rho_t_f * RHO_M
    w_rho = 0.5 * rho_t
    hist = {"t_gyr": [], "AinvPi": [], "Ainv_vgl": [],
            "balance": [], "rho_grad_ratio": [], "vrms": [],
            "rho_max": [], "Dl_Dt_well": [], "Dl_Dt_amb": [],
            "roll_frac_well": []}
    ell_prev = None
    pos_prev = None
    for step in range(NSTEPS):
        t = step * DT
        rho = cic_rho(pos)
        F, gF = landscape(rho, delta, rho_t, w_rho)
        gell = [delta * g for g in gF]          # grad lnA_spatial
        gvec = -grad3(poisson_phi(rho))
        a = cic3(gvec, pos)                     # Newtonian only
        if step > 0:
            vel += 0.5 * DT * a
        pos = (pos + DT * vel) % L
        rho = cic_rho(pos)
        F, gF = landscape(rho, delta, rho_t, w_rho)
        gell = [delta * g for g in gF]
        gvec = -grad3(poisson_phi(rho))
        a = cic3(gvec, pos)
        vel += 0.5 * DT * a
        # ---- diagnostics ----
        ell_p = Pi_bar * t + delta * cic_scalar(F, pos)   # lnA on worldlines
        Ainv_p = np.exp(-ell_p)
        gl_p = cic3(gell, pos)
        AinvPi = float(np.mean(Ainv_p) * Pi_bar)
        Ainv_vgl = float(np.mean(Ainv_p * np.sum(vel * gl_p, 1)))
        # gradient energy audit (grid rms)
        g2 = float(np.mean(np.sum(np.array(gell)**2, 0)))  # m^-2
        rho_grad_ratio = g2 / (6.0 * (H0 / C_SI)**2)
        # worldline derivative: Dl/Dt via stored previous step
        if ell_prev is not None:
            Dl = (ell_p - ell_prev) / DT
            well = cic_scalar(F, pos) < -0.5     # deep inside wells
            amb = cic_scalar(F, pos) > 0.5       # ambient
            hist["Dl_Dt_well"].append(float(np.mean(Dl[well])))
            hist["Dl_Dt_amb"].append(float(np.mean(Dl[amb])))
            hist["roll_frac_well"].append(
                float(np.mean(Dl[well]) / Pi_bar))
        ell_prev = ell_p
        hist["t_gyr"].append(t / GYR)
        hist["AinvPi"].append(AinvPi)
        hist["Ainv_vgl"].append(Ainv_vgl)
        hist["balance"].append(-(4.0 / 3.0) * Ainv_vgl / AinvPi
                               if AinvPi else 0.0)
        hist["rho_grad_ratio"].append(rho_grad_ratio)
        hist["vrms"].append(float(np.sqrt(np.mean(vel**2))))
        hist["rho_max"].append(float(rho.max() / RHO_M))
    return pos, vel, hist


def endpoint_redshift(pos, vel, delta, rho_t_f, Pi_bar, t_obs):
    """Pairwise endpoint z: observer particle i at t_obs, emitter j
    at t_e = t_obs - r_ij/c.  lnA(t) = Pi_bar t + delta F(x,t);
    evaluate F at current positions (slow-moving wells approx).
    Return binned z vs separation."""
    rho = cic_rho(pos)
    rho_t = rho_t_f * RHO_M
    F, _ = landscape(rho, delta, rho_t, 0.5 * rho_t)
    ell_now = Pi_bar * t_obs + delta * cic_scalar(F, pos)
    # subsample pairs
    idx = rng.choice(len(pos), 4000, replace=False)
    P = pos[idx]; E = ell_now[idx]
    # pairwise on a subset: 2000 obs x 2000 em
    io = rng.choice(len(P), 2000, replace=False)
    ie = rng.choice(len(P), 2000, replace=False)
    d = P[io][:, None, :] - P[ie][None, :, :]
    d = (d + L / 2) % L - L / 2
    r = np.sqrt(np.sum(d**2, -1))               # m
    dt_lt = r / C_SI                            # lookback time
    z = np.exp(E[io][:, None]
               - (Pi_bar * (t_obs - dt_lt) + E[ie][None, :]
                  - Pi_bar * t_obs)) - 1.0
    # = exp(Pi_bar dt_lt + delta(F_i - F_j)) - 1
    r_mpc = r / MPC
    bins = [2, 5, 10, 20, 40, 70, 100]
    out = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = (r_mpc.flatten() > lo) & (r_mpc.flatten() <= hi)
        if m.sum() > 0:
            out.append({"r_lo": lo, "r_hi": hi,
                        "z_mean": float(np.mean(z.flatten()[m])),
                        "z_std": float(np.std(z.flatten()[m])),
                        "n": int(m.sum())})
    return out


def main():
    res = {"units": "SI", "box_Mpc": L_MPC, "grid": N,
           "model": "lnA = Pi_bar t - delta tanh((rho_s-rho_t)/w); "
                    "Newtonian PM; passive shallow landscape; "
                    "A^-1-weighted balance measure"}
    configs = [
        ("delta=1e-6,Pi=H0",   dict(delta=1e-6, rho_t_f=1.0,
                                   Pi_bar=H0)),
        ("delta=1e-6,Pi=0",    dict(delta=1e-6, rho_t_f=1.0,
                                   Pi_bar=0.0)),
        ("delta=1e-5,Pi=H0",   dict(delta=1e-5, rho_t_f=1.0,
                                   Pi_bar=H0)),
        ("delta=1e-6,Pi=H0,rho_t=3", dict(delta=1e-6, rho_t_f=3.0,
                                         Pi_bar=H0)),
    ]
    out = {}
    for name, cfg in configs:
        pos, vel, h = run(cfg["delta"], cfg["rho_t_f"],
                          cfg["Pi_bar"])
        sl = slice(int(0.6 * NSTEPS), NSTEPS)
        z_bins = endpoint_redshift(pos, vel, cfg["delta"],
                                   cfg["rho_t_f"], cfg["Pi_bar"],
                                   NSTEPS * DT)
        s = {
            "AinvPi_s-1": float(np.mean(h["AinvPi"][sl])),
            "Ainv_vgl_s-1": float(np.mean(h["Ainv_vgl"][sl])),
            "balance_frac": float(np.mean(h["balance"][sl])),
            "rho_grad_over_rho_c": float(
                np.mean(h["rho_grad_ratio"][sl])),
            "vrms_over_c": float(np.mean(h["vrms"][sl]) / C_SI),
            "rho_max_over_mean": float(np.mean(h["rho_max"][sl])),
            "Dl_Dt_over_Pi_well": float(
                np.mean(h["roll_frac_well"][sl])),
            "Dl_Dt_amb_over_Pi": float(
                np.mean(np.array(h["Dl_Dt_amb"][sl]) /
                        cfg["Pi_bar"])) if cfg["Pi_bar"] else 0.0,
            "endpoint_z_bins": z_bins,
        }
        out[name] = s
        print(name, json.dumps(s)[:600], flush=True)
    res["scan"] = out
    res["verdict"] = {
        "statement": ("redshift survives iff Dl/Dt~Pi on well "
                      "worldlines; balance_frac measures the share "
                      "the flow term supplies; rho_grad/rho_c is the "
                      "energetic admissibility audit")}
    outdir = os.path.join(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "step_27_partition_audit.json"),
              "w") as f:
        json.dump(res, f, indent=2)
    print(json.dumps(res["verdict"], indent=1))


if __name__ == "__main__":
    main()
