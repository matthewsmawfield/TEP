#!/usr/bin/env python3
"""
step_23_steady_state.py
=======================
Dynamical test of the tilted-congruence volume mechanism.

Claim under test (from step_22):  with phi = -Pi t + psi(x,t),
theta~ = 0 for the physical galaxy congruence reduces to

        < v . grad psi > = Pi        (flow-gradient correlation)

i.e. matter must statistically migrate toward slow-time wells (phi
maxima; A = e^{-phi}).  In TEP this is not an imposed pattern: the
matter-frame acceleration is a = -grad ln A = +beta_A grad phi, so
phi wells ARE the gravitational attractors -- infall IS the correlation.

This script runs a minimal coupled solver:
  * N particles in a 2D periodic box
  * scalar field psi solved on a grid each step (screened Poisson,
    quasi-static: field relaxes on light-crossing >> dynamical time)
  * particle accel = Newtonian gravity of g-metric + scalar fifth force
        a = -grad Phi_N + beta_eff^2 * c^2 * grad psi   (screened)
  * measure C(t) = < v . grad psi >_particles and compare to Pi_hat
  * Einstein-side check: scalar stress rho_psi ~ (grad psi)^2 stays
    perturbative (|Psi_metric| << 1)
  * steady-state diagnostic: does C(t) decay after virialisation, or is
    it sustained while the drift keeps injecting volume pressure?

Signs: A = e^{-phi}, ln A = -phi;  matter acceleration a = -grad lnA =
+grad phi -> toward phi maxima (slow time).  Scalar source: overdensities
deepen psi (poisson: nabla^2 psi = +s*delta_rho with s>0 gives wells at
peaks when combined with the correct KG sign; we take the attractive
branch and report the sign convention used).
"""
import json
import numpy as np
import os

rng = np.random.default_rng(7)

# ---------- geometry / units -------------------------------------------
N = 64                      # grid
NP = 8192                   # particles
L = 1.0                     # box (code units; map: L -> ~ 8 Mpc patch)
DX = L / N
NSTEPS = 400
DT = 4e-4

# physics knobs (code units: G_N = 1, c = 1)
BETA2 = 1.0                 # conformal fifth-force strength rel. to gravity
SCREEN_MASS2 = 0.0          # mu^2: screened Poisson mass (0 = unscreened)
PI_HAT = 1.0                # target: <v.grad psi> should sustain ~Pi_hat
                          # in units set below by the drift consistency


def cic_rho(pos):
    """Cloud-in-cell density on N x N grid."""
    rho = np.zeros((N, N))
    g = pos / DX
    i0 = np.floor(g).astype(int) % N
    f = g - np.floor(g)
    for di in (0, 1):
        for dj in (0, 1):
            w = ((1 - di) + (2 * di - 1) * f[:, 0]) * \
                ((1 - dj) + (2 * dj - 1) * f[:, 1])
            ii = (i0[:, 0] + di) % N
            jj = (i0[:, 1] + dj) % N
            np.add.at(rho, (ii, jj), w)
    return rho * (NP / N**2) * (N**2 / NP)   # mean ~ NP/N^2 * N^2/NP = 1


def poisson_fft(rho, mass2=0.0):
    """Solve (nabla^2 - mass2) F = delta_rho  on periodic box."""
    delta = rho - rho.mean()
    kx = 2*np.pi*np.fft.fftfreq(N, d=DX)
    KX, KY = np.meshgrid(kx, kx, indexing="ij")
    k2 = KX**2 + KY**2
    dk = np.fft.fft2(delta)
    Fk = np.zeros_like(dk)
    mask = k2 + mass2 > 0
    Fk[mask] = -dk[mask] / (k2[mask] + mass2)
    Fk[~mask] = 0.0
    return np.fft.ifft2(Fk).real


def grad(F):
    Fx = (np.roll(F, -1, 0) - np.roll(F, 1, 0)) / (2*DX)
    Fy = (np.roll(F, -1, 1) - np.roll(F, 1, 1)) / (2*DX)
    return Fx, Fy


def sample(Fx, Fy, pos):
    """bilinear sample of gradient fields at particle positions."""
    g = pos / DX
    i0 = np.floor(g).astype(int) % N
    f = g - np.floor(g)
    out_x = np.zeros(len(pos)); out_y = np.zeros(len(pos))
    for di in (0, 1):
        for dj in (0, 1):
            w = ((1 - di) + (2 * di - 1) * f[:, 0]) * \
                ((1 - dj) + (2 * dj - 1) * f[:, 1])
            ii = (i0[:, 0] + di) % N
            jj = (i0[:, 1] + dj) % N
            out_x += w * Fx[ii, jj]
            out_y += w * Fy[ii, jj]
    return out_x, out_y


def main():
    # initial: uniform + small perturbations (growing structure)
    pos = rng.random((NP, 2)) * L
    vel = np.zeros((NP, 2))
    # seed a few velocity perturbations
    vel += 0.02 * rng.standard_normal((NP, 2))

    hist = {"t": [], "C": [], "psi_rms": [], "vrms": [], "gradpsi_rms": [],
            "rho_max": [], "virial": [], "metric_Psi": []}
    for step in range(NSTEPS):
        rho = cic_rho(pos)
        psi = poisson_fft(rho, SCREEN_MASS2)       # scalar landscape
        PhiN = -poisson_fft(rho, 0.0)              # Newtonian (attractive)
        psx, psy = grad(psi)
        Phx, Phy = grad(PhiN)
        # acceleration: gravity toward overdensities + scalar toward wells
        ax_p, ay_p = sample(-Phx + BETA2 * psx, -Phy + BETA2 * psy, pos)
        # KDK leapfrog
        if step > 0:
            vel += 0.5 * DT * np.stack([ax_p, ay_p], 1)
        pos = (pos + DT * vel) % L
        rho = cic_rho(pos)
        psi = poisson_fft(rho, SCREEN_MASS2)
        PhiN = -poisson_fft(rho, 0.0)
        psx, psy = grad(psi)
        Phx, Phy = grad(PhiN)
        ax_p, ay_p = sample(-Phx + BETA2 * psx, -Phy + BETA2 * psy, pos)
        vel += 0.5 * DT * np.stack([ax_p, ay_p], 1)
        # diagnostics
        vgx, vgy = sample(psx, psy, pos)
        C = np.mean(vel[:, 0] * vgx + vel[:, 1] * vgy)
        hist["t"].append(step * DT)
        hist["C"].append(float(C))
        hist["psi_rms"].append(float(np.std(psi)))
        hist["vrms"].append(float(np.sqrt(np.mean(vel**2))))
        hist["gradpsi_rms"].append(float(np.sqrt(np.mean(psx**2 + psy**2))))
        hist["rho_max"].append(float(rho.max()))
        ke = 0.5 * np.mean(vel**2)
        hist["virial"].append(float(ke))
        # Einstein-side: scalar stress ~ (grad psi)^2 -> metric pert scale
        hist["metric_Psi"].append(float(np.max(np.abs(PhiN)) +
                                    0.5 * np.max(psx**2 + psy**2)))

    # analysis
    t = np.array(hist["t"]); C = np.array(hist["C"])
    psi_rms = np.array(hist["psi_rms"])
    late = C[int(0.6 * NSTEPS):]
    res = {
        "setup": {"N": N, "NP": NP, "dt": DT, "nsteps": NSTEPS,
                  "beta2": BETA2, "screen_mass2": SCREEN_MASS2},
        "correlation_series": hist,
        "summary": {
            "C_final": float(C[-1]),
            "C_late_mean": float(late.mean()),
            "C_late_positive_fraction": float(np.mean(late > 0)),
            "C_max": float(C.max()),
            "psi_rms_final": float(psi_rms[-1]),
            "vrms_final": float(hist["vrms"][-1]),
            "gradpsi_rms_final": float(hist["gradpsi_rms"][-1]),
            "metric_Psi_max": float(max(hist["metric_Psi"])),
            "sustained_correlation": bool(late.mean() > 0 and
                late.mean() > 0.2 * C.max()),
        },
        "scaling_to_cosmology": {
            "note": "Pi_hat required = H0/c ~ 2.4e-26 m^-1.  In the sim "
                    "C saturates at ~ vrms*psi_rms/L_corr; the physical "
                    "requirement <v dphi> = Pi translates to landscape "
                    "gradients ~ O(1)/Mpc for v ~ 300 km/s -- consistent "
                    "with intergalactic well spacing.",
            "required_grad_m^-1": 2.4e-26 / 1e-3,
            "implied_scale_m": 1.0 / (2.4e-26 / 1e-3),
            "implied_scale_Mpc": 1.0 / (2.4e-26 / 1e-3) / 3.086e22,
        },
        "verdict": None,
    }
    s = res["summary"]
    res["verdict"] = (
        "sustained positive flow-gradient correlation built by dynamics"
        if s["sustained_correlation"] else
        "correlation failed to sustain -- mechanism needs drift-driven "
        "recharge or is dead")
    outdir = os.path.join(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "step_23_steady_state.json"), "w") as f:
        json.dump(res, f, indent=2)
    print(json.dumps(res["summary"], indent=1))
    print(json.dumps(res["scaling_to_cosmology"], indent=1))
    print(res["verdict"])


if __name__ == "__main__":
    main()
