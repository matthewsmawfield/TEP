#!/usr/bin/env python3
"""
step_24_physical_magnitude.py
=============================
Physical-units dynamical test of the statistical-flow non-expansion
mechanism (step_22 kinematics, step_23 toy dynamics).

Claim under test
----------------
In the inhomogeneous eternal cosmology the macroscopic expansion scalar
vanishes statistically,

    <theta~> = <div v> + 3 Pi + 3 <v . grad lnA> = 0 .

On a periodic (or statistically homogeneous) box <div v> is a boundary
term and vanishes, so the drift rate is carried by the flow-gradient
correlation:

    Pi = -<v . grad lnA>.

The scalar landscape is matter-sourced:  nabla^2 phi_hat = -4 pi G beta_A
delta_rho / c^2  (screened Poisson, quasi-static).  With lnA = beta_A
phi_hat this gives the DEF-type identity

    grad lnA = beta_A^2 * g_vec / c^2 ,        a_scalar = beta_A^2 g_N

where g_vec = -grad Phi_N is the Newtonian acceleration field.  Hence

    Pi_inferred = beta_A^2 * G1 ,   G1 = <v . g> / c^2  ,

and the sim measures G1 directly:  beta_A^2 = Pi / G1 is then a *required*
coupling, not an input.  Screening enters only through the velocities:
the physical force is (1 + S(rho) beta_A^2) g, so self-consistency demands
a (beta_A, S) pair that keeps v_rms ~ 10^-3 c while Pi_inferred = H0.

Equivalent screening form:  a_eff = S c^2 grad lnA = S beta_A^2 g, so

    Pi = <v . a_eff> / (S c^2)

-- the drift magnitude is fixed by observed kinematics x screening.

Phase A: pure Newtonian PM collapse in physical SI units -> G1(t),
         required beta_A, endpoint lnA scatter (redshift-scatter bound).
Phase B: force = (1 + S(rho) beta_A^2) g at the required beta_A, scanned
         over the screening density rho_scr; checks whether an (S, beta_A)
         regime sustains Pi_inferred ~ H0 with v_rms ~ 10^-3 c.

Weightings:  the proper-volume average of theta~ weights each mass
element by its proper volume (1/rho); we report mass-, volume- and
proper-volume-weighted G1.
"""
import json
import numpy as np
import os

rng = np.random.default_rng(11)

# ---------- physical constants (SI) ------------------------------------
G_SI = 6.674e-11
C_SI = 2.998e8
H0 = 2.27e-18                       # s^-1  (70 km/s/Mpc)
RHO_C = 3.0 * H0**2 / (8.0 * np.pi * G_SI)
OMEGA_M = 0.3
RHO_M = OMEGA_M * RHO_C             # ~2.6e-27 kg/m^3
MPC = 3.086e22
GYR = 3.156e16

# ---------- geometry ----------------------------------------------------
N = 64                              # grid^3
NP = 64**3                          # particles (one per cell lattice)
L_MPC = 100.0
L = L_MPC * MPC
DX = L / N
DT = 0.08 * GYR                     # 80 Myr
NSTEPS_A = 640                      # ~51 Gyr: static-universe collapse
                                    # time 1/sqrt(4 pi G rho) ~ 21 Gyr,
                                    # so this reaches developed structure
NSTEPS_B = 250                      # ~20 Gyr

FOUR_PI_G = 4.0 * np.pi * G_SI
PARTICLE_MASS = RHO_M * L**3 / NP


def kgrid():
    k = 2.0 * np.pi * np.fft.fftfreq(N, d=DX)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    k2 = KX**2 + KY**2 + KZ**2
    k2[0, 0, 0] = 1.0
    return k2


K2 = kgrid()


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
    return rho * PARTICLE_MASS / DX**3      # kg/m^3


def cic_field(field3, pos):
    """trilinear sample of a grid vector/scalar field at positions."""
    g = pos / DX
    i0 = np.floor(g).astype(int) % N
    f = g - np.floor(g)
    if field3.ndim == 3:
        out = np.zeros(len(pos))
        for di in (0, 1):
            for dj in (0, 1):
                for dk in (0, 1):
                    w = ((1 - di) + (2 * di - 1) * f[:, 0]) * \
                        ((1 - dj) + (2 * dj - 1) * f[:, 1]) * \
                        ((1 - dk) + (2 * dk - 1) * f[:, 2])
                    out += w * field3[(i0[:, 0] + di) % N,
                                      (i0[:, 1] + dj) % N,
                                      (i0[:, 2] + dk) % N]
        return out
    outs = []
    for comp in range(3):
        outs.append(cic_field(field3[comp], pos))
    return np.stack(outs, 1)


def poisson_phi(rho):
    """Solve nabla^2 Phi = 4 pi G (rho - mean). Returns Phi (m^2/s^2)."""
    delta = rho - RHO_M
    dk = np.fft.fftn(delta)
    pk = -FOUR_PI_G * dk / K2
    pk[0, 0, 0] = 0.0
    return np.fft.ifftn(pk).real


def grad3(F):
    gx = (np.roll(F, -1, 0) - np.roll(F, 1, 0)) / (2 * DX)
    gy = (np.roll(F, -1, 1) - np.roll(F, 1, 1)) / (2 * DX)
    gz = (np.roll(F, -1, 2) - np.roll(F, 1, 2)) / (2 * DX)
    return np.stack([gx, gy, gz])


def momentum_field(pos, vel):
    """mass-weighted velocity field on the grid (volume-weightable)."""
    mom = np.zeros((3, N, N, N)); mss = np.zeros((N, N, N))
    g = pos / DX
    i0 = np.floor(g).astype(int) % N
    f = g - np.floor(g)
    for di in (0, 1):
        for dj in (0, 1):
            for dk in (0, 1):
                w = ((1 - di) + (2 * di - 1) * f[:, 0]) * \
                    ((1 - dj) + (2 * dj - 1) * f[:, 1]) * \
                    ((1 - dk) + (2 * dk - 1) * f[:, 2])
                ii, jj, kk = (i0[:, 0] + di) % N, (i0[:, 1] + dj) % N, \
                    (i0[:, 2] + dk) % N
                for c in range(3):
                    np.add.at(mom[c], (ii, jj, kk),
                              w * vel[:, c])
                np.add.at(mss, (ii, jj, kk), w)
    mask = mss > 0
    vfield = np.zeros_like(mom)
    for c in range(3):
        vfield[c][mask] = mom[c][mask] / mss[mask]
    return vfield


def accel(rho, pos, beta2_eff):
    """g_vec = -grad Phi; total a = (1 + beta2_eff(rho)) g."""
    Phi = poisson_phi(rho)
    gvec = -grad3(Phi)
    if beta2_eff is not None:
        rho_p = cic_field(rho, pos)
        g_p = cic_field(gvec, pos)
        a = g_p * (1.0 + beta2_eff(rho_p))[:, None]
        return a, Phi, gvec
    a = cic_field(gvec, pos)
    return a, Phi, gvec


def initial_conditions():
    """Zel'dovich-style ICs: lattice + displacement from a smoothed
    Gaussian density field (delta_rms ~ 0.5 at grid scale)."""
    q = np.indices((int(round(NP ** (1/3))),) * 3).reshape(3, -1).T
    q = q.astype(float) * (L / NP ** (1/3))
    delta = rng.standard_normal((N, N, N))
    dk = np.fft.fftn(delta)
    k = np.sqrt(K2)
    k0 = 8.0 * 2 * np.pi / L            # coherence scale ~ L/8 = 12.5 Mpc
    pk_shape = 1.0 / (1.0 + (k / k0) ** 4)
    dk *= np.sqrt(pk_shape)
    d = np.fft.ifftn(dk).real
    d /= d.std()
    d *= 0.8                            # delta_rms = 0.8: partly developed
    # displacement potential: nabla^2 D = -delta -> D_k = delta_k/k2
    Dk = np.fft.fftn(d) / K2
    Dk[0, 0, 0] = 0.0
    D = np.fft.ifftn(Dk).real
    disp = grad3(D)                     # zel'dovich displacement ~ grad D
    # scale displacements to ~ DX/4 rms so particles perturb mildly
    disp *= (0.25 * DX) / np.sqrt(np.mean(np.sum(disp**2, 0)))
    pos = (q + cic_field(disp, q)) % L
    # Zel'dovich-consistent growing-mode velocities:  v = Psi / tau_grow,
    # tau_grow = 1/sqrt(4 pi G rho_bar) for a static (non-expanding) background
    tau_grow = 1.0 / np.sqrt(FOUR_PI_G * RHO_M)
    vel = cic_field(disp, q) / tau_grow
    return pos, vel


def run(pos, vel, nsteps, beta2_eff, collect):
    hist = {k: [] for k in
            ["t", "G1_mass", "G1_vol", "vrms", "rho_max",
             "phi_rms_c2", "vdiv_vol"]}
    for step in range(nsteps):
        rho = cic_rho(pos)
        a, Phi, gvec = accel(rho, pos, beta2_eff)
        if step > 0:
            vel += 0.5 * DT * a
        pos = (pos + DT * vel) % L
        rho = cic_rho(pos)
        a, Phi, gvec = accel(rho, pos, beta2_eff)
        vel += 0.5 * DT * a
        # ---- diagnostics ----
        g_p = cic_field(gvec, pos)
        vdotg = np.sum(vel * g_p, 1)
        # infall: v parallel to g (both toward wells) -> <v.g> > 0 and
        # Pi_inferred = beta_A^2 <v.g>/c^2  (grad lnA = -beta_A^2 g/c^2
        # since lnA = -phi_hat and phi_hat peaks at overdensities)
        G1_mass = np.mean(vdotg) / C_SI**2           # s^-1
        # proper-volume weighting: each particle weighs 1/rho_local
        rho_p = np.maximum(cic_field(rho, pos), 1e-3 * RHO_M)
        w = 1.0 / rho_p
        G1_vol = np.sum(vdotg * w) / np.sum(w) / C_SI**2
        # volume-weighted via gridded velocity field x grad f
        vfield = momentum_field(pos, vel)
        vdotg_vol = np.mean(np.sum(vfield * gvec, 0)) / C_SI**2
        hist["t"].append(step * DT / GYR)
        hist["G1_mass"].append(float(G1_mass))
        hist["G1_vol"].append(float(G1_vol))
        hist["vrms"].append(float(np.sqrt(np.mean(vel**2))))
        hist["rho_max"].append(float(rho.max() / RHO_M))
        hist["phi_rms_c2"].append(float(np.std(Phi) / C_SI**2))
        hist["vdiv_vol"].append(float(vdotg_vol))
    return hist


def summary_stats(hist):
    n = len(hist["t"])
    sl = slice(int(0.6 * n), n)
    return {
        "G1_mass_late_s-1": float(np.mean(hist["G1_mass"][sl])),
        "G1_vol_late_s-1": float(np.mean(hist["G1_vol"][sl])),
        "vrms_late_m_s": float(np.mean(hist["vrms"][sl])),
        "vrms_late_over_c": float(np.mean(hist["vrms"][sl]) / C_SI),
        "rho_max_late_over_mean": float(np.mean(hist["rho_max"][sl])),
        "phi_rms_c2_late": float(np.mean(hist["phi_rms_c2"][sl])),
    }


def main():
    res = {
        "units": "SI", "box_Mpc": L_MPC, "grid": N, "particles": NP,
        "rho_m_kg_m3": RHO_M, "H0_s-1": H0, "dt_Gyr": DT / GYR,
        "identity": "Pi_inferred = beta_A^2 * G1, G1 = -<v.g>/c^2",
    }

    # ---------------- Phase A: Newtonian only --------------------------
    pos, vel = initial_conditions()
    histA = run(pos.copy(), vel.copy(), NSTEPS_A, None, None)
    sA = summary_stats(histA)
    G1 = sA["G1_mass_late_s-1"]
    beta2_req = H0 / G1 if G1 > 0 else float("nan")
    beta_req = np.sqrt(beta2_req) if beta2_req > 0 else float("nan")
    # lnA landscape:  lnA = beta_A^2 * f with f = -Phi/c^2
    sigma_lnA = beta2_req * sA["phi_rms_c2_late"]
    res["phase_A_newtonian"] = {
        "stats": sA,
        "G1_series": histA,
        "required_beta_A2": float(beta2_req),
        "required_beta_A": float(beta_req),
        "sigma_lnA_endpoint_scatter": float(sigma_lnA),
        "grad_lnA_Mpc-1": float(H0 / sA["vrms_late_m_s"] * MPC),
        "a_bare_m_s2": float(C_SI**2 * H0 / sA["vrms_late_m_s"]),
        "S_for_a_eff_1e-10": float(1e-10 /
            (C_SI**2 * H0 / sA["vrms_late_m_s"])),
    }

    # ---------------- Phase B: screened scalar force -------------------
    # force multiplier (1 + S(rho) beta_A^2); S = min(1, rho_scr/rho)
    phaseB = {}
    for rho_scr_factor in [1e-4, 1e-3, 1e-2, 1e-1, 1.0]:
        rho_scr = rho_scr_factor * RHO_M
        if not np.isfinite(beta2_req):
            break
        def sfun(rho_p, _rs=rho_scr):
            return np.minimum(1.0, _rs / np.maximum(rho_p, 1e-6 * RHO_M))
        def b2(rho_p, _sf=sfun):
            return beta2_req * _sf(rho_p)
        posB, velB = initial_conditions()
        histB = run(posB, velB, NSTEPS_B, b2, None)
        sB = summary_stats(histB)
        Pi_inf = beta2_req * sB["G1_mass_late_s-1"]
        # effective scalar acceleration at mean density
        S_at_mean = min(1.0, rho_scr / RHO_M)
        a_eff = S_at_mean * beta2_req * 3e-10     # typ g ~ 3e-10 m/s^2
        phaseB[f"rho_scr/{rho_scr_factor}rho_mean"] = {
            "stats": sB,
            "Pi_inferred_s-1": float(Pi_inf),
            "Pi_over_H0": float(Pi_inf / H0),
            "S_at_mean_density": float(S_at_mean),
            "a_eff_at_mean_density_m_s2": float(a_eff),
            "S_beta2_at_mean": float(S_at_mean * beta2_req),
        }
    res["phase_B_screened"] = phaseB

    # ---------------- verdict ------------------------------------------
    okA = np.isfinite(beta_req) and 1.0 < beta_req < 1e4
    best = None
    for key, v in phaseB.items():
        score = abs(np.log10(max(v["Pi_over_H0"], 1e-30)))
        score += abs(np.log10(max(v["stats"]["vrms_late_over_c"], 1e-30)
                     / 1e-3))
        if best is None or score < best[0]:
            best = (score, key)
    res["verdict"] = {
        "G1_measured": "flow-gradient correlation builds dynamically",
        "beta_A_required": float(beta_req),
        "best_screening_window": best[1] if best else None,
        "consistency": (
            "self-consistent regime exists: Pi ~ H0 at v_rms ~ 10^-3 c"
            if best and abs(np.log10(
                phaseB[best[1]]["Pi_over_H0"])) < 0.7 and
                3e-4 < phaseB[best[1]]["stats"]["vrms_late_over_c"] < 3e-3
            else "no fully self-consistent (beta,S) window found -- "
                 "magnitude gap quantified in phase_A"),
    }
    outdir = os.path.join(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "step_24_physical_magnitude.json"),
              "w") as f:
        json.dump(res, f, indent=2)
    print(json.dumps(res["phase_A_newtonian"]["stats"], indent=1))
    print("beta_A required:", res["phase_A_newtonian"]["required_beta_A"])
    print(json.dumps({k: {"Pi/H0": v["Pi_over_H0"],
                          "vrms/c": v["stats"]["vrms_late_over_c"],
                          "S_beta2": v["S_beta2_at_mean"]}
                      for k, v in phaseB.items()}, indent=1))
    print(json.dumps(res["verdict"], indent=1))


if __name__ == "__main__":
    main()
