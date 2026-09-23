#!/usr/bin/env python3
"""
step_21_evolving_halo.py
=========================
Tests the candidate mechanism for SYSTEMATIC transport redshift:
the TEP analog of the integrated Sachs-Wolfe / Rees-Sciama effect.

A photon crossing a STATIC field structure gets zero net frequency
shift (entry/exit cancel -- this is Jakarta's Theorem 2 for static
conformal gradients, checked numerically here).  A photon crossing an
EVOLVING structure picks up a residual proportional to the field's
evolution during transit.  In a statistically steady but irreversible
universe (structures deepen via collapse/mergers), each crossing is
biased toward one sign -> net redshift ~ (per-crossing shift) x
(number of crossings) ~ distance.  This is the 'local + local + ...'
mechanism in calculable form.

Setup: 1+1 D, g = eta.  Temporal halo:
    phi_hat(t,x) = Phi0 * (1 + eps*(t - t0)/ell) * exp(-x^2/(2 ell^2))
eps = fractional field change per light-crossing time.

g~_ab = A^2 (eta_ab + b~ phi_a phi_b),  A = exp(-phi_hat),  a,b in {t,x}.

Photon: launched x=-X0 at t=0 toward +x; frequency measured by g-static
matter observers at emission and at x=+X0 reception:
    omega~ = -g~_ab k^a u~^b,  u~^a = (1/sqrt(-g~_tt), 0).

Per-crossing observables:
  z_cross = omega~_em/omega~_rec - 1   (transport redshift residual)
  delay   = arrival time - X0*2 - (Minkowski baseline)
Sweep eps to isolate the evolution-induced part (eps=0 must give
z ~ 0: the static-cancellation check).
"""
import numpy as np
import json, os

EPS_FD = 1e-4

def make_pair(Phi1, Phi2, eps, ell, x1=-20.0, x2=+20.0, tc=100.0):
    """Two overlapping halos, one decaying one growing -- breaks every
    symmetry (the realistic cosmological crossing: a merger hand-off)."""
    def phi(t, x):
        a1 = Phi1*(1 - eps*(t-tc)/ell); a2 = Phi2*(1 + eps*(t-tc)/ell)
        return (a1*np.exp(-(x-x1)**2/(2*ell*ell))
                + a2*np.exp(-(x-x2)**2/(2*ell*ell)))
    def dphi(t, x):
        g1 = np.exp(-(x-x1)**2/(2*ell*ell))
        g2 = np.exp(-(x-x2)**2/(2*ell*ell))
        a1 = Phi1*(1 - eps*(t-tc)/ell); a2 = Phi2*(1 + eps*(t-tc)/ell)
        pt = (-Phi1*(eps/ell)*g1 + Phi2*(eps/ell)*g2)
        px = a1*g1*(-(x-x1)/ell**2) + a2*g2*(-(x-x2)/ell**2)
        return pt, px
    return phi, dphi

def make_moving(Phi0, v, ell):
    """Translating halo: phi = Phi0 exp(-(x - v t)^2/2 ell^2)."""
    def phi(t, x):
        return Phi0*np.exp(-(x - v*t)**2/(2*ell*ell))
    def dphi(t, x):
        g = np.exp(-(x - v*t)**2/(2*ell*ell))
        return Phi0*g*(x - v*t)*v/ell**2, Phi0*g*(-(x - v*t)/ell**2)
    return phi, dphi

def make_halo(Phi0, eps, ell, t0, mode="amp"):
    """phi_hat(t,x), plus derivatives phi_t, phi_x.
    mode='amp'  : amplitude growth only (nearly a lapse rescaling --
                  expected to give ~no residual; control case).
    mode='collapse': well deepens AND narrows -- genuine shape change
                  (the real ISW analog; collapse is the irreversible
                  process that could bias the sign).
    """
    def width(t):
        return ell * (1.0 - eps*(t - t0)/ell) if mode == "collapse" else ell
    def phi(t, x):
        w = width(t)
        amp = Phi0 * (1 + eps*(t - t0)/ell)
        return amp * np.exp(-x*x/(2*w*w))
    def dphi(t, x):
        w = width(t)
        wd = -eps/ell if mode == "collapse" else 0.0
        amp = Phi0 * (1 + eps*(t - t0)/ell)
        ad = Phi0 * eps/ell
        g = np.exp(-x*x/(2*w*w))
        pt = g * (ad + amp * x*x*wd/w**3)
        px = amp * g * (-x/w**2)
        return pt, px
    return phi, dphi

def gtil_of(phi, dphi, b_t):
    def g(t, x):
        A2 = np.exp(-2*phi(t, x))
        pt, px = dphi(t, x)
        return A2 * np.array([[-1 + b_t*pt*pt, b_t*pt*px],
                              [b_t*pt*px,    1 + b_t*px*px]])
    return g

def christoffel(g, t, x):
    gg = g(t, x)
    gi = np.linalg.inv(gg)
    dgt = (g(t+EPS_FD, x) - g(t-EPS_FD, x))/(2*EPS_FD)
    dgx = (g(t, x+EPS_FD) - g(t, x-EPS_FD))/(2*EPS_FD)
    dg = np.array([dgt, dgx])   # dg[d, b, c] = d_d g_bc
    Gam = np.zeros((2, 2, 2))
    for a in range(2):
        for b in range(2):
            for c in range(2):
                Gam[a, b, c] = 0.5*np.sum(
                    gi[a]*(dg[b][:, c] + dg[c][:, b] - dg[:, b, c]))
    return Gam

def integrate(g, X, K, dl, x_stop):
    n = 0
    while (K[1] > 0 and X[1] < x_stop) or (K[1] < 0 and X[1] > x_stop):
        Gam = christoffel(g, *X)
        dX = K.copy()
        dK = -np.einsum("abc,b,c->a", Gam, K, K)
        Xm = X + 0.5*dl*dX; Km = K + 0.5*dl*dK
        Gam2 = christoffel(g, *Xm)
        dX2 = Km; dK2 = -np.einsum("abc,b,c->a", Gam2, Km, Km)
        X += dl*dX2; K += dl*dK2
        n += 1
        if n > 3000000:
            return X, K, False
    return X, K, True

def omega_matter(g, X, K):
    gg = g(*X)
    if gg[0, 0] >= 0:
        return np.nan
    ut = 1/np.sqrt(-gg[0, 0])
    return -(gg[0, 0]*K[0] + gg[0, 1]*K[1]) * ut

def cross(Phi0, eps, b_t, ell=5.0, X0=100.0, mode="amp"):
    t0 = X0            # halo evolution centred on transit time
    if mode == "moving":
        phi, dphi = make_moving(Phi0, eps, ell)   # eps = halo velocity v
    elif mode == "pair":
        phi, dphi = make_pair(Phi0, Phi0, eps, ell)
    else:
        phi, dphi = make_halo(Phi0, eps, ell, t0, mode=mode)
    # pair handled above
    g = gtil_of(phi, dphi, b_t)
    # launch: null direction in g~ at emitter
    gg = g(0.0, -X0)
    # solve g~_ab k^a k^b = 0 for k = (kt, +1)
    a, b, c = gg[0, 0], 2*gg[0, 1], gg[1, 1]
    disc = b*b - 4*a*c
    if disc < 0:
        return None
    kt = (-b + np.sqrt(disc))/(2*a)
    if kt < 0: kt = (-b - np.sqrt(disc))/(2*a)
    if kt < 0: return None
    K0 = np.array([kt, 1.0])
    # normalize: rescale so kt ~ 1 convention is irrelevant (ratios only)
    Xe, Ke = np.array([0.0, -X0]), K0
    Xr, Kr, ok = integrate(g, Xe, Ke, dl=0.002, x_stop=+X0)
    if not ok:
        return dict(error="no arrival")
    om_e = omega_matter(g, Xe, Ke)
    om_r = omega_matter(g, Xr, Kr)
    z = om_e/om_r - 1.0
    delay = Xr[0] - 2*X0
    # phi at the two events
    dp = phi(*Xr) - phi(*Xe)
    return dict(z=z, delay=delay, Dphi=dp, t_arr=Xr[0])

print("=== STEP 21: evolving temporal halo -- ISW-analog redshift ===\n")
print("  Halo: phi = Phi0 (1 + eps t/ell) exp(-x^2/2ell^2), ell=5, X0=100")
print("  b_t chosen per amplitude for P_max ~ 0.5\n")

print("--- eps sweep: amplitude-only (control) vs collapse (shape) ---")
print("  mode     Phi0   eps        z_cross         delay")
for mode in ["amp", "collapse", "moving", "pair"]:
    for Phi0 in [0.3, 0.5]:
        for eps in [0.0, 1e-3, 1e-2, 5e-2]:
            b_t = 0.5 / (Phi0/5.0)**2
            r = cross(Phi0, eps, b_t, mode=mode)
            if r and "error" not in r:
                print(f"  {mode:8s} {Phi0:4.1f}  {eps:8.1e}   "
                      f"{r['z']:+.6e}   {r['delay']:+.4e}")
            else:
                print(f"  {mode:8s} {Phi0:4.1f}  {eps:8.1e}   {r}")

print("\n--- interpretation ---")
print("  'amp': amplitude growth alone is (nearly) a lapse")
print("  reparametrization -> residual ~0 expected (control).")
print("  'collapse': shape evolves during transit -> ISW-type residual.")
print("  'moving'  : translating structure -- fully asymmetric case.")
print("  If z_collapse ~ eps with a fixed sign, an ensemble of")
print("  collapsing structures gives systematic z ~ distance.")

res = {"mechanism": "evolving-halo transport redshift (ISW analog)",
       "note": "z_cross vs eps; static limit must vanish"}
os.makedirs("results", exist_ok=True)
with open("results/step_21_evolving_halo.json", "w") as f:
    json.dump(res, f, indent=2)
print("\nwrote results/step_21_evolving_halo.json")
