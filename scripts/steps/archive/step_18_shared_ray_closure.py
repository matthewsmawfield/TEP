#!/usr/bin/env python3
"""
step_18_shared_ray_closure.py

Single-candidate closure test for the zero-volume-expansion cosmology.

Candidate geometry (geometric units c=1, phi_hat = phi/M_Pl dimensionless):

    phi_hat(t,r) = f(t+r) + D t          -- converging spherical wave + secular drift

    The pure wave part f(u), u = t+r, is a NULL gradient (Pi_hat = phi_hat_r = f'),
    which is the unique family in which:
      (i)  inward photons see a shared photon/GW null ray (k.phi = f'(1-v) ~ 0),
      (ii) all incoming directions are equivalent (isotropy for the central
           observer),
      (iii) the volume identity becomes an exact ODE:

            1+Q = 1/(1 - b f'^2) = e^{6f}   =>   f'^2 = (1 - e^{-6f})/b ,

      so phi_hat = f(t+r) with this f is an EXACT zero-volume field:
      A^6(1+Q) = 1 on every worldline (C = 1).  It carries zero redshift.
      The secular drift D is what must supply Delta phi = ln(1+z); it breaks
      exact nullity and therefore produces a photon-GW delay -- the two are
      traded through the same quantity.

Tests / gates on each candidate (b, D, f-amplitude):
  1. Volume-identity residual  d/dt ln[A^3 sqrt(1+Q)] on worldlines.
  2. Regularity: W = 1 - P + R > 0,  P < 1,  D_lapse = A^2 - B Pi^2 > 0.
  3. Inward-ray solve: exact quadratic for v(r), accumulate
        Delta phi_hat = int s dr ,  Delta T = int (1/v - 1) dr
        s = Pi + phi_r*(-v) = D + f'(1-v).
  4. Redshift match: |Delta phi_hat - ln(1+z_s)|.
  5. Multimessenger: Delta T <~ 1.7 s.
  6. Gravitational cost diagnostic: T_00 = Pi^2 + V,  T_0r = Pi*phi_r
     (directed flux, uncancellable by V) -- reported as field-strength scale
     kappa = sqrt(<Y>), not judged against rho_crit.
  7. Single-function check: b is const by construction -> B(phi) = b A^2/M_Pl^2
     passes the isosurface test identically (the question is whether the
     required b is compatible with the OTHER phenomenology constraints).
  8. Isotropy: displaced-observer ray test.

Sweep: (b, D) grid -> report the gate map and any viable island.
Output: results/step_18_shared_ray_closure.json
"""
import json
import numpy as np

C_LIGHT = 2.998e8      # m/s
YR = 365.25 * 86400.0
MPC = 3.086e22         # m

# ----------------------------------------------------------------------
# General exact ray speed (angle theta to grad phi)  -- as in step_16/17
# ----------------------------------------------------------------------

def v_ray(theta, P, R):
    c = np.cos(theta)
    return (-np.sqrt(np.abs(P * R)) * c
            + np.sqrt(np.clip(1.0 + R * c**2 - P, 0.0, None))) \
           / (1.0 + R * c**2)

# ----------------------------------------------------------------------
# Exact inward-ray speed for phi = f(t+r) + D t   (c=1)
# ----------------------------------------------------------------------
# s = D + f'(1-v);  null cone:  v^2 = 1 - b_t s^2  with  b_t = b/A^2 folded in.
# Let u = 1-v.  (1-u)^2 = 1 - b_t (D + f' u)^2
# => u^2(1 + b_t f'^2) + u(2 b_t D f' - 2) + b_t D^2 = 0
# physical root u ~ 0 for small D.

def solve_ray(b_t, D, fp):
    """Return (v, s, u) for inward ray; None if no real positive root."""
    a = 1.0 + b_t * fp**2
    bb = 2.0 * b_t * D * fp - 2.0
    cc = b_t * D**2
    disc = bb**2 - 4 * a * cc
    if disc < 0:
        return None
    u = (-bb - np.sqrt(disc)) / (2 * a)   # small root
    if not (0.0 <= u <= 1.0):
        return None
    v = 1.0 - u
    s = D + fp * u
    return v, s, u

# quick check: D=0 -> u=0 -> v=1 exactly (shared ray)
_v = solve_ray(1e6, 0.0, 1e-3)
assert abs(_v[0] - 1.0) < 1e-12

# ----------------------------------------------------------------------
# Volume-identity-compatible wave profile f(u)
# ----------------------------------------------------------------------
# f' = sqrt((1 - e^{-6f})/b), f(0)=0.  Integrate on a u-grid.

def build_f(b_t, ugrid):
    f = np.zeros_like(ugrid)
    du = ugrid[1] - ugrid[0]
    for i in range(1, len(ugrid)):
        rhs = np.sqrt(max(0.0, (1.0 - np.exp(-6 * f[i - 1])) / b_t))
        f[i] = f[i - 1] + du * rhs
    fp = np.gradient(f, ugrid)
    return f, fp

# ----------------------------------------------------------------------
# Gates evaluated on the full 2D field phi(t,r) = f(t+r) + D t
# ----------------------------------------------------------------------

def evaluate(b_t, D, fscale_u, label=""):
    """All quantities in SI-ish geometric units: lengths/times in metres (c=1
       via dividing seconds by c where needed).  u and r in metres; D, f' in 1/m."""
    # spatial domain: r in [0, L], time span a bit longer than L
    L = 40.0 * MPC
    nr, nt = 400, 400
    r = np.linspace(0, L, nr)
    t = np.linspace(-L, 0, nt)
    T, R = np.meshgrid(t, r, indexing="ij")
    u = T + R                                  # u = t + r (metres)
    # build f over the needed u-range
    umin = float(u.min())
    ug = np.linspace(max(0.0, umin - 10 * fscale_u), umin + u.max() - umin + 10 * fscale_u, 20000)
    f, fp = build_f(b_t, ug)
    fp_u = np.interp(u, ug, fp)                # f'(u) at each point, units 1/m
    f_u = np.interp(u, ug, f)
    phi = f_u + D * T
    Pi = D + fp_u                              # d_t phi
    Gr = fp_u                                  # d_r phi
    P = b_t * Pi**2
    Rr = b_t * Gr**2
    W = 1.0 - P + Rr
    Q = Rr / np.clip(1 - P, 1e-30, None)
    # volume identity residual: d/dt ln[A^3 sqrt(1+Q)] along worldlines (fixed r)
    A = np.exp(-phi)
    lnV = np.log(A**3 * np.sqrt(np.clip(1 + Q, 1e-30, None)))
    resid = np.gradient(lnV, t, axis=0)        # should be ~0 if identity holds
    # regularity
    reg_ok = bool(np.all(W > 0) and np.all(P < 1))
    # ---- inward ray: receiver at r=0,t=0; source at r=L; ray at fixed direction
    # ray crosses points (t(r) ~ -L + r, r).  Evaluate s, v along it.
    rr = np.linspace(0, L, 2000)
    tr = rr - L                                # arrival at t=0
    uu = tr + rr                               # = 2rr - L
    fp_r = np.interp(uu, ug, fp)
    vs = np.array([solve_ray(b_t, D, fpi) for fpi in fp_r])
    ok = np.array([x is not None for x in vs])
    if ok.all():
        vv = np.array([x[0] for x in vs]); ss = np.array([x[1] for x in vs])
        uu_ = np.array([x[2] for x in vs])
        dphi = float(np.trapezoid(ss, rr))     # Delta phi_hat along ray
        dT_m = float(np.trapezoid(1.0 / vv - 1.0, rr))  # metres of extra path
        dT_s = dT_m / C_LIGHT
    else:
        dphi, dT_s = np.nan, np.nan
    # energy diagnostics (field-strength scale, NOT rho_crit comparison)
    kappa = float(np.sqrt(np.mean(Gr**2)))     # mean |grad phi_hat| in 1/m
    flux_scale = float(np.mean(np.abs(Pi * Gr)))  # T_0r scale in 1/m^2 (x M_Pl^2)
    return dict(label=label, b_t=b_t, D=D,
                ident_resid_rms=float(np.sqrt(np.mean(resid**2))),
                regular=reg_ok,
                dphi_ray=dphi, delay_s=dT_s,
                kappa_per_m=kappa, flux_scale=flux_scale)

# ----------------------------------------------------------------------
# What D is required for the redshift?  D ~ ln(1+z)/L (since s ~ D on shared rays)
# ----------------------------------------------------------------------

z_s = 0.0098
L = 40.0 * MPC
D_req = np.log(1 + z_s) / L                  # 1/m
H0 = 2.2e-18
print(f"required secular drift D ~ ln(1+z)/L = {D_req:.3e} /m")
print(f"   (for reference H0/c = {H0/C_LIGHT:.3e} /m;  D/(H0/c) = {D_req*C_LIGHT/H0:.3f})")

# ----------------------------------------------------------------------
# SWEEP over (b_t, D)
# ----------------------------------------------------------------------

print("\n=== SWEEP (b_t, D): gates on the same field ===")
print(f"{'b_t [m^2]':>11} {'D/D_req':>8} {'ident.resid':>11} {'reg':>4} "
      f"{'dphi/ln(1+z)':>13} {'delay [s]':>10}")
grid = []
for b_t in [1e30, 1e34, 1e37, 1e40, 1e44]:
    for dfrac in [1.0, 10.0, 100.0]:
        r = evaluate(b_t, D_req * dfrac, fscale_u=L / 4)
        grid.append(r)
        dph = r["dphi_ray"] / np.log(1 + z_s) if np.isfinite(r["dphi_ray"]) else np.nan
        print(f"{r['b_t']:11.1e} {dfrac:8.0f} {r['ident_resid_rms']:11.3e} "
              f"{'ok' if r['regular'] else 'X':>4} {dph:13.3f} "
              f"{r['delay_s']:10.3e}")

# ----------------------------------------------------------------------
# The analytic frontier (fast, for interpretation)
# ----------------------------------------------------------------------

print("\n=== analytic frontier ===")
print("shared ray: u = 1-v ~ b_t D^2/(2(1+b_t f'^2))  (small-D limit)")
print("for b_t f'^2 ~ Q/(1+Q) ~ O(1):  delay ~ (b_t/2) D^2 L / c")
for b_t in [1e34, 1e37, 1e40]:
    dT = 0.5 * b_t * D_req**2 * L / C_LIGHT
    print(f"  b_t={b_t:9.1e} m^2:  DeltaT ~ {dT:9.2e} s   "
          f"(b_t^{1/2} ~ {np.sqrt(b_t)/3.086e16:8.1e} pc)")
print("kappa implied by Q=63 (b_t f'^2 ~ 1): kappa ~ 1/sqrt(b_t)")
for b_t in [1e34, 1e37, 1e40]:
    k = 1 / np.sqrt(b_t)
    print(f"  b_t={b_t:9.1e}: kappa={k:8.2e}/m ~ ({1/k/3.086e16:8.1e} pc)^-1")

# ----------------------------------------------------------------------
# displaced-observer isotropy test
# ----------------------------------------------------------------------

print("\n=== isotropy: displaced observer ===")
# field converges on r=0; an observer at r = d sees gradient still pointing at
# r=0, so only their rays along -r_hat get cancellation; O(1) anisotropy.
rng = np.random.default_rng(1)
n = rng.normal(size=(500, 3)); n /= np.linalg.norm(n, axis=1)[:, None]
for th in [0.0, np.pi / 3, np.pi / 2]:
    # cos angle between inward ray direction and -r_hat
    mu = n[:, 2] * np.cos(th)
    s = np.abs(1 + mu)  # s ~ f'(1 + cos angle) for transverse-ish geometry
    print(f"  misalignment {np.degrees(th):5.0f} deg: <s>/f' = {s.mean():.2f}, "
          f"spread {s.std():.2f}")

# ----------------------------------------------------------------------
# TEST H — the user's cancellation question, made precise
# ----------------------------------------------------------------------

print("\n=== TEST H: can spatial variation make the delay cancel? ===")
# For B>0:  v^2 - 1 = -(B/A^2) s^2 <= 0 at EVERY point and direction.
# Verify numerically over the full (P,R,theta) domain:
rng2 = np.random.default_rng(3)
Pv = rng2.uniform(0, 0.999, 20000)
Rv = rng2.uniform(0, 50, 20000)
thv = rng2.uniform(0, np.pi, 20000)
v_all = v_ray(thv, Pv, Rv)
print(f"  max v over random (P,R,theta) sample, B>0: {v_all.max():.10f}")
print("  -> B>0 photons are NEVER superluminal: delay is positive-definite,")
print("     no point-to-point cancellation of the delay integral is possible.")
print("     What DOES cancel is s itself (it integrates to the net Delta phi),")
print("     while the delay integrand ~ s^2 always accumulates.")
# Cauchy-Schwarz optimality: for fixed int s = Delta phi, int s^2 is minimized
# by UNIFORM s.  Oscillating field structure only increases the delay cost.
Ls = 40.0 * MPC
lg = np.linspace(0, Ls, 2000)
base = np.log(1 + z_s) / Ls
profiles = {
    "uniform": np.full_like(lg, base),
    "sinusoidal": base + 5 * base * np.sin(2 * np.pi * lg / (Ls / 7)),
    "sparse fronts (10%)": np.where((lg % (Ls / 10)) < 0.1 * Ls / 10,
                                    base * 10, 0.0),
}
for name, pr in profiles.items():
    dph = np.trapezoid(pr, lg)
    cost = np.trapezoid(pr**2, lg)
    print(f"    {name:22s}: int s = {dph:9.3e} (target {np.log(1+z_s):.3e}), "
          f"int s^2 = {cost:9.3e} (min {np.log(1+z_s)**2/Ls:.3e})")
print("  -> intermittency cannot beat the uniform bound: for fixed redshift,")
print("     lumpier s-profiles give LARGER delay.  The helpful kind of")
print("     structure is the near-null wave part (supplies Q at ~zero s),")
print("     not sign-alternating delay.")

# ----------------------------------------------------------------------
# TEST I — could B < 0 branches advance the photon and cancel net delay?
# ----------------------------------------------------------------------

print("\n=== TEST I: sign-changing B(phi)? ===")
print("  The volume identity SLAVES the sign of B to the field value:")
print("    Q = C e^{6 phi} - 1 = R/(1-P);  for P<1,  sign(B) = sign(C e^{6phi}-1).")
print("  Since B is one function of phi, its zero must lie at the same field")
print("  value on EVERY worldline -> forces uniform C = e^{-6 phi*}:")
print("    B(phi) > 0 for phi > phi*  (past: photons delayed)")
print("    B(phi) < 0 for phi < phi*  (post-crossing: photons ADVANCE)")
print("  The zero-crossing epoch is universal.  Net-delay cancellation is then")
print("  mathematically possible: a path crossing phi* samples both signs.")
print("  Costs: (a) tuned zero of B(phi); (b) B<0 gives superluminal photon")
print("  branches (causality questions; if B<0 TODAY, locally testable);")
print("  (c) the cancellation must hold on all source directions.")
# illustration: crossing at mid-path, |B| symmetric
bpos, bneg = 1e37, -1e37
s_uni = np.log(1 + z_s) / Ls
dT_pos = 0.5 * bpos * s_uni**2 * (Ls / 2) / C_LIGHT
dT_neg = 0.5 * bneg * s_uni**2 * (Ls / 2) / C_LIGHT
print(f"     illustration: B halves cancel -> net dT = {dT_pos+dT_neg:.2e} s")
print("     but each half still carries |dT| = "
      f"{abs(dT_pos):.2e} s of one-signed propagation effect.")
print("     STATUS: a genuine but tightly constrained escape route —")
print("     B(phi) sign structure is forced, not free.")

# ----------------------------------------------------------------------
# TEST J — the real knob: gradient scale vs delay at fixed redshift
# ----------------------------------------------------------------------

print("\n=== TEST J: delay vs wave-gradient scale (the actual trade-off) ===")
print("  delta ~ (b/2) D^2/(1 + b f'^2);  identity: b f'^2 ~ Q/(1+Q)")
print(f"{'kappa=(wave scale)^-1':>22} {'delay [s]':>10}")
for lam_pc in [1e4, 1e3, 1e2, 10.0, 1.0, 0.1]:
    lam = lam_pc * 3.086e16
    fp = 1.0 / lam
    b_req = 63.0 / (64.0 * fp**2)          # bf'^2 = Q/(1+Q), Q=63
    u = b_req * D_req**2 / (2 * (1 + b_req * fp**2))
    dT = u * Ls / C_LIGHT
    print(f"{lam_pc:18.0f} pc {dT:10.3e}")
print("  -> smaller-scale gradients (lumpier field) suppress the residual")
print("     delay.  Sub-~100-pc wave scales satisfy the 1.7 s bound;")
print("     the gravitational cost of that lumpiness is the open question.")

# ----------------------------------------------------------------------
# JSON
# ----------------------------------------------------------------------

result = {
    "exact_solution": {
        "form": "phi = f(t+r); f' = sqrt((1-e^{-6f})/b)",
        "property": "satisfies A^6(1+Q)=1 identically on every worldline; "
                    "pure null gradient; carries zero redshift",
    },
    "required_drift_per_m": D_req,
    "sweep": [{k: v for k, v in g.items() if k != "label"} for g in grid],
    "interpretation": (
        "The volume identity admits an exact null-gradient wave solution, but "
        "it is redshift-free (s=0 on inward rays).  Adding the drift D needed "
        "for redshift breaks nullity; the resulting delay scales as "
        "DeltaT ~ (b/2) D^2 L/c while the volume identity fixes b f'^2 ~ Q/(1+Q)."
        "Combining: DeltaT ~ (Q/(2(1+Q))) (D/f')^2 L/c.  Satisfying 1.7 s at "
        "z=0.0098 needs f' ~> (100-500 pc)^-1 -- i.e. cosmological gradients on "
        "sub-kpc scales, whose stress-energy must then be carried by the solved "
        "metric.  Whether Einstein's equations admit such structure while "
        "remaining consistent with measured local geometry is THE open "
        "existence question; this script does not import rho_crit to decide it."
    ),
}
out = "results/step_18_shared_ray_closure.json"
with open(out, "w") as fj:
    json.dump(result, fj, indent=2)
print("\nwrote", out)
