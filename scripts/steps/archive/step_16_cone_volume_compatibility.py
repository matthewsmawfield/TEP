#!/usr/bin/env python3
"""
step_16_cone_volume_compatibility.py

Cone-volume compatibility for the zero-volume-expansion cosmology (Jakarta v0.14, Sec. 8).

IMPORTANT FRAMING (revised after review):
  - The multimessenger observable is the difference between two null ARRIVAL EVENTS
    on the receiver's worldline (photon on g~, GW on g), expressed in the receiver's
    matter-frame proper time.  It is NOT an angular average of coordinate speeds.
  - Local matter-frame observers always measure c (local Lorentz invariance).  Nothing
    here claims otherwise.
  - The disformal term vanishes on rays with k^m phi_m = 0: a single ray can be null
    in BOTH metrics even when the full cones differ enormously.  Hence large Q does
    not force a large delay along every ray; it forces a difference somewhere.

Setup (local gravitational inertial frame, N=1, M_Pl=1):
    Pi = d_t phi,  Y = |grad phi|^2,  P = B Pi^2/A^2,  R = B Y/A^2,
    W = 1 - P + R > 0 (invertibility),  P < 1 (reference observer timelike).
    Q = R/(1-P);  zero volume expansion => 1+Q = C(x)[A/A_0(x)]^{-6}, C(x) = 1+Q_0.

Exact ray speed at angle theta to grad phi (g~_mn k^m k^n = 0, k^m=(w,k n^)):
    v = ( -sqrt(PR) cos th + sqrt(1 + R cos^2 th - P) ) / (1 + R cos^2 th)

Key ray quantity:  s = k^m phi_m / w = Pi + sqrt(Y) v cos th
    - speed deficit:      1 - v ~ (B/2A^2) s^2        (for s small)
    - field change/length: d phi_hat / dl = s / v ~ s  (phi_hat = phi/M_Pl)
    - so DeltaT = int (B/2A^2c) s^2 dl  and  Delta phi_hat = int s dl,
      giving the exact joint bound (Cauchy-Schwarz):
          DeltaT >= (B/2A^2c) (Delta phi_hat)^2 / L ,
      saturated iff s(l) is constant along the ray.

Tests:
  V.  closed-form checks incl. the P=R=63/64 null-gradient counterexample:
      Q=63, W=1, v(pi)=c exactly (shared null direction).
  A.  for required Q, per-direction deficit map and the minimum-direction deficit;
      existence of a zero-deficit direction iff R >= P (Y >= Pi^2, spacelike/null grad).
  B.  the joint frontier: minimum possible delay vs required endpoint redshift
      Delta phi_hat = ln(1+z); numerical bound vs GW170817; and the field condition
      needed to saturate it (near-null gradient, ray-aligned, all along the path).
  C.  transport-dominated redshift: shown to face the same joint integral (model
      dependent, not an exclusion).
  D.  tilted congruence: status flag only - the simple divergence argument is
      INCOMPLETE (needs Lorentz factors, acceleration terms); not demonstrated.
  E.  nonuniform C(x): conditional impossibility under B>=0 and positive drift sign.
  F.  single-candidate isosurface test - one trial field fails; not an exclusion.
  G.  B<0: signature alone does NOT exclude it (W>0 possible); same joint frontier
      applies with sign-flipped deformation (superluminal branch, |DeltaT| bounded
      the same way).

Output: results/step_16_cone_volume_compatibility.json
"""
import json
import numpy as np

C_LIGHT = 2.998e8          # m/s
YR = 365.25 * 24 * 3600.0
MPC = 3.086e22             # m

# ----------------------------------------------------------------------
# Exact ray speed + verification
# ----------------------------------------------------------------------

def v_ray(theta, P, R):
    c = np.cos(theta)
    return (-np.sqrt(P * R) * c + np.sqrt(np.clip(1.0 + R * c**2 - P, 0.0, None))) \
           / (1.0 + R * c**2)

def verify():
    ok = True
    ok &= abs(v_ray(np.pi / 2, 0.4, 2.0) - np.sqrt(0.6)) < 1e-12      # v_perp = sqrt(1-P)
    ok &= abs(v_ray(0.0, 0.0, 63.0) - 1.0 / 8.0) < 1e-12              # v_par = 1/sqrt(1+R)
    ok &= abs(v_ray(0.7, 0.0, 0.0) - 1.0) < 1e-12                     # B=0 -> c
    # reviewer counterexample: P = R = 63/64, Q = 63, v(pi) = c exactly
    P0 = R0 = 63.0 / 64.0
    Q0 = R0 / (1.0 - P0)
    ok &= abs(Q0 - 63.0) < 1e-12
    ok &= abs(v_ray(np.pi, P0, R0) - 1.0) < 1e-10                     # shared null ray
    ok &= abs(v_ray(0.0, P0, R0) - 1.0 / 127.0) < 1e-10               # v+ = c/127
    ok &= abs(v_ray(np.pi / 2, P0, R0) - np.sqrt(1 - P0)) < 1e-10     # v_perp = c/8
    W0 = 1.0 - P0 + R0
    ok &= abs(W0 - 1.0) < 1e-12                                       # invertible
    return ok

print("verification incl. null-gradient counterexample:",
      "PASS" if verify() else "FAIL")

TH = np.linspace(0.0, np.pi, 2001)
WGT = 0.5 * np.sin(TH)

# ----------------------------------------------------------------------
# TEST A — per-direction deficit: minimum vs average
# ----------------------------------------------------------------------

print("\n=== TEST A: per-direction deficit at required Q (minimum is the observable) ===")
print(f"{'drift z_d':>10} {'Q_min':>10} {'P*':>6} {'min|d|':>10} {'<d>':>10} {'worst':>8}  shared-ray?")
rowsA = []
for zd in [0.0098, 0.1, 0.5, 1.0, 3.0]:
    Q = (1 + zd) ** 6 - 1
    # scan splits; for each compute min-direction and average deficit
    best = None
    for P in np.linspace(0.0, 0.9995, 400):
        R = Q * (1 - P)
        v = v_ray(TH, P, R)
        d_min = float(np.min(1.0 - v))
        d_avg = float(np.trapezoid((1 - v) * WGT, TH))
        d_max = float(np.max(1.0 - v))
        cand = (d_min, d_avg, d_max, P)
        if best is None or cand[0] < best[0]:
            best = cand
    shared = "yes" if best[0] < 1e-8 else "no"
    rowsA.append(dict(z=zd, Q=Q, P_star=best[3], min_dir=best[0],
                      avg=best[1], worst=best[2]))
    print(f"{zd:10.4f} {Q:10.3f} {best[3]:6.3f} {best[0]:10.2e} "
          f"{best[1]:10.2e} {best[2]:8.4f}  {shared}")
print("  -> for R >= P (spatial/null gradient) a zero-deficit ray direction exists:")
print("     theta0 = arccos(-Pi/sqrtY).  Single multimessenger rays can evade")
print("     the cone deformation entirely; angular averages do not apply.")

# ----------------------------------------------------------------------
# TEST B — the joint delay-redshift frontier
# ----------------------------------------------------------------------

print("\n=== TEST B: joint frontier  DeltaT >= (B/2A^2 c) (Delta phi)^2 / L ===")
z_s = 0.0098
L = 40.0 * MPC
dphi = np.log(1.0 + z_s)                    # endpoint-law phi difference
Lt = L / C_LIGHT                            # 130 Myr
for boverA2 in [1.0, 1e10, 1e20, 1e23]:
    dT_min = boverA2 * dphi**2 / (2.0 * Lt)
    print(f"  B/A^2 = {boverA2:8.1e}:  DeltaT_min = {dT_min:9.2e} s   "
          f"(bound ~1.7 s)  {'OK' if dT_min < 1.7 else 'EXCLUDED'}")
print("  Saturation condition: s(l) = d phi/dl = const = ln(1+z)/L along the ray.")
print("  This requires R >= P (spacelike/null gradient) AND local grad-phi")
print("  oriented near theta0 = arccos(-Pi/sqrtY) along the WHOLE path - i.e. a")
print("  field whose gradient tracks the ray.  Existence for an isotropic sky of")
print("  sources is an open geometric constraint, not a demonstrated obstruction.")

# generic-orientation estimate for comparison (what a random field would give)
Q01 = (1 + z_s) ** 6 - 1
v_gen = v_ray(TH, 0.0, Q01)
dT_gen = np.trapezoid((1 / v_gen - 1) * WGT, TH) * Lt
print(f"  generic-orientation delay (for contrast): {dT_gen:9.2e} s - excluded.")

# ----------------------------------------------------------------------
# TEST C — transport-dominated redshift (status, not exclusion)
# ----------------------------------------------------------------------

print("\n=== TEST C: transport-dominated redshift ===")
print("  frequency shift accumulated on the ray and the delay integral share the")
print("  same deformation functional s(l); the joint frontier applies identically.")
print("  STATUS: model-dependent; NOT an exclusion.")

# ----------------------------------------------------------------------
# TEST D — tilted congruence: honest status
# ----------------------------------------------------------------------

print("\n=== TEST D: tilted congruence ===")
print("  theta~(u) = theta~(n) + (spatial divergence) + Lorentz/acceleration terms:")
print("  the simplified divergence argument is INCOMPLETE - a tilted physical")
print("  network is not ruled out by this script.  Requires the full congruence")
print("  expansion computed on a candidate solution.  STATUS: open, not closed.")

# ----------------------------------------------------------------------
# TEST E — nonuniform C(x): conditional check (B>=0, positive drift)
# ----------------------------------------------------------------------

print("\n=== TEST E: nonuniform C(x), under B>=0 and A_0/A(t') > 1 ===")
ok_E = False
for A0x in np.linspace(0.3, 1.0, 8):
    for drift in np.linspace(1.001, 4.0, 20):
        At = A0x / drift
        Q0_needed = (At / A0x) ** 6 - 1.0        # = drift^{-6} - 1 < 0
        ok_E |= (Q0_needed >= -1e-12)
print(f"  can a drifting region have Q=0 in the past AND Q_0>=0 today?  {ok_E}")
print("  -> under these sign assumptions: no.  For B<0 or reversed drift the")
print("     argument does not apply (see Test G).")

# ----------------------------------------------------------------------
# TEST F — single-candidate isosurface test
# ----------------------------------------------------------------------

print("\n=== TEST F: single-candidate B(phi) reconstruction ===")
Nt, Nxg = 600, 400
tg = np.linspace(0, 10, Nt); xg = np.linspace(0, 2 * np.pi, Nxg)
T, X = np.meshgrid(tg, xg, indexing="ij")
phi = (12.0 - 0.5 * T) + 0.3 * np.sin(2 * X - 0.4 * T)
Pi = np.gradient(phi, tg, axis=0); Y = np.gradient(phi, xg, axis=1) ** 2
Af = np.exp(-phi)
Qf = np.clip(Af ** -6 - 1.0, 1e-12, None)
B_req = Af ** 2 * Qf / (np.clip(Y, 1e-8, None) + Qf * Pi ** 2)
order = np.argsort(phi.ravel()); ph = phi.ravel()[order]; Br = B_req.ravel()[order]
nb = 60; bins = np.linspace(ph.min(), ph.max(), nb + 1); spread = []
for i in range(nb):
    m = (ph >= bins[i]) & (ph < bins[i + 1])
    if m.sum() > 5:
        spread.append(float(np.std(Br[m]) / (abs(np.mean(Br[m])) + 1e-30)))
print(f"  median fractional spread of B_req on phi-isosurfaces: {np.median(spread):.3f}")
print("  STATUS: this one trial field fails; does not exclude other configurations.")

# ----------------------------------------------------------------------
# TEST G — B < 0 branch
# ----------------------------------------------------------------------

print("\n=== TEST G: B < 0 ===")
print("  W = 1 - P + R > 0 can hold with B < 0 (P,R < 0): signature not violated")
print("  by the sign alone.  The deformation flips sign (v > 1, superluminal")
print("  branch); the SAME joint bound applies with |s|.  Exclusion must come")
print("  from stability/observation, not the sign.  STATUS: open.")

# ----------------------------------------------------------------------
# Verdict
# ----------------------------------------------------------------------

verdict = {
    "formulas_verified": bool(verify()),
    "key_correction": (
        "The multimessenger observable is a single-ray arrival-time difference, "
        "not an angular average.  For R >= P there exists a direction "
        "theta0 = arccos(-Pi/sqrtY) on which k^m phi_m = 0: the ray is null in "
        "both metrics.  Verified: P = R = 63/64 gives Q = 63 with v(pi) = c "
        "exactly.  A universal delay bound from Q alone is FALSE."
    ),
    "tradeoff": (
        "On an exactly shared ray, Delta phi = int k.phi = 0: no conformal "
        "endpoint redshift.  The joint constraint is "
        "DeltaT >= (B/2A^2c)(ln(1+z))^2/L - numerically ~1e-20 (B/A^2) s for "
        "z=0.0098, L=40 Mpc - compatible with the 1.7 s bound for a wide range "
        "of B.  Saturation needs s(l) uniform along the ray, i.e. a near-null "
        "field gradient oriented near theta0 for the whole path."
    ),
    "open_conditions": [
        "a field with R >= P (spacelike/null gradient) along observed rays",
        "local grad-phi near the cancellation direction along each observed path",
        "the same for rays from all sky directions (isotropy constraint on the "
        "gradient field - the hardest condition)",
        "the volume identity 1+Q = C(x)[A/A_0]^{-6} simultaneously satisfied",
        "single B(phi): isosurface test on the actual solution",
    ],
    "tests": {
        "A_min_direction_deficit": rowsA,
        "B_joint_bound_s_for_B_over_A2_1": 1.16e-20,
        "C_transport": "same joint functional; not excluded",
        "D_tilted_congruence": "argument incomplete; open",
        "E_nonuniform_C": "impossible under B>=0 and positive drift only",
        "F_isosurface_single_candidate": float(np.median(spread)),
        "G_negative_B": "not excluded by signature; same joint bound",
    },
    "conclusion": (
        "REVISED: the earlier universal-exclusion claim is withdrawn.  Large Q "
        "does not force large delay along every ray; for spatial/null-gradient "
        "fields a shared photon-GW null direction exists pointwise.  The real "
        "question is the joint one: can a single field configuration deliver "
        "Delta phi = ln(1+z) along each observed ray while keeping s(l) = "
        "k.phi/w small everywhere - i.e. near-null gradients aligned with all "
        "observed rays?  The bound DeltaT_min ~ 1e-20 (B/A^2) s shows the "
        "constraint is parametrically satisfiable; viability hinges on whether "
        "a real Einstein-scalar solution supplies the required gradient "
        "organisation isotropically.  The decisive next step is the coupled "
        "calculation on one actual field configuration: both channels "
        "propagated, delay AND redshift measured on the same rays, for "
        "several source directions."
    ),
}

out = "results/step_16_cone_volume_compatibility.json"
with open(out, "w") as f:
    json.dump(verdict, f, indent=2)
print("\nwrote", out)
