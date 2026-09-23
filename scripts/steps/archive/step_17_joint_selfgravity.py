#!/usr/bin/env python3
"""
step_17_joint_selfgravity.py

Joint self-gravity / redshift / multimessenger closure test (Jakarta v0.14, Sec. 8).

After step_16's correction: a single multimessenger event constrains ONE ray, and
for R >= P (spacelike/null gradient) a shared-null direction exists on which the
photon is undeformed.  This step tests whether that mechanism can also produce
redshift, given the field's gravitational self-energy.

Setup (units explicit; phi_hat = phi/M_Pl dimensionless, b = B M_Pl^2/A^2 has
dimensions of length^2):

  Along a ray, s = k.phi_hat / w = Pi_hat + sqrt(Y_hat) v cos(theta_hat):
      fractional speed offset  |1 - v| ~ (b/2) s^2      (leading order, weak def.)
      field change per length  d phi_hat / dl = s / v ~ s

      =>  DeltaT = int (b(l)/2c) s(l)^2 dl,
          DeltaT >= (Delta phi_hat)^2 / (2c int dl/b(l))   (Cauchy-Schwarz)
          equality iff s = const = Delta phi_hat / L.

  Endpoint-law redshift needs  Delta phi_hat = ln(1+z),  hence mean drift rate
      s_drift = ln(1+z)/L = z/L ~ H0/c   (low-z Hubble law: this IS H0/c).

  Null-gradient branch (P = R = Q/(1+Q), the minimum-energy family with a shared
  direction):  Y_hat = Pi_hat^2 = Q / ((1+Q) b).
  Scalar kinetic density   rho_kin = M_Pl^2 Pi_hat^2 = M_Pl^2 Q / ((1+Q) b).

  HOW LARGE CAN THE GRADIENT BE?  The permissible scalar energy density is set
  by the curvature the SOLVED geometry can carry while remaining consistent
  with measured weak-field gravity -- it is NOT set by rho_crit = 3 M_Pl^2 H0^2,
  which imports the expansion rate under test (circular).  Parameterize by the
  allowed cosmological gradient scale kappa:  sqrt(Y_hat) <= kappa.
      b >= b_min ~ (Q/(1+Q)) / kappa^2.

  The delay bound becomes conditional:
      DeltaT_min = b_min (Delta phi_hat)^2 / (2cL)
               ~ (Q/(2(1+Q))) (L/c) (s_drift / kappa)^2
  with s_drift = ln(1+z)/L ~ z/L.  Two regimes:
    - kappa ~ H0/c (gradients at the Hubble scale): s_drift/kappa ~ 1, no
      cancellation available, DeltaT_min ~ (Q/2(1+Q)) L/c ~ 1e14 s -- excluded.
    - kappa >> H0/c (cosmological gradients on much shorter scales, e.g. a
      lumpy intergalactic field): cancellation opens, DeltaT_min falls as
      1/kappa^2.  Under 1.7 s requires kappa ~> 1e7 H0/c ~ (sub-pc)^-1 scale
      gradients -- only a solved field can say whether such structure carries
      acceptable stress-energy.

  Combining:
      DeltaT_min = b_min (Delta phi_hat)^2 / (2cL)
               ~ (Q/(2(1+Q))) * (L/c) * (s_drift^2 / (H0/c)^2)^{-1}-correction...

  Cleaner form (used below):
      delta = (b/2A^2) s^2,   s ~ s_drift,  b Y_hat = Q/(1+Q) A^2
      =>  delta_min ~ (Q/(2(1+Q))) * (s_drift^2 / Y_hat_max)
      with  Y_hat_max ~ (H0/c)^2  and  s_drift ~ H0/c:
      =>  delta_min ~ Q/(2(1+Q))  -- O(Q), unavoidable.

  Equivalently: the delay bound times the field energy is a fixed product:
      DeltaT * rho_kin >= (Delta phi_hat)^2 M_Pl^2 Q / (2cL(1+Q)).

Tests:
  V.  cone-limit checks incl. null-gradient counterexample (P=R=63/64, v(pi)=c).
  A.  drift-gradient feasibility ratio r = s_drift / sqrt(Y_hat_max).
  B.  the b-independent product bound and DeltaT_min at rho_kin = rho_crit.
  C.  momentum-flux obstruction: T_0i = Pi d_i phi ~ directed flux; single-valued
      smooth field cannot cancel it pointwise; G_0i must respond -> static
      zero-shift ansatz inconsistent for a coherent directed field.
  D.  isotropy: spherical converging null wave phi = f(t+r) gives a shared ray
      for ALL incoming directions -- but only for the central observer.
      Superposition for many centres breaks the null-gradient property.
  E.  remaining open loopholes, stated precisely.

Output: results/step_17_joint_selfgravity.json
"""
import json
import numpy as np

C_LIGHT = 2.998e8          # m/s
YR = 365.25 * 24 * 3600.0
MPC = 3.086e22             # m
H0 = 2.2e-18               # s^-1 (~70 km/s/Mpc)
HUBBLE_LEN = C_LIGHT / H0  # c/H0 ~ 1.36e26 m

# ----------------------------------------------------------------------
# Cone verification (as step_16)
# ----------------------------------------------------------------------

def v_ray(theta, P, R):
    c = np.cos(theta)
    return (-np.sqrt(P * R) * c + np.sqrt(np.clip(1.0 + R * c**2 - P, 0.0, None))) \
           / (1.0 + R * c**2)

def verify():
    P0 = R0 = 63.0 / 64.0
    ok = abs(R0 / (1.0 - P0) - 63.0) < 1e-12
    ok &= abs(v_ray(np.pi, P0, R0) - 1.0) < 1e-10
    ok &= abs(v_ray(0.0, P0, R0) - 1.0 / 127.0) < 1e-10
    ok &= abs(v_ray(np.pi / 2, P0, R0) - np.sqrt(1 - P0)) < 1e-10
    ok &= abs(v_ray(0.0, 0.0, 63.0) - 1.0 / 8.0) < 1e-12
    return ok

print("cone verification incl. null-gradient counterexample:",
      "PASS" if verify() else "FAIL")

# ----------------------------------------------------------------------
# TEST A — drift rate vs maximum gradient (the feasibility ratio)
# ----------------------------------------------------------------------

print("\n=== TEST A: s_drift vs allowed gradient scale kappa ===")
# The allowed gradient is NOT rho_crit (that would import the expansion rate
# under test).  Parameterize by kappa = sqrt(Y_hat)_max, to be fixed by a
# solved geometry.  Reference values for orientation only:
#   Hubble-scale gradient:      kappa ~ H0/c ~ 7.3e-27 /m
#   galaxy-scale lumpiness:     kappa ~ (Mpc)^-1 ~ 3.2e-23 /m
#   sub-pc-scale lumpiness:     kappa ~ (0.01 pc)^-1 ~ 3.2e-18 /m
results = []
for z_s, L_Mpc in [(0.0098, 40.0), (0.1, 430.0), (1.0, 3400.0)]:
    L = L_Mpc * MPC
    s_drift = np.log(1 + z_s) / L           # m^-1
    Q = (1 + z_s) ** 6 - 1
    results.append(dict(z=z_s, L_Mpc=L_Mpc, s_drift=s_drift, Q=Q, L=L))
    print(f"  z={z_s:6.4f}: s_drift={s_drift:.2e}/m")
print("  Cancellation requires s_drift << kappa; redshift needs s ~ s_drift")
print("  integrated along the ray.  Whether kappa >> H0/c is consistent with")
print("  the field's stress-energy is a question for a solved geometry.")

# ----------------------------------------------------------------------
# TEST B — the b-independent delay bound
# ----------------------------------------------------------------------

print("\n=== TEST B: delay bound as a function of the allowed gradient scale ===")
# DeltaT_min ~ (Q/(2(1+Q))) (L/c) (s_drift/kappa)^2   [leading order, endpoint law]
bound_s = 1.7
kappas = {"Hubble scale (H0/c)": H0 / C_LIGHT,
          "galaxy scale (1/Mpc)": 1.0 / MPC,
          "0.01-pc scale": 1.0 / (0.01 * 3.086e16)}
for row in results:
    z_s, L, Q = row["z"], row["L"], row["Q"]
    line = f"  z={z_s:6.4f}: "
    for name, k in kappas.items():
        dT = 0.5 * Q / (1 + Q) * (L / C_LIGHT) * (row["s_drift"] / k) ** 2
        line += f"{name}: {dT:9.2e} s  "
    print(line)
print("\n  The multimessenger bound is satisfied only if the cosmological field")
print("  carries gradients at kappa ~> 1e7 H0/c (sub-pc-scale lumpiness) AND")
print("  the gradient field stays ray-aligned.  Whether Einstein's equations")
print("  permit that structure is the open existence question — NOT assumed here.")

for row in results:
    z_s, L, Q = row["z"], row["L"], row["Q"]
    # kappa needed to reach the bound
    k_req = row["s_drift"] * np.sqrt(0.5 * Q / (1 + Q) * (L / C_LIGHT) / bound_s)
    row["kappa_required"] = k_req
    print(f"  z={z_s:6.4f}: kappa needed for DeltaT<1.7 s: {k_req:.2e}/m "
          f"(~{1/k_req/3.086e16:.1e} pc wavelength)")

# ----------------------------------------------------------------------
# TEST C — momentum-flux obstruction
# ----------------------------------------------------------------------

print("\n=== TEST C: directed momentum flux (cannot be canceled by V) ===")
print("  Null-gradient field: T_mn = phi_m phi_n - g_mn V.")
print("  T_0i = Pi * d_i phi ~ Pi^2 * direction — a directed energy flux.")
print("  - V cannot cancel it (V enters only the trace part).")
print("  - A single-valued smooth field has ONE gradient direction per point:")
print("    pointwise cancellation impossible; G_0i must carry the flux")
print("    (gravitomagnetic/moving geometry) — inconsistent with the static")
print("    zero-shift ansatz used to derive the volume identity.")
print("  - Only an isotropic lump-superposition averages the flux to zero,")
print("    leaving radiation-like stress (p = rho/3) — and destroys the")
print("    ray-alignment the shared-null mechanism needs (each point cancels")
print("    only one ray direction).")

# ----------------------------------------------------------------------
# TEST D — isotropy: spherical converging-wave ansatz
# ----------------------------------------------------------------------

print("\n=== TEST D: isotropy and the Copernican problem ===")
# phi = f(t+r):  Pi = f', grad = f' r_hat -> null gradient, shared ray for ALL
# inward photons — but only for the observer at r=0.
# For a second observer displaced by vector d, the local gradient is still
# r_hat (pointing to us): their incoming rays see s ~ f'(1 - r_hat.v_hat'),
# direction-dependent -> anisotropic redshift of order the gradient itself.
# Numeric check: two observers, 100 random sky directions each.
f = 1.0
rng = np.random.default_rng(0)
dirs = rng.normal(size=(200, 3)); dirs /= np.linalg.norm(dirs, axis=1)[:, None]
# central observer: all rays inward -> s = f'(1 - v) ~ 0 for v~1
# displaced observer at d: ray from sky dir n arrives along v = -n; gradient
# still r_hat -> s = f'(1 - n.r_hat) + ... ; choose r_hat = z for displaced obs
r_hat = np.array([0.0, 0.0, 1.0])
s_disp = f * (1 + dirs @ r_hat)         # O(f) spread, mean f — anisotropic
print(f"  displaced observer: <s>/f' = {s_disp.mean():.2f}, "
      f"spread = {s_disp.std():.2f} — O(1) anisotropy, not a wash-out.")
print("  Only the central observer sees isotropic transport.  A universal")
print("  Hubble law from this mechanism requires a converging-wave structure")
print("  centred on every observer — impossible for one smooth field; the")
print("  'dynamic inhomogeneous network' version is an open existence question.")

# ----------------------------------------------------------------------
# TEST E — remaining open loopholes (precise statements)
# ----------------------------------------------------------------------

print("\n=== TEST E: what remains genuinely open ===")
opens = [
    "tilted congruence: if the physical galaxy network is not the zero-shift "
    "volume-normal congruence, the identity 1+Q = C A^{-6} changes form — the "
    "Q requirement itself is congruence-dependent (full expansion formula not "
    "yet derived on a candidate solution)",
    "B < 0: signature-safe (W>0 possible); flips deformation sign (superluminal "
    "branch) but |DeltaT| bound identical in magnitude — no escape shown",
    "non-endpoint redshift: transport redshift shares the same s(l) functional "
    "— delay and shift are the same integral; moving-endpoint Doppler is "
    "bounded by peculiar velocities and cannot produce a smooth Hubble law",
    "kappa >> H0/c: the needed gradient scale is set by the solved field, not "
    "by rho_crit; whether sub-pc-scale cosmological gradients are consistent "
    "with the Einstein equations is the open existence question",
]
for o in opens:
    print(f"  - {o}")

# ----------------------------------------------------------------------
# Verdict + JSON
# ----------------------------------------------------------------------

verdict = {
    "verified": bool(verify()),
    "test_A_feasibility": [
        {k: v for k, v in r.items()} for r in results
    ],
    "main_result": {
        "statement": (
            "Conditional bound (leading order, endpoint law assumed): "
            "DeltaT_min ~ (Q/(2(1+Q)))(L/c)(s_drift/kappa)^2 where kappa is the "
            "largest cosmological field gradient the solved geometry permits. "
            "Hubble-scale gradients give DeltaT ~ 1e14 s (excluded); "
            "kappa ~> 1e7 H0/c (sub-pc lumpiness) satisfies the 1.7 s bound "
            "parametrically — IF such gradients are gravitationally consistent "
            "and ray-aligned.  The previous unconditional ~1e13 exclusion is "
            "WITHDRAWN: it imported rho_crit = 3 M_Pl^2 H0^2, which is defined "
            "from the expansion rate under test."
        ),
        "assumptions_flagged": [
            "endpoint-law redshift Delta phi_hat = ln(1+z) (alternatives share "
            "the same integral)",
            "volume identity on the zero-shift volume-normal congruence "
            "(tilted congruence open)",
            "leading-order weak-deformation delay formula",
            "B >= 0 (B<0 same magnitude bound)",
        ],
        "kappa_required_per_z": [
            {"z": r["z"], "kappa_per_m": r["kappa_required"],
             "wavelength_pc": 1.0 / r["kappa_required"] / 3.086e16}
            for r in results
        ],
    },
    "test_C_momentum_flux": (
        "T_0i ~ Pi^2 direction is a directed flux V cannot cancel; needs "
        "gravitomagnetic metric response inconsistent with the static ansatz, "
        "or an isotropic lump average that destroys ray alignment."
    ),
    "test_D_isotropy": (
        "Spherical converging null wave solves ray alignment for one observer "
        "only — O(1) anisotropy for displaced observers.  Universal version "
        "requires per-observer converging structure — open existence question."
    ),
    "open_loopholes": opens,
    "conclusion": (
        "REVISED (v3): the shared-null-direction mechanism evades delay only "
        "for rays carrying no field change; redshift requires s ~ s_drift ~ "
        "z/L.  Whether cancellation is available depends on the allowed "
        "cosmological gradient scale kappa — which must come from a solved "
        "geometry, not an imported critical density.  If the field is lumpy "
        "on sub-pc scales with ray-aligned gradients, the bound is satisfied "
        "parametrically; if gradients are Hubble-scale, it fails by ~1e13. "
        "The multimessenger question is therefore NOT decided by inequalities "
        "— it requires the single coupled solution test (Einstein + scalar "
        "equations, both ray channels, redshift and delay on the same rays)."
    ),
}

out = "results/step_17_joint_selfgravity.json"
with open(out, "w") as fjson:
    json.dump(verdict, fjson, indent=2)
print("\nwrote", out)
