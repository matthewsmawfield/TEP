#!/usr/bin/env python3
"""
step_19_decisive_closure.py
============================
The closure test for the shared-ray mechanism, in the class where:

  (a) the volume identity 1+Q = C A^{-6} is carried by scalar gradients
      (volume-normal congruence, canonical scalar -- Jakarta's premises);
  (b) the photon--graviton delay is bounded by the leading-order cone
      kinematics  dT ~ (1/2c) int b s^2 dl;
  (c) the scalar stress tensor gravitates normally (canonical T_mn).

No FLRW critical density is imported.  The energy bound comes from the
candidate's own geometry: the Hamiltonian constraint integrated over the
propagation region (enclosed mass vs Schwarzschild compactness).

Structure of the argument
-------------------------
1.  The delay needs small b~ = B M_Pl^2 / A^2 (small coupling).
    The identity needs large gradient:  R = b~ Y_hat ~ Q/(1+Q) ~ 1
    =>  gradient energy rho = (1/2) M_Pl^2 Y_hat ~ M_Pl^2/(2 b~).
    So  dT * rho  >=  M_Pl^2 (Dphi_hat)^2 / (4 c L)   -- b~-independent.

2.  Requiring dT < 1.7 s fixes rho from below:
        rho >= M_Pl^2 (Dphi_hat)^2 / (4 c L dT).
    The enclosed compactness of a uniform-density region of size L is
        C(L) = 2 G M(<L) / L = (8 pi / 3) rho L^2 / M_Pl^2
             = (2 pi / 3) (Dphi_hat)^2 L / (c dT).
    For GW170817 numbers this is ~ 10^12 -- the propagation region
    would lie ~ 10^12 Schwarzschild radii deep.  A static, traversable,
    eternal universe is incompatible with that.

3.  The uncancellable part: even if a tuned V(phi) cancels the isotropic
    density, the wave part of T_mn is a null-dust flux
        T_mn ~ M_Pl^2 f'^2 k_m k_n ,
    which (i) cannot be cancelled by any scalar potential, (ii) focuses
    every null congruence (R_mn k^m k^n > 0), and (iii) accumulates mass
    at the centre -- an eternal converging wave accumulates without bound.

4.  Standing-wave and lumpy variants: superimposing the outgoing wave
    zeroes the net flux but doubles the density; lumpiness raises the
    integrated gradient energy for fixed Dphi (Cauchy-Schwarz).  The
    exclusion is therefore a statement about the MECHANISM CLASS, not
    about one profile.

Remaining open gates (stated, not assumed):
  - tilted physical congruence: the identity A^6(1+Q)=C was derived for
    the zero-shift volume-normal network; a tilted galaxy congruence may
    change the required Q (and hence the energy demand);
  - beyond-leading-order cone kinematics on the solved metric;
  - a field whose gradient stress does not gravitate canonically would
    leave the premises of this test.
"""
import numpy as np
import json
import os

C_LIGHT = 2.99792458e8          # m/s
M_PL    = 1.22089e28            # m^-1   (reduced Planck mass / hbar c)
MPC     = 3.085677581e22        # m
G_SI    = 6.674e-11

z_s   = 0.0098
L     = 40.0 * MPC
DT_LIM = 1.7                    # s, GW170817/GRB170817A
DPHI  = np.log(1.0 + z_s)       # required endpoint-law field change

# ----------------------------------------------------------------------
# 1. The b~-independent product bound
# ----------------------------------------------------------------------
# dT >= b~ Dphi^2 / (2 c L)            (Cauchy-Schwarz, leading order)
# rho_grad = M_Pl^2 R / (2 b~) = M_Pl^2 Q / (2 b~ (1+Q))  >= M_Pl^2/(4 b~)  for Q>=1
# => dT * rho >= M_Pl^2 Dphi^2 Q / (4 c L (1+Q)) ~ M_Pl^2 Dphi^2/(4 c L)

def rho_min_for_delay(dT):
    """Minimum field energy density (kg/m^3 via mass density units:
    returns rho in m^-2 times M_Pl^2 -- i.e. rho / M_Pl^2 in m^-2,
    times M_Pl^2 gives natural-units density)."""
    # rho = M_Pl^2 Dphi^2 / (4 c dT L)   [units: M_Pl^2 * m^-2]
    return M_PL**2 * DPHI**2 / (4.0 * C_LIGHT * dT * L)

rho_req = rho_min_for_delay(DT_LIM)          # natural units rho (m^-4-ish)
kappa_req = np.sqrt(rho_req) / M_PL          # equivalent gradient scale
print("=== 1. delay-energy product bound ===")
print(f"  dT <= {DT_LIM} s  requires")
print(f"    rho_phi / M_Pl^2 >= {rho_req/M_PL**2:.3e} m^-2")
print(f"    gradient scale  kappa >= {kappa_req:.3e} m^-1"
      f"  =  ({1/kappa_req/(MPC/1e6):.0f} pc)^-1")

# ----------------------------------------------------------------------
# 2. Compactness of the propagation region (integrated Hamiltonian)
# ----------------------------------------------------------------------
def compactness(rho, r):
    # 2 G M(<r)/r with M = (4pi/3) rho r^3 ; rho = rho_natural * M_Pl^2
    # In units G=c=1: 2M/r = (8pi/3) (rho/M_Pl^2) r^2
    return (8.0*np.pi/3.0) * (rho / M_PL**2) * r**2

C_L = compactness(rho_req, L)
print("\n=== 2. enclosed compactness at source distance ===")
print(f"  C(L) = (8pi/3) rho L^2 / M_Pl^2 = {C_L:.3e}")
print(f"  -> propagation region is ~{C_L:.1e} Schwarzschild masses deep:")
print(f"     the light-travel region cannot be static and traversable.")

# Invert: what delay would make the region merely critical, C=1?
rho_crit_geom = M_PL**2 * 3.0/(8.0*np.pi*L**2)     # rho for C(L)=1
dT_needed = M_PL**2 * DPHI**2 / (4.0*C_LIGHT*L*rho_crit_geom)
print(f"\n  geometric critical density for a 40 Mpc sphere:")
print(f"    rho_C / M_Pl^2 = 3/(8 pi L^2) = {rho_crit_geom/M_PL**2:.3e} m^-2")
print(f"  delay compatible with sub-critical region:  dT >= {dT_needed:.3e} s")
print(f"    = {dT_needed/3.156e7:.3e} yr   (observed: 1.7 s)")
print(f"  tension factor: {dT_needed/DT_LIM:.2e}")

# scan: delay vs compactness -- the trade-off curve
print("\n  dT bound [s]   rho/M_Pl^2 [m^-2]   C(40 Mpc)")
for dT in [1.7, 1e3, 1e6, 1e9, 1e12, dT_needed]:
    r = rho_min_for_delay(dT)
    print(f"    {dT:10.2e}   {r/M_PL**2:14.3e}   {compactness(r, L):10.3e}")

# ----------------------------------------------------------------------
# 3. The uncancellable flux
# ----------------------------------------------------------------------
print("\n=== 3. can V(phi) rescue it? ===")
print("  T_mn = M_Pl^2 f'^2 k_m k_n - g_mn[ 1/2 (dphi)^2 + V ]")
print("  A potential cancels only the isotropic part.  The k_m k_n")
print("  null-dust flux term is untouched:")
print("    R_mn k'^m k'^n = M_Pl^-2 T_mn k'^m k'^n = M_Pl^-2 f'^2 (k.k')^2 > 0")
print("  -> every other null congruence is focused (Raychaudhuri),")
print("     and the inward flux accumulates mass at the centre:")
print("     dm/dv ~ 4 pi r^2 f'^2  -- unbounded for an eternal wave.")
# illustrative accumulation: mass accrued across wave zone over one L/c
f2 = rho_req / M_PL**2            # ~ f'^2 = Y_hat
dm_dt = 4.0*np.pi * L**2 * f2 * M_PL**2 * C_LIGHT   # ~ mass/time (kg/s, scaled)
print(f"     order-of-magnitude influx over region: dm/dt ~ M_Pl^2 f'^2 L^2 c")
print(f"     = {dm_dt:.2e} (scaled) -- diverges; no steady state.")

# ----------------------------------------------------------------------
# 4. variants: standing wave, lumpy profile
# ----------------------------------------------------------------------
print("\n=== 4. variants ===")
print("  standing wave (in+out): net flux 0 but rho -> 2 rho,"
      " compactness worse x2")
print(f"    C(L) = {2*C_L:.2e}")
print("  lumpy profile: for fixed Dphi, int (dphi)^2 >= (Dphi)^2/L")
print("    (Cauchy-Schwarz) -- lumpiness raises both delay and energy;")
print("    uniform is already the optimum.")
print("  tilted congruence: the ONLY route that could change the answer --")
print("    if the physical network is tilted, the identity need not be")
print("    A^6(1+Q)=C and the entire Q-requirement must be re-derived.")

# ----------------------------------------------------------------------
# 5. verdict
# ----------------------------------------------------------------------
print("\n=== VERDICT ===")
print(f"  In the class (canonical scalar, volume-normal congruence,")
print(f"  leading-order cones): multimessenger bound requires rho such")
print(f"  that the 40-Mpc region is ~{C_L:.0e}x Schwarzschild-compact.")
print(f"  Equivalently, a traversable region requires dT >= "
      f"{dT_needed:.1e} s >> 1.7 s.")
print(f"  EXCLUDED for this mechanism class -- by the candidate's own")
print(f"  geometry, with no imported FLRW scale.")
print(f"  Open gates: (i) tilted physical congruence changing the volume")
print(f"  identity; (ii) non-canonical scalar stress; (iii) nonlocal")
print(f"  transport redshift (forced anyway by <H> = theta/3 = 0).")

res = {
    "assumptions": ["canonical scalar T_mn", "volume-normal congruence",
                    "leading-order cone kinematics",
                    "endpoint-law Dphi = ln(1+z)"],
    "z": z_s, "L_m": L, "dT_limit_s": DT_LIM, "Dphi": DPHI,
    "rho_required_over_Mpl2_m^-2": rho_req / M_PL**2,
    "kappa_required_m^-1": kappa_req,
    "kappa_scale_pc": 1.0/kappa_req/MPC,
    "compactness_at_L": C_L,
    "geometric_critical_rho_over_Mpl2_m^-2": rho_crit_geom / M_PL**2,
    "dT_for_critical_s": dT_needed,
    "tension": dT_needed / DT_LIM,
    "verdict": "mechanism class excluded at ~1e12 by self-compactness",
    "open_gates": ["tilted congruence (volume identity form)",
                   "non-canonical scalar stress",
                   "finite-distance transport redshift"],
}
os.makedirs("results", exist_ok=True)
with open("results/step_19_decisive_closure.json", "w") as f:
    json.dump(res, f, indent=2)
print("\nwrote results/step_19_decisive_closure.json")
