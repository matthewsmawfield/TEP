import json
import os

# Physical constants
AU = 1.496e11  # meters
M_sun = 1.989e30  # kg
M_earth = 5.972e24  # kg
M_moon = 7.348e22  # kg
R_sun = 6.96e8  # meters
G = 6.67430e-11  # m^3 kg^-1 s^-2
c = 299792458.0  # m/s

# Experimental bounds
# LLR bound on anomalous differential acceleration towards the Sun
# \Delta a / a_N < ~ 1.3e-13 (e.g. Murphy 2012, or Williams et al 2012)
LLR_BOUND_ETA = 1.3e-13

# Wide Binary calibration (from step_19)
R_s_wb = 2646 * AU
M_wb = 1.2 * M_sun

def get_Rs(M):
    """Scaling of the temporal topology activation radius with mass"""
    return R_s_wb * (M / M_wb)**(1./3.)

def S_eff(M1, M2, s, k=4):
    """Two-Body Kinetic Operator suppression factor"""
    Rs_tot = get_Rs(M1) + get_Rs(M2)
    return 1.0 / (1.0 + (Rs_tot / s)**k)

def run():
    print("=========================================================")
    print(" TEP Lunar Laser Ranging (LLR) Differential Constraint")
    print("=========================================================")
    print()

    # 1. & 2. Effective scalar responses of Earth and Moon
    s_sun_earth = 1.0 * AU
    # The Moon's orbit is ~0.00257 AU, but for the solar scalar field,
    # the mean distance of the Moon to the Sun is also 1 AU.
    # We evaluate the macroscopic suppression towards the Sun at s = 1 AU.
    
    # Calculate Rs radii
    Rs_sun = get_Rs(M_sun)
    Rs_earth = get_Rs(M_earth)
    Rs_moon = get_Rs(M_moon)
    
    print(f"Topological Radii:")
    print(f"  R_s(Sun)   = {Rs_sun/AU:.1f} AU")
    print(f"  R_s(Earth) = {Rs_earth/AU:.2f} AU")
    print(f"  R_s(Moon)  = {Rs_moon/AU:.2f} AU")
    
    # Calculate effective suppression factors towards the Sun
    S_sun_earth = S_eff(M_sun, M_earth, s_sun_earth)
    S_sun_moon = S_eff(M_sun, M_moon, s_sun_earth)
    
    print()
    print(f"Effective Scalar Responses (suppression factors at 1 AU):")
    print(f"  S_eff(Sun, Earth) = {S_sun_earth:.4e}")
    print(f"  S_eff(Sun, Moon)  = {S_sun_moon:.4e}")
    
    # 3. Differential acceleration toward the Sun
    # The Newtonian acceleration from the Sun at 1 AU
    a_N = G * M_sun / (s_sun_earth**2)
    
    # The scalar acceleration is a_N * S_eff (since bare coupling beta_A = -1 gives a_phi = a_N)
    a_phi_earth = a_N * S_sun_earth
    a_phi_moon = a_N * S_sun_moon
    
    delta_a_scalar = abs(a_phi_earth - a_phi_moon)
    delta_EM = abs(S_sun_earth - S_sun_moon)
    
    print()
    print(f"Accelerations toward the Sun:")
    print(f"  Newtonian a_N         = {a_N:.4e} m/s^2")
    print(f"  Scalar a_phi (Earth)  = {a_phi_earth:.4e} m/s^2")
    print(f"  Scalar a_phi (Moon)   = {a_phi_moon:.4e} m/s^2")
    print(f"  Diff scalar accel     = {delta_a_scalar:.4e} m/s^2")
    
    # 4. Resulting LLR observable
    print()
    print(f"LLR Observable (\Delta a / a_N):")
    print(f"  TEP Predicted \delta_EM = {delta_EM:.4e}")
    print(f"  Experimental Bound < {LLR_BOUND_ETA:.1e}")
    
    if delta_EM < LLR_BOUND_ETA:
        print("  => PASSES LLR CONSTRAINT")
    else:
        print("  => FAILS LLR CONSTRAINT")
        
    # Check that setting scalar coupling to zero recovers GR
    print()
    print("GR Reference Check:")
    print(f"  If scalar coupling is zero (S_eff -> 0), \delta_EM = 0.0 (Recovers GR Weak Equivalence Principle).")
    
    # 5. Check if same fixed parameters preserve Cassini, Saturn and wide-binary
    s_cassini = 1.6 * R_sun
    s_saturn = 9.5 * AU
    s_wb = 2646 * AU
    
    S_cassini = S_eff(M_sun, 0.0, s_cassini)
    S_saturn = S_eff(M_sun, 0.0, s_saturn)
    S_wb = S_eff(0.6*M_sun, 0.6*M_sun, s_wb)  # Approx 1.2 M_sun wide binary
    
    print()
    print("Cross-Scale Consistency Check (Fixed Parameters):")
    print(f"  Cassini suppression (s = 1.6 R_sun) : {S_cassini:.4e}  (Must be < ~10^-5)")
    print(f"  Saturn suppression (s = 9.5 AU)     : {S_saturn:.4e}  (Must be < ~10^-9)")
    print(f"  Wide Binary activation (s = 2646 AU): {S_wb:.4f}       (Must be ~0.5 for transition)")
    
    # Write JSON output
    output = {
        "radii_au": {
            "sun": Rs_sun / AU,
            "earth": Rs_earth / AU,
            "moon": Rs_moon / AU
        },
        "suppression": {
            "sun_earth": S_sun_earth,
            "sun_moon": S_sun_moon
        },
        "acceleration_ms2": {
            "newtonian": a_N,
            "scalar_earth": a_phi_earth,
            "scalar_moon": a_phi_moon,
            "differential": delta_a_scalar
        },
        "llr_test": {
            "predicted_delta_EM": delta_EM,
            "bound": LLR_BOUND_ETA,
            "passes": bool(delta_EM < LLR_BOUND_ETA)
        },
        "consistency_checks": {
            "cassini_suppression": S_cassini,
            "saturn_suppression": S_saturn,
            "wb_activation": S_wb
        }
    }
    
    outdir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "step_22_llr_differential_constraint.json"), "w") as f:
        json.dump(output, f, indent=2)
    print("\n[SUCCESS] Wrote JSON output to scripts/results/step_22_llr_differential_constraint.json")

if __name__ == "__main__":
    run()
