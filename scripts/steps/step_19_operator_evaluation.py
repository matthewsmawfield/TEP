import json
import os

# Physical constants
AU = 1.496e11  # meters
M_sun = 1.989e30  # kg
M_earth = 5.972e24  # kg
M_moon = 7.348e22  # kg
R_sun = 6.96e8  # meters
R_earth = 6.371e6  # meters

# Wide Binary calibration: R_s = 2646 AU measured for Paper 13's
# sample median total mass M ~ 1.2 M_sun (TEP-WB §2.2)
R_s_wb = 2646 * AU
M_wb = 1.2 * M_sun

def get_Rs(M):
    """Scaling of the temporal topology activation radius with mass"""
    # R_s(M) prop M^(1/3)
    return R_s_wb * (M / M_wb)**(1./3.)

def S_eff(M1, M2, s, k=4):
    """Two-Body Kinetic Operator suppression factor"""
    Rs_tot = get_Rs(M1) + get_Rs(M2)
    return 1.0 / (1.0 + (Rs_tot / s)**k)

def run():
    results = {}
    results["conventions"] = {
        "R_s_scaling": "R_s(M) = 2646 AU * (M / 1.2 M_sun)^(1/3), k = 4",
        "note": "R_s = 2646 AU is the Paper 13 wide-binary measurement at sample median total mass ~1.2 M_sun; R_s(1 M_sun) = 2490 AU under p = 1/3."
    }
    
    # 1. Cassini bound (s = 1.6 R_sun)
    # The spacecraft acts as a test mass, so we use M_sun and 0.
    s_cassini = 1.6 * R_sun
    S_cassini = S_eff(M_sun, 0.0, s_cassini)
    results["cassini_suppression"] = S_cassini
    results["cassini_suppression_Rs2646_at_1Msun"] = 1.0 / (1.0 + (2646 * AU / s_cassini)**4)
    
    # 2. Saturn ephemeris (s = 9.5 AU)
    s_saturn = 9.5 * AU
    S_saturn = S_eff(M_sun, 0.0, s_saturn) # Using test mass for generic S_eff bound, though M_saturn could be used
    results["saturn_suppression"] = S_saturn
    
    # 3. LLR (Earth-Moon differential)
    # Earth-Sun and Moon-Sun suppressions
    s_llr = 1.0 * AU
    S_earth_sun = S_eff(M_sun, M_earth, s_llr)
    S_moon_sun = S_eff(M_sun, M_moon, s_llr)
    
    # The absolute coupling is suppressed by this factor.
    # Differential acceleration is bounded by the max of these.
    results["llr_earth_sun_suppression"] = S_earth_sun
    results["llr_moon_sun_suppression"] = S_moon_sun
    results["llr_differential_bound"] = max(S_earth_sun, S_moon_sun)
    
    # 4. GP-B (s = 7013 km from Earth center)
    s_gpb = 7013e3 # meters
    S_gpb = S_eff(M_earth, 0.0, s_gpb)
    results["gpb_suppression"] = S_gpb
    
    # 5. Nested Hierarchy Decomposition at Earth surface
    # Two bases: (a) field-amplitude decomposition from the step_01
    # nonlinear nested solve — the basis quoted in manuscript §3;
    # (b) crude Newtonian-potential proxy, kept labeled for contrast.
    step01_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))), "results", "step_01_radial_ode.json")
    with open(step01_path) as f:
        nh = json.load(f)["nested_hierarchy"]["decomposition_at_earth_surface"]
    results["hierarchy_decomposition_field_amplitude_basis"] = {
        "galactic_percent": nh["galactic_fraction"] * 100,
        "solar_percent": nh["solar_fraction"] * 100,
        "terrestrial_percent": nh["earth_fraction"] * 100,
        "source": "results/step_01_radial_ode.json nested_hierarchy"
    }
    
    G = 6.674e-11
    M_gal = 1e11 * M_sun
    R_gal = 8e3 * 3.086e16
    Phi_gal = G * M_gal / R_gal
    
    Phi_sun = G * M_sun / AU
    Phi_earth = G * M_earth / R_earth
    
    Phi_tot = Phi_gal + Phi_sun + Phi_earth
    
    results["hierarchy_decomposition_newtonian_potential_basis"] = {
        "galactic_percent": (Phi_gal / Phi_tot) * 100,
        "solar_percent": (Phi_sun / Phi_tot) * 100,
        "terrestrial_percent": (Phi_earth / Phi_tot) * 100
    }
    
    # Print and save
    print(json.dumps(results, indent=2))
    
    outdir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "step_19_operator_evaluation.json"), 'w') as f:
        json.dump(results, f, indent=2)

if __name__ == "__main__":
    run()
