import json
import os
from scipy.optimize import brentq

# Physical constants
AU = 1.496e11  # meters
M_sun = 1.989e30  # kg
M_earth = 5.972e24  # kg
M_moon = 7.348e22  # kg
R_sun = 6.96e8  # meters
R_earth = 6.371e6  # meters
G = 6.674e-11
c = 299792458.0
H0 = 70e3 / 3.085677581e22
BETA_A = -1.0
g_t = c * H0 / (2 * BETA_A ** 2)   # canonical transition acceleration 3.4e-10

# Galactic embedding ambient (calibrated to the wide-binary plateau:
# observed plateau -> q_env^2 ~ 0.43 -> q_env ~ 0.66 -> X_gal ~ 0.52)
# PROVENANCE: this value was reverse-fit to the WB plateau under
# pre-propagator bookkeeping (vertices only). The data-derived
# alternative is the direct solar-circle estimate u_sun = a_sun/g_t =
# 0.570 -> X_sun = 0.32 (v_c=220 km/s, R_0=8.1 kpc; step_32
# direct_solar_circle_ambient), under which the fully-coupled
# propagator+vertex reading reproduces the plateau at ~3% without
# calibration (TEP-WB step_017 consistent_vertex_sweep). The ambient
# is a solved quantity in the theory; the physical value lies in the
# bracket u0 ~ 0.55-0.72. Keep X_GAL=0.52 for continuity with the
# published benchmark chain; the direct value is the canonical
# data-derived estimate going forward (issues 0-2, 13-4).
X_GAL = 0.52


def r_star_M(M):
    """Derived Temporal-Topology shell radius r* = sqrt(GM/g_t) (M^1/2 law)."""
    return (G * M / g_t) ** 0.5


def y_profile(x):
    """Exact flux-conserving profile y: y(1 + y^2 x^-4) = 1, x = r/r*.
    Asymptotes: y ~ x^(4/3) deep inside, y -> 1 outside."""
    if x >= 50.0:
        return 1.0
    return brentq(lambda yy: yy * (1.0 + yy * yy / x ** 4) - 1.0,
                  1e-30, 1.0, xtol=1e-14)


def S_sigma(X):
    """Inverse kinetic stiffness in ambient X (units of Lambda^4/2)."""
    return 1.0 / (1.0 + X)


def R_hierarchical(M_dom, s, X_extra=0.0):
    """Pair response for a test member embedded in a dominant source's
    pre-existing nonlinear field: vertex factors q = S_Sigma(X_env) with
    X_env the dominant field at the pair's scale plus external ambients,
    giving R = q_env^2 * y(s) -> y^3 ~ s^4 deep inside."""
    y = y_profile(s / r_star_M(M_dom))
    q_env = S_sigma((1.0 - y) / y + X_GAL + X_extra)
    return q_env ** 2 * y


def run():
    results = {}
    results["conventions"] = {
        "operator": "R(pair) = S_Sigma(X_env)^2 * y(s); "
                    "y(1 + y^2 (r*/r)^4) = 1; r* = sqrt(GM/g_t) (M^1/2)",
        "regime": "hierarchical channels: X_env = dominant-member field at "
                  "pair scale -> R -> y^3 ~ s^4; comparable-mass pairs: "
                  "X_env = galactic embedding -> S_env^2 * y(s) (4/3)",
        "g_t": g_t,
        "X_gal": X_GAL,
        "note": "Resolved nested operator (step_30, issue 0-27). The former "
                "bare-quartic S_eff = [1+(R_s_tot/s)^4]^-1 with M^(1/3) "
                "scaling is superseded; hierarchical benchmark values are "
                "numerically preserved through the y^3 deep-interior limit."
    }

    # 1. Cassini (s = 1.6 R_sun): the shear-channel path profile is y; the
    #    bound itself applies to the V-sector source charge (separate).
    s_cassini = 1.6 * R_sun
    results["cassini_path_profile"] = y_profile(s_cassini / r_star_M(M_sun))
    results["cassini_shear_response"] = R_hierarchical(M_sun, s_cassini)

    # 2. Saturn ephemeris (s = 9.5 AU)
    s_saturn = 9.5 * AU
    R_sat = R_hierarchical(M_sun, s_saturn)
    results["saturn_suppression"] = R_sat
    aN_sat = G * M_sun / s_saturn ** 2
    results["saturn_anomalous_accel_ms2"] = 2.0 * BETA_A ** 2 * R_sat * aN_sat

    # 3. LLR (Earth-Sun vs Moon-Sun differential at 1 AU)
    s_llr = 1.0 * AU
    X_em = (1.0 - y_profile(2.57e-3 * AU / r_star_M(M_earth))) \
        / y_profile(2.57e-3 * AU / r_star_M(M_earth))
    X_sun_1au = (1.0 - y_profile(s_llr / r_star_M(M_sun))) \
        / y_profile(s_llr / r_star_M(M_sun))
    S_earth_sun = S_sigma(X_sun_1au + X_GAL) ** 2 \
        * y_profile(s_llr / r_star_M(M_sun))
    S_moon_sun = S_sigma(X_sun_1au + X_em + X_GAL) ** 2 \
        * y_profile(s_llr / r_star_M(M_sun))
    results["llr_earth_sun_suppression"] = S_earth_sun
    results["llr_moon_sun_suppression"] = S_moon_sun
    results["llr_differential_bound"] = max(S_earth_sun, S_moon_sun)
    aN_1au = G * M_sun / s_llr ** 2
    results["llr_differential_accel_ms2"] = \
        abs(S_moon_sun - S_earth_sun) * 2.0 * aN_1au

    # 4. Earth-Moon pair response
    R_em = S_sigma(X_em + X_sun_1au + X_GAL) ** 2 \
        * y_profile(2.57e-3 * AU / r_star_M(M_earth))
    results["earth_moon_response"] = R_em

    # 5. GP-B (s = 7013 km from Earth centre)
    s_gpb = 7013e3
    results["gpb_suppression"] = R_hierarchical(M_earth, s_gpb,
                                                X_extra=X_sun_1au)

    # 6. Nested Hierarchy Decomposition at Earth surface
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
