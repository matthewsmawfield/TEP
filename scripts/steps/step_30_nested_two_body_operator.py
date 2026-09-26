#!/usr/bin/env python3
"""Nested two-body evaluation of the master-action screening operator.

The pairwise suppression originally quoted as
    S_eff = [1 + (R_s/s)^4]^-1
evaluated the environmental operator S_Sigma(E) at the *linear* Newtonian
profile g = GM/s^2. The flux-conserving radial solution of the stated
kinetic completion, P = X - V + X|X|/Lambda^4 with P_X phi' = C/r^2, gives
phi' ~ r^(-2/3) deep inside the shell and the exact implicit profile

    y (1 + y^2 (r*/r)^4) = 1,   y = phi'/phi'_lin,  r* = sqrt(GM/g_t).

The resolved two-body structure (bi-directional, nested) is

    R(pair) = S_Sigma(X_env)^2 * y(s),

i.e. the pair's mutual shear y(s) — the perturbation delta_phi on the
embedding ambient — propagated with the corpus's incremental-response
factor S_Sigma = 1/P_X evaluated at the ambient at the two vertices
(emission and response). X_env is the embedding ambient of the pair: the
field of the hierarchy level above the pair, excluding the pair's own
mutual shear field (suppressing a response by the response itself is the
double-counting the old s-quartic committed in the transition region).

Two regimes of the same operator:

  * Comparable-mass pairs (field wide binaries): the mutual field is the
    co-generated response; the embedding is the Galactic field and the
    bridge region between the members is Galactic-dominated by symmetry
    (equal contributions oppose and cancel at the saddle). The vertex
    factors saturate at S_env ~ 0.5, leaving R ~ S_env^2 * y(s) — the
    4/3 transition profile selected by the real 341k-binary forward
    model (TEP-WB step_015, chi2_red ~ 4.3-7.4 vs ~18-25 for steeper
    readings).

  * Hierarchical pairs (test member embedded in a dominant member's
    pre-existing nonlinear shell): the dominant member's field IS the
    embedding ambient — it exists independent of the pair interaction.
    By the flux-conservation identity 2X_dom/Lambda^4 = (1-y)/y the
    ambient factor equals the dominant member's own profile value,
    giving R ~ y^3 ~ s^4 deep inside — the quartic recovered as the
    nested deep-interior asymptote, and the reason the corpus's original
    Solar-System quartic bookkeeping was quantitatively correct.

Channels evaluated under the single rule:
  - Sun-Saturn, Sun-Mercury orbital perturbation (real ephemeris bound
    ~5e-13 m s^-2 supplementary acceleration, INPOP10e/Cassini tracking);
  - Earth-Moon pair (vs corpus-quoted ~2.5e-14);
  - Earth-Sun vs Moon-Sun differential (Nordtvedt channel);
  - Cassini photon path (single-body exterior charge — the separate
    V-sector channel, bound 5.8e-6);
  - wide-binary transition response under the saturated-ambient reading.
"""
import json
import numpy as np
from scipy.optimize import brentq

G = 6.674e-11
AU = 1.496e11
M_SUN = 1.989e30
M_EARTH = 5.972e24
M_MOON = 7.3477e22
M_SATURN = 5.683e26
M_MERCURY = 3.301e23
BETA = 1.0            # |beta_A|
G_T = 3.4e-10         # cH0/(2 beta_A^2), m s^-2
X_GAL = 0.52          # 2 X_galactic / Lambda^4 — calibrated so the vertex
                      # product S_Sigma^2 reproduces the observed plateau
                      # amplitude (q_env = 0.66 -> q_env^2 = 0.43 ~ the
                      # fitted S_eff(infty) = 0.434 of the alpha_sat = 0.366
                      # plateau); consistent with the corpus's S_Sigma(MW)
                      # ~ 0.74 ambient suppression


def Rs_AU(mass):
    return np.sqrt(G * mass / G_T) / AU


def y_profile(x):
    """Exact flux-conserving profile y = phi'/phi'_lin at x = r/r*."""
    if x <= 0:
        return 0.0
    return brentq(lambda y: y * (1 + y * y / x**4) - 1, 1e-30, 1.0,
                  xtol=1e-14)


def S_sigma(x_amb):
    """Environmental suppression at ambient kinetic value 2X/Lambda^4."""
    return 1.0 / (1.0 + x_amb)


def ambient_x_of_partner(y_partner, x_extra=X_GAL):
    """Ambient kinetic value inside a partner's nonlinear field.
    Identity: 2X_partner/Lambda^4 = (1 - y)/y."""
    return (1.0 - y_partner) / y_partner + x_extra


def run():
    out = {"conventions": {
        "operator": "R(pair) = S_Sigma(X_env)^2 * y(s); vertex ambient = "
                    "embedding field of the pair (hierarchy above), "
                    "excluding the pair's own mutual shear field",
        "y_eq": "y(1 + y^2 (r*/r)^4) = 1, flux-conserving P_X phi' = C/r^2 "
                "for P = X - V + X|X|/Lambda^4",
        "identity": "2X_dom/Lambda^4 = (1 - y)/y, so a pair embedded in a "
                    "dominant member's nonlinear field has vertex factors "
                    "= y(s), giving R ~ y^3 ~ s^4 (deep interior)",
        "regimes": "comparable-mass pair -> co-generated mutual field, "
                   "ambient saturates external (Galactic) -> R ~ "
                   "S_env^2*y ~ s^(4/3); hierarchical pair -> dominant "
                   "member's field is the ambient -> R ~ y^3 ~ s^4",
        "x_gal": X_GAL,
        "g_t": G_T,
    }}

    # --- shell radii ---
    shells = {k: Rs_AU(m) for k, m in
              [("Sun", M_SUN), ("Earth", M_EARTH), ("Moon", M_MOON),
               ("Saturn", M_SATURN), ("Mercury", M_MERCURY)]}
    out["shell_radii_AU"] = shells

    # --- exact profile sample ---
    xs = np.geomspace(1e-4, 100, 200)
    out["profile"] = {
        "x_over_Rs": xs.tolist(),
        "y": [y_profile(x) for x in xs],
    }
    kk = np.gradient(np.log([y_profile(x) for x in xs]), np.log(xs))
    m = (xs > 0.2) & (xs < 5)
    out["transition_effective_exponent"] = {
        "min": float(kk[m].min()), "median": float(np.median(kk[m])),
        "max": float(kk[m].max()),
        "note": "effective recovery exponent d ln y / d ln s of the exact "
                "profile over 0.2-5 r*; asymptotic interior value is 4/3",
    }

    q_gal = S_sigma(X_GAL)
    out["galactic_charge_factor"] = q_gal

    # --- Cassini: photon path through the Sun's suppressed field ---
    s_cass = 1.6 * 0.00465047   # 1.6 R_sun in AU
    y_cass = y_profile(s_cass / shells["Sun"])
    y_1AU = y_profile(1.0 / shells["Sun"])
    out["cassini"] = {
        "y_at_1p6_Rsun": y_cass,
        "y_at_1AU": y_1AU,
        "single_body_charge_bound": 5.8e-6,
        "note": "PPN probe reads the exterior charge (V-sector channel, "
                "corpus adopted ~2e-6); the kinetic profile values are the "
                "shear-channel suppression along/near the ray path, not "
                "the charge itself.",
    }

    # --- hierarchical orbital channels: X_env = dominant member's field
    # at the pair separation (the pair is embedded inside its shell).
    # q_env = S_Sigma(X_env) = y_dom(s) by the identity, so R = y^3.

    # Sun-Saturn
    y_sat = y_profile(9.5 / shells["Sun"])
    q_env_sat = S_sigma(ambient_x_of_partner(y_sat))
    R_sat = q_env_sat**2 * y_sat
    gN_sat = G * M_SUN / (9.5 * AU)**2
    out["saturn"] = {
        "s_AU": 9.5, "y_at_orbit": y_sat, "q_env": q_env_sat,
        "R_pair": R_sat,
        "delta_a_m_s2": 2 * BETA**2 * R_sat * gN_sat,
        "ephemeris_bound_m_s2": 5e-13,
        "ephemeris_bound_source": "INPOP10e/EPM Cassini-tracking "
            "supplementary acceleration; Pitjeva perihelion ~9e-15",
        "margin_vs_bound": 5e-13 / (2 * BETA**2 * R_sat * gN_sat),
        "old_quartic_corpus_value": 2.7e-11,
        "note": "R = y^3 = s^4 exactly equals the corpus's original "
                "quartic suppression at Saturn — the quartic was the "
                "correct hierarchical-nested response all along",
    }

    # Sun-Mercury
    y_merc = y_profile(0.387 / shells["Sun"])
    q_env_merc = S_sigma(ambient_x_of_partner(y_merc))
    R_merc = q_env_merc**2 * y_merc
    gN_merc = G * M_SUN / (0.387 * AU)**2
    out["mercury"] = {
        "s_AU": 0.387, "y_at_orbit": y_merc, "q_env": q_env_merc,
        "R_pair": R_merc,
        "delta_a_m_s2": 2 * BETA**2 * R_merc * gN_merc,
    }

    # Earth-Moon: embedded in Earth's nonlinear field at the lunar orbit
    # (X_env = Earth's field); the solar field at 1 AU is subdominant
    # (X_sun(1AU) ~ 6.7e4 vs X_earth(2.57e-3AU) ~ 4e4 — both deep; the
    # dominant-member reading uses Earth's shell).
    y_em = y_profile(2.57e-3 / shells["Earth"])
    q_env_em = S_sigma(ambient_x_of_partner(y_em,
                                          x_extra=ambient_x_of_partner(y_1AU)))
    R_em = q_env_em**2 * y_em
    out["earth_moon"] = {
        "s_AU": 2.57e-3, "y_Earth_at_moon": y_em,
        "q_env": q_env_em,
        "R_pair": R_em,
        "corpus_quoted": 2.5e-14,
        "ratio_to_corpus": R_em / 2.5e-14,
        "note": "embedding ambient = Earth's field at lunar orbit plus "
                "solar field at 1 AU; R ~ y_em^3",
    }

    # Earth-Sun and Moon-Sun (Nordtvedt differential)
    y_es = y_profile(1.0 / shells["Sun"])
    q_env_es = S_sigma(ambient_x_of_partner(y_es))
    R_earth_sun = q_env_es**2 * y_es
    # Moon's embedding additionally carries Earth's field at its orbit
    q_env_ms = S_sigma(ambient_x_of_partner(y_es)
                       + ambient_x_of_partner(y_em, x_extra=0.0))
    R_moon_sun = q_env_ms**2 * y_es
    gN_es = G * M_SUN / AU**2
    out["earth_sun"] = {"y_at_1AU": y_es, "q_env": q_env_es,
                        "R_pair": R_earth_sun}
    out["llr_differential"] = {
        "q_env_Earth": q_env_es, "q_env_Moon": q_env_ms,
        "R_EarthSun": R_earth_sun, "R_MoonSun": R_moon_sun,
        "differential_response": abs(R_moon_sun - R_earth_sun),
        "differential_delta_a_m_s2":
            abs(R_moon_sun - R_earth_sun) * 2 * BETA**2 * gN_es,
        "nordtvedt_bound_response": "~1e-13-level on the effective "
            "differential acceleration (eta < 4.4e-4)",
    }

    # --- wide-binary transition-region response ---
    wb_rows = []
    for s_over_rs in (0.3, 0.5, 1.0, 2.0, 4.0):
        yv = y_profile(s_over_rs)
        q_mut = S_sigma(ambient_x_of_partner(yv))
        wb_rows.append({
            "s_over_Rs": s_over_rs, "y": yv,
            "R_saturated_external_ambient": q_gal**2 * yv,
            "R_if_mutual_counted": q_mut**2 * yv,
        })
    out["wide_binary_responses"] = wb_rows
    out["wide_binary_forward_model"] = {
        "source": "TEP-WB step_015-derived-operator + nested-response scan",
        "chi2_red_single_y": 7.4,
        "chi2_red_one_vertex": 12.1,
        "chi2_red_two_vertices": 25.3,
        "verdict": "the population statistic selects the saturated-ambient "
                   "reading: for comparable-mass pairs the mutual field is "
                   "the measured response (co-generated), not a vertex "
                   "ambient — the fitted S_env ~ 0.45 plateau is the "
                   "external ambient vertex product",
    }

    # --- resolution ledger ---
    out["resolution"] = {
        "quartic_status": "R = y^3 ~ s^4 deep inside a dominant member's "
            "shell: the corpus's quartic Solar-System bookkeeping was "
            "quantitatively correct for hierarchical pairs; its error was "
            "applying the same form to comparable-mass wide binaries, "
            "where the vertex ambient saturates externally and the "
            "transition profile is s^(4/3)",
        "saturn_verdict": "delta_a = 3.5e-15 m s^-2 vs bound 5e-13 — "
            "passes by ~140x under the unified nested rule",
        "charge_vertex": "the exterior V-sector charge does not enter the "
            "shear vertex: a suppressed emitted flux would shrink the "
            "kinetic shell (r*^4 = C^2/Lambda^4), collapsing stellar "
            "shells to ~10 AU and erasing the observed 2.6 kAU "
            "transition — the corpus's charge/shear channel split is "
            "required, not optional",
        "falsifiable_refinement": "the mass-ratio dependence of the "
            "vertex ambient predicts unequal-mass wide binaries should "
            "show a mildly steeper effective transition (partial bridge "
            "cancellation) — a testable population prediction",
    }

    return out


if __name__ == "__main__":
    res = run()
    path = "results/step_30_nested_two_body_operator.json"
    with open(path, "w") as f:
        json.dump(res, f, indent=2)
    print("Saturn delta_a: "
          f"{res['saturn']['delta_a_m_s2']:.2e} m/s^2 "
          f"(bound 5e-13, margin {res['saturn']['margin_vs_bound']:.0f}x)")
    print("Mercury delta_a: "
          f"{res['mercury']['delta_a_m_s2']:.2e} m/s^2")
    print("Earth-Moon R: "
          f"{res['earth_moon']['R_pair']:.2e} (corpus 2.5e-14, "
          f"ratio {res['earth_moon']['ratio_to_corpus']:.2f})")
    print("LLR differential delta_a: "
          f"{res['llr_differential']['differential_delta_a_m_s2']:.2e} m/s^2")
    print("Transition effective exponent: "
          f"{res['transition_effective_exponent']['median']:.2f}")
    print(f"Saved {path}")
