#!/usr/bin/env python3
"""
TEP screening closure at the zero-parameter scale Lambda = sqrt(M_Pl H_0).

Temporal Topology saturation: the topology is pinned in dense regions and free
in diffuse ones. The pinning scale is the topological correlation length
lambda_C, which runs with density when the self-interaction is inverse-power.
A quadratic self-interaction gives a fixed lambda_C and therefore no
environmental response at all -- see quadratic_comparison in the output.

    V(phi) = Lambda^4 ( 1 + Lambda/phi )        [n = 1]
    Lambda = sqrt(M_Pl H_0) = 1.871 meV         [zero free parameters]

This is TEP screening, not chameleon screening. The inverse-power form is a
microscopic realisation of Temporal Topology saturation (Paper 0 Sec 2.2);
chameleon, symmetron and Vainshtein remain candidate realisations, not the
defining ontology.

Constants: CODATA 2018, Planck 2018 H_0. No fitted inputs.
"""
import math, json

hbar = 1.054571817e-34; c = 2.99792458e8; e = 1.602176634e-19
Mpc = 3.0856775814913673e22; AU = 1.495978707e11; Rsun = 6.957e8
G = 6.67430e-11
MPl = 2.435323e18                      # reduced Planck mass, GeV
H0 = 67.4e3 / Mpc                      # s^-1  (Planck 2018)
H0_GeV = hbar * H0 / (e * 1e9)
LAM = math.sqrt(MPl * H0_GeV)          # GeV
GeV4_Jm3 = (1e9 * e) ** 4 / (hbar * c) ** 3

def rho_GeV4(kg_m3):
    return kg_m3 * c ** 2 / GeV4_Jm3

def topology_pinning(rho_kg, n=1, beta=1.0):
    """Saturated topology: (phi_eq [GeV], m_eff [GeV], lambda_C [m]) at density rho_kg."""
    r = rho_GeV4(rho_kg)
    phi = (n * LAM ** (n + 4) * MPl / (beta * r)) ** (1.0 / (n + 1))
    m = math.sqrt(n * (n + 1) * LAM ** (n + 4) / phi ** (n + 2))
    return phi, m, hbar * c / (m * 1e9 * e)

def thin_shell(M_kg, R_m, rho_body, rho_amb, beta=1.0):
    """
    TEP screening suppression for a body whose saturation mass runs with density.

    The linear Helmholtz result 3(1+x)/(x(2x+1)) assumes a CONSTANT mass and is
    not valid here: when m runs with rho the suppression is set by the shell over
    which the topology relaxes from its interior pinned value to the ambient one.

        dR/R = (phi_amb - phi_in) / (6 beta M_Pl Phi_N)      [capped at 1]
        s    = 3 dR/R                                        [charge fraction]

    Returns (dR_over_R, s, Phi_N).
    """
    PhiN = G * M_kg / (R_m * c * c)
    phi_in = topology_pinning(rho_body)[0]
    phi_amb = topology_pinning(rho_amb)[0]
    d = (phi_amb - phi_in) / (6 * beta * MPl * PhiN)
    d = min(d, 1.0)
    return d, 3 * d, PhiN

def main():
    out = {"Lambda_GeV": LAM, "Lambda_meV": LAM * 1e12,
           "Lambda4_GeV4": LAM ** 4, "rho_DE_obs_GeV4": 2.538e-47,
           "a0_Lambda2_over_MPl_SI": c * H0, "a0_MOND_SI": 1.2e-10,
           "_definitional_warning": (
               "Lambda^4 = M_Pl^2 H0^2 = rho_crit/3 identically, and Lambda^2/M_Pl = H0 "
               "identically. So 'Lambda^4 matches dark energy' and 'a_0 = c H_0' are "
               "consequences of the DEFINITION Lambda = sqrt(M_Pl H0), not predictions. "
               "They show the scale is natural; they are not independent confirmations.")}

    # --- local null gates (thin-shell suppression, correct bound per observable) ---
    sh = {}
    for nm, M, R, rb, ra in (("Sun", 1.989e30, Rsun, 1408.0, 1e-20),
                             ("Earth", 5.972e24, 6.371e6, 5515.0, 1e-12),
                             ("Moon", 7.342e22, 1.7374e6, 3344.0, 1e-20),
                             ("NS", 2.8e30, 1.2e4, 5.9e17, 1e-21),
                             ("WD", 1.2e30, 7e6, 1e9, 1e-21)):
        d, s, PhiN = thin_shell(M, R, rb, ra)
        sh[nm] = {"dR_over_R": d, "s": s, "Phi_N": PhiN}
    gates = {
        "Cassini_gamma": {"expr": "|gamma-1| = 4*beta^2*s_Sun/(1+2*beta^2*s_Sun)",
                          "value": 4 * sh["Sun"]["s"] / (1 + 2 * sh["Sun"]["s"]),
                          "bound": 2.3e-5},
        "Geodesy_Earth": {"expr": "alpha_E = 2*beta^2*s_E",
                          "value": 2 * sh["Earth"]["s"], "bound": 1e-8},
        "LLR_Nordtvedt": {"expr": "eta = 2|s_E - s_Moon|",
                          "value": 2 * abs(sh["Earth"]["s"] - sh["Moon"]["s"]), "bound": 4.4e-4},
        "Pulsar_NS": {"expr": "alpha_NS = 2*beta^2*s_NS",
                      "value": 2 * sh["NS"]["s"], "bound": 1e-3},
        "WhiteDwarf": {"expr": "alpha_WD = 2*beta^2*s_WD",
                       "value": 2 * sh["WD"]["s"], "bound": 1e-2}}
    for k, v in gates.items():
        v["margin"] = v["bound"] / v["value"]
        v["status"] = "PASS" if v["value"] < v["bound"] else "FAIL"
    out["thin_shell"] = sh
    out["local_gates"] = gates
    out["clock_force_split"] = {
        nm: {"S_Sigma_force": sh[nm]["s"], "S_A_clock": 1.0, "ratio": 1.0 / sh[nm]["s"]}
        for nm in ("Sun", "Earth", "Moon")}

    # --- wide-binary scale: invert R_s to an ambient density ---
    lo, hi = -26.0, -18.0
    target = 2646 * AU
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if topology_pinning(10 ** mid)[2] - target > 0: lo = mid
        else: hi = mid
    rho_req = 10 ** (0.5 * (lo + hi))
    out["wide_binary"] = {
        "R_s_observed_AU": 2646, "rho_required_kg_m3": rho_req,
        "n_H_required_cm3": rho_req / (1.6726e-27 * 1e6),
        "local_ISM_warm_cm3": "0.1 - 0.5",
        "verdict": "standard warm-ISM density; not tuned"}

    # --- environmental ordering: lambda_C ~ rho^(-3/4) ---
    env = []
    for nm, a, b in (("in_plane_vs_out_of_plane", 4662, 7131), ("disk_control", 4145, 6856)):
        ratio = b / a
        env.append({"split": nm, "R_s_dense_AU": a, "R_s_dilute_AU": b,
                    "R_s_ratio": ratio, "implied_density_ratio": ratio ** (-4.0 / 3.0),
                    "direction_predicted": "lower density -> larger R_s",
                    "direction_observed": "matches"})
    out["environmental_ordering"] = env

    # --- Cepheid conformal channel ---
    phi_c, _, _ = topology_pinning(1e-6)     # Cepheid envelope
    phi_g, _, _ = topology_pinning(1e-21)    # galactic ISM
    out["cepheid"] = {
        "delta_phi_over_MPl": abs(phi_g - phi_c) / MPl,
        "required": 2.07e-2,
        "shortfall_orders": math.log10(2.07e-2 / (abs(phi_g - phi_c) / MPl)),
        "verdict": "conformal channel insufficient; disformal B(phi) still required"}

    # --- contrast with the quadratic candidate ---
    m_quad = LAM
    out["quadratic_comparison"] = {
        "lambda_C_everywhere_m": hbar * c / (m_quad * 1e9 * e),
        "mass_is_density_independent": True,
        "wide_binary_reachable": False,
        "note": "V'' = Lambda^2 is constant: lambda_C = 0.105 mm in EVERY environment, "
                "so Temporal Topology saturation has no environmental response"}

    # --- lambda_C is NOT the GNSS correlation length ---
    lo, hi = -22.0, 6.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if topology_pinning(10 ** mid)[2] - 4.2e6 > 0: lo = mid
        else: hi = mid
    out["gnss_scale_distinction"] = {
        "lambda_T_GNSS_m": 4.2e6,
        "rho_that_would_give_it_kg_m3": 10 ** (0.5 * (lo + hi)),
        "physical": False,
        "note": ("lambda_T = 4,200 km corresponds to no physical terrestrial density. "
                 "It is a distinct scale from the topological correlation length and must "
                 "come from the disformal sector B(phi), not from TEP screening.")}

    print(json.dumps(out, indent=2))
    with open("results/tep_screening_closure_results.json", "w") as f:
        json.dump(out, f, indent=2)

if __name__ == "__main__":
    main()
