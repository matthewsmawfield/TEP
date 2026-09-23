"""Legacy potential-scan and nested-hierarchy computations (v0.13 restored).

These sections were computed by the pre-refactor step_01 and remain cited by
the manuscript: the three-potential comparison (V = 0, quartic, quadratic),
the 3-zone density profile, the stellar-population scan, and the nested
screening hierarchy (galaxy -> Sun -> Earth decomposition).

Conventions: normalized field Psi = varphi/psi_uns (O(1) unscreened),
mu0 the effective dimensionless potential strength. The
gamma_1au_screened_source fields evaluate the screened-source PPN formula
gamma - 1 = -4 beta_A^2 S_Sigma/(1 + 2 beta_A^2 S_Sigma), linear in S_Sigma
because the photon probe is unscreened (Burrage & Sakstein 2018, Eq. 3.31).
This far-field PPN check at 1 AU is the binding Cassini
test. The conjunction-point (1.6 R_sun) diagnostics in step_01 are
near-field quantities reported separately.
"""
import numpy as np
from scipy.integrate import solve_bvp
import warnings
warnings.filterwarnings('ignore')
from tep_model import (HBAR_C as HBAR_C_GeV_m, M_PL as M_PL_GEV, G as G_SI,
                       C as C_SI, M_SUN as M_SUN_SI, R_SUN as R_SUN_SI,
                       R_EARTH as R_EARTH_SI, M_EARTH as M_EARTH_SI, AU as AU_SI)

CASSINI_BOUND = 2.3e-5


def compute_source_params(M_gev, R_gev):
    M_natural = M_gev / M_PL_GEV
    R_natural = R_gev * M_PL_GEV
    psi_uns = M_natural / (4.0 * np.pi * R_natural)
    rho = 3.0 * M_natural / (4.0 * np.pi * R_natural**3)
    return psi_uns, rho, M_natural, R_natural


def si_to_gev_mass(M_si_kg):
    return M_si_kg / 1.7827e-27


def si_to_gev_length(L_si_m):
    return L_si_m / HBAR_C_GeV_m


M_SUN_G = si_to_gev_mass(M_SUN_SI)
R_SUN_G = si_to_gev_length(R_SUN_SI)
M_EARTH_G = si_to_gev_mass(M_EARTH_SI)
R_EARTH_G = si_to_gev_length(R_EARTH_SI)
M_GAL_G = 1e11 * M_SUN_G
R_GAL_G = si_to_gev_length(50 * 3.086e19)

PSI_SUN = compute_source_params(M_SUN_G, R_SUN_G)[0]
PSI_EARTH = compute_source_params(M_EARTH_G, R_EARTH_G)[0]
M_NS_G = 1.4 * M_SUN_G
R_NS_G = si_to_gev_length(1e4)
PSI_NS = compute_source_params(M_NS_G, R_NS_G)[0]
M_BH_G = 1e6 * M_SUN_G
R_BH_SI = 2 * G_SI * 1e6 * M_SUN_SI / C_SI**2
R_BH_G = si_to_gev_length(R_BH_SI)
PSI_BH = compute_source_params(M_BH_G, R_BH_G)[0]
PSI_GAL = compute_source_params(M_GAL_G, R_GAL_G)[0]


def solve_bvp_conformal(psi_uns, potential_type="none", lam=0.0, x_max=1e6,
                        rho_bg_ratio=0.0, n_points=500, guess=None):
    """Normalized-field BVP: Psi = varphi/psi_uns, x = r/R.

    Interior: nabla^2 Psi = -3 exp(-psi_uns Psi) + V_term; exterior with
    rho_bg_ratio. V_term: 0 / mu0 Psi^3 (quartic) / mu0 Psi (quadratic).
    """
    mu0 = lam if potential_type in ("quartic", "quadratic") else 0.0

    def ode(x, y):
        psi, dpsi = y
        x_safe = np.maximum(x, 1e-10)
        if potential_type == "quartic":
            v_term = mu0 * psi**3
        elif potential_type == "quadratic":
            v_term = mu0 * psi
        else:
            v_term = 0.0
        delta = 0.01
        theta = 0.5 * (1.0 - np.tanh((x - 1.0) / delta))
        rho_eff = 3.0 * (theta + rho_bg_ratio * (1.0 - theta))
        source = -rho_eff * np.exp(-psi_uns * psi) + v_term
        return np.vstack([dpsi, source - 2.0 * dpsi / x_safe])

    def bc(ya, yb):
        return np.array([ya[1], yb[0]])

    if guess is not None:
        if hasattr(guess, 'x') and hasattr(guess, 'y'):
            x_init, y_init = guess.x, guess.y
        else:
            x_init, y_init = guess
    else:
        x_init = np.logspace(-4, np.log10(x_max), n_points)
        x_init[0] = 1e-4
        x_init[-1] = x_max
        psi_init = np.where(x_init < 1, 1.5 - x_init**2 / 2.0, 1.0 / x_init)
        psi_init = np.maximum(psi_init, 0)
        dpsi_init = np.where(x_init < 1, -x_init, -1.0 / x_init**2)
        y_init = np.vstack([psi_init, dpsi_init])

    sol = solve_bvp(ode, bc, x_init, y_init, tol=1e-6, max_nodes=200000, verbose=0)
    if not sol.success:
        return None, None, None, {"success": False, "message": sol.message}

    x_dense = np.logspace(-4, np.log10(x_max), 1000)
    x_dense[0] = 1e-4
    x_dense[-1] = x_max
    info = {"success": True, "message": sol.message,
            "psi_center": float(sol.sol(1e-4)[0]),
            "psi_surface": float(sol.sol(1.0)[0]),
            "n_nodes": len(sol.x), "sol": sol}
    return x_dense, sol.sol(x_dense)[0], sol.sol(x_dense)[1], info


def compute_screening(x, psi, dpsi, psi_uns):
    """S_Sigma(x) = x^2 |Psi'| (1 unscreened); S_A from A = exp(-psi_uns Psi)."""
    S_sigma = x**2 * np.abs(dpsi)
    A = np.exp(-psi_uns * psi)
    A_uns = np.exp(-psi_uns / x)
    with np.errstate(divide='ignore', invalid='ignore'):
        S_A = np.where(np.abs(A_uns - 1) > 1e-30, (A - 1) / (A_uns - 1), 1.0)
    return S_sigma, S_A


def compute_gamma(S_sigma):
    """Screened PPN formula (Burrage & Sakstein 2018, Living Rev. Relativity
    21:1, Eq. 3.31): gamma - 1 ~ -4 beta_A^2 (1-M(r_s)/M), LINEAR in the
    thin-shell/unscreened mass fraction, which is S_sigma here. A squared-S
    formula double-counts the screening (it was used in an earlier version
    of this file and has been corrected). The Pade form below reduces to
    this linear result for S<<1 and to the standard unscreened DEF result
    gamma-1=-2*alpha_DEF^2/(1+alpha_DEF^2) at S=1 (alpha_DEF^2=2*beta_A^2,
    beta_A=-1). S_sigma is the screened source charge evaluated in the far
    field (1 AU), where the profile has joined its Coulomb exterior; this
    is the binding Cassini test. The conjunction-point (1.6 R_sun)
    near-field diagnostics are reported separately in step_01.
    """
    return -4.0 * S_sigma / (1.0 + 2.0 * S_sigma)


def interp_log(x, vals, xt):
    if xt < x[0] or xt > x[-1]:
        return None
    return float(np.interp(np.log10(xt), np.log10(x), vals))


def _scan_row(x, psi, dpsi, psi_uns, info, x_eval, x_wb=None):
    S_sig, _ = compute_screening(x, psi, dpsi, psi_uns)
    s_eval = interp_log(x, S_sig, x_eval)
    g_eval = compute_gamma(np.array([s_eval]))[0] if s_eval is not None else None
    row = {
        "psi_center": float(info["psi_center"]),
        "s_sigma_surface": float(interp_log(x, S_sig, 1.0)),
        "s_sigma_1au": float(s_eval) if s_eval is not None else None,
        "gamma_1au_screened_source": float(g_eval) if g_eval is not None else None,
        "proxy_below_bound_1au": bool(abs(g_eval) < CASSINI_BOUND) if g_eval is not None else None,
    }
    if x_wb is not None:
        s_wb = interp_log(x, S_sig, x_wb)
        row["s_sigma_wb"] = float(s_wb) if s_wb is not None else None
        row["wb_unscreened"] = bool(s_wb is not None and s_wb > 0.1)
    return row


def run():
    """Run the legacy potential scans and nested-hierarchy decomposition."""
    out = {"note": ("Normalized-field solver (Psi = varphi/psi_uns); "
                    "gamma_1au_screened_source is the screened-source PPN "
                    "gamma -4 beta_A^2 S/(1+2 beta_A^2 S) (linear in S; "
                    "name is historical) evaluated with the 1 AU source "
                    "charge — the binding Cassini test (see step_01 "
                    "solar_scan for the near-field conjunction "
                    "diagnostics)."),
           "potentials": {}, "nested_hierarchy": {}}

    x_1au = AU_SI / R_SUN_SI
    x_wb = 2646.0 * x_1au

    # ---- V = 0 across sources ----
    pot_none = {}
    for name, psi_uns in [("Earth", PSI_EARTH), ("Sun", PSI_SUN),
                          ("NeutronStar", PSI_NS), ("BH_1e6", PSI_BH),
                          ("Galaxy", PSI_GAL)]:
        x, psi, dpsi, info = solve_bvp_conformal(psi_uns, "none", x_max=1e6)
        if x is None:
            continue
        row = {"psi_uns": float(psi_uns)}
        row.update(_scan_row(x, psi, dpsi, psi_uns, info, x_1au))
        pot_none[name] = row
    out["potentials"]["none"] = pot_none

    # ---- Quartic mu0 scan (Sun) ----
    quartic = {}
    prev = None
    for mu0 in [0, 0.001, 0.01, 0.1, 1.0, 3.0, 10.0, 30.0, 100.0,
                300.0, 1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6, 3e6]:
        x, psi, dpsi, info = solve_bvp_conformal(
            PSI_SUN, "quartic", lam=mu0, x_max=1e6, guess=prev)
        if x is None and prev is not None:
            x, psi, dpsi, info = solve_bvp_conformal(
                PSI_SUN, "quartic", lam=mu0, x_max=1e6, guess=None)
        if x is None:
            continue
        prev = info["sol"]
        row = {"mu0": float(mu0)}
        row.update(_scan_row(x, psi, dpsi, PSI_SUN, info, x_1au))
        quartic[f"mu0_{mu0}"] = row
    out["potentials"]["quartic"] = quartic

    # ---- Quadratic mu0 scan (Sun), with WB column ----
    quadratic = {}
    for mu0 in [0, 1e-6, 1e-4, 1e-2, 0.1, 1.0, 10.0, 100.0, 1e3, 1e4, 1e6]:
        x, psi, dpsi, info = solve_bvp_conformal(
            PSI_SUN, "quadratic", lam=mu0, x_max=1e7)
        if x is None:
            continue
        row = {"mu0": float(mu0)}
        row.update(_scan_row(x, psi, dpsi, PSI_SUN, info, x_1au, x_wb))
        quadratic[f"mu0_{mu0}"] = row
    out["potentials"]["quadratic"] = quadratic

    # ---- 3-zone density profile (quadratic + quartic) ----
    x_helio = 2e4

    def solve_3zone(psi_uns, mu0, rho_ip=1e-25, rho_ism=1e-27, x_max=1e7,
                    potential="quadratic", guess=None):
        def ode(x, y):
            psi, dpsi = y
            x_safe = np.maximum(x, 1e-10)
            if potential == "quadratic":
                v_term = mu0 * psi
            elif potential == "quartic":
                v_term = mu0 * psi**3
            else:
                v_term = 0.0
            theta1 = 0.5 * (1.0 - np.tanh((x - 1.0) / 0.01))
            theta2 = 0.5 * (1.0 - np.tanh((x - x_helio) / 0.5))
            rho = theta1 + rho_ip * (1 - theta1) * theta2 + \
                rho_ism * (1 - theta1) * (1 - theta2)
            source = -3.0 * rho * np.exp(-psi_uns * psi) + v_term
            return np.vstack([dpsi, source - 2.0 * dpsi / x_safe])

        def bc(ya, yb):
            return np.array([ya[1], yb[0]])

        if guess is not None and hasattr(guess, 'x'):
            x_init, y_init = guess.x, guess.y
        else:
            x_init = np.logspace(-4, np.log10(x_max), 1000)
            x_init[0] = 1e-4
            x_init[-1] = x_max
            psi_init = np.where(x_init < 1, 1.5 - x_init**2 / 2, 1.0 / x_init)
            psi_init = np.maximum(psi_init, 0)
            dpsi_init = np.where(x_init < 1, -x_init, -1.0 / x_init**2)
            y_init = np.vstack([psi_init, dpsi_init])
        sol = solve_bvp(ode, bc, x_init, y_init, tol=1e-6,
                        max_nodes=200000, verbose=0)
        if not sol.success:
            return None
        x_d = np.logspace(-4, np.log10(x_max), 2000)
        x_d[0] = 1e-4
        x_d[-1] = x_max
        return x_d, sol.sol(x_d)[0], sol.sol(x_d)[1], sol

    zone3 = {}
    for mu0 in [0.01, 0.1, 1.0, 10.0, 100.0]:
        r = solve_3zone(PSI_SUN, mu0)
        if r is None:
            continue
        x, psi, dpsi, _ = r
        row = {"mu0": float(mu0)}
        row.update(_scan_row(x, psi, dpsi, PSI_SUN,
                             {"psi_center": float(np.interp(np.log10(1e-4), np.log10(x), psi))},
                             x_1au, x_wb))
        zone3[f"mu0_{mu0}"] = row
    # ISM-density scan at mu0 = 10
    for rho_ism_exp in [-30, -27, -25, -22, -20, -18, -15]:
        r = solve_3zone(PSI_SUN, 10.0, rho_ism=10.0**rho_ism_exp)
        if r is None:
            continue
        x, psi, dpsi, _ = r
        row = {"mu0": 10.0, "rho_ism_over_rho_sun": float(10.0**rho_ism_exp)}
        row.update(_scan_row(x, psi, dpsi, PSI_SUN,
                             {"psi_center": float(np.interp(np.log10(1e-4), np.log10(x), psi))},
                             x_1au, x_wb))
        zone3[f"ism_{rho_ism_exp}"] = row
    out["potentials"]["quadratic_3zone"] = zone3

    # Quartic 3-zone with physical lambda conversion
    FACTOR_SUN_NAT = (M_SUN_G / (4.0 * np.pi * M_PL_GEV)) ** 2
    quartic_3zone = {}
    prev3 = None
    for mu0 in [1e5, 3e5, 1e6]:
        r = solve_3zone(PSI_SUN, mu0, potential="quartic", guess=prev3)
        if r is None:
            prev3 = None
            continue
        x, psi, dpsi, sol3 = r
        prev3 = sol3
        row = {"mu0": float(mu0), "lambda_physical": float(mu0 / FACTOR_SUN_NAT)}
        row.update(_scan_row(x, psi, dpsi, PSI_SUN,
                             {"psi_center": float(np.interp(np.log10(1e-4), np.log10(x), psi))},
                             x_1au, x_wb))
        psi_min_sun = (3.0 / mu0) ** (1.0 / 3.0)
        psi_min_ism = (3.0e-27 / mu0) ** (1.0 / 3.0)
        m2_sun = 3.0 * mu0 * psi_min_sun**2
        m2_ism = 3.0 * mu0 * psi_min_ism**2
        row["compton_sun_au"] = float((1.0 / np.sqrt(m2_sun)) * R_SUN_SI / AU_SI)
        row["compton_ism_pc"] = float((1.0 / np.sqrt(m2_ism)) * R_SUN_SI / 3.086e16)
        quartic_3zone[f"mu0_{mu0:g}"] = row
    out["potentials"]["quartic_3zone"] = quartic_3zone

    # ---- Stellar-population scan (fixed physical lambda, mu0 = 1e5) ----
    lam_phys = 1e5 / FACTOR_SUN_NAT
    stellar_pop = {}
    prev_pop = None
    for name, m_msun, r_msun in [("M dwarf", 0.15, 0.20), ("K dwarf", 0.45, 0.55),
                                 ("G dwarf", 0.80, 0.80), ("Sun", 1.00, 1.00),
                                 ("F dwarf", 1.25, 1.15), ("A dwarf", 2.00, 1.70),
                                 ("Subgiant", 1.50, 3.50), ("Red giant", 1.20, 12.0),
                                 ("Massive", 5.00, 3.50), ("Very massive", 8.0, 5.0)]:
        M_g = m_msun * M_SUN_G
        R_g = r_msun * R_SUN_G
        R_si = r_msun * R_SUN_SI
        psi_uns_i = compute_source_params(M_g, R_g)[0]
        x_wb_i = 2646 * AU_SI / R_si
        x_max_i = max(1e6, 10.0 * x_wb_i)
        x, psi, dpsi, info = solve_bvp_conformal(
            psi_uns_i, "quartic", lam=1e5, x_max=x_max_i, guess=prev_pop)
        if x is None:
            x, psi, dpsi, info = solve_bvp_conformal(
                psi_uns_i, "quartic", lam=1e5, x_max=x_max_i, guess=None)
        if x is None:
            stellar_pop[name] = {"M_msun": m_msun, "R_msun": r_msun,
                                 "psi_uns": float(psi_uns_i), "failed": True}
            prev_pop = None
            continue
        prev_pop = info["sol"]
        S_sig, _ = compute_screening(x, psi, dpsi, psi_uns_i)
        s_1au = interp_log(x, S_sig, AU_SI / R_si)
        s_wb = interp_log(x, S_sig, x_wb_i)
        rs_pred_m13 = 2646.0 * (m_msun / 1.3) ** (1.0 / 3.0)
        rs_pred_tep = np.sqrt(G_SI * (m_msun * M_SUN_SI) / (2.0 * 5.0e-10)) / AU_SI
        screening_factor_analytic = (3.0 * (4.0 * np.pi)**2 * M_PL_GEV**2 /
                                     (lam_phys * (m_msun * M_SUN_G)**2)) ** (1.0 / 3.0)
        phi_star_wb = psi_uns_i * float(np.interp(np.log10(x_wb_i), np.log10(x), psi))
        amb_ratio = 2.85e-7 / abs(phi_star_wb) if phi_star_wb != 0 else float('inf')
        stellar_pop[name] = {
            "M_msun": float(m_msun), "R_msun": float(r_msun),
            "psi_uns": float(psi_uns_i),
            "s_sigma_surface": float(interp_log(x, S_sig, 1.0)),
            "s_sigma_1au": float(s_1au) if s_1au else None,
            "s_sigma_wb": float(s_wb) if s_wb else None,
            "rs_predicted_m13_au": float(rs_pred_m13),
            "rs_predicted_tep_au": float(rs_pred_tep),
            "gradient_screening_factor_analytic": float(screening_factor_analytic),
            "ambient_to_star_ratio_wb": float(amb_ratio),
            "failed": False,
        }
    out["potentials"]["stellar_population_scan"] = stellar_pop

    # ---- Nested screening hierarchy ----
    def solve_nested(psi_uns, phi_env, mu0=0.0, potential="linear",
                     x_max=1e6, n_init=800, guess=None):
        """Perturbation dPsi on environmental field phi_env; total field
        phi_env + psi_uns dPsi; dPsi -> 0 asymptotically."""
        def ode(x, y):
            d, dp = y
            x_safe = np.maximum(x, 1e-10)
            theta = 0.5 * (1.0 - np.tanh((x - 1.0) / 0.01))
            rho_eff = 3.0 * theta
            v_term = mu0 * d**3 if potential == "quartic" else mu0 * d
            source = -rho_eff * np.exp(-phi_env - psi_uns * d) + v_term
            return np.vstack([dp, source - 2.0 * dp / x_safe])

        def bc(ya, yb):
            return np.array([ya[1], yb[0]])

        if guess is not None and hasattr(guess, 'x'):
            sol = solve_bvp(ode, bc, guess.x, guess.y, tol=1e-6,
                            max_nodes=200000, verbose=0)
        else:
            x_init = np.logspace(-4, np.log10(x_max), n_init)
            x_init[0] = 1e-4
            x_init[-1] = x_max
            d_init = np.where(x_init < 1, 1.5 - x_init**2 / 2, 1.0 / x_init)
            d_init = np.maximum(d_init, 0)
            dp_init = np.where(x_init < 1, -x_init, -1.0 / x_init**2)
            sol = solve_bvp(ode, bc, x_init, np.vstack([d_init, dp_init]),
                            tol=1e-6, max_nodes=200000, verbose=0)
        if not sol.success:
            return None
        x_d = np.logspace(-4, np.log10(x_max), 2000)
        x_d[0] = 1e-4
        x_d[-1] = x_max
        return x_d, sol.sol(x_d)[0], sol.sol(x_d)[1], sol

    nested = {"description": ("phi_total = phi_env + delta_phi: one continuous "
                              "multi-scale profile, no boundaries between levels"),
              "chain": "cosmological -> galaxy(8 kpc) -> sun(1 AU) -> earth"}

    X_SUN_GAL = 8.0 / 50.0
    gal_variants = {}
    phi_gal_sun = None
    xg = pg = dg = None
    PSI_GAL_USED = None
    for lab, M_fac, R_kpc in [("fiducial", 1.0, 50.0), ("light", 0.5, 50.0),
                              ("heavy", 2.0, 50.0), ("compact", 1.0, 30.0),
                              ("extended", 1.0, 70.0)]:
        m_g = M_fac * M_GAL_G
        r_g = si_to_gev_length(R_kpc * 3.086e19)
        psi_g = compute_source_params(m_g, r_g)[0]
        x_sp = 8.0 / R_kpc
        r_v = solve_nested(psi_g, 0.0, x_max=1e7)
        if r_v is None:
            continue
        xv, pv, dv, _ = r_v
        pv_sun = interp_log(xv, pv, x_sp)
        phi_v = psi_g * pv_sun
        gal_variants[lab] = {"M_factor": M_fac, "R_kpc": R_kpc,
                             "psi_uns": float(psi_g), "phi_at_8kpc": float(phi_v)}
        if lab == "fiducial":
            xg, pg, dg = xv, pv, dv
            phi_gal_sun = phi_v
            psi_gal_sun = pv_sun
            PSI_GAL_USED = psi_g

    if phi_gal_sun is not None:
        v_c = 220e3
        gal_vals = [v["phi_at_8kpc"] for v in gal_variants.values()]
        nested["galaxy"] = {
            "phi_at_sun_pos": float(phi_gal_sun),
            "psi_at_sun_pos": float(psi_gal_sun),
            "model_variants": gal_variants,
            "phi_range": [float(min(gal_vals)), float(max(gal_vals))],
            "vc_crosscheck": float((v_c / C_SI) ** 2),
        }
        gal_field_vs_x = {}
        for xc, lab in [(0.02, "inner_1kpc"), (0.16, "sun_8kpc"),
                        (0.5, "25kpc"), (1.0, "edge_50kpc"), (3.0, "halo_150kpc")]:
            v = interp_log(xg, pg, xc)
            if v is not None:
                gal_field_vs_x[lab] = float(PSI_GAL_USED * v)
        nested["galaxy"]["field_vs_radius"] = gal_field_vs_x
        nested["galaxy"]["dlnA_sun_vs_inner"] = float(
            gal_field_vs_x["sun_8kpc"] - gal_field_vs_x["inner_1kpc"])
        nested["galaxy"]["dlnA_sun_vs_halo"] = float(
            gal_field_vs_x["sun_8kpc"] - gal_field_vs_x["halo_150kpc"])

        R_GAL_SI = 50.0 * 3.086e19
        dphi_dr = PSI_GAL_USED * abs(interp_log(xg, dg, X_SUN_GAL)) / R_GAL_SI
        nested["ambient_shear"] = {
            "a_amb_uniform_sphere_ms2": float(C_SI**2 * dphi_dr),
            "a_amb_isothermal_ms2": float(2.0 * v_c**2 / (8.0e3 * 3.086e16)),
            "g_TEP_ms2": 5.0e-10,
            "a_internal_wb_2646au_ms2": float(
                G_SI * 1.24 * M_SUN_SI / (2646.0 * AU_SI) ** 2),
            "interpretation": ("Ambient galactic shear at the solar circle "
                               "sets the acceleration scale matching g_TEP."),
        }

        # Sun perturbation on galactic background vs isolated
        sun_nested = {}
        phi_sun_1au = phi_sun_wb = None
        for lab, mu0_v, pot_t in [("V0", 0.0, "linear"),
                                  ("quad_0.01", 0.01, "linear"),
                                  ("quartic_1e5", 1e5, "quartic")]:
            if pot_t == "quartic":
                g_n = g_i = None
                r_s = r_i = None
                for m_sub in [1.0, 100.0, 1e4, mu0_v]:
                    r_s = solve_nested(PSI_SUN, phi_gal_sun, mu0=m_sub,
                                       potential="quartic", x_max=1e6, guess=g_n)
                    if r_s:
                        g_n = r_s[3]
                    r_i = solve_nested(PSI_SUN, 0.0, mu0=m_sub,
                                       potential="quartic", x_max=1e6, guess=g_i)
                    if r_i:
                        g_i = r_i[3]
            else:
                r_s = solve_nested(PSI_SUN, phi_gal_sun, mu0=mu0_v,
                                   potential=pot_t, x_max=1e6)
                r_i = solve_nested(PSI_SUN, 0.0, mu0=mu0_v,
                                   potential=pot_t, x_max=1e6)
            if r_s is None or r_i is None:
                continue
            xs, ps, ds, _ = r_s
            xi, pi_, di, _ = r_i
            s_n = x_1au**2 * abs(interp_log(xs, ds, x_1au))
            s_i = x_1au**2 * abs(interp_log(xi, di, x_1au))
            p_1au = interp_log(xs, ps, x_1au)
            p_wb = interp_log(xs, ps, x_wb)
            sun_nested[lab] = {
                "mu0": mu0_v, "potential": pot_t,
                "s_sigma_1au_nested": float(s_n),
                "s_sigma_1au_isolated": float(s_i),
                "gamma_1au_screened_source_nested": float(
                    compute_gamma(np.array([s_n]))[0]),
                "phi_perturb_1au": float(PSI_SUN * p_1au),
                "phi_perturb_wb": float(PSI_SUN * p_wb) if p_wb else None,
            }
            if lab == "V0":
                phi_sun_1au = PSI_SUN * p_1au
                phi_sun_wb = PSI_SUN * p_wb
        nested["sun"] = sun_nested

        r_s0 = solve_nested(PSI_SUN, phi_gal_sun, x_max=1e6)
        if r_s0 is not None:
            xs0, ps0, ds0, _ = r_s0
            psi_env_sun = phi_gal_sun / PSI_SUN
            idx_c = np.where(ps0 < psi_env_sun)[0]
            x_cross = float(xs0[idx_c[0]]) if len(idx_c) > 0 else None
            nested["sun"]["x_crossover_Rsun"] = x_cross
            nested["sun"]["crossover_AU"] = (x_cross * R_SUN_SI / AU_SI
                                             if x_cross else None)
            if phi_sun_1au:
                nested["sun"]["ambient_ratio_at_1au"] = float(
                    phi_gal_sun / phi_sun_1au)
            if phi_sun_wb:
                nested["sun"]["ambient_ratio_at_wb"] = float(
                    phi_gal_sun / phi_sun_wb)

        # Earth on galactic + solar background
        phi_env_earth = phi_gal_sun + (phi_sun_1au or 0.0)
        r_e = solve_nested(PSI_EARTH, phi_env_earth, x_max=1e6)
        if r_e is not None:
            xe, pe, de, _ = r_e
            phi_e_surf = PSI_EARTH * interp_log(xe, pe, 1.0)
            nested["earth"] = {
                "phi_env": float(phi_env_earth),
                "phi_perturb_surface": float(phi_e_surf),
                "s_sigma_surface": float(abs(interp_log(xe, de, 1.0))),
                "ambient_ratio_at_surface": float(phi_env_earth / phi_e_surf),
            }

            # Common-mode structure
            x_grid = np.logspace(np.log10(x_1au), np.log10(1e6), 200)
            dphi_local = PSI_SUN * np.interp(
                np.log10(x_grid), np.log10(xs0), ps0)
            local_variation = float(dphi_local[0] - dphi_local[-1])
            nested["common_mode"] = {
                "local_field_variation_1au_to_edge": local_variation,
                "ambient_field": float(phi_gal_sun),
                "variation_over_ambient": float(local_variation / phi_gal_sun),
            }

            # Decomposition at Earth surface
            phi_total_e = phi_env_earth + phi_e_surf
            nested["decomposition_at_earth_surface"] = {
                "phi_galactic": float(phi_gal_sun),
                "phi_solar_1au": float(phi_sun_1au or 0.0),
                "phi_earth_surface": float(phi_e_surf),
                "phi_total": float(phi_total_e),
                "galactic_fraction": float(phi_gal_sun / phi_total_e),
                "solar_fraction": float((phi_sun_1au or 0.0) / phi_total_e),
                "earth_fraction": float(phi_e_surf / phi_total_e),
            }

        # Cosmological baseline sensitivity
        for pc in [1e-7, 1e-6]:
            r_c = solve_nested(PSI_SUN, phi_gal_sun + pc, x_max=1e6)
            if r_c is None:
                continue
            xc_, pc_, dc_, _ = r_c
            s_c = x_1au**2 * abs(interp_log(xc_, dc_, x_1au))
            nested.setdefault("cosmo_baseline", {})[f"phi_cosmo_{pc:.0e}"] = {
                "s_sigma_1au": float(s_c)}

        nested["interpretation"] = (
            "Cassini bounds the Sun's local gradient perturbation — the "
            "constraint is local, not global. The field value at any point "
            "is dominated by the containing environment; local fluctuations "
            "are small common-mode perturbations on the ambient baseline "
            "and cancel in local comparisons. The galactic baseline varies "
            "with galactocentric radius, so the ambient clock rate is "
            "position-dependent across the galaxy.")

    out["nested_hierarchy"] = nested
    return out
