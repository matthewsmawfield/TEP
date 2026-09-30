#!/usr/bin/env python3
"""Binary-pulsar scalar-dipole constraint under the resolved nested operator.

Gate-B kill-shot (b1-pulsar-dipole): the standard executioner of
|alpha_0| ~ 1 scalar-tensor theories. The canonical kinetic-screening
operator of the corpus is

    P(X) = X - V + X|X|/Lambda^4,   g_t = c H_0/(2 beta_A^2),
    r* = sqrt(GM/g_t),   y(1 + y^2 (r*/r)^4) = 1,
    vertex factor q = S_Sigma(X_amb) = 1/(1 + X_amb),
    X_amb = (1 - y_partner)/y_partner + X_GAL.

Charge structure (corpus charge_vertex resolution, step_30):
the emitted scalar flux is conserved through the nonlinear shell —
suppressing the charge would shrink r* and erase the observed
wide-binary transition — so the bare DEF charge is NOT suppressed
as a flux; suppression enters through the fluctuation propagator.
The exponential coupling gives beta_0 = 0 — no spontaneous
scalarization — but does NOT remove the compactness/sensitivity
asymmetry: self-gravitating bodies carry alpha_A = alpha_0(1-2 s_A)
with s_NS ~ 0.1-0.3 and s_WD ~ 0, the same strong-EP mechanism that
bounds Brans-Dicke. The resolved radiation calculation (solved l=1
transmission on the pair-shell background) is step_69, which
supersedes the three readings below; this step retains them as the
ambient/background bookkeeping and as the historical bracket.

Three charge readings are retained for comparison only:

  1. 'static_vertex' — alpha_i = alpha_0 * S_Sigma(X_amb,i); static
     response suppression misapplied to a charge.
  2. 'canonical_sqrtZ' — alpha_i = alpha_0/sqrt(Z_i) at the pair
     scale; a point evaluation, not the solved propagator.
  3. 'vainshtein_emission' — (a/r*)^(9/2) impedance borrowed from
     cubic-Galileon screening; the specified P(X) operator contains
     no such term. Not a TEP reading.

The observational bound is the measured P_b residual versus the GR
quadrupole prediction, converted to a bound on |alpha_A - alpha_B|
through the standard dipole luminosity

    Pdot_dip = -(4 pi^2 G / (P_b c^3)) (mA mB/M) (Delta_alpha)^2 K_D,

with K_D an O(1) theory factor evaluated on [0.5, 4.0].

Systems: J1738+0333, J0348+0432, J0737-3039 (parameters from the
published timing solutions).
"""
import json
import math
from scipy.optimize import brentq

G = 6.674e-11
C = 299792458.0
AU = 1.496e11
M_SUN = 1.989e30
T_SUN = G * M_SUN / C**3          # 4.925e-6 s
H0 = 70e3 / 3.085677581e22
BETA_A = -1.0
ALPHA_0 = math.sqrt(2.0) * abs(BETA_A)   # DEF normalization of the bare charge
G_T = C * H0 / (2 * BETA_A**2)
X_GAL = 0.52


def r_star(M_kg):
    return math.sqrt(G * M_kg / G_T)


def y_profile(x):
    if x <= 0:
        return 0.0
    if x >= 50.0:
        return 1.0
    return brentq(lambda y: y * (1.0 + y * y / x**4) - 1.0,
                  1e-30, 1.0, xtol=1e-14)


def ambient_x(M_partner_kg, s):
    """Kinetic ambient at the pair scale inside the partner's shell."""
    y = y_profile(s / r_star(M_partner_kg))
    return (1.0 - y) / y + X_GAL


def semimajor(mA, mB, P_b):
    return (G * (mA + mB) * P_b**2 / (4 * math.pi**2)) ** (1.0 / 3.0)


def pdot_gr(mA, mB, P_b, e):
    """Peters quadrupole prediction, masses in kg."""
    f = (1 + 73 * e**2 / 24 + 37 * e**4 / 96) / (1 - e**2) ** 3.5
    return (-(192 * math.pi / 5) * (2 * math.pi / P_b) ** (5.0 / 3.0)
            * G ** (5.0 / 3.0) / C**5 * mA * mB / (mA + mB) ** (1.0 / 3.0) * f)


def evaluate(name, mA_msun, mB_msun, P_b_s, e, pdot_resid, pdot_resid_err):
    mA, mB = mA_msun * M_SUN, mB_msun * M_SUN
    a = semimajor(mA, mB, P_b_s)
    v = 2 * math.pi * a / P_b_s

    # partner-shell ambients at the pair scale (self-field excluded by rule)
    XA = ambient_x(mB, a)          # ambient at A: B's shell + Galactic
    XB = ambient_x(mA, a)          # ambient at B: A's shell + Galactic

    qA, qB = 1.0 / (1.0 + XA), 1.0 / (1.0 + XB)          # static vertices
    ZA, ZB = 1.0 + 3.0 * XA, 1.0 + 3.0 * XB               # perturbation norm
    aA = {r: ALPHA_0 * f for r, f in
          [("static_vertex", qA), ("canonical_sqrtZ", 1.0 / math.sqrt(ZA)),
           ("vainshtein_emission", 1.0)]}
    aB = {r: ALPHA_0 * f for r, f in
          [("static_vertex", qB), ("canonical_sqrtZ", 1.0 / math.sqrt(ZB)),
           ("vainshtein_emission", 1.0)]}

    # emission-impedance suppression of the power (cubic-Galileon class)
    rstar_pair = r_star((mA + mB))
    vains_power_supp = (a / rstar_pair) ** 4.5

    pdGR = pdot_gr(mA, mB, P_b_s, e)
    resid_allow = abs(pdot_resid) + 2 * abs(pdot_resid_err)

    rows = {}
    for r in aA:
        dAlpha = abs(aA[r] - aB[r])
        if r == "vainshtein_emission":
            # bare asymmetry is zero (beta_0 = 0); the emitted power is
            # additionally shell-suppressed relative to an unsuppressed
            # fiducial dipole at Delta_alpha = ALPHA_0.
            dAlpha = ALPHA_0
        coeff = (4 * math.pi**2 * G / (P_b_s * C**3)) * (mA * mB / (mA + mB))
        pd_dip = coeff * dAlpha**2
        if r == "vainshtein_emission":
            pd_dip *= vains_power_supp
        # bound on Delta_alpha for K_D in [0.5, 4]
        bound = math.sqrt(resid_allow / (coeff * 0.5))
        bound_hi = math.sqrt(resid_allow / (coeff * 4.0))
        rows[r] = {
            "alpha_A": aA[r], "alpha_B": aB[r],
            "delta_alpha_pred": dAlpha,
            "pdot_dip_pred_s_per_s": pd_dip,
            "pdot_dip_over_gr": pd_dip / abs(pdGR),
            "allowed_residual_s_per_s_2sigma": resid_allow,
            "delta_alpha_bound_KD_0p5": bound,
            "delta_alpha_bound_KD_4p0": bound_hi,
            "margin_vs_bound_KD_0p5": bound / max(dAlpha, 1e-300),
            "passes": pd_dip <= resid_allow,
        }

    return {
        "system": name,
        "masses_msun": [mA_msun, mB_msun], "P_b_s": P_b_s, "ecc": e,
        "semimajor_m": a, "v_over_c": v / C,
        "pair_shell_radius_m": rstar_pair,
        "a_over_rstar": a / rstar_pair,
        "ambient_X_at_A": XA, "ambient_X_at_B": XB,
        "vertex_qA": qA, "vertex_qB": qB,
        "pdot_gr_s_per_s": pdGR,
        "readings": rows,
    }


def run():
    systems = [
        # name, mA, mB (Msun), P_b (s), e, observed residual vs GR, err
        ("J1738+0333", 1.46, 0.181, 8.5 * 3600, 3e-7, 0.18e-14, 0.32e-14),
        ("J0348+0432", 2.01, 0.172, 2.46 * 3600, 2.4e-6, 0.05e-12, 0.20e-12),
        ("J0737-3039", 1.337, 1.250, 2.45 * 3600, 0.0877, 0.0e-12, 0.02e-12),
    ]
    out = {
        "conventions": {
            "operator": "P = X - V + X|X|/Lambda^4; r* = sqrt(GM/g_t); "
                        "y(1+y^2(r*/r)^4)=1; X_amb=(1-y)/y+X_GAL",
            "g_t": G_T, "X_gal": X_GAL, "alpha_0_DEF": ALPHA_0,
            "bare_charge_asymmetry": "superseded",
            "bare_charge_note": "beta_0 = 0 removes spontaneous "
                "scalarization only; sensitivities s_A (strong-EP "
                "binding response) still make NS and WD charges differ. "
                "The resolved radiation calculation is step_69, which "
                "conserves the charge and suppresses the l=1 mode "
                "through the shell impedance.",
            "charge_vertex": "emitted flux conserved through the shell "
                "(step_30 resolution); propagator solved in step_69 "
                "(Z_par = 1+6u, Z_perp = 1+2u; T_amp ~ a/r*)",
            "dipole_formula": "Pdot_dip = -(4 pi^2 G/(P_b c^3))(mA mB/M)"
                "(Delta_alpha)^2 K_D, K_D scanned on [0.5,4]",
        },
        "systems": [evaluate(*s) for s in systems],
    }

    # --- stellar fifth force (b2): interior scalar-gravity suppression ---
    # A fluid element inside a star is embedded in the star's own nonlinear
    # shell: the hierarchical reading gives R_local = q_env^2 y ~ y^3, so
    # G_eff/G - 1 = 2 beta_A^2 R_local. Evaluated at representative radii.
    R_SUN = 6.96e8
    stellar = []
    for name, M_ms, R_rs, frac in [
        ("solar_core", 1.0, 1.0, 0.25),
        ("solar_surface", 1.0, 1.0, 1.0),
        ("rgb_core_0p5Msun_12Rsun_envelope", 0.5, 12.0, 0.001),
        ("rgb_envelope_midpoint", 0.5, 12.0, 0.5),
    ]:
        x = frac * R_rs * R_SUN / r_star(M_ms * M_SUN)
        y = y_profile(x)
        R_loc = y**3
        stellar.append({
            "location": name, "x_over_rstar": x, "y": y,
            "G_eff_over_G_minus_1": 2 * BETA_A**2 * R_loc,
        })
    out["stellar_fifth_force"] = {
        "rule": "hierarchical interior reading: R_local ~ y^3 with y the "
                "star's own flux-conserving profile at the element radius",
        "rgb_tip_bound_G_eff_over_G_minus_1": 0.02,
        "points": stellar,
        "verdict": "passes by ~20+ orders of magnitude at every interior "
                   "point; no stellar-structure constraint is active",
    }

    # corpus-level verdict
    worst = min(
        row["margin_vs_bound_KD_0p5"]
        for s in out["systems"] for row in [s["readings"]["canonical_sqrtZ"]]
    )
    out["verdict"] = {
        "tightest_margin": worst,
        "tightest_reading": "canonical_sqrtZ",
        "superseded_by": "step_69_dipole_transmission",
        "result": (
            "Historical bracket retained for comparison. The resolved "
            "calculation is step_69: charge conserved through the shell, "
            "asymmetry sourced by sensitivities alpha_i = alpha_0(1-2s_i), "
            "emission suppressed by the solved l=1 transmission "
            "T_amp ~ a/r*; all systems pass by 2-4 orders of magnitude. "
            "The beta_0 = 0 'identical charges' premise is removed: it "
            "confused the absence of scalarization with the absence of "
            "sensitivities."
        ),
    }
    return out


if __name__ == "__main__":
    res = run()
    import os
    outdir = os.path.join(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))), "results")
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, "step_31_binary_pulsar_operator.json")
    with open(path, "w") as f:
        json.dump(res, f, indent=2)
    for s in res["systems"]:
        cz = s["readings"]["canonical_sqrtZ"]
        sv = s["readings"]["static_vertex"]
        print(f"{s['system']}: a/r* = {s['a_over_rstar']:.1e}, "
              f"X_amb(A,B) = {s['ambient_X_at_A']:.1e},{s['ambient_X_at_B']:.1e}")
        print(f"   sqrtZ:  dAlpha={cz['delta_alpha_pred']:.2e} "
              f"vs bound {cz['delta_alpha_bound_KD_0p5']:.2e} "
              f"-> margin {cz['margin_vs_bound_KD_0p5']:.0f}x "
              f"(dip/GR={cz['pdot_dip_over_gr']:.1e})")
        print(f"   static: dAlpha={sv['delta_alpha_pred']:.2e} "
              f"-> margin {sv['margin_vs_bound_KD_0p5']:.0f}x")
    print(f"\nVERDICT: {res['verdict']['result'][:120]}...")
    print(f"Saved {path}")
