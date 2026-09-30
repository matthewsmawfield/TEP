#!/usr/bin/env python3
"""Master-sector admissibility ledger (plan AUD-0 / T0.1-T0.3, AUD-1 / T1.0).

Evaluates candidate functional forms for the master sector
{P(X,phi), V(phi), B(phi)} against the constraint ledger assembled from
Paper 0 (Sections 2.2, 4, 7, 8; Appendix E R1-R14), the recovery
constraints F1-F6, and tep-rules.md Rules 4, 9-13, 16, 20, 22, 23.

Nothing here edits manuscripts.  All candidate forms are evaluated;
each rejection carries its failure mode (plan doctrine D7).

Sections:
  A. Kinetic sector P_X  -- admissible-space scan, deep-MOND tail
     exponent uniqueness, stability, benchmark propagation, the
     cosmological branch sequence, and calibration of the single new
     coefficient k against a_0.
  B. Potential sector V(phi) -- minimum existence on the beta_A = -1
     branch, density-scaling of the effective mass (environmental
     amplitude screening), the deep-well floor V_0, and cosmological
     drive compatibility.
  C. Disformal sector B(phi) -- the admissibility gate table.
  D. Cosmological-branch screening sequence under each P_X candidate.

Conventions (step_27 / step_32):
  xi  = |X| / Lambda^4 = (g / g_t)^2  on the static branch,
  Lambda^4 = M_Pl^2 H0^2,  g_t = c H0 / (2 beta_A^2) = 3.4e-10 m/s^2.
The flux form of the static spherical field equation
    r^2 P_X(xi) phi' = C(r)  = |beta_A| M(r) / (4 pi M_Pl)
is exact (no BVP needed for the spherical RAR); the deep-MOND tail is
read off it algebraically.
"""
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

# ---------------------------------------------------------------------------
# constants (SI where applicable)
# ---------------------------------------------------------------------------
c    = 2.998e8                       # m/s
G    = 6.674e-11                     # m^3 kg^-1 s^-2
H0   = 70.0e3 / 3.086e22             # s^-1
M_Pl = 2.435e18                      # GeV (reduced)
beta_A = -1.0

g_t  = c * H0 / (2.0 * beta_A**2)    # shear threshold, m/s^2
A0_MOND = 1.2e-10                    # empirical deep-MOND scale, m/s^2

OUT = {"step": "master_sector_admissibility",
       "constants": {"H0_s-1": H0, "g_t": g_t, "beta_A": beta_A,
                     "Lambda4": "M_Pl^2 H0^2",
                     "a0_empirical": A0_MOND}}


# ===========================================================================
# A. Kinetic sector: admissible P_X space
# ===========================================================================
#
# Candidate list (xi = |X|/Lambda^4;  P_X = dP/dX is the stiffness):
#   baseline     1 + 2 xi                          (v6.0 incumbent)
#   two_branch   k sqrt(xi) + 2 xi                 (plan AUD-1 candidate)
#   twobranch_u  1 + k sqrt(xi) + 2 xi             (two-branch + unit term)
#   exp_interp   (1 - e^{-k sqrt(xi)}) + 2 xi      (smooth interpolant)
#   aqual_only   k sqrt(xi)                        (screening removed)
#   power-law family k xi^p + 2 xi for p in {1/3,1/2,2/3,1}
#
# Analytic anchors:
#   Deep tail (P_X ~ k xi^p):  flux k (phi'^2/2L^4)^p phi' = C/r^2
#     => phi' ~ r^{-2/(2p+1)},  a_phi ~ r^{-2/(2p+1)}.
#   Flat rotation asymptote (a_phi ~ 1/r)  <=>  p = 1/2.  UNIQUE power.
#   For p = 1/2:  a_phi^2 = a_eff g_N,
#     a_eff = 2 sqrt(2) |beta_A|^3 Lambda^2 / (k M_Pl)
#           = 2 sqrt(2) |beta_A|^3 H0 / k          (natural units)
#           = 4 sqrt(2) |beta_A| g_t / k           (SI).
#   Radial sound speed (static branch X<0): c_s^2 = P_X/(P_X + 2X P_XX)
#     two-branch: P_X + 2X P_XX = 2 k sqrt(xi) + 6 xi > 0,
#     c_s^2 -> 1/2 as xi -> 0, -> 1/3 as xi -> inf.
#   Static-branch energy rho_phi = 2X P_X - P:
#     small-xi term gives -(4k/3)|X|^{3/2}/Lambda^2 < 0
#     (AQUAL-type negative static energy; cosmological branch positive).

def a_eff_from_k(k):
    """SI deep-MOND acceleration scale for small-xi coefficient k."""
    return 4.0 * np.sqrt(2.0) * abs(beta_A) * g_t / k


def k_from_a_eff(a_eff):
    return 4.0 * np.sqrt(2.0) * abs(beta_A) * g_t / a_eff


K_STAR = k_from_a_eff(A0_MOND)         # calibrated coefficient

CANDS = {
    "baseline":    lambda xi: 1.0 + 2.0 * xi,
    "two_branch":  lambda xi: K_STAR * np.sqrt(xi) + 2.0 * xi,
    "twobranch_u": lambda xi: 1.0 + K_STAR * np.sqrt(xi) + 2.0 * xi,
    "exp_interp":  lambda xi: (1.0 - np.exp(-K_STAR * np.sqrt(xi))) + 2.0 * xi,
    "aqual_only":  lambda xi: K_STAR * np.sqrt(xi),
}
POWER_FAMILIES = {p: (lambda xi, p=p: K_STAR * xi**p + 2.0 * xi)
                  for p in (1.0/3.0, 0.5, 2.0/3.0, 1.0)}


def legendre_factor(fun, xi):
    """Radial stability/Legendre factor P_X + 2 X P_{,XX}.

    With P_X = f(xi), xi = |X|/Lambda^4:
      static branch (X<0):  P_{,XX} = f' * dxi/dX = -f'/Lambda^4
        => P_X + 2 X P_{,XX} = f + 2 xi f'
      timelike branch (X>0): dxi/dX = +1/Lambda^4 gives the same
        expression f + 2 xi f'  (sign-symmetric)."""
    xi = np.maximum(np.asarray(xi, dtype=float), 1e-30)
    d = (fun(xi * 1.001) - fun(xi * 0.999)) / (0.002 * xi)
    return fun(xi) + 2.0 * xi * d


def tail_exponent(fun, C=1.0, r_lo=1e2, r_hi=1e6):
    """Solve P_X(xi) phi' = C/r^2 on a radius grid; return log-slope of
    a_phi = phi' (normalisation constants drop out of the slope) and of
    a_phi^2 vs the Newtonian g_N = C/r^2 proxy."""
    r = np.logspace(np.log10(r_lo), np.log10(r_hi), 400)
    phi = np.empty_like(r)
    for i, ri in enumerate(r):
        rhs = C / ri**2
        # solve P_X(w^2/(2)) * w = rhs  for w = phi'/sqrt(2)L^2-scaled:
        # work in xi directly: P_X(xi)*sqrt(2 xi) = rhs/... -- use w-form
        f = lambda w: fun(max(w * w / 2.0, 1e-300)) * w - rhs
        hi = 1.0
        while f(hi) < 0:
            hi *= 10.0
            if hi > 1e30:
                break
        try:
            phi[i] = brentq(f, 1e-30, hi, xtol=1e-40, rtol=1e-12)
        except Exception:
            phi[i] = np.nan
    ok = np.isfinite(phi)
    slope = np.polyfit(np.log(r[ok]), np.log(phi[ok]), 1)[0]
    return slope, r[ok], phi[ok]


secA = {"k_star_from_a0": float(K_STAR),
        "a_eff_formula": "a_eff = 4 sqrt(2) |beta_A| g_t / k",
        "deep_tail_exponent_rule": "a_phi ~ r^{-2/(2p+1)}; flat tail requires p = 1/2 uniquely",
        "candidates": {}}

XI_GRID = np.logspace(-12, 12, 200)

for name, fun in CANDS.items():
    rec = {}
    px = np.array([fun(x) for x in XI_GRID])
    rec["P_X_positive"] = bool(np.all(px > 0))
    leg = np.array([legendre_factor(fun, x) for x in XI_GRID])
    rec["legendre_positive"] = bool(np.all(np.array(leg) > 0))
    cs2 = px / np.maximum(leg, 1e-300)
    rec["c_s2_range"] = [float(np.min(cs2)), float(np.max(cs2))]
    slope, _, _ = tail_exponent(fun)
    rec["tail_slope_phi"] = float(slope)
    rec["flat_tail"] = bool(abs(slope + 1.0) < 0.05)
    secA["candidates"][name] = rec

# power-law family: exponent uniqueness demonstration
secA["power_law_scan"] = {}
for p, fun in POWER_FAMILIES.items():
    slope, _, _ = tail_exponent(fun)
    secA["power_law_scan"][f"p={p:.4f}"] = {
        "analytic_tail": -2.0 / (2.0 * p + 1.0),
        "numeric_tail": float(slope),
        "flat": bool(abs(slope + 1.0) < 0.05)}

# branch crossover xi* where k sqrt(xi) = 2 xi  ->  xi* = (k/2)^2
secA["branch_crossover_xi"] = float((K_STAR / 2.0) ** 2)

# which environments sit on which branch (xi = (g/g_t)^2)
env = {
    "Sun conjunction limb 1.6 Rsun": 107.0,
    "Earth surface":                 9.82,
    "Sun field at Saturn 9.5 AU":    6.57e-5,
    "Earth-Moon (LLR)":              2.70e-3,
    "Sun field at 1 AU":             5.93e-3,
    "Wide binary 2646 AU":           1.05e-9,
    "MW solar circle (ambient)":     2.0e-10,
    "GC core":                       7.0e-8,
    "Void / cosmic web":             1.0e-12,
}
secA["environment_branches"] = {}
xi_star = secA["branch_crossover_xi"]
for name, g in env.items():
    xi = (g / g_t) ** 2
    pxb, pxt = float(CANDS["baseline"](xi)), float(CANDS["two_branch"](xi))
    secA["environment_branches"][name] = {
        "g": g, "xi": xi,
        "P_X_baseline": pxb, "P_X_two_branch": pxt,
        "S_Sigma_baseline": 1.0 / pxb, "S_Sigma_two_branch": 1.0 / pxt,
        "benchmark_shift": pxt / pxb,
        "branch": "large-X" if xi > xi_star else "small-X"}

OUT["kinetic_sector"] = secA


# ===========================================================================
# B. Potential sector V(phi) family scan
# ===========================================================================
#
# Requirements on V (ledger):
#   V-1 (F1):  admissible minimum on the beta_A = -1 matter-coupled
#              branch  V_eff(u) = V(u) + rho A(u), A = e^{-u}.
#   V-2:       density-dependent m_eff -- fixed-mass forms kill the
#              environmental amplitude screening (quadratic excluded).
#   V-3 (R3):  locally flat floor at deep u -- constant V_,u + Q_m for
#              the rolling-floor self-consistency (closure 7e-14);
#              steep/exponential forms fail secularly.
#   V-4 (R4):  closed-static benchmark reconstruction V -> 2 M_Pl^2/a^2
#              floor, V_0 identified with the curvature scale.
#   V-5:       cosmological range Delta u ~ ln(1+z) while flattening to
#              the floor (master-family compatibility).
#   V-6 (Rule 20):  floor V_0 supplies the scale-invariant well core
#              density; N_min > 0 (mechanism: source starvation --
#              matter term ~ rho A_,phi -> 0 as u -> inf).
#   V-7:       unsourced void sector must not run secularly (|du| < 1e-2).
#   V-8:       V''>0 at the matter-domain minimum (stability, Sec. 4).

LAM_REF = 7.5e-71          # reference quartic coupling (Appendix E, step_03)
U_S = 10.0               # transition scale of the master family (R3)


def V_master(u, lam=LAM_REF, V0=1.0, u_s=U_S):
    """Master family  V(u) = lam u^4/4 e^{-(u/u_s)^4} + V0 e^{-(u_s/u)^4}
    in M_Pl^4 units (V0 in units of the floor value)."""
    u = np.asarray(u, dtype=float)
    matter = lam * u**4 / 4.0 * np.exp(-(u / u_s)**4)
    floor = V0 * np.exp(-(u_s / np.maximum(u, 1e-12))**4)
    return matter + floor


from scipy.optimize import minimize_scalar


def effective_minimum(Vfun, rho):
    """Locate min of V_eff(u) = V(u) + rho e^{-u} on u>0 by bounded
    minimisation in t = ln u (robust where V_eff is numerically flat).
    Returns (u_min, status): status is 'interior', 'ambient_pinned'
    (u -> 0 boundary; the dense-side drive stalls at the ambient field)
    or 'runaway' (u -> inf boundary; no admissible minimum)."""
    T_LO, T_HI = -140.0, 7.0
    f = lambda t: Vfun(np.exp(t)) + rho * np.exp(-np.exp(t))
    res = minimize_scalar(f, bounds=(T_LO, T_HI), method="bounded")
    umin = float(np.exp(res.x))
    if res.x > T_HI - 0.5:
        return umin, "runaway"
    if res.x < T_LO + 0.5:
        return umin, "ambient_pinned"
    return umin, "interior"


def meff_scaling(Vfun, rhos):
    """m_eff^2 = V''(u_min) + rho e^{-u_min} (the source curvature is
    added analytically -- differencing the full V_eff loses it to
    float noise at void densities).  Returns density exponent of m_eff,
    du_min/d rho sign, and the interior-minimum fraction."""
    us, ms, statuses = [], [], []
    for rho in rhos:
        umin, status = effective_minimum(Vfun, rho)
        statuses.append(status)
        if status != "interior":
            continue
        h = max(umin * 1e-4, 1e-30)
        d2 = (Vfun(umin + h) - 2.0 * Vfun(umin) + Vfun(umin - h)) / h**2 \
            + rho * np.exp(-umin)
        us.append(umin); ms.append(d2)
    if len(us) < 2:
        return None, None, 0.0, statuses
    us = np.asarray(us); ms = np.asarray(ms)
    ok_rhos = np.array([r for r, s in zip(rhos, statuses)
                        if s == "interior"])
    p = np.polyfit(np.log(ok_rhos), np.log(np.sqrt(np.abs(ms))), 1)[0]
    du_sign = float(np.sign(us[-1] - us[0]))
    frac = float(len(us) / len(rhos))
    return float(p), du_sign, frac, statuses


RHOS = np.array([1e-120, 1e-95, 1e-80, 1e-60])  # M_Pl^4 units: void -> compact

VCANDS = {
    "quadratic m^2 u^2/2":
        (lambda u: 0.5e-84 * u**2,
         "fixed mass -> density-independent lambda_C; kills environmental screening"),
    "inverse_power L^4(1+1/u)":
        (lambda u: 1e-84 * (1.0 + 1.0 / np.maximum(u, 1e-30)),
         "beta_A=-1 branch: V' < 0 and source < 0 -> no admissible minimum (F1)"),
    "quartic lam u^4/4":
        (lambda u: LAM_REF * u**4 / 4.0,
         "reference amplitude-sector form (v6.0)"),
    "exponential L_V^4 e^u":
        (lambda u: 1e-84 * np.exp(np.minimum(u, 300.0)),
         "runaway + coupling; min exists at e^{2u} = rho/L_V^4"),
    "symmetron -mu^2 u^2/2 + lam u^4/4":
        (lambda u: -1e-84 * u**2 / 2.0 + LAM_REF * u**4 / 4.0,
         "wrong environmental direction: symmetric phase in dense media"),
    "master_family quartic+floor":
        (lambda u: V_master(u),
         "V = lam u^4/4 e^{-(u/u_s)^4} + V0 e^{-(u_s/u)^4}; u_s ~ 9-10 (R3)"),
    "constant floor V0":
        (lambda u: np.ones_like(np.asarray(u, dtype=float)),
         "u->inf limit of the master family; no matter-domain minimum alone"),
}

secB = {"candidates": {}}
for name, (fun, note) in VCANDS.items():
    rec = {"note": note}
    umin, status = effective_minimum(fun, 1e-80)
    rec["minimum_at_rho_1e-80"] = {"u_min": umin, "status": status}
    p, sgn, frac, statuses = meff_scaling(fun, RHOS)
    rec["m_eff_density_exponent"] = p if p is not None else \
        "no density range admits interior minima"
    rec["du_min_drho_sign"] = sgn
    rec["density_fraction_with_minima"] = frac
    rec["minimum_statuses"] = statuses
    # floor / flatness diagnostic: relative slope at deep u
    u_deep = np.array([50.0, 100.0, 200.0])
    vv = fun(u_deep)
    rec["deep_u_flatness_dV_over_V"] = (
        float(np.abs(np.gradient(vv, u_deep)[-1] / max(vv[-1], 1e-300)))
        if np.all(np.isfinite(vv)) else None)
    secB["candidates"][name] = rec

# master-family structural checks: the quartic term governs only while
# lam u^4/4 e^{-(u/u_s)^4} > V0 e^{-(u_s/u)^4}; the crossover u_c marks
# where the floor takes over.  For V0 ~ M_Pl^4 the floor dominates above
# u ~ few -- the matter-domain minima (u_min << 1) are safely quartic,
# but the well-interior field (u ~ u_s ~ 9-10, R3) lives on the floor,
# which is exactly the Rule-20 role.
def u_crossover(V0, u_s=U_S, lam=LAM_REF):
    f = lambda u: (lam * u**4 / 4.0 * np.exp(-(u / u_s)**4)
                   - V0 * np.exp(-(u_s / u)**4))
    grid = np.logspace(-3, 2, 4000)
    fv = f(grid)
    idx = np.where(np.diff(np.sign(fv)) != 0)[0]
    if len(idx) == 0:
        return None
    i = int(idx[0])
    return float(brentq(f, grid[i], grid[i + 1]))


secB["master_family"] = {
    "u_crossover_V0_eq_1": u_crossover(1.0),
    "u_crossover_V0_eq_1e-30": u_crossover(1e-30),
    "matter_domain": "u << u_c: quartic minima (u_min ~ 5e-4 at Earth "
                     "density) sit deep inside the quartic regime",
    "floor_domain": "u > u_c: V -> V0, locally flat (R3 rolling-floor "
                    "condition); well interiors land here",
    "V_at_u_100_over_V0": float(V_master(100.0)),
    "V_at_u_200_over_V0": float(V_master(200.0)),
    "u_transition_R3": "u ~ 9-10 at 1.01 r_h (step_16_interior_roll)",
    "V0_from_closed_static": "V_0 = 2 M_Pl^2 / a^2 (step_12, R4)",
    "source_starvation_note": ("matter coupling rho A_,phi ~ -rho e^{-u} -> 0 "
                               "at u >> 1: the deep-well field self-pins; "
                               "candidate mechanism for N_min > 0 (Rule 20), "
                               "to be verified by the T-B1 well-family solve")}
secB["global_floor_question"] = (
    "In the master family the floor V(u -> inf) = V0 is the SAME limit the "
    "ambient power-law branch approaches (A_clock -> 0 <-> u -> inf): the "
    "floor is global, not merely effective.  On a flat floor V_,u -> 0 the "
    "roll asymptotes -- A -> 0 is approached over infinite coordinate time "
    "and never completed, while well cores pin at finite u through source "
    "starvation.  Whether this satisfies Rule 20's stated N_min > 0 "
    "verbatim (finite minimum rather than asymptotic approach) is the open "
    "structural question for T0.2/T-B1.")

OUT["potential_sector"] = secB


# ===========================================================================
# C. Disformal sector B(phi) admissibility gate table
# ===========================================================================
#
# Gates (Paper 0 SS4; step_52 margins; GW170817; Rule 13):
#   B-1 signature:   B (dphi)^2 > -A^2  on realized profiles
#   B-2 hyperbolic:  B >= 0 -> Z_t, Z_s >= 1 unconditionally;
#                    B < 0 -> pointwise admissibility (all realized
#                    negative-branch reconstructions fail: Q -> -1, R5)
#   B-3 GW170817:    B -> 0 at phi ~ 0; path integral of B(dphi)^2 below
#                    the |c_g - c_gamma| ~ few e-15 bound
#   B-4 holonomy:    admissible residual circulation requires spatially
#                    varying (B/A^2)(n.dphi) around a loop (Rule 13)
#   B-5 Rule 16:     g-sector wave remains a pure tensor ripple; B lives
#                    in the matter metric and must not inject scalar
#                    polarisation at leading order
secC = {
    "candidates": {
        "B ~ phi^2 (pure)": {
            "B(0)=0": True,
            "verdict": "REJECTED",
            "failure": "nonzero wherever phi != 0; fails the Solar-System "
                       "cone-tightness requirement (Paper 29 exclusion)"},
        "B ~ phi^5": {
            "B(0)=0": True,
            "verdict": "REJECTED",
            "failure": "Lorentzian signature of the matter metric lost "
                       "along high-z propagation paths (z >~ 3.4); "
                       "fails GW170817 window (Paper 29 exclusion)"},
        "bump B0 phi^2 exp(-phi^2/phi_c^2)": {
            "B(0)=0": True,
            "B_galactic_nonzero": True,
            "signature": "margins: GW-host wall 1e15 above multimessenger "
                         "bound (step_52)",
            "verdict": "ADMISSIBLE_CANDIDATE",
            "open": "B0 window [2.1e21, 6.1e164] dimensionless; fixed-action "
                    "normalisation unresolved (R6)"},
        "Paper-28 envelope B0 u^2/(1+u^2) exp(-u^4/2 sigma_B^4)": {
            "B(0)=0": True,
            "verdict": "ADMISSIBLE_CANDIDATE",
            "open": "prescribed strong-field form; weak-field limit "
                    "~ B0 u^2; unique normalisation part of the closure "
                    "problem (Sec. 2.2)"},
        "B < 0 branch": {
            "verdict": "CONDITIONALLY_EXCLUDED",
            "failure": "pointwise admissibility fails on every realized "
                       "reconstruction to date (volume-balance b<0 drives "
                       "Q -> -1: degenerate matter cone, superluminal "
                       "longitudinal photons; R5, step_52)"},
    },
    "structural_conditions": {
        "signature": "B (dphi)^2 > -A^2",
        "hyperbolicity": "Z_t = 1 + B(1+w) rho/A^2 > 0; "
                         "Z_s = 1 + B(1-w) rho/(3 A^2) > 0",
        "causality": "matter cone inside gravitational cone on the "
                     "realized deformation",
        "holonomy_kernel": "d(delta sigma) = -(B/A^2) d(phi_dot/N^2) ^ d phi"},
    "note": ("B(0) = 0 and a mid-field bump are both enforced by data: "
             "GW170817 ties the cosmological ambient to the conformal "
             "limit while Rule 13 requires non-exact transport to remain "
             "possible somewhere.  The surviving forms are precisely the "
             "bump/envelope family -- constrained, not free.")}

OUT["disformal_sector"] = secC


# ===========================================================================
# D. Cosmological branch: xi_cosmo(z) and S_Sigma sequence
# ===========================================================================
#
# phi = -M_Pl ln A, ambient map A_clock = (1+z)^{-1}  =>
#   phi_dot = -M_Pl H  =>  X = M_Pl^2 H^2 / 2  =>  xi_cosmo = H^2/(2 H0^2)
# with the conformal-image bookkeeping H^2/H0^2 = Om (1+z)^3 + OL.
OM, OL = 0.3, 0.7


def H2_over_H02(z):
    return OM * (1.0 + z) ** 3 + OL


secD = {"xi_cosmo": "H^2(z) / (2 H0^2)",
        "sequence": {}}
for z in (0.0, 0.5, 1.0, 2.0, 3.0, 10.0, 100.0):
    xi = H2_over_H02(z) / 2.0
    row = {"xi": float(xi)}
    for name in ("baseline", "two_branch", "exp_interp"):
        px = CANDS[name](xi)
        row[f"S_Sigma_{name}"] = float(1.0 / px)
        row[f"G_eff_over_G_{name}"] = float(1.0 + 2.0 * beta_A**2 / px)
    secD["sequence"][f"z={z}"] = row

secD["note"] = (
    "Under the two-branch completion the late-time boost is much weaker "
    "than the incumbent's G_eff ~ 2G at z=0 -- the small-X branch makes "
    "the cosmological branch stiffer (P_X ~ 12.3 vs 2).  The sign of the "
    "corrected growth prediction relative to the below-LCDM clustering "
    "data is the AUD-3 discriminant; X_env on linear scales must come "
    "from the constraint-slice landscape, not from H(z) alone (plan R3.2).")

OUT["cosmological_branch"] = secD


# ===========================================================================
# E. Summary verdicts
# ===========================================================================
OUT["verdicts"] = {
    "P_X": ("two-branch  k sqrt(xi) + 2 xi  is the unique minimal "
            "completion satisfying: flat-tail rotation asymptote "
            "(p=1/2 unique among power laws), large-X screening "
            "unchanged (benchmarks pass identically), stability "
            "P_X>0 and Legendre>0 on both branches, and the WB/"
            "Solar-System crossover placed correctly.  Unit-term "
            "variants fail the RAR (linear small-xi limit).  "
            "k = 4 sqrt(2) g_t/a0 ~ "
            f"{K_STAR:.1f} (supersedes the plan's k^2~11 estimate)."),
    "V":   ("master family V = lam u^4/4 e^{-(u/u_s)^4} + V0 e^{-(u_s/u)^4} "
            "remains the only candidate meeting V-1..V-8 jointly; "
            "exponential V passes the amplitude-sector gates and is the "
            "runner-up family to test for the ambient drift; quadratic "
            "and inverse-power fail F-adjacent/F1; symmetron runs the "
            "wrong direction."),
    "B":   ("bump/envelope family (B(0)=0, mid-field support, "
            "large-field damping) is the surviving class; pure powers "
            "are excluded; B<0 is conditionally excluded by its own "
            "null cone."),
    "open": ("global-floor question (V0 global vs effective) and the "
             "N_min mechanism (source starvation) go to T0.2/T-B1; "
             "hyperbolicity extension to the noncanonical branch is "
             "T1.4; uniqueness across the wider K(X) class remains "
             "open beyond the power-law result.")}

# ---------------------------------------------------------------------------
dest = Path(__file__).resolve().parents[2] / "results" / \
    "step_54_master_sector_admissibility.json"
dest.write_text(json.dumps(OUT, indent=2))
print(json.dumps(OUT["verdicts"], indent=2))
print("wrote", dest)
