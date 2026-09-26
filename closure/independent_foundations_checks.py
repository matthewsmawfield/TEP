#!/usr/bin/env python3
"""
Independent foundations checks of the TEP v5.x architecture (2026-09-13).

Every check uses only ingredients the corpus itself states (Paper 0 v0.11 section 2.2 and 8;
CORPUS_PLAN.md section 7; superseded/tep_architecture.py):

  S = int sqrt(-g) [ M_Pl^2 R/2 + P(X,phi) ] + S_m[g~],   P = K X + mu^2 sqrt(2X) - V
  g~ = A^2 g,   A = exp(beta phi/M_Pl),   beta = -1
  static flat Einstein frame,  1+z = A_0/A_em,  matter frame a~ = A,  H~(z) ~ H_LCDM(z)

  C1  static background with the full P(X): K X cannot be dropped -> P_X < 0 (gradient instability)
  C2  matter-clock age and Einstein-frame time back to A -> 0 are both finite
  C3  G m^2/(hbar c) carries conformal weight A^2 -> Gdot/G = 2 H~_0
  C4  linear growth from recombination with gravity weight A^2
  C5  redshift-acceleration identity: clock-rate structure implies a gravity anomaly
  C6  cuscuton local branch: sqrt(2X) undefined for static gradients; floor acceleration
  C7  Paper 17 eta_resid mapped to Earth-Moon differential acceleration
  C8  GNSS shared-satellite fraction vs baseline (a mundane common-mode model; honest negative)
  C9  Paper 27 geodesically complete branch A = C eta^-p: implied H~(z) and distances

Run from the TEP/closure directory:  python3 independent_foundations_checks.py
Units: reduced M_Pl = 1 and H0 = 1 for cosmology unless stated.
"""
import json, math, os
import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp
from scipy.optimize import curve_fit

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = {}
c = 2.99792458e8
Mpc = 3.0856775814913673e22
AU = 1.495978707e11
GMsun = 1.32712440018e20
yr = 365.25 * 86400.0
Gyr = 1e9 * yr
H0 = 67.4e3 / Mpc
Om, Or = 0.3153, 9.1e-5
OL = 1 - Om - Or
E = lambda z: math.sqrt(Om * (1 + z) ** 3 + Or * (1 + z) ** 4 + OL)

# ---------------------------------------------------------------- C1
phid = sp.symbols('phidot', positive=True)
K, rm, rr = sp.symbols('K rho_m rho_r', positive=True)
m2, V = sp.symbols('mu2 V', real=True)
X = sp.symbols('X', positive=True)
P = K * X + m2 * sp.sqrt(2 * X) - V
PX, PXX = sp.diff(P, X), sp.diff(P, X, 2)
rho_phi = sp.simplify(2 * X * PX - P)
bg = {X: phid ** 2 / 2}
sol_full = sp.solve([sp.Eq((rho_phi + rm + rr).subs(bg), 0),
                     sp.Eq((P + rr / 3).subs(bg), 0)], [V, m2], dict=True)[0]
sol_corp = sp.solve([sp.Eq(V + rm + rr, 0),
                     sp.Eq((m2 * sp.sqrt(2 * X) - V + rr / 3).subs(bg), 0)], [V, m2], dict=True)[0]
eps, a_, b_ = sp.symbols('epsilon dphidot gradphi', real=True)
Xp = (phid + eps * a_) ** 2 / 2 - (eps * b_) ** 2 / 2
c2 = sp.diff(K * Xp + m2 * sp.sqrt(2 * Xp), eps, 2).subs(eps, 0) / 2
T_coef = sp.simplify(sp.diff(c2, a_, 2).subs(sol_full))    # L2 = T/2 dphidot^2 - G/2 (grad dphi)^2
G_coef = sp.simplify(-sp.diff(c2, b_, 2).subs(sol_full))
OUT["C1"] = {"rho_phi": str(rho_phi), "V_full": str(sp.simplify(sol_full[V])),
             "mu2_full": str(sp.simplify(sol_full[m2])),
             "P_X_full": str(sp.simplify(PX.subs(bg).subs(sol_full))),
             "P_X_corpus_KX_dropped": str(sp.simplify(PX.subs(bg).subs(sol_corp))),
             "fluct_time_coeff": str(T_coef), "fluct_gradient_coeff": str(G_coef)}
rows = []
for z in (0.0, 1.0, 10.0, 1100.0, 1e5, 4e8):
    A = 1 / (1 + z); pd = A * E(z); rho_m, rho_r = 3 * Om * A, 3 * Or; KX = 0.5 * pd ** 2
    rows.append({"z": z, "KX_over_rho_m_plus_rho_r": KX / (rho_m + rho_r),
                 "V_full": -(rho_m + rho_r) - KX,
                 "cs2_full_K1": -(rho_m + 4 / 3 * rho_r) / pd ** 2,
                 "cs2_corpus_K1": 1 - (rho_m + 4 / 3 * rho_r) / pd ** 2})
OUT["C1"]["rows_K1"] = rows
OUT["C1"]["efold_time_s"] = {lab: 1 / (math.sqrt(abs(rows[0]["cs2_full_K1"])) * c * 2 * math.pi / L)
                             for lab, L in (("1_AU_mode", AU), ("1_km_mode", 1e3))}

# ---------------------------------------------------------------- C2
OUT["C2"] = {
    "matter_clock_age_Gyr": quad(lambda z: 1 / ((1 + z) * E(z)), 0, np.inf, limit=400)[0] / H0 / Gyr,
    "einstein_frame_time_Gyr": quad(lambda z: 1 / E(z), 0, np.inf, limit=400)[0] / H0 / Gyr,
    "coasting_constant_phidot_age_Gyr": 1 / H0 / Gyr}

# ---------------------------------------------------------------- C3
OUT["C3"] = {"Gdot_over_G_per_yr": 2 * H0 * yr,
             "ratio_to_LLR_1sigma_7.6e-14": 2 * H0 * yr / 7.6e-14,
             "G_rec_over_G0": (1 / 1101) ** 2, "G_BBN_over_G0": (1 / 4e8) ** 2}

# ---------------------------------------------------------------- C4
def growth(model, gs=1.0):
    def rhs(N, y):
        a = math.exp(N); e2 = Om / a ** 3 + Or / a ** 4 + OL
        dlnE = (-3 * Om / a ** 3 - 4 * Or / a ** 4) / (2 * e2)
        src = 1.5 * Om / (a ** 3 * e2) * (gs * a ** 2 if model == "tep" else 1.0)
        return [y[1], -(2 + dlnE) * y[1] + src * y[0]]
    a0 = 1 / 1101
    s = solve_ivp(rhs, [math.log(a0), 0.0], [a0, a0], rtol=1e-9, atol=1e-14)
    return float(s.y[0, -1] / a0)
OUT["C4"] = {"LCDM": growth("lcdm"), "TEP_screened": growth("tep", 1.0),
             "TEP_unscreened_1plus2beta2": growth("tep", 3.0)}

# ---------------------------------------------------------------- C5
lamT = 4.2e6
OUT["C5"] = {f"dlnA={d:g}": c ** 2 * d / lamT for d in (1e-13, 1e-14, 1e-15, 1e-18)}
OUT["C5"]["max_time_variable_dlnA_gravity_noise_1e-10"] = 1e-10 * lamT / c ** 2

# ---------------------------------------------------------------- C6
cH0 = c * H0
OUT["C6"] = {"X_sign_flip_accel_m_s2": cH0 / 2,
             "X_sign_flip_radius_solar_mass_AU": math.sqrt(GMsun / (cH0 / 2)) / AU,
             "floor_accel_m_s2": 3 * Om * cH0,
             "floor_accel_corpus_factor_half_m_s2": 3 * Om * cH0 / 2,
             "corpus_T2_saturn_bound_m_s2": 1e-14,
             "floor_over_bound": 3 * Om * cH0 / 2 / 1e-14,
             "kappa_times_X": 3 * Om / 2}

# ---------------------------------------------------------------- C7
eta, s_eta = -3.91e-4, 5.63e-5
OUT["C7"] = {"Delta_mg_mi": eta * (-4.64e-10 + 0.19e-10),
             "cosD_amplitude_mm": 13.1e3 * abs(eta),
             "eta_tension_with_HM2018_sigma": (eta + 0.2e-4) / math.hypot(s_eta, 1.1e-4)}

# ---------------------------------------------------------------- C8
R = 6371.0
def psi(h, mask=10.0):
    e = math.radians(mask); return math.acos(R * math.cos(e) / (R + h)) - e
rng = np.random.default_rng(1)
v = rng.normal(size=(600000, 3)); v /= np.linalg.norm(v, axis=1)[:, None]
def shared(d_km, ps):
    th = d_km / R
    va = v[:, 2] > math.cos(ps)
    vb = (v[:, 0] * math.sin(th) + v[:, 2] * math.cos(th)) > math.cos(ps)
    return (va & vb).sum() / va.sum()
ds = np.linspace(0, 13000, 53)
OUT["C8"] = {}
for name, h in (("GPS", 20200), ("Galileo", 23222), ("GLONASS", 19130)):
    ps = psi(h); J = np.array([shared(d, ps) for d in ds]); sel = ds <= 8000
    lam1 = curve_fit(lambda d, lam: np.exp(-d / lam), ds[sel], J[sel], p0=[4000])[0][0]
    OUT["C8"][name] = {"lambda_pure_exp_0_8000km": float(lam1), "shared_fraction_at_4200km": float(shared(4200, ps))}

# ---------------------------------------------------------------- C9
chi_l = lambda z: quad(lambda x: 1 / E(x), 0, z, limit=400)[0]
OUT["C9"] = {}
for p in (0.25, 0.5):
    chi_p = lambda z: p * ((1 + z) ** (1 / p) - 1)
    OUT["C9"][f"p={p}"] = {"H_over_H0_z1": 2 ** ((p - 1) / p), "H_over_H0_z1100": 1101 ** ((p - 1) / p),
                           "delta_mu_mag_z1": 5 * math.log10(chi_p(1) / chi_l(1)),
                           "acoustic_angle_ratio_z1090": chi_l(1090) / chi_p(1090)}
OUT["C9"]["LCDM"] = {"H_over_H0_z1": E(1), "H_over_H0_z1100": E(1100)}

os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
json.dump(OUT, open(os.path.join(HERE, "results", "independent_foundations_checks.json"), "w"), indent=2)
print(json.dumps({k: v for k, v in OUT.items() if k != "C1"}, indent=1))
print("C1 P_X full   :", OUT["C1"]["P_X_full"])
print("C1 P_X corpus :", OUT["C1"]["P_X_corpus_KX_dropped"])
print("C1 fluct coeffs: time", OUT["C1"]["fluct_time_coeff"], " gradient", OUT["C1"]["fluct_gradient_coeff"])
for r in rows:
    print("   z=%-8g KX/(rho_m+rho_r)=%-9.3g cs2 full=%-10.3g cs2 corpus=%-7.3g V=%.3g"
          % (r["z"], r["KX_over_rho_m_plus_rho_r"], r["cs2_full_K1"], r["cs2_corpus_K1"], r["V_full"]))
