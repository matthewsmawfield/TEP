#!/usr/bin/env python3
"""
TEP gravity-weight completion (2026-09-13).

Premise: G carries s^-2 (and m^3 kg^-1). If the time field rescales matter's seconds, rods and
masses by A(phi), the dimensionless strength of gravity G m^2/(hbar c) carries weight A^2 unless
the gravitational term carries the same factor. The completion puts the time field into gravity:

    S = int d^4x sqrt(-g) [ (M^2/2) F(phi) R + K X + mu^2(phi) sqrt(2X) - V(phi) ] + S_m[A^2 g]
    F(phi) = A^2(phi) Psi(phi),   A = exp(beta phi/M),  beta = -1

Psi = 1 is exact gravity-clock weight matching; psi_1 = d ln Psi / d(phi/M) is the residual.

Checks
  W1  Euler-Lagrange background equations (lapse N, scale a, phi) derived symbolically for general F.
  W2  Static frame a = 1 with Psi = 1, K = 0: V(phi), mu^2(phi) built from an exact LCDM clock history;
      all three static-frame equations are satisfied (this verifies the frame map).
  W3  Same with K > 0: V grows without bound toward the temporal horizon  -> K = 0 is required.
  W4  Fluctuation coefficients in the clock frame: no ghost (6 beta^2 + K/A^2 > 0), c_s^2 = s/(6 beta^2 + K/A^2) >= 0.
      Contrast: F = 1 (Paper 0 v0.11) gives c_s^2 = -(rho_m + 4 rho_r/3)/(K phidot^2) < 0.
  W5  Bounds on the residual weight mismatch psi_1: Gdot/G (LLR) and Cassini; matter coupling alpha.
  W6  Rate/gradient transition of the cuscuton around a mass: g_x and r_x for the allowed alpha.
Run from closure/:  python3 tep_gravity_weight_completion.py
"""
import json, math, os
import sympy as sp
from sympy.calculus.euler import euler_equations

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = {}

# ------------------------------------------------------------------ W1 symbolic background equations
t = sp.symbols('t')
M = sp.symbols('M', positive=True)
b, K, sig = sp.symbols('beta K sigma', real=True)
rm0, rr0 = sp.symbols('rho_m0 rho_r0', positive=True)
N, a, ph = [sp.Function(n)(t) for n in ('N', 'a', 'phi')]
mu2, V, Psi = sp.Function('mu2'), sp.Function('V'), sp.Function('Psi')

def lagrangian(F):
    A = sp.exp(b * ph / M)
    phd = ph.diff(t)
    X = phd ** 2 / (2 * N ** 2)
    sqrt2X = sig * phd / N                      # sigma = sign(phidot); TEP branch sigma = -1
    R = 6 * (a.diff(t, 2) / (a * N ** 2) - a.diff(t) * N.diff(t) / (a * N ** 3) + a.diff(t) ** 2 / (a ** 2 * N ** 2))
    return a ** 3 * N * (M ** 2 / 2 * F * R + K * X + mu2(ph) * sqrt2X - V(ph)) - N * rm0 * A - N * rr0 / a

Vs, Vp, m2s, m2p, Ps, Pp, Ppp = sp.symbols('V Vp mu2 mu2p Psi Psip Psipp')
p0, p1, p2, a0s, a1, a2 = sp.symbols('phi0 phidot phiddot a0 adot addot')

def clean(e, static=True):
    e = e.doit()
    e = e.subs({N.diff(t, 2): 0, N.diff(t): 0}).subs(N, 1)
    if static:
        e = e.subs({a.diff(t, 2): 0, a.diff(t): 0}).subs(a, 1)
    else:
        e = e.subs({a.diff(t, 2): a2, a.diff(t): a1}).subs(a, a0s)
    rep = {}
    for d in e.atoms(sp.Derivative):
        f = d.expr
        if getattr(f, 'func', None) == V: rep[d] = Vp
        elif getattr(f, 'func', None) == mu2: rep[d] = m2p
        elif getattr(f, 'func', None) == Psi: rep[d] = Pp if d.derivative_count == 1 else Ppp
    e = e.xreplace(rep)
    e = e.subs({ph.diff(t, 2): p2, ph.diff(t): p1})
    e = e.xreplace({V(ph): Vs, mu2(ph): m2s, Psi(ph): Ps}).subs(ph, p0)
    return sp.simplify(e)

A_s = sp.exp(b * p0 / M)
eqs_general = [clean(e.lhs - e.rhs, static=False)
               for e in euler_equations(lagrangian(sp.exp(2 * b * ph / M) * Psi(ph)), [N, a, ph], t)]
OUT["W1_FRW_equations_general_F"] = {"E_N (energy constraint)": str(eqs_general[0]),
                                      "E_a (pressure equation)": str(eqs_general[1]),
                                      "E_phi (scalar equation)": str(eqs_general[2])}
eqs_static_gw = [clean(e.lhs - e.rhs) for e in euler_equations(lagrangian(sp.exp(2 * b * ph / M) * Psi(ph)), [N, a, ph], t)]
eqs_static_F1 = [clean(e.lhs - e.rhs) for e in euler_equations(lagrangian(sp.Integer(1)), [N, a, ph], t)]
OUT["W1_static_frame_equations"] = {
    "gravity_weight_F=A^2 Psi": [str(e) for e in eqs_static_gw],
    "Paper0_v0.11_F=1": [str(e) for e in eqs_static_F1]}

# ------------------------------------------------------------------ W2/W3 exact LCDM clock history
Om, Or = 0.3153, 9.1e-5
OL = 1 - Om - Or
z = sp.symbols('z', nonnegative=True)
H = sp.sqrt(Om * (1 + z) ** 3 + Or * (1 + z) ** 4 + OL)          # clock-frame H~(z), units H0 = M = 1
Az = 1 / (1 + z)
Hdot = -(1 + z) * H * sp.diff(H, z)                                 # dH~/dt~
phidot_static = -Az * H                                             # dphi/deta, phi = M ln(1+z)
phiddot_static = -Az ** 2 * (H ** 2 + Hdot)

def design(Kval):
    """Clock-frame design: P~_X = 0 (scalar has rho+p = 0), rho~_phi = 3 Omega_L. Returns V(z), mu2(z)."""
    kin = 6 + Kval / Az ** 2                                         # beta^2 = 1
    mu2z = -kin * Az ** 3 * H
    Vz = Az ** 4 * (3 * OL - kin * H ** 2 / 2)
    return sp.simplify(Vz), sp.simplify(mu2z)

def residuals(Kval, weight=True):
    Vz, mu2z = design(Kval)
    subsd = {b: -1, M: 1, sig: -1, K: Kval, Ps: 1, Pp: 0, Ppp: 0, rm0: 3 * Om, rr0: 3 * Or}
    eqs = eqs_static_gw if weight else eqs_static_F1
    rows = []
    for zv in (0.0, 0.5, 2.0, 10.0, 1100.0, 1e5, 4e8):
        vals = {p0: math.log(1 + zv), p1: float(phidot_static.subs(z, zv)), p2: float(phiddot_static.subs(z, zv)),
                Vs: float(Vz.subs(z, zv)), m2s: float(mu2z.subs(z, zv)),
                Vp: float((sp.diff(Vz, z) * (1 + z)).subs(z, zv)), m2p: float((sp.diff(mu2z, z) * (1 + z)).subs(z, zv))}
        scale = 3 * Om / (1 + zv) + 3 * Or
        res = [abs(float(e.subs(subsd).subs(vals))) / scale for e in eqs]
        rows.append({"z": zv, "V": vals[Vs], "residual_E_N": res[0], "residual_E_a": res[1], "residual_E_phi": res[2]})
    return rows

OUT["W2_static_frame_K0_Psi1"] = residuals(0.0)
OUT["W3_static_frame_K1_Psi1"] = residuals(1.0)
OUT["W3_note"] = "K>0: V ~ -K Omega_r (1+z)^2/2 at high z (unbounded); K=0: V -> -3 Omega_r (bounded)"

# ------------------------------------------------------------------ W4 stability
Xt, s_, Kk, Aa = sp.symbols('Xtilde s K A', positive=True)
m2t = sp.symbols('mu2t', real=True)
Pt = (6 + Kk / Aa ** 2) * Xt + m2t * Aa ** -3 * sp.sqrt(2 * Xt)
PtX = sp.diff(Pt, Xt); PtXX = sp.diff(Pt, Xt, 2)
m2_for_s = sp.solve(sp.Eq(PtX, s_), m2t)[0]
OUT["W4_clock_frame_fluctuations"] = {
    "time_coefficient P_X+2X P_XX": str(sp.simplify((PtX + 2 * Xt * PtXX).subs(m2t, m2_for_s))),
    "gradient_coefficient P_X": str(sp.simplify(PtX.subs(m2t, m2_for_s))),
    "c_s^2": "s/(6 beta^2 + K/A^2)  >= 0 for s >= 0",
    "Paper0_v0.11_F=1_static": "c_s^2 = -(rho_m + 4 rho_r/3)/(K phidot^2) < 0 (see independent_foundations_checks.py C1)",
    "dark_energy_equation_of_state": "1 + w_phi = s H~^2 / rho~_phi"}

# ------------------------------------------------------------------ W5 bounds on psi_1
H0_yr = 67.4e3 / 3.0856775814913673e22 * 365.25 * 86400
LLR_1sig = 7.6e-14
psi1_LLR = LLR_1sig / H0_yr                     # Gdot/G = psi_1 H~  (phi/M = ln(1+z), d(phi/M)/dt~ = -H~)
kin_hat = 6.0
alpha = lambda p: p / (2 * math.sqrt(kin_hat))
psi1_cassini = 2 * math.sqrt(kin_hat) * math.sqrt(2.3e-5 / 2)
OUT["W5_bounds"] = {
    "Gdot_over_G": "psi_1 * H~_0",
    "alpha_matter": "-psi_1 / (2 sqrt(6 beta^2 + 1.5 psi_1^2))",
    "psi1_max_LLR": psi1_LLR, "alpha_max_LLR": alpha(psi1_LLR),
    "gamma_minus_1_at_LLR_limit": 2 * alpha(psi1_LLR) ** 2,
    "psi1_max_Cassini": psi1_cassini,
    "tightest": "LLR"}

# ------------------------------------------------------------------ W6 cuscuton rate/gradient transition
c = 2.99792458e8; H0 = 67.4e3 / 3.0856775814913673e22; GMsun = 1.32712440018e20; AU = 1.495978707e11
al = alpha(psi1_LLR)
g_x = math.sqrt(kin_hat) * c * H0 / (2 * al)
OUT["W6_transition"] = {"g_x_formula": "sqrt(6) c H0 / (2|alpha|)", "g_x_at_LLR_alpha_m_s2": g_x,
                        "r_x_Sun_AU": math.sqrt(GMsun / g_x) / AU,
                        "r_x_scaling": "r_x = sqrt(2|alpha| G M / (sqrt(6) c H0))  ~ M^(1/2)",
                        "Psi_equal_1": "no matter-sourced gradient; X stays timelike wherever |Phi|/c^2 << 1"}

# ------------------------------------------------------------------ W7 eternal temporal-horizon branch (static frame)
# K = 0, Psi = 1. For any static-frame history A(t): V from E_N, mu^2 from E_a, then
# s = 6 beta^2 + mu^2/(A^2 |phidot|)  (clock-weighted gradient coefficient; c_s^2 = s/6).
def s_static(A_, Ad, Add):
    pd = -Ad / A_                                   # beta = -1, M = 1: phi = -ln A
    pdd = -Add / A_ + Ad ** 2 / A_ ** 2
    subsd = {b: -1, M: 1, sig: (-1 if pd < 0 else 1), K: 0, Ps: 1, Pp: 0, Ppp: 0, rm0: 3 * Om, rr0: 3 * Or}
    Vval = -(3 * Om * A_ + 3 * Or)
    m2sol = float(sp.solve(eqs_static_gw[1].subs(subsd).subs({p0: -math.log(A_), p1: pd, p2: pdd, Vs: Vval}), m2s)[0])
    s_eq = 6 + m2sol / (A_ ** 2 * abs(pd))
    s_formula = 4 - 2 * A_ * Add / Ad ** 2 - (3 * Om * A_ + 4 * Or) / Ad ** 2
    return {"A": A_, "s_from_field_equations": s_eq, "s_closed_form": s_formula}
W7 = {}
for p in (0.25, 0.5):
    W7[f"power_law_p={p}"] = [dict(t=tv, **s_static((-tv) ** -p, p * (-tv) ** (-p - 1), p * (p + 1) * (-tv) ** (-p - 2)))
                              for tv in (-10.0, -1e3, -1e5)]
W7["exponential_lambda=1"] = [dict(t=tv, **s_static(math.exp(tv), math.exp(tv), math.exp(tv))) for tv in (-5.0, -10.0, -20.0)]
zz = 0.3
Hv, Hp = float(H.subs(z, zz)), float(sp.diff(H, z).subs(z, zz)); Av = 1 / (1 + zz)
Hdv = -(1 + zz) * Hv * Hp
W7["LCDM_check_z0.3"] = s_static(Av, Av ** 2 * Hv, Av ** 3 * (2 * Hv ** 2 + Hdv))
W7["closed_form"] = "s = 4 - 2 A Addot/Adot^2 - (rho_m0 A + 4 rho_r0/3)/Adot^2   (static-frame time t)"
OUT["W7_horizon_branch_stability"] = W7

# ------------------------------------------------------------------ W8 gradient-rate transition vs Paper 13
alpha_sat, Rs, Rs_stat, Rs_tot, Mwb = 0.366, 2646.0, 182.0, 609.0, 1.2
alpha_wb = math.sqrt(((1 + alpha_sat) ** 2 - 1) / 2)          # force boost 2 alpha^2 = vtilde^2 - 1
rx_wb = math.sqrt(2 * alpha_wb * GMsun * Mwb / (math.sqrt(kin_hat) * c * H0)) / AU
OUT["W8_rx_vs_Rs"] = {"alpha_from_alpha_sat": alpha_wb, "r_x_AU_M1.2": rx_wb, "R_s_AU": Rs,
                      "difference_in_total_sigma": (rx_wb - Rs) / Rs_tot,
                      "open": "which side of r_x carries the boost requires the nonlinear static solution; alpha ~ 0.66 here vs ambient bound 2.3e-4"}
print("W7", json.dumps(W7, indent=1))
print("W8", json.dumps(OUT["W8_rx_vs_Rs"], indent=1))

os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
json.dump(OUT, open(os.path.join(HERE, "results", "tep_gravity_weight_completion.json"), "w"), indent=2, default=str)
print("W1 static-frame equations, F = A^2 Psi:")
for lab, e in zip(("E_N", "E_a", "E_phi"), eqs_static_gw): print("  ", lab, "=", e)
print("W1 static-frame equations, F = 1 (v0.11):")
for lab, e in zip(("E_N", "E_a", "E_phi"), eqs_static_F1): print("  ", lab, "=", e)
for key in ("W2_static_frame_K0_Psi1", "W3_static_frame_K1_Psi1"):
    print(key)
    for r in OUT[key]:
        print("   z=%-8g V=%-12.4g res E_N=%.1e E_a=%.1e E_phi=%.1e" % (r["z"], r["V"], r["residual_E_N"], r["residual_E_a"], r["residual_E_phi"]))
print("W4", json.dumps(OUT["W4_clock_frame_fluctuations"], indent=1))
print("W5", json.dumps(OUT["W5_bounds"], indent=1))
print("W6", json.dumps(OUT["W6_transition"], indent=1))
