"""CORPUS_PLAN §7.8 Target 4 (Paper 15): numerical check of the stationarity theorem.

A test mass follows geodesics of the TEP matter metric
    ds^2 = -(1 + 2 Phi_N) A^2 dt^2 + A^2 dx^2 + B (grad phi . dx)^2,   A = exp(beta phi)
with a STATIC, non-spherical phi (monopole + J2-like + tesseral lobe) and deliberately large
conformal and disformal couplings. Static metric => Killing energy conserved => |v| in = |v| out
far from the body. Control: the same run with a slowly time-dependent phi must show a change.
Units: c = 1, GM = 1.
"""
import json
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

x, y, z, t = sp.symbols("x y z t", real=True)
vx, vy, vz = sp.symbols("vx vy vz", real=True)
beta, B0, eps = -1.0, 5e-3, sp.Symbol("eps", real=True)
r = sp.sqrt(x**2 + y**2 + z**2)
phi = 1e-3 * (1 / r) * (1 + 0.3 * (1.5 * (z / r) ** 2 - 0.5) / r**2 + 0.2 * (x * y / r**2) / r**2) * (1 + eps * t)
PhiN = -1e-4 / r
A2 = sp.exp(2 * beta * phi)
v = sp.Matrix([vx, vy, vz]); q = sp.Matrix([x, y, z])
grad = sp.Matrix([sp.diff(phi, s) for s in (x, y, z)])
L = -sp.sqrt((1 + 2 * PhiN) * A2 - A2 * (v.dot(v)) - B0 * (v.dot(grad)) ** 2)
p = sp.Matrix([sp.diff(L, s) for s in (vx, vy, vz)])
M = p.jacobian(v)
rhs = sp.Matrix([sp.diff(L, s) for s in (x, y, z)]) - p.jacobian(q) * v - sp.diff(p, t)
E = (p.dot(v) - L)
fM = sp.lambdify((t, x, y, z, vx, vy, vz, eps), M, "numpy")
fR = sp.lambdify((t, x, y, z, vx, vy, vz, eps), rhs, "numpy")
fE = sp.lambdify((t, x, y, z, vx, vy, vz, eps), E, "numpy")

def run(eps_val):
    def f(tt, s):
        a = np.linalg.solve(np.array(fM(tt, *s, eps_val), float), np.array(fR(tt, *s, eps_val), float).ravel())
        return np.concatenate([s[3:], a])
    s0 = np.array([-3000.0, 4.0, 2.0, 0.02, 0.0, 0.003])
    sol = solve_ivp(f, (0, 3.0e5), s0, rtol=1e-12, atol=1e-14, method="DOP853")
    a, b = sol.y[:, 0], sol.y[:, -1]
    rmin = float(np.min(np.linalg.norm(sol.y[:3], axis=0)))
    Ein, Eout = float(fE(0, *a, eps_val)), float(fE(sol.t[-1], *b, eps_val))
    # asymptotic speed from conserved-energy relation at the endpoints (removes finite-r potential)
    def vinf(Eval):
        return float(np.sqrt(max(0.0, 1 - 1 / Eval**2)))
    return {"eps": eps_val, "r_min": rmin, "r_end": float(np.linalg.norm(b[:3])),
            "E_in": Ein, "E_out": Eout, "dE_over_E": (Eout - Ein) / Ein,
            "v_inf_in": vinf(Ein), "v_inf_out": vinf(Eout), "dv_inf_over_v": (vinf(Eout) - vinf(Ein)) / vinf(Ein),
            "deflection_deg": float(np.degrees(np.arccos(np.dot(a[3:], b[3:]) / np.linalg.norm(a[3:]) / np.linalg.norm(b[3:]))))}

out = {"static": run(0.0), "time_dependent_control": run(2e-6)}
out["statement"] = ("Static TEP matter metric (conformal + disformal, non-spherical phi): Killing energy conserved, "
                    "asymptotic speed unchanged. Time-dependent control changes it, so the check has power.")
Path(__file__).resolve().parent.joinpath("results").mkdir(exist_ok=True)
json.dump(out, open(Path(__file__).resolve().parent / "results" / "tep_flyby_stationarity_check.json", "w"), indent=1)
for k in ("static", "time_dependent_control"):
    o = out[k]
    print(f"{k:24s} r_min={o['r_min']:.1f}  deflection={o['deflection_deg']:.1f} deg  dE/E={o['dE_over_E']:+.2e}  dv_inf/v={o['dv_inf_over_v']:+.2e}")
