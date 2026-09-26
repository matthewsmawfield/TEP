#!/usr/bin/env python3
"""
Does TEP's static background need the cuscuton at all? (2026-09-13)

Claim to check: with the gravity weight F=A^2*Psi in place, an ORDINARY field-dependent-normalization
scalar P(X,phi) = K(phi)*X - V(phi)  [NO sqrt(2X), well-defined for X<0 too, i.e. fine for local static
gradients] can support the SAME static background, with K(phi) coming out POSITIVE (no ghost) at every
epoch -- making the cuscuton's non-analytic sqrt(2X) term dispensable.

Why this can work where Paper 0 v0.11 (F=1) cannot: for P=K(phi)X-V(phi), P_X=K(phi) exactly (X-
independent, since K depends only on phi), so P_XX=0 identically and
    c_s^2 = P_X/(P_X+2X P_XX) = K(phi)/K(phi) = 1   whenever K(phi) > 0.
Stability reduces PURELY to the sign of K(phi) -- no cancellation games, no non-analytic term needed.
The only question is whether the SAME two background equations (rho_tot=0, p_tot=0) that fixed V(phi)
and mu^2(phi) before now fix V(phi) and K(phi) with K(phi) > 0 throughout, once F=A^2 is included.

Reuses the exact symbolic field equations already derived and verified in
tep_gravity_weight_completion.py (eqs_static_gw), just with mu2 -> 0 and K promoted to a function
evaluated pointwise in z (i.e. solved at each z rather than held fixed).
"""
import math, json
import sympy as sp
from sympy.calculus.euler import euler_equations

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
    sqrt2X = sig * phd / N
    R = 6 * (a.diff(t, 2) / (a * N ** 2) - a.diff(t) * N.diff(t) / (a * N ** 3) + a.diff(t) ** 2 / (a ** 2 * N ** 2))
    return a ** 3 * N * (M ** 2 / 2 * F * R + K * X + mu2(ph) * sqrt2X - V(ph)) - N * rm0 * A - N * rr0 / a

Vs, Vp, m2s, Ps, Pp, Ppp = sp.symbols('V Vp mu2 Psi Psip Psipp')
p0, p1, p2 = sp.symbols('phi0 phidot phiddot')

def clean(e):
    e = e.doit().subs({N.diff(t, 2): 0, N.diff(t): 0}).subs(N, 1)
    e = e.subs({a.diff(t, 2): 0, a.diff(t): 0}).subs(a, 1)
    rep = {}
    for d in e.atoms(sp.Derivative):
        f = d.expr
        if getattr(f, 'func', None) == V: rep[d] = Vp
    e = e.xreplace(rep).subs({ph.diff(t, 2): p2, ph.diff(t): p1})
    e = e.xreplace({V(ph): Vs, mu2(ph): m2s, Psi(ph): Ps}).subs(ph, p0)
    return sp.simplify(e)

eqs_gw = [clean(e.lhs - e.rhs) for e in euler_equations(lagrangian(sp.exp(2 * b * ph / M) * Psi(ph)), [N, a, ph], t)]
E_N, E_a = eqs_gw[0], eqs_gw[1]
print("E_N (rho_tot=0):", E_N)
print("E_a (p_tot+rho_r/3=0, times 3... check form):", E_a)

# Set mu2 -> 0 (no cuscuton). Solve E_N for V (algebraic, same as before). Substitute into E_a and
# solve for K -- now interpreted as the LOCAL value needed at this z, i.e. K(phi(z)).
V_sol = sp.solve(sp.Eq(E_N.subs(m2s, 0), 0), Vs)[0]
K_sol = sp.solve(sp.Eq(E_a.subs(m2s, 0).subs(Vs, V_sol), 0), K)
print("\nV(z) solution (mu2=0):", V_sol)
print("K(z) solution(s) (mu2=0, solving E_a):", K_sol)

# Numeric scan, Psi=1 (Psip=Psipp=0), beta=-1, M=1, LCDM-consistent background phidot, phiddot.
Om, Or = 0.3153, 9.1e-5
OL = 1 - Om - Or
Hf = lambda z: math.sqrt(Om * (1 + z) ** 3 + Or * (1 + z) ** 4 + OL)
K_expr = (K_sol[0] if len(K_sol) == 1 else K_sol)
K_num = sp.lambdify((p0, p1, p2, b, M, Ps, Pp, Ppp, rm0, rr0), K_expr, 'math')   # K_sol has no V or K left in it

rows = []
for z in (0.0, 0.5, 1.0, 3.0, 10.0, 100.0, 1100.0, 1e5, 4e8, 1e12, 1e20):
    A = 1 / (1 + z)
    phi0 = math.log(1 + z)
    phidot = -A * Hf(z)
    Hd = -(1 + z) * Hf(z) * (Hf(z + 1e-6 * max(z, 1)) - Hf(z)) / (1e-6 * max(z, 1)) if z < 1e6 else None
    # phiddot via finite difference in z is unreliable at huge z; use analytic d(phidot)/dt = A^2*(H^2+Hdot)
    h = max(z * 1e-6, 1e-9)
    Hp = (Hf(z + h) - Hf(z - h)) / (2 * h)
    Hdot_t = -(1 + z) * Hf(z) * Hp
    phiddot = A ** 2 * (Hf(z) ** 2 + Hdot_t)
    Kval = K_num(phi0, phidot, phiddot, -1.0, 1.0, 1.0, 0.0, 0.0, 3 * Om, 3 * Or)
    Vval = -(Kval * phidot ** 2 / 2 + 3 * Om * A + 3 * Or)
    rows.append({"z": z, "K": Kval, "V": Vval, "K_positive": Kval > 0})
    print(f"z={z:<10g}  K(z)={Kval:<14.6g}  V(z)={Vval:<12.5g}  {'OK (no ghost)' if Kval>0 else 'GHOST (K<0)'}")

json.dump(rows, open("results/tep_no_cuscuton.json", "w"), indent=2)
allpos = all(r["K_positive"] for r in rows)
print(f"\nK(z) > 0 at every sampled epoch: {allpos}")
