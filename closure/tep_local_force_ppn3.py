#!/usr/bin/env python3
"""
TEP local force / PPN, take 3: bake the bookkeeping parameter eps into the ansatz from the start
(so no post-hoc substitution into Derivative objects is needed), vary the exact nonlinear Lagrangian,
then extract the O(eps) linear equation via a Gateaux derivative d/deps|_{eps=0} -- robust, avoids
sympy's series()/degree() machinery choking on a complex expression (what broke take 2).
"""
import sympy as sp
from sympy.calculus.euler import euler_equations

r, t, th = sp.symbols('r t theta')
M, beta, Kk, p_exp, eps = sp.symbols('M beta K p epsilon', real=True)
Phi1, Lam1, dphi1 = [sp.Function(s)(r) for s in ('Phi1', 'Lam1', 'dphi1')]
Phi, Lam, dphi = eps * Phi1, eps * Lam1, eps * dphi1     # eps baked in from the start

g = sp.diag(-(1 + 2 * Phi), 1 / (1 - 2 * Lam), r ** 2, r ** 2 * sp.sin(th) ** 2)
coords = [t, r, th, sp.Symbol('varphi')]
ginv = g.inv()
n = 4
Gamma = [[[sp.together(sp.Rational(1, 2) * sum(ginv[a, d] * (sp.diff(g[d, b], coords[c]) + sp.diff(g[d, c], coords[b]) - sp.diff(g[b, c], coords[d])) for d in range(n)))
           for c in range(n)] for b in range(n)] for a in range(n)]
def Ric(b, c):
    s1 = sum(sp.diff(Gamma[a][b][c], coords[a]) for a in range(n))
    s2 = sum(sp.diff(Gamma[a][b][a], coords[c]) for a in range(n))
    s3 = sum(Gamma[a][a][d] * Gamma[d][b][c] for a in range(n) for d in range(n))
    s4 = sum(Gamma[a][c][d] * Gamma[d][b][a] for a in range(n) for d in range(n))
    return s1 - s2 + s3 - s4
Rmn = sp.Matrix(4, 4, lambda i, j: Ric(i, j) if i <= j else 0)
for i in range(4):
    for j in range(i):
        Rmn[i, j] = Rmn[j, i]
Rs = sum(ginv[i, j] * Rmn[i, j] for i in range(4) for j in range(4))
print("Built exact R[g] (unsimplified, eps-exact). Building full Lagrangian...")

A = sp.exp(beta * dphi / M)
F = A ** p_exp
sqrtg = sp.sqrt(-g.det())
X = -sp.Rational(1, 2) * ginv[1, 1] * sp.diff(dphi, r) ** 2

Ltot = sqrtg * (F * Rs + Kk * X)
print("Lagrangian built. Computing Euler-Lagrange via sympy's euler_equations "
      "(handles the Phi''/Lam'' terms in R[g] correctly; this is the expensive step)...")

eqs = euler_equations(Ltot, [Phi1, Lam1, dphi1], r)
eqPhi_full = sp.Eq(eqs[0].lhs - eqs[0].rhs, 0).lhs
eqLam_full = sp.Eq(eqs[1].lhs - eqs[1].rhs, 0).lhs
eqDphi_full = sp.Eq(eqs[2].lhs - eqs[2].rhs, 0).lhs
print("  Euler-Lagrange done. Now linearizing via d/deps|eps=0 ...")

def linear_part(expr):
    # eqPhi_full = eps * [functional derivative of L w.r.t. Phi, evaluated at Phi=eps*Phi1],
    # and that functional derivative itself vanishes at Phi=0 (flat space trivially solves vacuum
    # eqns), so eqPhi_full = O(eps^2), not O(eps^1). Extract the eps^2 Taylor coefficient instead.
    d2 = sp.diff(expr, eps, 2).subs(eps, 0)
    return sp.simplify(d2)

eqPhi_lin = linear_part(eqPhi_full)
eqLam_lin = linear_part(eqLam_full)
eqDphi_lin = linear_part(eqDphi_full)

print("\nLinear field equations (vacuum, r>0), general F=A^p (overall factor of sin(theta) stripped below):")
sinth = sp.sin(th)
eqPhi_lin = sp.simplify(eqPhi_lin / sinth)
eqLam_lin = sp.simplify(eqLam_lin / sinth)
eqDphi_lin = sp.simplify(eqDphi_lin / sinth)
print("  eqPhi  = 0 :", eqPhi_lin)
print("  eqLam  = 0 :", eqLam_lin)
print("  eqDphi = 0 :", eqDphi_lin)

print("\n--- Specializing p = +2, -2, 0 (F=1) ---")
for pv in (sp.Integer(2), sp.Integer(-2), sp.Integer(0)):
    print(f"\np={pv}:")
    print("  eqPhi  :", sp.simplify(eqPhi_lin.subs(p_exp, pv)))
    print("  eqLam  :", sp.simplify(eqLam_lin.subs(p_exp, pv)))
    print("  eqDphi :", sp.simplify(eqDphi_lin.subs(p_exp, pv)))
