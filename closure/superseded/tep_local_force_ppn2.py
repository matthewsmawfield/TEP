#!/usr/bin/env python3
"""
TEP local force / PPN, take 2: vary the EXACT nonlinear Lagrangian, THEN linearize the resulting
differential equations. This avoids the truncated-action pitfall (varying an action already cut to
linear order in a field gives a spurious constant equation, not a field equation for it).

ds^2 = -(1+2Phi(r)) dt^2 + (1-2Lam(r)) dr^2 + r^2 dOmega^2,   phi(r) = phi_c + dphi(r), phi_c ~ 0 (ambient)
S = int d^4x sqrt(-g) [ (MPl2/2) F(phi) R[g] + K X ] + S_m[point mass, A^2 g]
F(phi) = A(phi)^p,  A = exp(beta*phi/M),  X = -(1/2) g^{rr} dphi'(r)^2 = -(1-2Lam)^{-1} dphi'(r)^2/2 (<0: static
gradient -- the well-known cuscuton-branch issue; here K is a small ORDINARY kinetic coefficient, K>0,
exactly the "K small and positive, not exactly zero" fix already flagged as needed for the local sector).
Run from closure/:  python3 tep_local_force_ppn2.py
"""
import sympy as sp

r, t, th = sp.symbols('r t theta')
M, beta, Kk, p_exp = sp.symbols('M beta K p', real=True)
Phi, Lam, dphi = [sp.Function(s)(r) for s in ('Phi', 'Lam', 'dphi')]

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
    return sp.simplify(s1 - s2 + s3 - s4)
Rmn = sp.Matrix(4, 4, lambda i, j: Ric(i, j) if i <= j else 0)
for i in range(4):
    for j in range(i):
        Rmn[i, j] = Rmn[j, i]
Rs = sp.simplify(sum(ginv[i, j] * Rmn[i, j] for i in range(4) for j in range(4)))
print("Exact R[g] for this ansatz:")
print("  ", Rs)

A = sp.exp(beta * dphi / M)
F = A ** p_exp
sqrtg = sp.sqrt(-g.det())
X = -sp.Rational(1, 2) * ginv[1, 1] * sp.diff(dphi, r) ** 2       # ginv[1,1] = g^{rr} = (1-2Lam)

Ltot = sp.expand(sqrtg * (F * Rs + Kk * X))
print("\nBuilt full nonlinear Lagrangian density (this may take a moment to simplify)...")

def EL(field):
    dL_df = sp.diff(Ltot, field)
    dL_ddf = sp.diff(Ltot, sp.diff(field, r))
    return sp.simplify(dL_df - sp.diff(dL_ddf, r))

eqPhi_full = EL(Phi)
eqLam_full = EL(Lam)
eqDphi_full = EL(dphi)
print("Done. Now linearizing about Phi=Lam=dphi=0 ...")

# Linearize: substitute field -> eps*field, series to O(eps^2), take the O(eps^1) coefficient.
eps = sp.symbols('epsilon')
sub = {Phi: eps * Phi, Lam: eps * Lam, dphi: eps * dphi}
def linearize(expr):
    e = expr.subs({sp.diff(Phi, r): eps * sp.diff(Phi, r), sp.diff(Lam, r): eps * sp.diff(Lam, r),
                   sp.diff(dphi, r): eps * sp.diff(dphi, r), sp.diff(dphi, r, 2): eps * sp.diff(dphi, r, 2),
                   sp.diff(Phi, r, 2): eps * sp.diff(Phi, r, 2)}).subs(sub)
    ser = sp.series(e, eps, 0, 2).removeO()
    return sp.simplify(sp.expand(ser).coeff(eps, 1))

eqPhi_lin = linearize(eqPhi_full)
eqLam_lin = linearize(eqLam_full)
eqDphi_lin = linearize(eqDphi_full)

print("\nLinear field equations (vacuum, r>0), general F=A^p:")
print("  eqPhi  = 0 :", eqPhi_lin)
print("  eqLam  = 0 :", eqLam_lin)
print("  eqDphi = 0 :", eqDphi_lin)

print("\n--- Specializing p = +2, -2, 0 (F=1) ---")
for pv in (sp.Integer(2), sp.Integer(-2), sp.Integer(0)):
    print(f"\np={pv}:")
    print("  eqPhi  :", sp.simplify(eqPhi_lin.subs(p_exp, pv)))
    print("  eqLam  :", sp.simplify(eqLam_lin.subs(p_exp, pv)))
    print("  eqDphi :", sp.simplify(eqDphi_lin.subs(p_exp, pv)))
