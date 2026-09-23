#!/usr/bin/env python3
"""Minisuperspace derivation of the scalar-Gauss-Bonnet interior system.

Metric:  ds^2 = -N(T)^2 dT^2 + a(T)^2 dx^2 + b(T)^2 dOmega^2
Scalar:  u = U(T) + sigma x   (roll gradient along the spatialized t)

Action (kappa^2 = 8 pi):
    S = ∫ N a b^2 [ R/(2 kappa^2) + U'^2/(2 N^2) - sigma^2/(2 a^2)
                    - V(u) ]  +  ∫ N a b^2 lambda f(u) G

The GB piece reduces (exactly, up to a boundary term) to
    L_GB = -8 lambda f'(u) U' a' (1 + b'^2/N^2)/N
so the field equations stay second order.  Euler-Lagrange equations
w.r.t. (N, a, b, u) are assembled into

    M(y) . [a'', b'', u'']^T  =  S(y)

printed for embedding in the numerical solver.

Checks:
  * G on Schwarzschild interior (N=1/a, b=T) -> 48 M^2/T^6
  * lambda -> 0 limit reproduces the canonical-scalar KS equations
    already validated in step_33.
"""
import sympy as sp

t, th = sp.symbols('t th')
Nf, af, bf, uf, ff = (sp.Function(s) for s in ('N', 'a', 'b', 'u', 'f'))
N, a, b, u = Nf(t), af(t), bf(t), uf(t)
fp = sp.Function('fp')          # f'(u)
kap2, lam, sig = sp.symbols('kap2 lam sig')

coords = [t, sp.symbols('x'), th, sp.symbols('ph')]
x, ph = coords[1], coords[3]
g = sp.diag(-N**2, a**2, b**2, b**2*sp.sin(th)**2)
gi = g.inv(); n = 4

Gam = [[[0]*n for _ in range(n)] for _ in range(n)]
for l in range(n):
    for m in range(n):
        for nn in range(n):
            Gam[l][m][nn] = sum(
                gi[l, s]*(sp.diff(g[s, nn], coords[m])
                          + sp.diff(g[m, s], coords[nn])
                          - sp.diff(g[m, nn], coords[s]))/2
                for s in range(n))
Riem = [[[[0]*n for _ in range(n)] for _ in range(n)] for _ in range(n)]
for l in range(n):
    for m in range(n):
        for s in range(n):
            for nn in range(n):
                Riem[l][m][s][nn] = (
                    sp.diff(Gam[l][nn][m], coords[s])
                    - sp.diff(Gam[l][s][m], coords[nn])
                    + sum(Gam[l][s][k]*Gam[k][nn][m]
                          - Gam[l][nn][k]*Gam[k][s][m] for k in range(n)))
Riem_dn = [[[[sum(g[m, k]*Riem[k][nn][s][r] for k in range(n))
              for r in range(n)] for s in range(n)]
            for nn in range(n)] for m in range(n)]
Ric = [[sum(Riem[l][m][l][nn] for l in range(n))
        for nn in range(n)] for m in range(n)]
Rsc = sum(gi[m, nn]*Ric[m][nn] for m in range(n) for nn in range(n))
R2 = sum(Riem_dn[m][nn][s][r]*Riem_dn[A][B][C][D]
         * gi[m, A]*gi[nn, B]*gi[s, C]*gi[r, D]
         for m in range(n) for nn in range(n)
         for s in range(n) for r in range(n)
         for A in range(n) for B in range(n)
         for C in range(n) for D in range(n))
Ric2 = sum(Ric[m][nn]*Ric[A][B]*gi[m, A]*gi[nn, B]
           for m in range(n) for nn in range(n)
           for A in range(n) for B in range(n))
GB = sp.simplify(Rsc**2 - 4*Ric2 + R2)

M = sp.symbols('M', positive=True)
asch = sp.sqrt(2*M/t - 1)
GB_schw = sp.simplify(GB.subs({bf(t): t, af(t): asch,
                               Nf(t): 1/asch}).doit())
print("CHECK G(Schw interior) =", GB_schw, " (expect 48 M^2/t^6)")

# ---------- reduced Lagrangians ----------
ap, bp, up = sp.diff(a, t), sp.diff(b, t), sp.diff(u, t)
Vf = sp.Function('V')
L_EH = (-(2*b*ap*bp + a*bp**2)/N + N*a) / kap2
L_u = N*a*b**2*(up**2/(2*N**2) - sig**2/(2*a**2) - Vf(u))
# GB: Nab^2 f G  ==  8 f d/dt[ a'(1+b'^2/N^2)/N ]  ->  -8 f' u' Q
L_GB = -8*lam*fp(u)*up*ap*(1 + bp**2/N**2)/N

subsN = {N: 1, sp.diff(N, t): 0}


def el(L, q):
    qp = sp.diff(q, t)
    return sp.simplify(
        sp.diff(sp.diff(L, qp), t) - sp.diff(L, q)).subs(subsN)


def el_N(L):
    return sp.simplify(sp.diff(L, N)).subs(subsN)


Ltot = L_EH + L_u + L_GB
E_N = el_N(Ltot)
E_a = el(Ltot, a)
E_b = el(Ltot, b)
E_u = el(Ltot, u)

print("\n--- lambda=0 checks vs step_33 equations ---")
EN0 = sp.simplify((E_N).subs(lam, 0))
print("E_N(lam=0) =", EN0)
Ea0 = sp.simplify((E_a).subs(lam, 0))
print("E_a(lam=0) =", Ea0)
Eb0 = sp.simplify((E_b).subs(lam, 0))
print("E_b(lam=0) =", Eb0)
Eu0 = sp.simplify((E_u).subs(lam, 0))
print("E_u(lam=0) =", Eu0)

# ---------- mass matrix ----------
app, bpp, upp = sp.symbols('app bpp upp')
acc = [sp.diff(a, t, 2), sp.diff(b, t, 2), sp.diff(u, t, 2)]
eqs = [E_a, E_b, E_u]
Mmat = sp.Matrix(3, 3, lambda i, j:
                 sp.expand(eqs[i]).coeff(acc[j]))
Svec = sp.Matrix(3, 1, lambda i, j:
                 sp.simplify(-(eqs[i] - sum(
                     Mmat[i, k]*acc[k] for k in range(3)))))
print("\nM11 =", sp.simplify(Mmat[0, 0]))
print("M12 =", sp.simplify(Mmat[0, 1]))
print("M13 =", sp.simplify(Mmat[0, 2]))
print("M21 =", sp.simplify(Mmat[1, 0]))
print("M22 =", sp.simplify(Mmat[1, 1]))
print("M23 =", sp.simplify(Mmat[1, 2]))
print("M31 =", sp.simplify(Mmat[2, 0]))
print("M32 =", sp.simplify(Mmat[2, 1]))
print("M33 =", sp.simplify(Mmat[2, 2]))
print("\nS1 =", sp.simplify(Svec[0]))
print("S2 =", sp.simplify(Svec[1]))
print("S3 =", sp.simplify(Svec[2]))
print("\nCONSTRAINT E_N =", sp.simplify(E_N))
