#!/usr/bin/env python3
"""
TEP local force / PPN derivation from first principles (2026-09-13).

Question at stake: for gravity weighted by F(phi) = A^p(phi) (matter coupled to g~=A^2 g as always),
what is the actual weak-field force on a test particle near a point mass, and how does it depend on p?
No formula is assumed from memory -- the Ricci scalar is built from the metric via explicit Christoffel/
Ricci contraction, and the coupled (metric, scalar) linear field equations are derived by varying the
quadratic action, exactly as tep_gravity_weight_completion.py did for the FRW background.

Ansatz: static, spherically symmetric, weak field.
    ds^2 = -(1+2*Phi(r)) dt^2 + (1-2*Lam(r)) dr^2 + r^2 dOmega^2      (c=1, small Phi, Lam)
    phi(r) = phi_c + dphi(r)                                          (phi_c = ambient value, ~0 locally)
Action:  S = int d^4x sqrt(-g) [ (M^2/2) F(phi) R  +  K X ] + S_m[A^2 g]      (K small, local test-mass limit)
Point mass at origin sources the MATTER action via g~=A^2 g, i.e. rho enters as A(phi_c)^... rho at
leading order (matter minimally coupled to g~, so its stress tensor in the g~-frame is ordinary; we read
off the coupling to phi and to (Phi,Lam) directly from expanding S_m[A^2 g] to the needed order).

Steps
  P1  Build R[g] exactly (Christoffel + Ricci contraction) for this ansatz, expand to O(Phi,Lam,dphi).
  P2  Quadratic action for (Phi, Lam, dphi) with a point-mass source; vary to get the three linear
      field equations (Phi eq, Lam eq, dphi eq) -- this is where the "F-derivative term" the user
      described lives, and its sign/coefficient come out of the algebra, not memory.
  P3  Solve the linear system for Phi(r), dphi(r) outside a point mass.
  P4  Geodesic equation for a test particle in g~ = A^2 g: extract the *physical* acceleration.
  P5  Compare F = A^2 (my original proposal) vs F = A^{-2} (the sign my hand check flagged) vs F = 1
      (Paper 0 as written), report gamma_PPN - 1 and the net coupling to Phi_total for each.
Run from closure/:  python3 tep_local_force_ppn.py
"""
import sympy as sp

r, t = sp.symbols('r t')
M, beta, Kk, p_exp = sp.symbols('M beta K p', real=True)   # p_exp = power in F = A^p
Phi, Lam, dphi = [sp.Function(s)(r) for s in ('Phi', 'Lam', 'dphi')]
eps = sp.symbols('epsilon')   # bookkeeping order parameter; drop O(eps^3)

# ---------------------------------------------------------------- P1: exact R for static spherical metric
g = sp.diag(-(1 + 2 * eps * Phi), (1 + 2 * eps * Lam), r ** 2, r ** 2 * sp.sin(sp.Symbol('theta')) ** 2)
coords = [t, r, sp.Symbol('theta'), sp.Symbol('varphi')]
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
Rs_series = sp.series(Rs, eps, 0, 2).removeO()
Rs_lin = sp.simplify(sp.expand(Rs_series).coeff(eps, 1)) * eps + Rs_series.coeff(eps, 0)
print("R[g] to linear order in (Phi,Lam):")
print("  ", sp.simplify(Rs_lin))

# ---------------------------------------------------------------- P2: quadratic action, point mass source
# sqrt(-g) to needed order:
sqrtg = sp.sqrt(-g.det())
sqrtg_lin = sp.series(sqrtg, eps, 0, 2).removeO()
Mm = sp.symbols('M_point', positive=True)
MPl2 = sp.symbols('MPl2', positive=True)          # M_Pl^2 prefactor on the gravity term
dphi_o = sp.Function('dphi')(r)
# Track ALL THREE perturbations (Phi, Lam, dphi) at the SAME formal order eps, then truncate consistently.
subs_order = {Phi: eps * Phi, Lam: eps * Lam, dphi_o: eps * dphi}
A = sp.exp(beta * (eps * dphi) / M)                # A(phi_c + dphi), phi measured from ambient (A(phi_c)=1)
F = A ** p_exp
Lgrav = sp.expand(sqrtg_lin.subs({Phi: eps * Phi, Lam: eps * Lam}) * F * Rs_lin.subs({Phi: eps * Phi, Lam: eps * Lam}))
Lgrav_quad = sp.expand(sp.series(Lgrav, eps, 0, 3).removeO())    # O(eps^1)+O(eps^2), consistent orders now
# scalar-gradient term K X = -(K/2) g^{rr} dphi'(r)^2 ~ -(K/2) dphi'(r)^2 at background g^rr=1, O(eps^2):
Lscalar = -sp.Rational(1, 2) * Kk * (eps * sp.diff(dphi, r)) ** 2

Lg2 = sp.expand(MPl2 * sp.Rational(1, 2) * Lgrav_quad + Lscalar)
# keep only the O(eps^2) [quadratic] part -> gives LINEAR field equations after variation, divide out eps^2
Lg2_quad = sp.Poly(Lg2, eps).coeff_monomial(eps ** 2) if sp.Poly(Lg2, eps).degree() >= 2 else 0
Lg2_lin_src = sp.Poly(Lg2, eps).coeff_monomial(eps ** 1) if sp.Poly(Lg2, eps).degree() >= 1 else 0
print("\nO(eps^1) piece of the gravity+scalar Lagrangian (should be a total derivative / vanish in vacuum):")
print("  ", sp.simplify(Lg2_lin_src))

# point-mass matter term (weak-field point particle coupled through g~=A^2 g, delta-function source in 3D
# represented via a boundary term: standard trick is to work with the EL vacuum equations for r>0 and fix
# the integration constant by matching total "charge" via Gauss's law, so no explicit delta needed below)

def EL(field, sym):
    dL_df = sp.diff(Lg2_quad, sym)
    dL_ddf = sp.diff(Lg2_quad, sp.diff(sym, r))
    return sp.simplify(sp.expand(r ** -2 * (dL_df - sp.diff(dL_ddf, r))))     # strip common r^2 sin(theta) weight later

eqPhi = EL(Phi, Phi)
eqLam = EL(Lam, Lam)
eqDphi = EL(dphi, dphi)

print("\nP2 Euler-Lagrange (vacuum r>0), quadratic-order fields, general F=A^p (times M^2 sin(theta)):")
print("  d/dPhi  :", eqPhi)
print("  d/dLam  :", eqLam)
print("  d/ddphi :", eqDphi)

# ---------------------------------------------------------------- P3: solve outside a point mass, for p = +2, -2, 0
G, cc = sp.symbols('G c', positive=True)
results = {}
for pv in (sp.Integer(2), sp.Integer(-2), sp.Integer(0)):
    eqP = sp.simplify(eqPhi.subs(p_exp, pv) / sp.Abs(sp.sin(sp.Symbol('theta'))))
    eqL = sp.simplify(eqLam.subs(p_exp, pv) / sp.Abs(sp.sin(sp.Symbol('theta'))))
    eqD = sp.simplify(eqDphi.subs(p_exp, pv) / sp.Abs(sp.sin(sp.Symbol('theta'))))
    results[str(pv)] = {"eqPhi": str(eqP), "eqLam": str(eqL), "eqDphi": str(eqD)}
    print(f"\n--- p={pv} ---")
    print("  eqPhi  :", eqP)
    print("  eqLam  :", eqL)
    print("  eqDphi :", eqD)
