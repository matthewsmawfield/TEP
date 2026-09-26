"""TEP V(phi) BVP programme v2: proper chameleon setup.

Key fixes from v1:
  1. Ambient density = cosmic mean (not ~0).  The chameleon minimum exists only
     when there's ambient matter to balance V'.  With rho_ambient = 0, the field
     has no equilibrium and no screening develops.
  2. InvPower potential regularized near phi=0 (V = Lam^(4+n) / (|phi|^n + phi0^n)).
  3. Initial guess uses the algebraic equilibrium values for interior and ambient.
  4. Charge extraction checks consistency at multiple radii.

The scalar equation (static, radial):
  phi'' + (2/r) phi' = V'(phi) + (beta_A / M_Pl) * rho(r) * exp(beta_A phi / M_Pl)

With beta_A = -1, the matter term is negative.  Equilibrium: V'(phi_min) = (rho/M_Pl) exp(-phi_min/M_Pl).
Requires V' > 0 at phi_min.

For InvPow V = Lam^(4+n)/|phi|^n (n even, phi < 0):
  V' = n Lam^(4+n) / |phi|^(n+1) > 0 for phi < 0  ✓
  Equilibrium: n Lam^(4+n) / |psi|^n+1) = (rho/M_Pl) exp(psi/M_Pl),  psi = |phi|
  m_eff^2 = n(n+1) Lam^(4+n) / |psi|^(n+2) ∝ rho^((n+2)/(n+1))
  p = (n+2)/2(n+1):  n=2 -> p=2/3, n=4 -> p=3/5, n->inf -> p=1/2
"""
import numpy as np
from scipy.integrate import solve_bvp
from scipy.optimize import brentq
import json, os

c = 2.998e8; G = 6.674e-11; hbar = 6.582e-25
M_Pl = 2.435e18; H0 = 70e3/3.086e22
Msun = 1.989e30; AU = 1.496e11; pc = 3.086e16; kpc = 1e3*pc; Mpc = 1e6*pc
BETA = -1.0
GeV_per_g = 5.61e23; GeV_inv_per_cm = 5.07e13
Lam_DE = np.sqrt(M_Pl * H0 * hbar)

def rho_to_GeV4(rho_g_cc): return rho_g_cc * GeV_per_g / GeV_inv_per_cm**3
def cm_to_GeVinv(cm): return cm * GeV_inv_per_cm
def GeVinv_to_cm(gi): return gi / GeV_inv_per_cm
def M_to_GeV(kg): return kg * GeV_per_g * 1e3

RHO_COSMIC = rho_to_GeV4(9.2e-30)  # cosmic mean density in GeV^4

# --- Potentials ---

class Potential:
    def dV(self, phi): raise NotImplementedError
    def d2V(self, phi): raise NotImplementedError
    def equilibrium(self, rho_GeV4): raise NotImplementedError
    def describe(self): raise NotImplementedError

class CoshPotential(Potential):
    """V = Lam^4 (cosh(phi/M_Pl) - 1).  Min at phi=0.  V' > 0 for phi > 0."""
    def __init__(self, Lam): self.Lam = Lam; self.L4 = Lam**4
    def dV(self, phi): return (self.L4 / M_Pl) * np.sinh(phi / M_Pl)
    def d2V(self, phi): return (self.L4 / M_Pl**2) * np.cosh(phi / M_Pl)
    def equilibrium(self, rho):
        if rho <= 0: return 0.0
        # (L4/M_Pl) sinh(u) = (rho/M_Pl) exp(-u),  u = phi/M_Pl
        # L4 sinh(u) = rho exp(-u)
        # L4 (e^u - e^-u)/2 = rho e^-u
        # L4 (e^(2u) - 1) = 2 rho
        # u = 0.5 ln(1 + 2 rho / L4)
        return 0.5 * M_Pl * np.log(1 + 2 * rho / self.L4)
    def m_eff(self, phi, rho):
        return np.sqrt(self.d2V(phi) + (BETA**2 / M_Pl**2) * rho * np.exp(BETA * phi / M_Pl))
    def describe(self): return f"Cosh V, Lambda = {self.Lam*1e12:.2f} meV"

class InvPowPotential(Potential):
    """V = Lam^(4+n) / (|phi|^n + phi0^n), n even, phi < 0.
    Regularized: phi0 prevents singularity at phi=0.
    V' = n Lam^(4+n) |phi|^(n-1) / (|phi|^n + phi0^n)^2  [for phi < 0, this is > 0]
    Equilibrium psi = |phi|:  n Lam^(4+n) psi^(n-1) / (psi^n + phi0^n)^2 = (rho/M_Pl) exp(psi/M_Pl)
    """
    def __init__(self, Lam, n, phi0_factor=1e-5):
        self.Lam = Lam; self.n = n; self.L4pn = Lam**(4 + n)
        self.phi0 = phi0_factor * M_Pl
    def _w(self, phi):
        psi = np.maximum(np.abs(phi), 1e-50)
        return psi**self.n + self.phi0**self.n
    def dV(self, phi):
        psi = np.maximum(np.abs(phi), 1e-50)
        w = self._w(phi)
        return self.n * self.L4pn * psi**(self.n - 1) / w**2
    def d2V(self, phi):
        psi = np.maximum(np.abs(phi), 1e-50)
        w = self._w(phi)
        n = self.n
        term1 = n * (n-1) * psi**(n-2) / w**2
        term2 = -2 * n * psi**(n-1) * n * psi**(n-1) / w**3
        return self.L4pn * (term1 + term2)
    def equilibrium(self, rho):
        if rho <= 0: return -10 * M_Pl  # no finite minimum in vacuum
        # Find psi > 0 such that n L4pn psi^(n-1) / (psi^n + phi0^n)^2 = (rho/M_Pl) exp(psi/M_Pl)
        f = lambda psi: self.n * self.L4pn * psi**(self.n-1) / (psi**self.n + self.phi0**self.n)**2 - (rho/M_Pl) * np.exp(psi/M_Pl)
        # Search range: from very small psi to large psi
        try:
            psi_min = brentq(f, 1e-10 * M_Pl, 10 * M_Pl, xtol=1e-10 * M_Pl)
            return -psi_min  # phi < 0
        except (ValueError, RuntimeError):
            return -0.1 * M_Pl
    def m_eff(self, phi, rho):
        return np.sqrt(np.maximum(self.d2V(phi) + (BETA**2/M_Pl**2) * rho * np.exp(BETA*phi/M_Pl), 1e-50))
    def describe(self): return f"InvPow V=Lam^{4+self.n}/(|phi|^{self.n}+phi0^{self.n}), phi<0, Lambda={self.Lam*1e12:.2f} meV"

# --- Source profiles ---

class Source:
    def __init__(self, r_cm, rho_GeV4, M_body_GeV, name=""):
        self.r = r_cm; self.rho = rho_GeV4; self.M = M_body_GeV; self.name = name
    def rho_at_cm(self, r_cm):
        r_cm = np.asarray(r_cm)
        scalar = r_cm.ndim == 0
        res = np.interp(np.atleast_1d(r_cm), self.r, self.rho)
        res = np.where(np.atleast_1d(r_cm) <= self.r[0], self.rho[0], res)
        res = np.where(np.atleast_1d(r_cm) >= self.r[-1], self.rho[-1], res)
        return float(res[0]) if scalar else res

def uniform_sphere(rho_g_cc, R_cm, M_kg, r_max_cm=None, n=500, rho_ambient_g_cc=9.2e-30):
    if r_max_cm is None: r_max_cm = max(100 * R_cm, 1e20)
    r = np.logspace(np.log10(R_cm * 1e-4), np.log10(r_max_cm), n)
    rho = np.where(r <= R_cm, rho_g_cc, rho_ambient_g_cc)
    return Source(r, rho_to_GeV4(rho), M_to_GeV(M_kg), "uniform")

# --- BVP solver ---

def solve_bvp_phi(potential, source, R_body_cm, phi_ambient=None,
                  tol=1e-6, max_nodes=100000, verbose=False):
    """Solve phi'' + (2/r) phi' = V'(phi) + (beta/M_Pl) rho(r) A(phi).
    BC: phi'(r_min) = 0, phi(r_max) = phi_ambient.
    """
    r_min = source.r[0]
    r_max = source.r[-1]

    # Compute equilibrium values
    rho_core = source.rho_at_cm(r_min)
    rho_amb = source.rho_at_cm(r_max)
    phi_core_eq = potential.equilibrium(rho_core)
    if phi_ambient is None:
        phi_ambient = potential.equilibrium(rho_amb)

    if verbose:
        print(f"  phi_core_eq = {phi_core_eq:.4e} GeV, phi_ambient = {phi_ambient:.4e} GeV")
        m_core = potential.m_eff(phi_core_eq, rho_core)
        m_amb = potential.m_eff(phi_ambient, rho_amb)
        lam_core = GeVinv_to_cm(1/m_core) if m_core > 0 else 1e30
        lam_amb = GeVinv_to_cm(1/m_amb) if m_amb > 0 else 1e30
        print(f"  m_core = {m_core:.4e} GeV, lambda_core = {lam_core:.4e} cm (R_body = {R_body_cm:.4e} cm)")
        print(f"  m_amb  = {m_amb:.4e} GeV, lambda_amb  = {lam_amb:.4e} cm")
        print(f"  Screening ratio R_body/lambda_core = {R_body_cm/lam_core:.2f} (need >>1 for thin shell)")

    # Initial guess: tanh transition from core equilibrium to ambient
    r_GeVinv = np.array([cm_to_GeVinv(ri) for ri in np.logspace(np.log10(r_min), np.log10(r_max), 500)])
    R_body_GeVinv = cm_to_GeVinv(R_body_cm)
    x_trans = np.log10(r_GeVinv / R_body_GeVinv)
    phi_init = phi_core_eq + (phi_ambient - phi_core_eq) * 0.5 * (1 + np.tanh(x_trans))

    def ode(r, y):
        phi = y[0]; pp = y[1]
        r_cm_arr = np.array([GeVinv_to_cm(ri) for ri in r])
        rho = source.rho_at_cm(r_cm_arr)
        matter = (BETA / M_Pl) * rho * np.exp(BETA * phi / M_Pl)
        vp = potential.dV(phi)
        r_safe = np.maximum(r, r_GeVinv[0] * 0.1)
        ppr = vp + matter - (2.0 / r_safe) * pp
        return np.vstack([pp, ppr])

    def bc(ya, yb):
        return np.array([ya[1], yb[0] - phi_ambient])

    sol = solve_bvp(ode, bc, r_GeVinv, np.vstack([phi_init, np.zeros_like(phi_init)]),
                    tol=tol, max_nodes=max_nodes, verbose=verbose)
    if not sol.success:
        return None, None, None, sol

    r_dense = np.logspace(np.log10(r_GeVinv[0]), np.log10(r_GeVinv[-1]), 3000)
    y_dense = sol.sol(r_dense)
    r_cm = np.array([GeVinv_to_cm(ri) for ri in r_dense])
    return r_cm, y_dense[0], y_dense[1], sol

# --- Charge extraction ---

def extract_charge(r_cm, phi, phi_prime, source, R_body_cm):
    """Q = -4 pi r^2 phi'(r) at exterior radii.  Q_0 = |beta| M / M_Pl."""
    radii_ext = [2 * R_body_cm, 5 * R_body_cm, 10 * R_body_cm, 50 * R_body_cm]
    charges = []
    for r_ext in radii_ext:
        idx = np.argmin(np.abs(r_cm - r_ext))
        if r_cm[idx] < R_body_cm * 1.5: continue
        Q = -4 * np.pi * (cm_to_GeVinv(r_cm[idx]))**2 * phi_prime[idx]
        charges.append(float(Q))
    if not charges:
        return dict(Q=None, Q_0=None, S_Sigma=None, consistency=None)
    Q = np.mean(charges)
    Q_0 = abs(BETA) * source.M / M_Pl
    S = Q / Q_0 if Q_0 > 0 else 0.0
    consistency = (max(charges) - min(charges)) / (abs(Q) + 1e-30) if len(charges) > 1 else 0.0
    return dict(Q=float(Q), Q_0=float(Q_0), S_Sigma=float(S),
                consistency=float(consistency), charges=charges)

# --- Bodies ---

BODY_DEFS = {
    'Sun':     dict(rho=1.41,    R=6.96e10, M=1.989e30),
    'Earth':   dict(rho=5.51,    R=6.371e8, M=5.972e24),
    'Moon':    dict(rho=3.34,    R=1.737e8, M=7.342e22),
    'Jupiter': dict(rho=1.33,    R=6.99e9,  M=1.898e27),
    'NS':      dict(rho=2.8e14,  R=1.2e6,   M=2.784e30),
    'WD':      dict(rho=1e6,     R=7e8,     M=1.2e30),
    'WB_star': dict(rho=1.41,    R=6.96e10, M=2.467e30),
}

def make_body(name, r_max_factor=100):
    d = BODY_DEFS[name]
    return uniform_sphere(d['rho'], d['R'], d['M'], r_max_cm=max(r_max_factor*d['R'], 1e20))

# --- Gates ---

def compute_gates(results):
    gates = {}
    def gate(name, val, bound, ok, corpus=None, gtype='null'):
        gates[name] = dict(value=float(val), bound=float(bound), pass_=bool(ok),
                           corpus_claimed=float(corpus) if corpus is not None else None, type=gtype)
    S = lambda b: results.get(b, {}).get('S_Sigma', None)
    if S('Sun') is not None:
        gate("Cassini S_Sun < 5.75e-6", S('Sun'), 5.75e-6, abs(S('Sun'))<5.75e-6)
    if S('Earth') is not None:
        gate("Geodesy alpha_Earth < 1e-8", 2*BETA**2*S('Earth'), 1e-8, abs(2*BETA**2*S('Earth'))<1e-8)
    if S('Earth') is not None and S('Moon') is not None:
        sE = (G*5.972e24/(6.371e8*c**2))*S('Earth')
        sM = (G*7.342e22/(1.737e8*c**2))*S('Moon')
        gate("LLR |eta_N| < 4.4e-4", abs(4*BETA**2*(sE-sM)), 4.4e-4,
             abs(4*(sE-sM))<4.4e-4, corpus=-3.91e-4)
    if S('NS') is not None:
        gate("Pulsar |alpha_NS| < 1e-3", abs(BETA)*S('NS'), 1e-3, abs(BETA)*S('NS')<1e-3)
    if S('WD') is not None:
        gate("WD |alpha_WD| < 1e-2", abs(BETA)*S('WD'), 1e-2, abs(BETA)*S('WD')<1e-2)
    return gates

# --- Main ---

def run_potential(potential, body_names, verbose=False):
    results = {}
    for name in body_names:
        d = BODY_DEFS[name]
        source = make_body(name)
        r_cm, phi, pp, sol = solve_bvp_phi(potential, source, d['R'], verbose=verbose)
        if sol.success:
            ch = extract_charge(r_cm, phi, pp, source, d['R'])
            ch['success'] = True
            ch['phi_core'] = float(phi[0])
            ch['phi_surface'] = float(phi[np.argmin(np.abs(r_cm - d['R']))])
            ch['phi_max'] = float(np.max(phi))
            ch['phi_min'] = float(np.min(phi))
            results[name] = ch
        else:
            results[name] = dict(success=False, status=sol.message, S_Sigma=None)
    return results

def main():
    body_names = ['Sun', 'Earth', 'Moon', 'Jupiter', 'NS', 'WD', 'WB_star']
    all_results = {}

    # Lambda scan for Cosh
    cosh_lams = [Lam_DE, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0]
    # Lambda scan for InvPow n=2
    invpow_lams = [Lam_DE, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1]

    potentials = []
    for lam in cosh_lams:
        potentials.append((f"Cosh_{lam*1e12:.1f}meV", CoshPotential(lam)))
    for lam in invpow_lams:
        potentials.append((f"InvPow_n2_{lam*1e12:.1f}meV", InvPowPotential(lam, 2)))
    potentials.append(("InvPow_n4_DE", InvPowPotential(Lam_DE, 4)))

    print("=" * 80)
    print("TEP V(phi) BVP PROGRAMME v2 — proper chameleon setup")
    print("=" * 80)
    print(f" Lambda_DE = {Lam_DE*1e12:.2f} meV,  rho_cosmic = {RHO_COSMIC:.3e} GeV^4")
    print(f" beta_A = {BETA},  ambient density = 9.2e-30 g/cm^3 (cosmic mean)")
    print()

    for pname, pot in potentials:
        print(f"\n{'='*70}")
        print(f"  {pname}: {pot.describe()}")
        print(f"{'='*70}")
        results = run_potential(pot, body_names, verbose=False)

        print(f"  {'Body':10s} {'S_Sigma':>10s} {'consist':>8s} {'phi_core':>12s} {'phi_surf':>12s}")
        for name in body_names:
            r = results.get(name, {})
            if r.get('success'):
                print(f"  {name:10s} {r['S_Sigma']:10.3e} {r['consistency']:8.1e} {r['phi_core']:12.3e} {r['phi_surface']:12.3e}")
            else:
                print(f"  {name:10s}  FAILED: {r.get('status','')}")

        gates = compute_gates(results)
        npass = sum(g['pass_'] for g in gates.values())
        print(f"\n  Gates: {npass}/{len(gates)} pass")
        for gn, g in gates.items():
            cc = f"  (corpus: {g['corpus_claimed']:.2e})" if g['corpus_claimed'] is not None else ""
            print(f"    {'PASS' if g['pass_'] else 'FAIL'}  {gn:40s} {g['value']:.2e}  (bound {g['bound']:.1e}){cc}")
        all_results[pname] = dict(potential=pot.describe(), results=results, gates=gates)

    # Summary
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")
    print(f"{'Potential':25s} {'Sun':>10s} {'Earth':>10s} {'Moon':>10s} {'Jup':>10s} {'NS':>10s} {'WD':>10s} {'Gates':>6s}")
    for pname, res in all_results.items():
        S = lambda b: res['results'].get(b, {}).get('S_Sigma', None)
        ng = sum(g['pass_'] for g in res['gates'].values())
        nt = len(res['gates'])
        vals = [S(b) for b in ['Sun','Earth','Moon','Jupiter','NS','WD']]
        vstr = [f"{v:.2e}" if v is not None else " FAIL " for v in vals]
        print(f"{pname:25s} {vstr[0]:>10s} {vstr[1]:>10s} {vstr[2]:>10s} {vstr[3]:>10s} {vstr[4]:>10s} {vstr[5]:>10s} {ng}/{nt}")

    os.makedirs(os.path.join(os.path.dirname(__file__), 'results'), exist_ok=True)
    with open(os.path.join(os.path.dirname(__file__), 'results', 'vphi_bvp_v2_results.json'), 'w') as f:
        json.dump(all_results, f, indent=1, default=float)
    print(" wrote results/vphi_bvp_v2_results.json")

if __name__ == '__main__':
    main()
