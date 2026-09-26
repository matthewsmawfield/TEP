"""TEP V(phi) BVP programme: flexible potentials, proper body-charge extraction.

Jakarta freezes A(phi) = exp(beta_A phi/M_Pl), beta_A = -1.  Screening is an
observable PROJECTION Sigma^obs = S_Sigma(E) grad ln A, NOT a modification of A.

The scalar equation (static, radial):
  phi'' + (2/r) phi' = V_{,phi}(phi) + (beta_A / M_Pl) * rho(r) * exp(beta_A phi / M_Pl)

With beta_A = -1, the matter term is NEGATIVE.  Equilibrium requires V_{,phi} > 0.

Potential classes tested:
  Cosh:     V = Lam^4 (cosh(phi/M_Pl) - 1),  V' = (Lam^4/M_Pl) sinh(phi/M_Pl) > 0 for phi > 0
  Power:    V = Lam^4 (phi/M_Pl)^n, n>1,     V' = n Lam^4 phi^(n-1) / M_Pl^n > 0 for phi > 0
  InvPow:   V = Lam^(4+n) / |phi|^n, n even, phi < 0:  V' = n Lam^(4+n) / |phi|^(n+1) > 0
            (chameleon-type with phi < 0; gives p = (n+2)/2(n+1), up to 2/3 for n=2)

Body charge extraction (Jakarta definition):
  Q = -4 pi r^2 phi'(r)  evaluated outside the source (Gauss's law, exact for any V)
  Q_0 = |beta_A| M_body / M_Pl  (unscreened charge, linear theory)
  S_Sigma = Q / Q_0  (integrated exterior source-charge projection)

This is NOT a local S(a) evaluation.  It requires solving the BVP and extracting
the exterior field gradient.
"""
import numpy as np
from scipy.integrate import solve_bvp
import json, os

# --- Constants (SI / GeV) ---
c = 2.998e8; G = 6.674e-11; hbar = 6.582e-25
M_Pl = 2.435e18       # GeV
H0 = 70e3/3.086e22    # s^-1
Msun = 1.989e30; AU = 1.496e11; pc = 3.086e16; kpc = 1e3*pc; Mpc = 1e6*pc
BETA = -1.0
GeV_per_g = 5.61e23; GeV_inv_per_cm = 5.07e13

# Derived scale
Lam_DE = np.sqrt(M_Pl * H0 * hbar)   # GeV, dark-energy scale ~1.9 meV

def rho_to_GeV4(rho_g_cc): return rho_g_cc * GeV_per_g / GeV_inv_per_cm**3
def cm_to_GeVinv(cm): return cm * GeV_inv_per_cm
def GeVinv_to_cm(gi): return gi / GeV_inv_per_cm
def M_to_GeV(kg): return kg * GeV_per_g * 1e3  # kg -> g -> GeV

# --- Potential classes ---

class Potential:
    """Base: provides V'(phi) and V''(phi) in GeV^3 (phi in GeV)."""
    def dV(self, phi): raise NotImplementedError
    def d2V(self, phi): raise NotImplementedError
    def describe(self): raise NotImplementedError

class CoshPotential(Potential):
    """V = Lam^4 (cosh(phi/M_Pl) - 1).  Minimum at phi=0.  V' > 0 for phi > 0."""
    def __init__(self, Lam): self.Lam = Lam; self.Lam4 = Lam**4
    def dV(self, phi): return (self.Lam4 / M_Pl) * np.sinh(phi / M_Pl)
    def d2V(self, phi): return (self.Lam4 / M_Pl**2) * np.cosh(phi / M_Pl)
    def describe(self): return f"Cosh V, Lambda = {self.Lam*1e12:.2f} meV"

class PowerPotential(Potential):
    """V = Lam^4 (phi/M_Pl)^n, n > 1.  V' > 0 for phi > 0.  p = (n-2)/2(n-1) <= 1/2."""
    def __init__(self, Lam, n): self.Lam = Lam; self.n = n; self.Lam4 = Lam**4
    def dV(self, phi):
        phi = np.asarray(phi, dtype=float)
        return self.n * self.Lam4 * np.maximum(phi, 1e-30)**(self.n - 1) / M_Pl**self.n
    def d2V(self, phi):
        phi = np.asarray(phi, dtype=float)
        if self.n <= 1: return np.zeros_like(phi)
        return self.n * (self.n - 1) * self.Lam4 * np.maximum(phi, 1e-30)**(self.n - 2) / M_Pl**self.n
    def describe(self): return f"Power V=(phi/M_Pl)^{self.n}, Lambda = {self.Lam*1e12:.2f} meV"

class InvPowerPotential(Potential):
    """V = Lam^(4+n) / |phi|^n, n even, phi < 0.  Chameleon-type with phi < 0.
    V' = n Lam^(4+n) / |phi|^(n+1) > 0 for phi < 0.  p = (n+2)/2(n+1), up to 2/3 for n=2."""
    def __init__(self, Lam, n):
        self.Lam = Lam; self.n = n; self.Lam4pn = Lam**(4 + n)
    def dV(self, phi):
        phi = np.asarray(phi, dtype=float)
        apsi = np.maximum(np.abs(phi), 1e-30 * M_Pl)
        return self.n * self.Lam4pn / apsi**(self.n + 1)
    def d2V(self, phi):
        phi = np.asarray(phi, dtype=float)
        apsi = np.maximum(np.abs(phi), 1e-30 * M_Pl)
        return -self.n * (self.n + 1) * self.Lam4pn / apsi**(self.n + 2)
    def describe(self): return f"InvPower V=Lam^{4+self.n}/|phi|^{self.n} (phi<0), Lambda = {self.Lam*1e12:.2f} meV"

# --- Source profiles ---

class Source:
    """Radial density profile rho(r) in GeV^4, with body mass M in GeV."""
    def __init__(self, r_cm, rho_GeV4, M_body_GeV, name=""):
        self.r = r_cm; self.rho = rho_GeV4; self.M = M_body_GeV; self.name = name
    def rho_at_cm(self, r_cm):
        r_cm = np.asarray(r_cm)
        scalar = r_cm.ndim == 0
        res = np.interp(np.atleast_1d(r_cm), self.r, self.rho)
        res = np.where(np.atleast_1d(r_cm) <= self.r[0], self.rho[0], res)
        res = np.where(np.atleast_1d(r_cm) >= self.r[-1], self.rho[-1], res)
        return float(res[0]) if scalar else res

def uniform_sphere(rho_g_cc, R_cm, M_kg, r_max_cm=None, n=400):
    """Uniform sphere with ambient density 0."""
    if r_max_cm is None: r_max_cm = 100 * R_cm
    r = np.logspace(np.log10(R_cm * 1e-4), np.log10(r_max_cm), n)
    rho = np.where(r <= R_cm, rho_g_cc, 1e-30)  # near-zero ambient
    return Source(r, rho_to_GeV4(rho), M_to_GeV(M_kg), "uniform")

def layered_sphere(layers_g_cc, R_cm, M_kg, r_max_cm=None, n=600):
    """Layered sphere: layers = [(R_boundary_cm, rho_g_cc), ...] from center outward."""
    if r_max_cm is None: r_max_cm = 100 * R_cm
    r = np.logspace(np.log10(R_cm * 1e-4), np.log10(r_max_cm), n)
    rho = np.full_like(r, layers_g_cc[0][1])
    for i, (Rb, rv) in enumerate(layers_g_cc):
        if i == 0: continue
        mask = r > layers_g_cc[i-1][0]
        rho[mask] = rv
    rho[r > layers_g_cc[-1][0]] = 1e-30
    return Source(r, rho_to_GeV4(rho), M_to_GeV(M_kg), "layered")

# --- BVP solver ---

def solve_bvp_phi(potential, source, phi_ambient=0.0, r_max_cm=None,
                  tol=1e-6, max_nodes=80000, verbose=False):
    """Solve phi'' + (2/r) phi' = V'(phi) + (beta/M_Pl) rho(r) A(phi).

    BC: phi'(r_min) = 0 (regular centre), phi(r_max) = phi_ambient.
    Returns r_cm, phi_GeV, phi_prime_GeV_per_cm, sol.
    """
    r_min = source.r[0]
    if r_max_cm is None: r_max_cm = source.r[-1]

    # Work in GeV^-1 units for the ODE
    r_GeVinv = np.array([cm_to_GeVinv(ri) for ri in np.logspace(np.log10(r_min), np.log10(r_max_cm), 400)])
    phi_init = np.full_like(r_GeVinv, phi_ambient)

    def ode(r, y):
        phi = y[0]; pp = y[1]
        r_cm = np.array([GeVinv_to_cm(ri) for ri in r])
        rho = source.rho_at_cm(r_cm)
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

    r_dense = np.logspace(np.log10(r_GeVinv[0]), np.log10(r_GeVinv[-1]), 2000)
    y_dense = sol.sol(r_dense)
    r_cm = np.array([GeVinv_to_cm(ri) for ri in r_dense])
    return r_cm, y_dense[0], y_dense[1], sol

# --- Charge extraction ---

def extract_charge(r_cm, phi, phi_prime, source, r_body_cm):
    """Extract body charge Q = -4 pi r^2 phi'(r) at multiple exterior radii.
    Q_0 = |beta| M_body / M_Pl (unscreened, linear theory).
    S_Sigma = Q / Q_0.
    """
    # Evaluate at several radii outside the body
    r_body_GeVinv = cm_to_GeVinv(r_body_cm)
    radii_ext = [2 * r_body_cm, 5 * r_body_cm, 10 * r_body_cm, 50 * r_body_cm]
    charges = []
    for r_ext in radii_ext:
        idx = np.argmin(np.abs(r_cm - r_ext))
        if r_cm[idx] < r_body_cm: continue
        Q = -4 * np.pi * (cm_to_GeVinv(r_cm[idx]))**2 * phi_prime[idx]
        charges.append(Q)

    if not charges:
        return dict(Q=None, Q_0=None, S_Sigma=None, consistency=None)

    Q = np.mean(charges)
    Q_0 = abs(BETA) * source.M / M_Pl
    S = Q / Q_0 if Q_0 > 0 else 0.0
    consistency = (max(charges) - min(charges)) / (abs(Q) + 1e-30) if len(charges) > 1 else 0.0

    return dict(Q=float(Q), Q_0=float(Q_0), S_Sigma=float(S),
                consistency=float(consistency), charges=[float(q) for q in charges])

# --- Body definitions ---

def make_bodies():
    """Standard body set for the V(phi) programme."""
    bodies = {}
    # Sun
    bodies['Sun'] = uniform_sphere(1.41, 6.96e10, 1.989e30, r_max_cm=100*6.96e10)
    # Earth
    bodies['Earth'] = uniform_sphere(5.51, 6.371e8, 5.972e24, r_max_cm=50*6.371e8)
    # Moon
    bodies['Moon'] = uniform_sphere(3.34, 1.737e8, 7.342e22, r_max_cm=50*1.737e8)
    # Jupiter
    bodies['Jupiter'] = uniform_sphere(1.33, 6.99e9, 1.898e27, r_max_cm=50*6.99e9)
    # NS
    bodies['NS'] = uniform_sphere(2.8e14, 1.2e6, 2.784e30, r_max_cm=1e8)
    # WD
    bodies['WD'] = uniform_sphere(1e6, 7e8, 1.2e30, r_max_cm=1e11)
    # WB star (1.24 Msun)
    bodies['WB_star'] = uniform_sphere(1.41, 6.96e10, 2.467e30, r_max_cm=100*6.96e10)
    return bodies

BODY_RADIUS_CM = {'Sun': 6.96e10, 'Earth': 6.371e8, 'Moon': 1.737e8,
                  'Jupiter': 6.99e9, 'NS': 1.2e6, 'WD': 7e8, 'WB_star': 6.96e10}

# --- Gate definitions ---

def compute_gates(results):
    """Compute gates from extracted body charges."""
    gates = {}
    def gate(name, val, bound, ok, corpus=None, gtype='quant'):
        gates[name] = dict(value=float(val), bound=float(bound), pass_=bool(ok),
                           corpus_claimed=float(corpus) if corpus is not None else None, type=gtype)

    S = lambda b: results.get(b, {}).get('S_Sigma', None)

    # Cassini: S_Sigma^Sun < 5.75e-6 (Jakarta l.544)
    if S('Sun') is not None:
        gate("Cassini S_Sigma^Sun < 5.75e-6", S('Sun'), 5.75e-6, abs(S('Sun')) < 5.75e-6, gtype='null')

    # Geodesy: alpha_Earth = 2 beta^2 S_Sigma^Earth < 1e-8 (at Earth-radius range)
    if S('Earth') is not None:
        alpha_E = 2 * BETA**2 * S('Earth')
        gate("Geodesy alpha_Earth < 1e-8", alpha_E, 1e-8, abs(alpha_E) < 1e-8, gtype='null')

    # LLR Nordtvedt: eta_N = 4 beta^2 (s_E - s_M), s_A = (GM_A/(R_A c^2)) S_Sigma^A
    if S('Earth') is not None and S('Moon') is not None:
        s_E = (G * 5.972e24 / (6.371e8 * c**2)) * S('Earth')
        s_M = (G * 7.342e22 / (1.737e8 * c**2)) * S('Moon')
        eta_N = abs(4 * BETA**2 * (s_E - s_M))
        gate("LLR |eta_N| < 4.4e-4", eta_N, 4.4e-4, abs(eta_N) < 4.4e-4,
             corpus=-3.91e-4, gtype='null')

    # Pulsar: |alpha_NS| = |beta| S_Sigma^NS < 1e-3
    if S('NS') is not None:
        gate("Pulsar |alpha_NS| < 1e-3", abs(BETA) * S('NS'), 1e-3,
             abs(BETA) * S('NS') < 1e-3, gtype='null')

    # WD: |beta| S_Sigma^WD < 1e-2
    if S('WD') is not None:
        gate("WD |alpha_WD| < 1e-2", abs(BETA) * S('WD'), 1e-2,
             abs(BETA) * S('WD') < 1e-2, gtype='null')

    return gates

# --- Main ---

def run_potential(potential, bodies, phi_ambient=0.0, verbose=False):
    """Solve BVP for all bodies with a given potential. Return results dict."""
    results = {}
    for name, source in bodies.items():
        R = BODY_RADIUS_CM[name]
        r_cm, phi, pp, sol = solve_bvp_phi(potential, source, phi_ambient=phi_ambient,
                                           r_max_cm=max(source.r[-1], 100*R),
                                           verbose=verbose)
        if sol.success:
            ch = extract_charge(r_cm, phi, pp, source, R)
            results[name] = ch
            results[name]['success'] = True
            results[name]['phi_max'] = float(np.max(np.abs(phi)))
            results[name]['phi_min'] = float(np.min(phi))
            results[name]['potential'] = potential.describe()
        else:
            results[name] = dict(success=False, status=sol.message, S_Sigma=None)
            if verbose: print(f"  {name}: BVP FAILED - {sol.message}")
    return results

def main():
    bodies = make_bodies()
    all_results = {}

    # --- Test multiple potentials ---
    potentials = [
        ("Cosh_LamDE", CoshPotential(Lam_DE)),
        ("Cosh_1eV", CoshPotential(1e-3)),       # 1 meV
        ("Cosh_1eV_10", CoshPotential(1e-2)),     # 10 meV
        ("InvPow_n2_LamDE", InvPowerPotential(Lam_DE, 2)),
        ("InvPow_n2_1eV", InvPowerPotential(1e-3, 2)),
        ("InvPow_n2_1eV_10", InvPowerPotential(1e-2, 2)),
        ("InvPow_n4_LamDE", InvPowerPotential(Lam_DE, 4)),
        ("Power_n3_LamDE", PowerPotential(Lam_DE, 3)),
        ("Power_n3_1eV", PowerPotential(1e-3, 3)),
    ]

    print("=" * 80)
    print("TEP V(phi) BVP PROGRAMME")
    print("=" * 80)
    print(f" Lambda_DE = {Lam_DE*1e12:.2f} meV")
    print(f" beta_A = {BETA}")
    print(f" Bodies: {list(bodies.keys())}")
    print()

    for pname, pot in potentials:
        print(f"\n{'='*60}")
        print(f"Potential: {pname} — {pot.describe()}")
        print(f"{'='*60}")

        # For InvPower, try phi_ambient < 0 (field sits at negative values)
        if 'InvPow' in pname:
            phi_amb = -0.1 * M_Pl  # start with small negative ambient
        else:
            phi_amb = 0.0

        results = run_potential(pot, bodies, phi_ambient=phi_amb, verbose=True)

        print(f"\n  {'Body':12s} {'S_Sigma':>12s} {'Q/Q_0':>12s} {'consistency':>12s} {'phi_range':>20s}")
        for name in bodies:
            r = results.get(name, {})
            if r.get('success'):
                print(f"  {name:12s} {r['S_Sigma']:12.4e} {r['S_Sigma']:12.4e} {r['consistency']:12.2e} [{r['phi_min']:.2e},{r['phi_max']:.2e}]")
            else:
                print(f"  {name:12s}  FAILED: {r.get('status','')}")

        gates = compute_gates(results)
        npass = sum(g['pass_'] for g in gates.values())
        print(f"\n  Gates: {npass}/{len(gates)} pass")
        for gname, g in gates.items():
            cc = f"  (corpus: {g['corpus_claimed']:.2e})" if g['corpus_claimed'] is not None else ""
            print(f"    {'PASS' if g['pass_'] else 'FAIL'}  {gname:45s} {g['value']:.2e}  (bound {g['bound']:.1e}){cc}")

        all_results[pname] = dict(potential=pot.describe(), results=results, gates=gates)

    # --- Summary ---
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")
    print(f"{'Potential':25s} {'Sun S':>10s} {'Earth S':>10s} {'Moon S':>10s} {'NS S':>10s} {'Gates':>8s}")
    for pname, res in all_results.items():
        S = lambda b: res['results'].get(b, {}).get('S_Sigma', None)
        ng = sum(g['pass_'] for g in res['gates'].values())
        nt = len(res['gates'])
        vals = [S(b) for b in ['Sun', 'Earth', 'Moon', 'NS']]
        vstr = [f"{v:.2e}" if v is not None else "  FAIL " for v in vals]
        print(f"{pname:25s} {vstr[0]:>10s} {vstr[1]:>10s} {vstr[2]:>10s} {vstr[3]:>10s} {ng}/{nt:>3d}")

    os.makedirs(os.path.join(os.path.dirname(__file__), 'results'), exist_ok=True)
    with open(os.path.join(os.path.dirname(__file__), 'results', 'vphi_bvp_results.json'), 'w') as f:
        json.dump(all_results, f, indent=1, default=float)
    print(f"\n wrote results/vphi_bvp_results.json")

if __name__ == '__main__':
    main()
