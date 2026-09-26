"""TEP V(phi) BVP programme v4: fixed mass floor + InvPow regularization + Lambda scan.

Fixes from v3:
  1. m_eff_at floor 1e-50 -> 1e-200 (was clipping m^2 ~ 1e-84 to 1e-50, showing
     all Compton wavelengths as 0.  Actual Cosh ambient lambda_C = 2.63 Gpc.)
  2. InvPow phi0 = 1e-6*Lambda (was 1.0*Lambda, which made potential quadratic
     near origin with mass ~ Lambda, too heavy for cosmology).
  3. InvPow Lambda scan includes 10-100 MeV (where background lambda_C ~ Gpc
     and Earth has thin-shell screening).
  4. InvPow equilibrium search uses correct bracketing for |phi| >> phi0 branch.
"""
import numpy as np
from scipy.integrate import solve_bvp
from scipy.optimize import brentq
import json, os, warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

c = 2.998e8; G = 6.674e-11; hbar = 6.582e-25
hbar_c = 1.973e-14  # GeV cm
M_Pl = 2.435e18; H0 = 70e3/3.086e22
Msun = 1.989e30; pc_cm = 3.086e18; Mpc = 1e6*pc_cm
BETA = -1.0
GeV_per_g = 5.61e23; GeV_inv_per_cm = 5.07e13
Lam_DE = np.sqrt(M_Pl * H0 * hbar)

def rho_to_GeV4(rho_g_cc): return rho_g_cc * GeV_per_g / GeV_inv_per_cm**3
def cm_to_GeVinv(cm): return cm * GeV_inv_per_cm
def GeVinv_to_cm(gi): return gi / GeV_inv_per_cm
def M_to_GeV(kg): return kg * GeV_per_g * 1e3
def compton_Mpc(m_GeV): return hbar_c / m_GeV / Mpc if m_GeV > 0 else float('inf')

RHO_COSMIC = rho_to_GeV4(9.2e-30)

# --- Potentials ---

class CoshPotential:
    def __init__(self, Lam): self.Lam = Lam; self.L4 = Lam**4
    def dV(self, phi): return (self.L4 / M_Pl) * np.sinh(phi / M_Pl)
    def d2V(self, phi): return (self.L4 / M_Pl**2) * np.cosh(phi / M_Pl)
    def equilibrium(self, rho):
        if rho <= 0: return 0.0
        return 0.5 * M_Pl * np.log(1 + 2 * rho / self.L4)
    def m_eff_at(self, phi, rho):
        m2 = self.d2V(phi) + (BETA**2/M_Pl**2)*rho*np.exp(BETA*phi/M_Pl)
        return np.sqrt(np.maximum(m2, 1e-200))  # FIXED: 1e-200 not 1e-50
    def describe(self): return f"Cosh, Lambda={self.Lam*1e12:.2f} meV"

class InvPowPotential:
    """V = Lam^(4+n) / (|phi|^n + phi0^n), n even, phi < 0.
    phi0 = 1e-6 * Lambda (FIXED: was 1.0*Lambda in v3, 1e-5*M_Pl in v2).
    """
    def __init__(self, Lam, n, phi0_factor=1e-6):
        self.Lam = Lam; self.n = n; self.L4pn = Lam**(4 + n)
        self.phi0 = phi0_factor * Lam
    def _psi(self, phi): return np.maximum(np.abs(phi), 1e-200)
    def _w(self, phi):
        psi = self._psi(phi)
        return psi**self.n + self.phi0**self.n
    def dV(self, phi):
        psi = self._psi(phi); w = self._w(phi)
        return self.n * self.L4pn * psi**(self.n - 1) / w**2
    def d2V(self, phi):
        psi = self._psi(phi); w = self._w(phi); n = self.n
        term1 = n * (n-1) * psi**(n-2) / w**2
        term2 = -2 * n**2 * psi**(2*n-2) / w**3
        return -self.L4pn * (term1 + term2)  # FIXED: chain rule dpsi/dphi=-1 for phi<0
    def equilibrium(self, rho):
        if rho <= 0: return -10 * M_Pl
        f = lambda psi: self.n * self.L4pn * psi**(self.n-1) / (psi**self.n + self.phi0**self.n)**2 - (rho/M_Pl) * np.exp(psi/M_Pl)
        # For |phi| >> phi0: f'(psi) ~ n/psi^(n+1).  Equilibrium: n*Lam^(4+n)/psi^(n+1) = rho/M_Pl
        # psi_eq ~ (n*Lam^(4+n)*M_Pl/rho)^(1/(n+1))
        psi_guess = (self.n * self.L4pn * M_Pl / rho)**(1.0/(self.n+1))
        # Search around this guess
        lo = max(self.phi0 * 10, psi_guess * 0.01)
        hi = min(M_Pl, psi_guess * 100)
        try:
            psi_eq = brentq(f, lo, hi, xtol=psi_guess*1e-10)
            return -psi_eq
        except (ValueError, RuntimeError):
            try:
                psi_eq = brentq(f, self.phi0*10, M_Pl, xtol=1e-30)
                return -psi_eq
            except:
                return -psi_guess
    def m_eff_at(self, phi, rho):
        m2 = self.d2V(phi) + (BETA**2/M_Pl**2)*rho*np.exp(BETA*phi/M_Pl)
        return np.sqrt(np.maximum(m2, 1e-200))  # FIXED
    def describe(self): return f"InvPow n={self.n}, Lambda={self.Lam*1e6:.1f} keV, phi0={self.phi0:.1e}"

# --- Source ---

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
    r_min = source.r[0]; r_max = source.r[-1]
    rho_core = source.rho_at_cm(r_min)
    rho_amb = source.rho_at_cm(r_max)
    phi_core_eq = potential.equilibrium(rho_core)
    if phi_ambient is None:
        phi_ambient = potential.equilibrium(rho_amb)

    m_core = potential.m_eff_at(phi_core_eq, rho_core)
    m_amb = potential.m_eff_at(phi_ambient, rho_amb)
    lam_core_cm = hbar_c / m_core if m_core > 0 else 1e30
    lam_amb_cm = hbar_c / m_amb if m_amb > 0 else 1e30
    lam_amb_Mpc = lam_amb_cm / Mpc

    if verbose:
        print(f"  phi_core_eq = {phi_core_eq:.4e} GeV ({phi_core_eq/M_Pl:.4e} M_Pl)")
        print(f"  phi_ambient = {phi_ambient:.4e} GeV ({phi_ambient/M_Pl:.4e} M_Pl)")
        print(f"  m_core = {m_core:.4e} GeV, lambda_core = {lam_core_cm:.4e} cm = {lam_core_cm/R_body_cm:.2f} R_body")
        print(f"  m_amb  = {m_amb:.4e} GeV, lambda_amb  = {lam_amb_Mpc:.2f} Mpc")
        print(f"  R_body/lambda_core = {R_body_cm/lam_core_cm:.2e} (need >>1 for thin shell)")

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
    def gate(name, val, bound, ok, corpus=None):
        gates[name] = dict(value=float(val), bound=float(bound), pass_=bool(ok),
                           corpus_claimed=float(corpus) if corpus is not None else None)
    S = lambda b: results.get(b, {}).get('S_Sigma', None)
    if S('Sun') is not None:
        gate("Cassini", S('Sun'), 5.75e-6, abs(S('Sun'))<5.75e-6)
    if S('Earth') is not None:
        gate("Geodesy", 2*BETA**2*abs(S('Earth')), 1e-8, 2*BETA**2*abs(S('Earth'))<1e-8)
    if S('Earth') is not None and S('Moon') is not None:
        sE = (G*5.972e24/(6.371e8*c**2))*S('Earth')
        sM = (G*7.342e22/(1.737e8*c**2))*S('Moon')
        gate("LLR", abs(4*BETA**2*(sE-sM)), 4.4e-4, abs(4*(sE-sM))<4.4e-4, corpus=-3.91e-4)
    if S('NS') is not None:
        gate("Pulsar", abs(BETA)*abs(S('NS')), 1e-3, abs(BETA)*abs(S('NS'))<1e-3)
    if S('WD') is not None:
        gate("WD", abs(BETA)*abs(S('WD')), 1e-2, abs(BETA)*abs(S('WD'))<1e-2)
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
            results[name] = ch
        else:
            results[name] = dict(success=False, status=sol.message, S_Sigma=None)
    return results

def main():
    body_names = ['Sun', 'Earth', 'Moon', 'Jupiter', 'NS', 'WD', 'WB_star']
    all_results = {}

    # Cosh: scan from dark-energy scale upward
    cosh_lams = [Lam_DE, 3e-12, 10e-12, 30e-12, 100e-12, 300e-12, 1e-9]

    # InvPow n=2: scan from meV to ~100 MeV (where background lambda_C ~ Gpc)
    invpow_lams = [Lam_DE, 10e-12, 100e-12, 1e-9, 10e-9, 30e-9, 66e-9, 100e-9]

    potentials = []
    for lam in cosh_lams:
        potentials.append((f"Cosh_{lam*1e12:.1f}meV", CoshPotential(lam)))
    for lam in invpow_lams:
        potentials.append((f"InvPow_n2_{lam*1e9:.0f}neV", InvPowPotential(lam, 2, phi0_factor=1e-6)))

    print("=" * 100)
    print("TEP V(phi) BVP PROGRAMME v4 — fixed mass floor + InvPow regularization + Lambda scan")
    print("=" * 100)
    print(f" Lambda_DE = {Lam_DE*1e12:.2f} meV,  rho_cosmic = {RHO_COSMIC:.3e} GeV^4")
    print(f" beta_A = {BETA},  ambient density = 9.2e-30 g/cm^3 (cosmic mean)")
    print(f" Cosh vacuum lambda_C at Lambda_DE: {compton_Mpc(Lam_DE**2/M_Pl):.0f} Mpc = {compton_Mpc(Lam_DE**2/M_Pl)/1e3:.1f} Gpc")
    print()

    for pname, pot in potentials:
        print(f"\n{'='*90}")
        print(f"  {pname}: {pot.describe()}")
        phi_bg = pot.equilibrium(RHO_COSMIC)
        m_bg = pot.m_eff_at(phi_bg, RHO_COSMIC)
        lam_bg_Mpc = compton_Mpc(m_bg) if m_bg > 0 else float('inf')
        print(f"  Background: phi_amb = {phi_bg:.3e} GeV ({phi_bg/M_Pl:.3e} M_Pl), m_amb = {m_bg:.3e} GeV, lambda_C = {lam_bg_Mpc:.2f} Mpc")

        # Also show vacuum Compton for reference
        if isinstance(pot, CoshPotential):
            m_vac = pot.Lam**2 / M_Pl
        else:
            m_vac = pot.m_eff_at(0, 0)
        lam_vac_Mpc = compton_Mpc(m_vac) if m_vac > 0 else float('inf')
        print(f"  Vacuum: m_vac = {m_vac:.3e} GeV, lambda_C_vac = {lam_vac_Mpc:.2f} Mpc")
        print(f"{'='*90}")

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
            print(f"    {'PASS' if g['pass_'] else 'FAIL'}  {gn:12s} {g['value']:.2e}  (bound {g['bound']:.1e}){cc}")
        all_results[pname] = dict(potential=pot.describe(), lambda_bg_Mpc=float(lam_bg_Mpc),
                                   lambda_vac_Mpc=float(lam_vac_Mpc),
                                   results=results, gates=gates)

    # Summary
    print(f"\n{'='*120}")
    print("SUMMARY")
    print(f"{'='*120}")
    print(f"{'Potential':25s} {'lam_C_bg':>10s} {'lam_C_vac':>10s} {'Sun S':>10s} {'Earth S':>10s} {'Moon S':>10s} {'Jup S':>10s} {'NS S':>10s} {'WD S':>10s} {'Gates':>6s}")
    for pname, res in all_results.items():
        S = lambda b: res['results'].get(b, {}).get('S_Sigma', None)
        ng = sum(g['pass_'] for g in res['gates'].values())
        nt = len(res['gates'])
        vals = [S(b) for b in ['Sun','Earth','Moon','Jupiter','NS','WD']]
        vstr = [f"{v:.2e}" if v is not None else " FAIL " for v in vals]
        lbg = f"{res['lambda_bg_Mpc']:.1f}M" if res['lambda_bg_Mpc'] < 1e4 else f"{res['lambda_bg_Mpc']/1e3:.1f}G"
        lvac = f"{res['lambda_vac_Mpc']:.1f}M" if res['lambda_vac_Mpc'] < 1e4 else f"{res['lambda_vac_Mpc']/1e3:.1f}G"
        print(f"{pname:25s} {lbg:>10s} {lvac:>10s} {vstr[0]:>10s} {vstr[1]:>10s} {vstr[2]:>10s} {vstr[3]:>10s} {vstr[4]:>10s} {vstr[5]:>10s} {ng}/{nt}")

    os.makedirs(os.path.join(os.path.dirname(__file__), 'results'), exist_ok=True)
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'vphi_bvp_v4_results.json')
    with open(outpath, 'w') as f:
        json.dump(all_results, f, indent=1, default=float)
    print(f"\n wrote {outpath}")

if __name__ == '__main__':
    main()
