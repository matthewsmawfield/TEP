"""TEP V(phi) BVP programme v5: shooting method for chameleon BVP.

The previous BVP (solve_bvp) converged to trivial solutions because the ODE is
stiff (mass changes by 10^19 across the body surface). This version uses a
shooting method with adaptive ODE integration (solve_ivp with Radau).

Key result from analytical analysis:
  InvPow n=2 at Lambda ~ 66 MeV:
    - Ambient lambda_C = 4.3 Gpc (cosmological coherence)
    - Earth lambda_C = 1862 km < R = 6371 km (screened)
    - Chameleon exponent p = 2/3 bridges the 19-order gap

The chameleon field transitions CONTINUOUSLY from phi_amb to phi_body over
a distance ~ lambda_C_body, which can be << R even when lambda_C_amb >> R.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import json, os, warnings
warnings.filterwarnings('ignore')

c = 2.998e8; G = 6.674e-11; hbar = 6.582e-25
hbar_c = 1.973e-14  # GeV cm
M_Pl = 2.435e18; H0 = 70e3/3.086e22
pc_cm = 3.086e18; Mpc = 1e6*pc_cm
BETA = -1.0
GeV_per_g = 5.61e23; GeV_inv_per_cm = 5.07e13
Lam_DE = np.sqrt(M_Pl * H0 * hbar)

def rho_to_GeV4(rho_g_cc): return rho_g_cc * GeV_per_g / GeV_inv_per_cm**3
def cm_to_GeVinv(cm): return cm * GeV_inv_per_cm
def GeVinv_to_cm(gi): return gi / GeV_inv_per_cm
def M_to_GeV(kg): return kg * GeV_per_g * 1e3
def compton_Mpc(m_GeV): return hbar_c / m_GeV / Mpc if m_GeV > 0 else float('inf')

RHO_COSMIC = rho_to_GeV4(9.2e-30)

# --- Potential ---

class InvPowPotential:
    """V = Lam^(4+n) / (|phi|^n + phi0^n), n even, phi < 0.
    phi0 = 1e-6 * Lambda (regularization).
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
        return -self.L4pn * (term1 + term2)  # chain rule dpsi/dphi=-1
    def equilibrium(self, rho):
        if rho <= 0: return -10 * M_Pl
        f = lambda psi: self.n * self.L4pn * psi**(self.n-1) / (psi**self.n + self.phi0**self.n)**2 - (rho/M_Pl) * np.exp(psi/M_Pl)
        psi_guess = (self.n * self.L4pn * M_Pl / rho)**(1.0/(self.n+1))
        lo = max(self.phi0 * 10, psi_guess * 0.01)
        hi = min(M_Pl, psi_guess * 100)
        try:
            psi_eq = brentq(f, lo, hi, xtol=psi_guess*1e-10)
            return -psi_eq
        except:
            return -psi_guess
    def m_eff_at(self, phi, rho):
        m2 = self.d2V(phi) + (BETA**2/M_Pl**2)*rho*np.exp(BETA*phi/M_Pl)
        return np.sqrt(np.maximum(m2, 1e-200))
    def describe(self): return f"InvPow n={self.n}, Lambda={self.Lam*1e3:.1f} MeV"

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
    if r_max_cm is None: r_max_cm = max(1000 * R_cm, 1e20)
    r = np.logspace(np.log10(R_cm * 1e-4), np.log10(r_max_cm), n)
    rho = np.where(r <= R_cm, rho_g_cc, rho_ambient_g_cc)
    return Source(r, rho_to_GeV4(rho), M_to_GeV(M_kg), "uniform")

# --- Shooting method BVP ---

def solve_shooting(potential, source, R_body_cm, phi_ambient=None, verbose=False):
    """Solve the chameleon BVP using a shooting method.

    Shoot from the center: phi(0) = phi_core, phi'(0) = 0.
    Adjust phi_core until phi(r_max) = phi_amb.
    """
    r_min_cm = R_body_cm * 1e-6  # start slightly off-center to avoid singularity
    r_max_cm = source.r[-1]
    rho_core = source.rho_at_cm(r_min_cm)
    rho_amb = source.rho_at_cm(r_max_cm)
    phi_body_eq = potential.equilibrium(rho_core)
    if phi_ambient is None:
        phi_ambient = potential.equilibrium(rho_amb)

    # Compton wavelengths for diagnostics
    m_body = potential.m_eff_at(phi_body_eq, rho_core)
    m_amb = potential.m_eff_at(phi_ambient, rho_amb)
    lam_body = hbar_c / m_body if m_body > 0 else 1e30
    lam_amb = hbar_c / m_amb if m_amb > 0 else 1e30

    if verbose:
        print(f"  phi_body_eq = {phi_body_eq:.4e} GeV ({phi_body_eq/M_Pl:.4e} M_Pl)")
        print(f"  phi_amb     = {phi_ambient:.4e} GeV ({phi_ambient/M_Pl:.4e} M_Pl)")
        print(f"  m_body = {m_body:.4e} GeV, lambda_C_body = {lam_body/1e5:.1f} km, R = {R_body_cm/1e5:.1f} km")
        print(f"  m_amb  = {m_amb:.4e} GeV, lambda_C_amb  = {lam_amb/Mpc:.1f} Mpc")
        print(f"  lambda_C_body / R = {lam_body/R_body_cm:.4f}")

    # Convert to GeV^-1
    r_min = cm_to_GeVinv(r_min_cm)
    r_max = cm_to_GeVinv(r_max_cm)
    R_body = cm_to_GeVinv(R_body_cm)

    def ode_interior(r, y):
        """ODE inside the body (r < R_body)."""
        phi, pp = y
        r_cm = GeVinv_to_cm(r)
        rho = source.rho_at_cm(r_cm)
        matter = (BETA / M_Pl) * rho * np.exp(BETA * phi / M_Pl)
        vp = potential.dV(phi)
        ppr = vp + matter - (2.0 / max(r, r_min*0.1)) * pp
        return [pp, ppr]

    def ode_exterior(r, y):
        """ODE outside the body (r > R_body)."""
        phi, pp = y
        r_cm = GeVinv_to_cm(r)
        rho = source.rho_at_cm(r_cm)
        matter = (BETA / M_Pl) * rho * np.exp(BETA * phi / M_Pl)
        vp = potential.dV(phi)
        ppr = vp + matter - (2.0 / max(r, r_min*0.1)) * pp
        return [pp, ppr]

    def shoot(phi_core):
        """Integrate from center with phi(0)=phi_core, phi'(0)=0. Return phi(r_max) - phi_amb."""
        # Start at r_min with regular center BC: phi(r_min) ≈ phi_core, phi'(r_min) ≈ 0
        y0 = [phi_core, 0.0]

        # Integrate interior (r_min to R_body)
        sol_in = solve_ivp(ode_interior, [r_min, R_body], y0,
                          method='Radau', rtol=1e-10, atol=1e-12,
                          max_step=(R_body - r_min)/100)
        if not sol_in.success:
            return 1e30

        # Integrate exterior (R_body to r_max)
        y_R = sol_in.y[:, -1]
        # Use logarithmic steps for exterior
        r_ext = np.logspace(np.log10(R_body), np.log10(r_max), 500)
        r_ext = np.clip(r_ext, R_body*1.0001, r_max*0.9999)
        sol_out = solve_ivp(ode_exterior, [R_body, r_max], y_R,
                           method='Radau', rtol=1e-10, atol=1e-12,
                           t_eval=r_ext)
        if not sol_out.success:
            return 1e30

        return sol_out.y[0, -1] - phi_ambient

    # Bracket: phi_core between phi_amb and phi_body_eq
    # For thin shell: phi_core ≈ phi_body_eq (deep interior at body equilibrium)
    # For thick shell: phi_core ≈ phi_amb (field doesn't transition)
    f_amb = shoot(phi_ambient)
    f_body = shoot(phi_body_eq)

    if verbose:
        print(f"  shoot(phi_amb) = {f_amb:.4e}")
        print(f"  shoot(phi_body) = {f_body:.4e}")

    if f_amb * f_body > 0:
        # No sign change - try a wider bracket
        if verbose:
            print(f"  WARNING: no sign change in bracket")
        # Try intermediate values
        for frac in [0.1, 0.5, 0.9, 0.99, 0.999]:
            phi_try = phi_ambient + frac * (phi_body_eq - phi_ambient)
            f_try = shoot(phi_try)
            if verbose:
                print(f"  shoot({frac:.3f} * (body-amb)) = {f_try:.4e}")
            if f_amb * f_try < 0:
                phi_core = brentq(shoot, phi_ambient, phi_try, xtol=abs(phi_body_eq-phi_ambient)*1e-8)
                break
            f_amb = f_try
        else:
            return None, None, None, dict(success=False, message="No sign change found")
    else:
        phi_core = brentq(shoot, phi_ambient, phi_body_eq,
                         xtol=abs(phi_body_eq - phi_ambient) * 1e-10)

    # Final integration with the solution
    y0 = [phi_core, 0.0]
    sol_in = solve_ivp(ode_interior, [r_min, R_body], y0,
                      method='Radau', rtol=1e-10, atol=1e-12,
                      dense_output=True)
    y_R = sol_in.y[:, -1]
    r_ext = np.logspace(np.log10(R_body), np.log10(r_max), 1000)
    r_ext = np.clip(r_ext, R_body*1.0001, r_max*0.9999)
    sol_out = solve_ivp(ode_exterior, [R_body, r_max], y_R,
                       method='Radau', rtol=1e-10, atol=1e-12,
                       t_eval=r_ext, dense_output=True)

    # Combine solutions
    r_in_cm = np.array([GeVinv_to_cm(r) for r in sol_in.t])
    r_out_cm = np.array([GeVinv_to_cm(r) for r in sol_out.t])
    r_all = np.concatenate([r_in_cm, r_out_cm])
    phi_all = np.concatenate([sol_in.y[0], sol_out.y[0]])
    pp_all = np.concatenate([sol_in.y[1], sol_out.y[1]])

    info = dict(success=True, phi_core=float(phi_core),
                m_body=float(m_body), m_amb=float(m_amb),
                lam_body_km=float(lam_body/1e5), lam_amb_Mpc=float(lam_amb/Mpc))
    return r_all, phi_all, pp_all, info

# --- Charge extraction ---

def extract_charge(r_cm, phi, phi_prime, source, R_body_cm):
    radii_ext = [2 * R_body_cm, 5 * R_body_cm, 10 * R_body_cm, 50 * R_body_cm, 100 * R_body_cm]
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

def make_body(name, r_max_factor=1000):
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

def main():
    body_names = ['Sun', 'Earth', 'Moon', 'Jupiter', 'NS', 'WD', 'WB_star']

    # Lambda values to scan around the analytical optimum (~66 MeV)
    lambdas = [30e-3, 50e-3, 66e-3, 80e-3, 100e-3, 150e-3, 200e-3]

    print("=" * 110)
    print("TEP V(phi) BVP PROGRAMME v5 — shooting method, InvPow n=2")
    print("The chameleon field transitions CONTINUOUSLY from phi_amb to phi_body")
    print("over a distance ~ lambda_C_body, which can be << R even when lambda_C_amb >> R.")
    print("=" * 110)
    print(f" rho_cosmic = {RHO_COSMIC:.3e} GeV^4, beta_A = {BETA}")
    print()

    all_results = {}

    for Lam in lambdas:
        pot = InvPowPotential(Lam, 2, phi0_factor=1e-6)
        phi_bg = pot.equilibrium(RHO_COSMIC)
        m_bg = pot.m_eff_at(phi_bg, RHO_COSMIC)
        lam_bg_Mpc = compton_Mpc(m_bg)

        print(f"\n{'='*100}")
        print(f"  Lambda = {Lam*1e3:.1f} MeV: {pot.describe()}")
        print(f"  Background: phi_amb = {phi_bg:.3e} GeV ({phi_bg/M_Pl:.2f} M_Pl), "
              f"m_amb = {m_bg:.3e} GeV, lambda_C_amb = {lam_bg_Mpc:.1f} Mpc = {lam_bg_Mpc/1e3:.2f} Gpc")
        print(f"{'='*100}")

        results = {}
        for name in body_names:
            d = BODY_DEFS[name]
            source = make_body(name)
            r_cm, phi, pp, info = solve_shooting(pot, source, d['R'], verbose=(name=='Earth'))

            if info.get('success'):
                ch = extract_charge(r_cm, phi, pp, source, d['R'])
                ch['success'] = True
                ch['phi_core'] = float(phi[0])
                ch['phi_surface'] = float(phi[np.argmin(np.abs(r_cm - d['R']))])
                ch['lam_body_km'] = info.get('lam_body_km', 0)
                ch['lam_amb_Mpc'] = info.get('lam_amb_Mpc', 0)
                results[name] = ch
                print(f"  {name:10s}  S_Sigma = {ch['S_Sigma']:+.4e}  consist = {ch['consistency']:.2e}  "
                      f"phi_core = {ch['phi_core']:.3e}  lam_C_body/R = {ch['lam_body_km']/(d['R']/1e5):.4f}")
            else:
                results[name] = dict(success=False, status=info.get('message',''), S_Sigma=None)
                print(f"  {name:10s}  FAILED: {info.get('message','')}")

        gates = compute_gates(results)
        npass = sum(g['pass_'] for g in gates.values())
        print(f"\n  Gates: {npass}/{len(gates)} pass")
        for gn, g in gates.items():
            cc = f"  (corpus: {g['corpus_claimed']:.2e})" if g['corpus_claimed'] is not None else ""
            print(f"    {'PASS' if g['pass_'] else 'FAIL'}  {gn:12s}  {g['value']:.3e}  (bound {g['bound']:.1e}){cc}")

        all_results[f"Lam_{Lam*1e3:.1f}MeV"] = dict(
            potential=pot.describe(), lambda_bg_Mpc=float(lam_bg_Mpc),
            results=results, gates=gates)

    # Summary
    print(f"\n{'='*130}")
    print("SUMMARY")
    print(f"{'='*130}")
    print(f"{'Lambda':12s} {'lam_C_bg':>10s} {'Sun S':>12s} {'Earth S':>12s} {'Moon S':>12s} {'Jup S':>12s} {'NS S':>12s} {'WD S':>12s} {'Gates':>6s}")
    for pname, res in all_results.items():
        S = lambda b: res['results'].get(b, {}).get('S_Sigma', None)
        ng = sum(g['pass_'] for g in res['gates'].values())
        nt = len(res['gates'])
        vals = [S(b) for b in ['Sun','Earth','Moon','Jupiter','NS','WD']]
        vstr = [f"{v:+.2e}" if v is not None else " FAIL " for v in vals]
        lbg = f"{res['lambda_bg_Mpc']:.0f}M" if res['lambda_bg_Mpc'] < 1e4 else f"{res['lambda_bg_Mpc']/1e3:.1f}G"
        print(f"{pname:12s} {lbg:>10s} {vstr[0]:>12s} {vstr[1]:>12s} {vstr[2]:>12s} {vstr[3]:>12s} {vstr[4]:>12s} {vstr[5]:>12s} {ng}/{nt}")

    os.makedirs(os.path.join(os.path.dirname(__file__), 'results'), exist_ok=True)
    outpath = os.path.join(os.path.dirname(__file__), 'results', 'vphi_bvp_v5_results.json')
    with open(outpath, 'w') as f:
        json.dump(all_results, f, indent=1, default=float)
    print(f"\n wrote {outpath}")

if __name__ == '__main__':
    main()
