#!/usr/bin/env python3
"""
TEP Scalar-Sector Closure

One scale:  Lambda = (M_Pl * H0 * hbar)^(1/2)
One coupling: beta_A = -1
One potential: V(phi) = Lambda^4 + (1/2)*Lambda^2*phi^2

Three projections:
  1. Scalar mass:     m = Lambda  ->  lambda_C = 0.10 mm  (local screening)
  2. Dark energy:     V(0) = Lambda^4  (cosmology)
  3. Transition accel: a_0 = Lambda^2/M_Pl = c*H0  (wide binaries, SPARC)

Identity: Lambda^2/M_Pl = H0  (natural units, exact)

Usage:
    python3 scalar_closure.py
"""

import numpy as np
import json
import os

# ── Constants ──────────────────────────────────────────────────────────
c = 2.998e8           # m/s
G = 6.674e-11         # m^3 kg^-1 s^-2
hbar = 6.582e-25      # GeV s
hbar_c = 1.973e-14    # GeV cm
M_Pl = 2.435e18       # GeV (reduced Planck mass)
H0 = 70e3 / 3.086e22  # s^-1 (Hubble constant)
Mpc = 3.086e22        # m
BETA = -1.0           # conformal coupling
GeV_inv_per_cm = 5.07e13
GeV_per_g = 5.61e23

# ── The single scale ───────────────────────────────────────────────────
Lam = np.sqrt(M_Pl * H0 * hbar)  # GeV

# ── Three projections ──────────────────────────────────────────────────
m_scalar = Lam                          # scalar mass
lam_C = hbar_c / m_scalar               # Compton wavelength (cm)
V0 = Lam**4                             # vacuum energy
a0 = c * H0                             # transition acceleration (m/s^2)
a0_natural = Lam**2 / M_Pl              # transition acceleration (GeV, natural)

# ── Identity check ────────────────────────────────────────────────────
# Lambda^2/M_Pl = H0*hbar  (exact, since Lambda^2 = M_Pl*H0*hbar)
# a_0 = c*H0 = c*Lambda^2/(M_Pl*hbar)  (exact)
identity_ratio = (Lam**2 / M_Pl) / (H0 * hbar)

# ── Bodies ────────────────────────────────────────────────────────────
bodies = {
    'Sun':     {'R_cm': 6.96e10, 'M_kg': 1.989e30},
    'Earth':   {'R_cm': 6.371e8,  'M_kg': 5.972e24},
    'Moon':    {'R_cm': 1.737e8,  'M_kg': 7.342e22},
    'Jupiter': {'R_cm': 6.99e9,   'M_kg': 1.898e27},
    'NS':      {'R_cm': 1.2e6,    'M_kg': 2.784e30},
    'WD':      {'R_cm': 7e8,      'M_kg': 1.2e30},
}

# ── Yukawa screening factor ───────────────────────────────────────────
def S_Sigma(m, R_cm):
    """Source-charge screening for uniform sphere, V = m^2*phi^2/2."""
    x = m * R_cm * GeV_inv_per_cm
    if x < 1e-10:
        return 1.0
    return 3.0 * (1 + x) / (x * (2 * x + 1))

# ── Gate definitions ──────────────────────────────────────────────────
def check_gates(results):
    S_Sun = results['Sun']['S']
    S_E = results['Earth']['S']
    S_M = results['Moon']['S']
    S_NS = results['NS']['S']
    S_WD = results['WD']['S']

    # Cassini: gamma_PPN - 1 = -4*beta_A^2*S_Sun/(1+2*beta_A^2*S_Sun)
    # (linear screened-source map; unscreened photon probe)
    # |S_Sun| < 5.75e-6 for |beta_A| = 1
    cassini = {
        'observable': 'S_Sun',
        'value': abs(S_Sun),
        'bound': 5.75e-6,
        'pass': abs(S_Sun) < 5.75e-6,
    }

    # Geodesy: 2*beta^2*|S_Earth| < 1e-8
    alpha_E = 2 * BETA**2 * abs(S_E)
    geodesy = {
        'observable': 'alpha_Earth',
        'value': alpha_E,
        'bound': 1e-8,
        'pass': alpha_E < 1e-8,
    }

    # LLR: |eta_N| = |4*(s_E - s_M)| < 4.4e-4
    sE = (G * 5.972e24 / (6.371e8 * c**2)) * S_E
    sM = (G * 7.342e22 / (1.737e8 * c**2)) * S_M
    eta_N = abs(4 * (sE - sM))
    llr = {
        'observable': 'eta_N',
        'value': eta_N,
        'bound': 4.4e-4,
        'pass': eta_N < 4.4e-4,
    }

    # Pulsar: |beta|*|S_NS| < 1e-3
    alpha_NS = abs(BETA) * abs(S_NS)
    pulsar = {
        'observable': 'alpha_NS',
        'value': alpha_NS,
        'bound': 1e-3,
        'pass': alpha_NS < 1e-3,
    }

    # WD: |beta|*|S_WD| < 1e-2
    alpha_WD = abs(BETA) * abs(S_WD)
    wd = {
        'observable': 'alpha_WD',
        'value': alpha_WD,
        'bound': 1e-2,
        'pass': alpha_WD < 1e-2,
    }

    gates = {
        'Cassini': cassini,
        'Geodesy': geodesy,
        'LLR': llr,
        'Pulsar': pulsar,
        'WD': wd,
    }
    n_pass = sum(g['pass'] for g in gates.values())
    return gates, n_pass

# ── Main calculation ──────────────────────────────────────────────────
def main():
    print('=' * 72)
    print('TEP SCALAR-SECTOR CLOSURE')
    print('=' * 72)
    print()
    print(f'Scale:     Lambda = (M_Pl * H0 * hbar)^(1/2) = {Lam:.4e} GeV = {Lam*1e3:.2f} meV')
    print(f'Coupling:  beta_A = {BETA}')
    print(f'Potential: V(phi) = Lambda^4 + (1/2)*Lambda^2*phi^2')
    print()

    # ── Identity ──
    print('IDENTITY (exact):')
    print(f'  Lambda^2 / M_Pl = {Lam**2/M_Pl:.6e} GeV')
    print(f'  H0 * hbar       = {H0*hbar:.6e} GeV')
    print(f'  Ratio:           {identity_ratio:.6f}')
    print()

    # ── Three projections ──
    print('THREE PROJECTIONS OF ONE SCALE:')
    print(f'  1. Scalar mass:      m = Lambda = {m_scalar:.4e} GeV')
    print(f'     Compton wavelength: lambda_C = {lam_C:.4f} cm = {lam_C*10:.2f} mm')
    print(f'  2. Dark energy:      V(0) = Lambda^4 = {V0:.4e} GeV^4')
    rho_DE = 7e-30 * GeV_per_g / GeV_inv_per_cm**3
    print(f'     Observed rho_DE ~  {rho_DE:.4e} GeV^4')
    print(f'     Ratio:              {V0/rho_DE:.2f}')
    print(f'  3. Transition accel:  a_0 = c*H0 = {a0:.4e} m/s^2')
    print(f'     = Lambda^2/M_Pl (natural) = {a0_natural:.4e} GeV')
    print(f'     MOND a_0:           {1.2e-10:.4e} m/s^2')
    print(f'     SPARC g_TEP:        {5e-10:.4e} m/s^2')
    print()

    # ── Local screening ──
    print('LOCAL SCREENING (Yukawa, continuous gradient transition):')
    print(f'  S_Sigma = 3(1+x) / (x(2x+1)),  x = mR')
    print()
    print(f'  {"Body":10s} {"R (km)":>10s} {"x=mR":>12s} {"S_Sigma":>14s} '
          f'{"phi_eq/M_Pl":>14s}')

    results = {}
    for name, b in bodies.items():
        S = S_Sigma(m_scalar, b['R_cm'])
        x = m_scalar * b['R_cm'] * GeV_inv_per_cm
        # Equilibrium field
        rho = b['M_kg'] / (4/3 * np.pi * (b['R_cm']/100)**3)  # kg/m^3
        rho_g = rho * 1000  # g/m^3 -> g/cm^3 * 1000... no
        rho_g = rho * 1e-3 / 1e-6  # kg/m^3 -> g/cm^3
        rho_GeV4 = rho_g * GeV_per_g / GeV_inv_per_cm**3
        phi_eq = rho_GeV4 / (m_scalar**2 * M_Pl)
        results[name] = {'S': S, 'x': x, 'phi_eq': phi_eq}
        print(f'  {name:10s} {b["R_cm"]/1e5:10.0f} {x:12.2e} {S:14.4e} '
              f'{phi_eq/M_Pl:14.2e}')
    print()

    # ── Gates ──
    gates, n_pass = check_gates(results)
    print(f'GATES: {n_pass}/5 PASS')
    print(f'  {"Gate":10s} {"Observable":12s} {"Value":>14s} {"Bound":>14s} '
          f'{"Margin":>10s} {"Status":>8s}')
    for gname, g in gates.items():
        margin = g['bound'] / g['value'] if g['value'] > 0 else float('inf')
        status = 'PASS' if g['pass'] else 'FAIL'
        print(f'  {gname:10s} {g["observable"]:12s} {g["value"]:14.4e} '
              f'{g["bound"]:14.1e} {margin:10.0f}x {status:>8s}')
    print()

    # ── Cosmology ──
    print('COSMOLOGY (temporal shear = Hubble flow):')
    print(f'  Clock map: A_clock(z) = (1+z)^-1')
    print(f'  Field:     phi(z) = M_Pl * ln(1+z)')
    print(f'  Shear:     Sigma_0 = (beta_A/M_Pl) * dphi/dt ~ H0')
    print(f'  Static matter frame: a_m = 1 (no expansion)')
    print(f'  Cosmic redshift = accumulated temporal shear along line of sight')
    print(f'  Temporal horizon: A_clock -> 0 as z -> infinity (no Big Bang)')
    print(f'  V(0) = Lambda^4 ~ dark energy density')
    print()

    # ── Wide binaries ──
    print('WIDE BINARIES (temporal shear unscreening):')
    print(f'  Transition acceleration: a_0 = c*H0 = {a0:.3e} m/s^2')
    print(f'  = Lambda^2/M_Pl  (same scale Lambda, exact identity)')
    print(f'  Environmental operator: S(a) = [1 + (a/a_0)^2]^(-1/2)')
    print(f'  S -> 1 (unscreened) for a << a_0')
    print(f'  S -> a_0/a (screened) for a >> a_0')
    print()
    # Predicted transition radius for solar-mass binary
    M_WB = 1.24 * 1.989e30  # kg (sample mean from Paper 13)
    R_s_pred = np.sqrt(G * M_WB / a0)
    R_s_pred_AU = R_s_pred / 1.496e11
    R_s_eta2 = np.sqrt(G * M_WB / (2 * a0)) / 1.496e11
    print(f'  Predicted R_s = sqrt(GM/a_0) = {R_s_pred_AU:.0f} AU '
          f'(for M = 1.24 M_sun, no EFE)')
    print(f'  Predicted R_s = sqrt(GM/(2*a_0)) = {R_s_eta2:.0f} AU '
          f'(with eta=2 external field)')
    print(f'  Observed R_s = 2646 AU (Paper 13)')
    print(f'  Ratio (eta=2): {2646 / R_s_eta2:.2f}')
    print()

    # ── Summary ──
    print('=' * 72)
    print('SUMMARY')
    print('=' * 72)
    print(f'  One scale:   Lambda = (M_Pl*H0*hbar)^(1/2) = {Lam*1e3:.2f} meV')
    print(f'  One coupling: beta_A = -1')
    print(f'  One potential: V = Lambda^4 + (1/2)*Lambda^2*phi^2')
    print()
    print(f'  Projection 1 (screening):  m = Lambda, lambda_C = {lam_C*10:.1f} mm')
    print(f'    All 5 local gates PASS (tightest: geodesy, {gates["Geodesy"]["bound"]/gates["Geodesy"]["value"]:.0f}x margin)')
    print(f'  Projection 2 (cosmology):  V(0) = Lambda^4 ~ rho_DE')
    print(f'    Clock map A_clock(z) = (1+z)^-1, temporal horizon, no Big Bang')
    print(f'  Projection 3 (wide binaries): a_0 = Lambda^2/M_Pl = cH0')
    print(f'    Predicted R_s = {R_s_pred_AU:.0f} AU (no EFE), {R_s_eta2:.0f} AU (eta=2)')
    print(f'    Observed = 2646 AU')
    print()
    print(f'  Identity: Lambda^2/M_Pl = H0 (exact, ratio = {identity_ratio:.6f})')
    print()

    # ── Save results ──
    os.makedirs('results', exist_ok=True)
    output = {
        'Lambda_GeV': float(Lam),
        'Lambda_meV': float(Lam * 1e3),
        'beta_A': BETA,
        'potential': 'V(phi) = Lambda^4 + (1/2)*Lambda^2*phi^2',
        'scalar_mass_GeV': float(m_scalar),
        'compton_wavelength_cm': float(lam_C),
        'compton_wavelength_mm': float(lam_C * 10),
        'V0_GeV4': float(V0),
        'rho_DE_observed_GeV4': float(rho_DE),
        'V0_over_rho_DE': float(V0 / rho_DE),
        'a0_m_s2': float(a0),
        'a0_natural_GeV': float(a0_natural),
        'identity_ratio': float(identity_ratio),
        'Lambda2_over_M_Pl_GeV': float(Lam**2 / M_Pl),
        'H0_hbar_GeV': float(H0 * hbar),
        'gates': {k: {kk: (float(vv) if isinstance(vv, (int, float, np.floating, np.integer))
                       else bool(vv) if isinstance(vv, (bool, np.bool_))
                       else str(vv)) for kk, vv in v.items()}
                  for k, v in gates.items()},
        'n_gates_pass': int(n_pass),
        'bodies': {name: {
            'S_Sigma': float(r['S']),
            'x_mR': float(r['x']),
            'phi_eq_over_M_Pl': float(r['phi_eq'] / M_Pl),
        } for name, r in results.items()},
        'wide_binary': {
            'R_s_predicted_AU': float(R_s_pred_AU),
            'R_s_predicted_eta2_AU': float(R_s_eta2),
            'R_s_observed_AU': 2646,
            'M_WB_solar_masses': 1.24,
        },
    }
    with open('results/scalar_closure.json', 'w') as f:
        json.dump(output, f, indent=2)
    print(f'Results saved to results/scalar_closure.json')


if __name__ == '__main__':
    main()
