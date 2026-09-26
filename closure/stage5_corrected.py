"""Stage 5 CORRECTED: Optimal eta with proper Nordtvedt and GNSS formulas.

Key corrections:
1. Nordtvedt: eta_N = 4*beta_A^2 * (s_Earth - s_Moon)
   where s_A = (GM_A/(R_A c^2)) * S_A  (sensitivity = binding energy * screening)
   NOT just S_A directly.

2. GNSS: The TEP residual from GR is:
   Delta_TEP - Delta_GR = S_Earth^2 * Delta_GR / 2
   (from the PPN gamma correction: gamma - 1 = 2*beta_A^2*S^2)

3. Cassini: gamma - 1 = -4*beta_A^2*S_Sun/(1+2*beta_A^2*S_Sun), |gamma-1| < 2.3e-5
   => S_Sun < 2.3e-5/(4*beta_A^2) ~ 5.75e-6 (linear screened-source map)
"""

import numpy as np
import sys
sys.path.insert(0, '/Users/matthewsmawfield/www/Temporal Equivalence Principle/closure')

from radial_bvp_dimless import (
    DimlessCosh, make_layered_sphere_dimless, make_uniform_sphere_dimless,
    solve_bvp_dimless, extract_observables_dimless,
    BETA_A, to_natural_density, to_natural_length, compute_lambda, M_Pl,
)

# Physical constants
c = 3e10  # cm/s
G = 6.674e-8  # cm^3/(g s^2)

# Body parameters
bodies = {
    'Earth': {'M_g': 5.972e27, 'R_cm': 6.371e8, 'rho_gcc': 5.51,
              'layers': [(0.55, 13.0/5.51), (0.90, 4.5/5.51), (1.0, 2.7/5.51)]},
    'Sun':   {'M_g': 1.989e33, 'R_cm': 6.96e10, 'rho_gcc': 1.41,
              'layers': [(0.2, 150.0/1.41), (1.0, 1.0)]},
    'Moon':  {'M_g': 7.342e25, 'R_cm': 1.737e8, 'rho_gcc': 3.34,
              'layers': None},
}

# Self-gravity binding energy fractions
for name, b in bodies.items():
    b['E_grav_frac'] = G * b['M_g'] / (b['R_cm'] * c**2)
    b['R_nat'] = to_natural_length(b['R_cm'])
    b['rho_nat'] = to_natural_density(b['rho_gcc'])
    b['lambda'] = compute_lambda(b['rho_nat'], b['R_nat'])

R_earth_nat = bodies['Earth']['R_nat']

print("=" * 80)
print("CORRECTED OPTIMAL ETA SCAN")
print("=" * 80)
print()
print("Self-gravity binding energy fractions:")
for name, b in bodies.items():
    print(f"  {name}: GM/(Rc^2) = {b['E_grav_frac']:.4e}")
print()

# GNSS parameters
Delta_Phi_GR = 5.28e-10  # GPS gravitational time dilation
GNSS_threshold = 1e-16

# Constraint thresholds
Cassini_S_Sun_max = 5.75e-6
LLR_eta_N_max = 4.4e-4  # APOLLO bound

print(f"Cassini: S_Sun < {Cassini_S_Sun_max}")
print(f"LLR: |eta_N| < {LLR_eta_N_max}")
print(f"GNSS: residual > {GNSS_threshold}")
print()

# Create sources
for name, b in bodies.items():
    if b['layers']:
        x_arr = np.logspace(-4, np.log10(50.0), 800)
        s_arr = np.zeros_like(x_arr)
        prev = 0
        for x_bound, s_val in b['layers']:
            mask = (x_arr > prev) & (x_arr <= x_bound)
            s_arr[mask] = s_val
            prev = x_bound
        s_arr[x_arr > prev] = 0.0
        from radial_bvp_dimless import DimlessSource
        b['source'] = DimlessSource(x_arr, s_arr, name)
    else:
        b['source'] = make_uniform_sphere_dimless(s_core=1.0, x_max=50.0)
    b['source'].lambda_param = b['lambda']

results = []

for log_eta in np.arange(-0.5, 2.5, 0.1):
    eta_earth = 10.0**log_eta

    # Solve each body
    solved = {}
    for name, b in bodies.items():
        eta = eta_earth * (b['R_nat'] / R_earth_nat)**2
        pot = DimlessCosh(eta)
        u_eq = pot.equilibrium(b['lambda'], 1.0)

        x_init = np.logspace(-4, 2, 300)
        u_guess = np.where(x_init <= 1.0, u_eq, u_eq * np.exp(-(x_init - 1.0)))

        x, u, u_p, sol = solve_bvp_dimless(pot, b['source'], 0.0,
                                           tol=1e-10, u_init_arr=u_guess,
                                           max_nodes=100000)
        if not sol.success:
            solved[name] = None
            continue

        obs = extract_observables_dimless(x, u, u_p, b['source'], pot, 0.0, b['lambda'])
        solved[name] = {'S_Sigma': abs(obs['S_Sigma']), 'obs': obs}

    if any(v is None for v in solved.values()):
        continue

    S_earth = solved['Earth']['S_Sigma']
    S_sun = solved['Sun']['S_Sigma']
    S_moon = solved['Moon']['S_Sigma']

    # CORRECTED Nordtvedt: eta_N = 4*beta_A^2*(s_Earth - s_Moon)
    # s_A = (GM_A/(R_A c^2)) * S_A
    s_earth = bodies['Earth']['E_grav_frac'] * S_earth
    s_moon = bodies['Moon']['E_grav_frac'] * S_moon
    eta_N = 4 * BETA_A**2 * (s_earth - s_moon)

    # CORRECTED GNSS: residual = |gamma-1|/2 * Delta_GR with the linear
    # screened-source map gamma-1 = -4*beta_A^2*S/(1+2*beta_A^2*S)
    gnss_residual = (2 * BETA_A**2 * S_earth / (1 + 2*BETA_A**2*S_earth)) * Delta_Phi_GR

    # CORRECTED Cassini: gamma-1 = -4*beta_A^2*S_Sun/(1+2*beta_A^2*S_Sun)
    # (linear screened-source cross-term, unscreened photon probe)
    cassini_gamma = abs(4 * BETA_A**2 * S_sun / (1 + 2*BETA_A**2*S_sun))

    cassini_ok = cassini_gamma < 2.3e-5
    llr_ok = abs(eta_N) < LLR_eta_N_max
    gnss_detectable = gnss_residual > GNSS_threshold

    results.append({
        'log_eta': log_eta,
        'eta': eta_earth,
        'S_earth': S_earth,
        'S_sun': S_sun,
        'S_moon': S_moon,
        's_earth': s_earth,
        's_moon': s_moon,
        'eta_N': eta_N,
        'gnss_residual': gnss_residual,
        'cassini_gamma': cassini_gamma,
        'cassini_ok': cassini_ok,
        'llr_ok': llr_ok,
        'gnss_detectable': gnss_detectable,
    })

# Print results
print(f"{'log_eta':>8s}  {'S_earth':>10s}  {'S_sun':>10s}  {'S_moon':>10s}  "
      f"{'eta_N':>12s}  {'GNSS_res':>12s}  {'gamma-1':>12s}  {'C':>2s} {'L':>2s} {'G':>2s}")
print("-" * 95)
for r in results:
    c = 'Y' if r['cassini_ok'] else 'N'
    l = 'Y' if r['llr_ok'] else 'N'
    g = 'Y' if r['gnss_detectable'] else 'N'
    print(f"{r['log_eta']:8.2f}  {r['S_earth']:10.3e}  {r['S_sun']:10.3e}  "
          f"{r['S_moon']:10.3e}  {r['eta_N']:12.4e}  {r['gnss_residual']:12.4e}  "
          f"{r['cassini_gamma']:12.4e}  {c} {l} {g}")

print()

# Find optimal: all constraints satisfied AND GNSS detectable
all_ok = [r for r in results if r['cassini_ok'] and r['llr_ok'] and r['gnss_detectable']]
if all_ok:
    print("ALL CONSTRAINTS SATISFIED WITH GNSS DETECTABLE:")
    for r in all_ok:
        print(f"  log(eta) = {r['log_eta']:.2f}: S_earth={r['S_earth']:.3e}, "
              f"eta_N={r['eta_N']:.3e}, GNSS={r['gnss_residual']:.3e}")
    best = max(all_ok, key=lambda r: r['gnss_residual'])
    print(f"\n  OPTIMAL: log(eta) = {best['log_eta']:.2f}")
    print(f"  S_earth = {best['S_earth']:.4e}")
    print(f"  S_sun = {best['S_sun']:.4e}")
    print(f"  S_moon = {best['S_moon']:.4e}")
    print(f"  eta_N = {best['eta_N']:.4e} (bound: {LLR_eta_N_max})")
    print(f"  gamma-1 = {best['cassini_gamma']:.4e} (bound: 2.3e-5)")
    print(f"  GNSS residual = {best['gnss_residual']:.4e} (threshold: {GNSS_threshold})")
else:
    llr_cassini = [r for r in results if r['cassini_ok'] and r['llr_ok']]
    if llr_cassini:
        best = max(llr_cassini, key=lambda r: r['gnss_residual'])
        print(f"Best (LLR+Cassini OK, GNSS marginal):")
        print(f"  log(eta) = {best['log_eta']:.2f}, GNSS = {best['gnss_residual']:.3e}")
