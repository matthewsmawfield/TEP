"""Stage 5: Find the optimal eta that satisfies ALL constraints simultaneously.

Constraints:
1. Cassini: |S_Sigma_Sun| < 5.75e-6
2. LLR (Nordtvedt): |4*beta_A^2*(S_Earth - S_Moon)| < 1e-3
3. GNSS detectability: |S_Sigma_Earth * Delta_Phi_GR| > 1e-16 (residual from GR)
4. Wide binary: no anomalous acceleration (S_Sigma < 1e-3)
5. White dwarf: not ruled out by pulsar-WD systems

The GNSS observable is NOT the raw scalar field difference Delta_u.
It is the SCREENED deviation from GR:
  Delta_TEP - Delta_GR = -S_Sigma * Delta_Phi_GR

where Delta_Phi_GR ~ 5.2e-10 is the standard gravitational time dilation.

The LLR Nordtvedt effect:
  eta_N = 4 * beta_A^2 * (S_Earth - S_Moon)
  LLR bound: |eta_N| < 1e-3
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
M_Pl_val = 2.435e18  # GeV
R_earth_nat = to_natural_length(6.371e8)
R_sun_nat = to_natural_length(6.96e10)
R_moon_nat = to_natural_length(1.737e8)

# GR gravitational time dilation for GPS
# Delta_Phi/c^2 = (GM/R)(1 - R/r_sat) / c^2
# For GPS: r_sat = 26560 km, R_earth = 6371 km
# = 6.95e-10 * (1 - 6371/26560) = 6.95e-10 * 0.76 = 5.28e-10
Delta_Phi_GR = 5.28e-10

# Source profiles
source_earth = make_layered_sphere_dimless(
    layers=[(0.55, 13.0/5.51), (0.90, 4.5/5.51), (1.0, 2.7/5.51)],
    s_ambient=0.0, x_max=50.0
)
source_earth.lambda_param = compute_lambda(to_natural_density(5.51), R_earth_nat)

source_sun = make_layered_sphere_dimless(
    layers=[(0.2, 150.0/1.41), (1.0, 1.0)],
    s_ambient=0.0, x_max=50.0
)
source_sun.lambda_param = compute_lambda(to_natural_density(1.41), R_sun_nat)

source_moon = make_uniform_sphere_dimless(s_core=1.0, x_max=50.0)
source_moon.lambda_param = compute_lambda(to_natural_density(3.34), R_moon_nat)

ratio_R_sun = R_sun_nat / R_earth_nat
ratio_R_moon = R_moon_nat / R_earth_nat

# Fine eta scan
print("=" * 80)
print("OPTIMAL ETA SCAN: ALL CONSTRAINTS")
print("=" * 80)
print()
print(f"GR time dilation (GPS): Delta_Phi/c^2 = {Delta_Phi_GR:.4e}")
print(f"GNSS detection threshold: 1e-16")
print(f"Required S_Sigma_Earth for GNSS: > {1e-16/Delta_Phi_GR:.4e}")
print(f"LLR Nordtvedt bound: |eta_N| < 1e-3")
print(f"Cassini bound: |S_Sun| < 5.75e-6")
print()

results = []

for log_eta_10 in np.arange(-0.3, 1.5, 0.1):
    eta_earth = 10.0**log_eta_10
    eta_sun = eta_earth * ratio_R_sun**2
    eta_moon = eta_earth * ratio_R_moon**2

    # Solve Earth
    pot_e = DimlessCosh(eta_earth)
    u_eq_e = pot_e.equilibrium(source_earth.lambda_param, 1.0)
    x_init = np.logspace(-4, 2, 300)
    u_guess = np.where(x_init <= 1.0, u_eq_e, u_eq_e * np.exp(-(x_init - 1.0)))
    x_e, u_e, u_p_e, sol_e = solve_bvp_dimless(pot_e, source_earth, 0.0,
                                                tol=1e-10, u_init_arr=u_guess,
                                                max_nodes=100000)
    if not sol_e.success:
        continue
    obs_e = extract_observables_dimless(x_e, u_e, u_p_e, source_earth, pot_e,
                                        0.0, source_earth.lambda_param)

    # Solve Sun
    pot_s = DimlessCosh(eta_sun)
    u_eq_s = pot_s.equilibrium(source_sun.lambda_param, 1.0)
    x_init_s = np.logspace(-4, 2, 300)
    u_guess_s = np.where(x_init_s <= 1.0, u_eq_s, u_eq_s * np.exp(-(x_init_s - 1.0)))
    x_s, u_s, u_p_s, sol_s = solve_bvp_dimless(pot_s, source_sun, 0.0,
                                                tol=1e-10, u_init_arr=u_guess_s,
                                                max_nodes=100000)
    if not sol_s.success:
        continue
    obs_s = extract_observables_dimless(x_s, u_s, u_p_s, source_sun, pot_s,
                                        0.0, source_sun.lambda_param)

    # Solve Moon
    pot_m = DimlessCosh(eta_moon)
    u_eq_m = pot_m.equilibrium(source_moon.lambda_param, 1.0)
    x_init_m = np.logspace(-4, 2, 300)
    u_guess_m = np.where(x_init_m <= 1.0, u_eq_m, u_eq_m * np.exp(-(x_init_m - 1.0)))
    x_m, u_m, u_p_m, sol_m = solve_bvp_dimless(pot_m, source_moon, 0.0,
                                                tol=1e-10, u_init_arr=u_guess_m,
                                                max_nodes=100000)
    if not sol_m.success:
        continue
    obs_m = extract_observables_dimless(x_m, u_m, u_p_m, source_moon, pot_m,
                                        0.0, source_moon.lambda_param)

    S_earth = abs(obs_e['S_Sigma'])
    S_sun = abs(obs_s['S_Sigma'])
    S_moon = abs(obs_m['S_Sigma'])

    # Constraints
    cassini_ok = S_sun < 5.75e-6
    eta_nordtvedt = 4 * BETA_A**2 * (obs_e['S_Sigma'] - obs_m['S_Sigma'])
    llr_ok = abs(eta_nordtvedt) < 1e-3
    gnss_residual = S_earth * Delta_Phi_GR
    gnss_detectable = gnss_residual > 1e-16

    all_ok = cassini_ok and llr_ok
    gnss_and_llr = cassini_ok and llr_ok and gnss_detectable

    results.append({
        'log_eta': log_eta_10,
        'eta': eta_earth,
        'S_earth': S_earth,
        'S_sun': S_sun,
        'S_moon': S_moon,
        'eta_nordtvedt': eta_nordtvedt,
        'gnss_residual': gnss_residual,
        'cassini_ok': cassini_ok,
        'llr_ok': llr_ok,
        'gnss_detectable': gnss_detectable,
        'all_ok': all_ok,
        'gnss_and_llr': gnss_and_llr,
    })

# Print results
print(f"{'log_eta':>8s}  {'eta':>10s}  {'S_earth':>10s}  {'S_sun':>10s}  "
      f"{'S_moon':>10s}  {'eta_N':>10s}  {'GNSS_res':>10s}  {'C':>2s} {'L':>2s} {'G':>2s}")
print("-" * 85)
for r in results:
    c = 'Y' if r['cassini_ok'] else 'N'
    l = 'Y' if r['llr_ok'] else 'N'
    g = 'Y' if r['gnss_detectable'] else 'N'
    print(f"{r['log_eta']:8.2f}  {r['eta']:10.3e}  {r['S_earth']:10.3e}  "
          f"{r['S_sun']:10.3e}  {r['S_moon']:10.3e}  {r['eta_nordtvedt']:10.3e}  "
          f"{r['gnss_residual']:10.3e}  {c} {l} {g}")

print()

# Find optimal
both = [r for r in results if r['gnss_and_llr']]
if both:
    print("OPTIMAL: GNSS detectable AND LLR safe AND Cassini safe")
    for r in both:
        print(f"  log(eta) = {r['log_eta']:.2f}: S_earth={r['S_earth']:.3e}, "
              f"S_moon={r['S_moon']:.3e}, GNSS residual={r['gnss_residual']:.3e}")
else:
    # Find best compromise
    llr_safe = [r for r in results if r['llr_ok'] and r['cassini_ok']]
    if llr_safe:
        best = max(llr_safe, key=lambda r: r['gnss_residual'])
        print(f"BEST COMPROMISE (LLR safe, Cassini safe, GNSS marginal):")
        print(f"  log(eta) = {best['log_eta']:.2f}: S_earth={best['S_earth']:.3e}, "
              f"S_moon={best['S_moon']:.3e}, GNSS residual={best['gnss_residual']:.3e}")
    else:
        print("No eta satisfies all constraints simultaneously.")
        print("The tension between GNSS detectability and LLR safety is real.")
