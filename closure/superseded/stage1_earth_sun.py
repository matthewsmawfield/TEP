"""Stage 1: Earth and Sun profiles with the cosh potential.

V = Lambda^4 (cosh(phi/M_Pl) - 1)

This is the first real TEP closure test. The potential has:
  - Minimum at phi=0 (vacuum)
  - V_{,phi} > 0 for phi > 0 (equilibrium condition with beta_A = -1)
  - Stable everywhere
  - Connects to phi -> +inf (temporal horizon)

The dimensionless equation is:
  u'' + (2/x) u' = eta * sinh(u) + lambda * s(x) * beta_A * exp(beta_A u)

where:
  u = phi/M_Pl, x = r/R, s = rho_*/rho_ref
  eta = Lambda^4 R^2 / M_Pl^2
  lambda = rho_ref R^2 / M_Pl^2

For the Earth and Sun, lambda is very small (~10^-9 and ~10^-5 respectively).
The field excursion is tiny: u ~ lambda/eta for small lambda.

The key observable is S_Sigma = Q/Q_0, extracted from the solved exterior profile.
Cassini requires S_Sigma^(Sun) < 5.75e-6.

NO thin-shell estimates. NO chameleon mechanism. The Temporal Topology is the
ontology; V(phi) generates phi(x); S_Sigma is read off the solved profile.
"""

import numpy as np
import sys
sys.path.insert(0, '/Users/matthewsmawfield/www/Temporal Equivalence Principle/closure')

from radial_bvp_dimless import (
    DimlessCosh, make_uniform_sphere_dimless, make_layered_sphere_dimless,
    solve_bvp_dimless, extract_observables_dimless,
    BETA_A, to_natural_density, to_natural_length, compute_lambda, M_Pl,
)


# ---------------------------------------------------------------------------
# Physical parameters
# ---------------------------------------------------------------------------

M_Pl_val = 2.435e18  # GeV

# Earth
rho_earth_mean_gcc = 5.51    # g/cm^3
rho_earth_core_gcc = 13.0
R_earth_cm = 6.371e8         # cm
M_earth_g = 5.972e27         # g

# Sun
rho_sun_mean_gcc = 1.41
rho_sun_core_gcc = 150.0
R_sun_cm = 6.96e10
M_sun_g = 1.989e33

# Convert to natural units
rho_earth_nat = to_natural_density(rho_earth_mean_gcc)
R_earth_nat = to_natural_length(R_earth_cm)
rho_sun_nat = to_natural_density(rho_sun_mean_gcc)
R_sun_nat = to_natural_length(R_sun_cm)

# Dimensionless groups
lambda_earth = compute_lambda(rho_earth_nat, R_earth_nat)
lambda_sun = compute_lambda(rho_sun_nat, R_sun_nat)

print("TEP Scalar Field Closure Test — Stage 1")
print("Potential: V = Lambda^4 (cosh(phi/M_Pl) - 1)")
print()
print("Physical parameters:")
print(f"  M_Pl = {M_Pl_val:.4e} GeV")
print(f"  Earth: rho = {rho_earth_mean_gcc} g/cm^3, R = {R_earth_cm:.3e} cm")
print(f"  Sun:   rho = {rho_sun_mean_gcc} g/cm^3, R = {R_sun_cm:.3e} cm")
print()
print("Dimensionless groups (lambda = rho_ref R^2 / M_Pl^2):")
print(f"  lambda_earth = {lambda_earth:.4e}")
print(f"  lambda_sun   = {lambda_sun:.4e}")
print(f"  lambda_sun / lambda_earth = {lambda_sun/lambda_earth:.4e}")
print()

# Cassini constraint
S_Sun_max = 5.75e-6
print(f"Cassini constraint: S_Sigma^(Sun) < {S_Sun_max}")
print()


def solve_body(pot, source, u_ambient, name, body_mass_g, R_body_cm):
    """Solve BVP for a spherical body and extract observables."""
    print(f"--- {name} ---")
    print(f"  Potential: {pot.describe()}")

    # Better initial guess: equilibrium inside, ambient outside
    u_eq = pot.equilibrium(source.lambda_param, 1.0)
    print(f"  Equilibrium u_min(rho_ref) = {u_eq:.6e}")

    x_init = np.logspace(-4, np.log10(source.x[-1]), 300)
    u_guess = np.where(x_init <= 1.0, u_eq, u_ambient + (u_eq - u_ambient) * np.exp(-(x_init - 1.0)))

    x, u, u_p, sol = solve_bvp_dimless(pot, source, u_ambient, tol=1e-10,
                                       u_init_arr=u_guess, max_nodes=100000)

    if not sol.success:
        print(f"  BVP FAILED: {sol.message}")
        # Try relaxed
        x, u, u_p, sol = solve_bvp_dimless(pot, source, u_ambient, tol=1e-6,
                                           u_init_arr=u_guess, max_nodes=200000)
        if not sol.success:
            print(f"  BVP FAILED (relaxed): {sol.message}")
            return None
        print(f"  BVP CONVERGED (relaxed tol=1e-6)")
    else:
        print(f"  BVP CONVERGED (tol=1e-10)")

    obs = extract_observables_dimless(x, u, u_p, source, pot, u_ambient,
                                      source.lambda_param)

    # Print key observables
    u_center = u[np.argmin(np.abs(x - 0.01))]
    u_surface = u[np.argmin(np.abs(x - 1.0))]
    u_far = u[np.argmin(np.abs(x - 10.0))]

    print(f"  u(0.01) = {u_center:.6e}  (equilibrium = {u_eq:.6e})")
    print(f"  u(1.0)  = {u_surface:.6e}")
    print(f"  u(10)   = {u_far:.6e}  (ambient = {u_ambient:.6e})")
    print(f"  q = {obs['q']:.6e}")
    print(f"  q_0 = {obs['q_0']:.6e}")
    print(f"  S_Sigma = {obs['S_Sigma']:.6e}")

    # Temporal Shear profile
    Sigma_center = obs['Sigma_r_times_R'][np.argmin(np.abs(x - 0.01))]
    Sigma_surface = obs['Sigma_r_times_R'][np.argmin(np.abs(x - 1.0))]
    print(f"  |Sigma_r * R| at center  = {abs(Sigma_center):.4e}")
    print(f"  |Sigma_r * R| at surface = {abs(Sigma_surface):.4e}")

    # Disformal deformation shape
    eps_peak = np.max(obs['eps_B_shape'])
    eps_at_surface = obs['eps_B_shape'][np.argmin(np.abs(x - 1.0))]
    print(f"  epsilon_B/B0 peak = {eps_peak:.4e}")
    print(f"  epsilon_B/B0 at surface = {eps_at_surface:.4e}")

    print()
    return obs


# ---------------------------------------------------------------------------
# Scan over eta (the single potential parameter)
# ---------------------------------------------------------------------------

print("=" * 70)
print("ETA SCAN: V = Lambda^4 (cosh(phi/M_Pl) - 1)")
print("=" * 70)
print()
print("eta = Lambda^4 R^2 / M_Pl^2 controls the potential strength.")
print("Small eta -> weak potential -> field barely responds -> large S_Sigma")
print("Large eta -> strong potential -> field pinned -> small S_Sigma")
print()
print("For each eta, solve Earth and Sun with the SAME potential parameter.")
print("Check: does S_Sigma^(Sun) < 5.75e-6?")
print()

# Use Earth as the reference body (R = R_earth, rho_ref = rho_earth_mean)
# For the Sun, rescale: x = r/R_sun, s = rho/rho_sun_mean
# eta_sun = Lambda^4 R_sun^2 / M_Pl^2 = eta_earth * (R_sun/R_earth)^2

ratio_R = R_sun_nat / R_earth_nat
print(f"R_sun/R_earth = {ratio_R:.4f}")
print(f"eta_sun = eta_earth * (R_sun/R_earth)^2 = eta_earth * {ratio_R**2:.4f}")
print()

# Source profiles
source_earth = make_layered_sphere_dimless(
    layers=[(0.55, 13.0/5.51), (0.90, 4.5/5.51), (1.0, 2.7/5.51)],
    s_ambient=0.0, x_max=50.0
)
source_earth.lambda_param = lambda_earth

source_sun = make_layered_sphere_dimless(
    layers=[(0.2, 150.0/1.41), (1.0, 1.0)],
    s_ambient=0.0, x_max=50.0
)
source_sun.lambda_param = lambda_sun

results = []

for log_eta in range(-4, 15):
    eta_earth = 10.0**log_eta
    eta_sun = eta_earth * ratio_R**2

    print(f"eta_earth = 10^{log_eta} = {eta_earth:.4e}, eta_sun = {eta_sun:.4e}")

    pot_earth = DimlessCosh(eta_earth)
    pot_sun = DimlessCosh(eta_sun)

    obs_earth = solve_body(pot_earth, source_earth, 0.0, "Earth",
                           M_earth_g, R_earth_cm)
    obs_sun = solve_body(pot_sun, source_sun, 0.0, "Sun",
                         M_sun_g, R_sun_cm)

    if obs_earth is not None and obs_sun is not None:
        S_sun = obs_sun['S_Sigma']
        S_earth = obs_earth['S_Sigma']
        cassini_ok = abs(S_sun) < S_Sun_max

        results.append({
            'log_eta': log_eta,
            'eta_earth': eta_earth,
            'eta_sun': eta_sun,
            'S_earth': S_earth,
            'S_sun': S_sun,
            'cassini_ok': cassini_ok,
            'u_earth_center': obs_earth['u'][np.argmin(np.abs(obs_earth['x']-0.01))],
            'u_sun_center': obs_sun['u'][np.argmin(np.abs(obs_sun['x']-0.01))],
        })

        status = "CASSINI OK" if cassini_ok else "CASSINI FAIL"
        print(f"  -> S_earth = {S_earth:.4e}, S_sun = {S_sun:.4e} [{status}]")
    else:
        results.append({
            'log_eta': log_eta,
            'eta_earth': eta_earth,
            'eta_sun': eta_sun,
            'S_earth': None,
            'S_sun': None,
            'cassini_ok': False,
        })
        print(f"  -> SOLVER FAILED")
    print()

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

print("=" * 70)
print("SCAN SUMMARY")
print("=" * 70)
print()
print(f"{'log_eta':>8s}  {'eta_earth':>12s}  {'S_earth':>12s}  {'S_sun':>12s}  {'Cassini':>8s}")
print("-" * 60)
for r in results:
    s_e = f"{r['S_earth']:.4e}" if r['S_earth'] is not None else "FAILED"
    s_s = f"{r['S_sun']:.4e}" if r['S_sun'] is not None else "FAILED"
    cassini = "OK" if r['cassini_ok'] else "FAIL"
    print(f"{r['log_eta']:8d}  {r['eta_earth']:12.4e}  {s_e:>12s}  {s_s:>12s}  {cassini:>8s}")

print()
passing = [r for r in results if r['cassini_ok']]
if passing:
    print(f"Cassini constraint satisfied for {len(passing)} eta values:")
    for r in passing:
        print(f"  log(eta_earth) = {r['log_eta']}: S_sun = {r['S_sun']:.4e}, S_earth = {r['S_earth']:.4e}")
else:
    print("No eta values satisfy Cassini. The cosh potential may be insufficient")
    print("for the minimal K=1 closure, or the parameter range needs extension.")
