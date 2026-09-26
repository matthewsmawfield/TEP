"""Stage 3: Multi-domain profile scan.

Solve the scalar BVP for every TEP environment using the SAME cosh V
parameter (eta). This tests whether one V generates correct screening
across all weak-field domains:

  - Earth (GNSS, clock networks)
  - Sun (Cassini, LLR)
  - Wide binary (Paper 14 held-out prediction)
  - Ultra-cool dwarf (Paper 15)
  - Diffuse ambient / cosmological void (Paper 13, Paper 31)
  - White dwarf (Paper 15 compact object)

For each environment, extract:
  - S_Sigma (source-charge screening)
  - S_A (clock-rate amplitude, from Delta ln A)
  - Sigma_r profile (Temporal Shear)
  - epsilon_B shape (disformal deformation)

No thin-shell. No chameleon. Direct BVP solve and observable extraction.
"""

import numpy as np
import sys
sys.path.insert(0, '/Users/matthewsmawfield/www/Temporal Equivalence Principle/closure')

from radial_bvp_dimless import (
    DimlessCosh, DimlessSource, solve_bvp_dimless,
    extract_observables_dimless, BETA_A,
    to_natural_density, to_natural_length, compute_lambda, M_Pl,
)
from scipy.integrate import trapezoid

# ---------------------------------------------------------------------------
# Environment definitions
# ---------------------------------------------------------------------------

environments = []

def add_env(name, rho_mean_gcc, R_cm, M_g, rho_core_gcc=None,
            rho_ambient_gcc=0.0, x_max=50.0, layers=None):
    """Add an environment to the scan."""
    rho_ref = to_natural_density(rho_mean_gcc)
    R_nat = to_natural_length(R_cm)
    lam = compute_lambda(rho_ref, R_nat)

    if layers:
        # Layered profile
        x_arr = np.logspace(-4, np.log10(x_max), 800)
        s_arr = np.zeros_like(x_arr)
        prev_bound = 0
        for x_bound, s_val in layers:
            mask = (x_arr > prev_bound) & (x_arr <= x_bound)
            s_arr[mask] = s_val
            prev_bound = x_bound
        s_arr[x_arr > prev_bound] = rho_ambient_gcc / rho_mean_gcc
        source = DimlessSource(x_arr, s_arr, name)
    else:
        # Uniform sphere
        x_arr = np.logspace(-4, np.log10(x_max), 500)
        s_arr = np.where(x_arr <= 1.0, 1.0, rho_ambient_gcc / rho_mean_gcc)
        source = DimlessSource(x_arr, s_arr, name)

    source.lambda_param = lam
    environments.append({
        'name': name,
        'rho_mean': rho_mean_gcc,
        'R_cm': R_cm,
        'M_g': M_g,
        'rho_ref_nat': rho_ref,
        'R_nat': R_nat,
        'lambda': lam,
        'source': source,
        'rho_ambient': rho_ambient_gcc,
    })

# Earth: mean 5.51 g/cm^3, R=6371 km, M=5.97e27 g
add_env('Earth', 5.51, 6.371e8, 5.972e27,
        layers=[(0.55, 13.0/5.51), (0.90, 4.5/5.51), (1.0, 2.7/5.51)],
        rho_ambient_gcc=1e-24, x_max=50.0)

# Sun: mean 1.41 g/cm^3, R=6.96e10 cm, M=1.989e33 g
add_env('Sun', 1.41, 6.96e10, 1.989e33,
        layers=[(0.2, 150.0/1.41), (1.0, 1.0)],
        rho_ambient_gcc=1e-24, x_max=50.0)

# Wide binary: two ~1 M_sun stars, separation ~7000 AU
# Treat as single body with effective radius = separation/2
# Ambient density: ~0.1 atoms/cm^3 ~ 1.7e-25 g/cm^3
# But for screening, the relevant "body" is each star in the ambient field
# Use M_sun parameters with very low ambient
add_env('WideBinary_star', 1.41, 6.96e10, 1.989e33,
        layers=[(0.2, 150.0/1.41), (1.0, 1.0)],
        rho_ambient_gcc=1.7e-25, x_max=500.0)  # lower ambient, wider domain

# Ultra-cool dwarf: M~0.08 M_sun, R~0.1 R_sun, rho~10 g/cm^3
add_env('UCD', 10.0, 0.1 * 6.96e10, 0.08 * 1.989e33,
        rho_ambient_gcc=1e-24, x_max=50.0)

# White dwarf: M~0.6 M_sun, R~0.01 R_sun, rho~10^6 g/cm^3
add_env('WhiteDwarf', 1e6, 0.01 * 6.96e10, 0.6 * 1.989e33,
        rho_ambient_gcc=1e-24, x_max=50.0)

# Neutron star: M~1.4 M_sun, R~10 km, rho~10^14 g/cm^3
add_env('NeutronStar', 1e14, 1e6, 1.4 * 1.989e33,
        rho_ambient_gcc=1e-24, x_max=50.0)

# Diffuse ambient (cosmological void): rho~10^-29 g/cm^3, no body
# Use a very large "body" with density = ambient (tests the unscreened branch)
add_env('DiffuseVoid', 1e-29, 1e26, 0,  # 1 Mpc scale, no real body
        rho_ambient_gcc=1e-29, x_max=10.0)

# Moon: M~7.3e25 g, R~1.7e8 cm, rho~3.34 g/cm^3
add_env('Moon', 3.34, 1.737e8, 7.342e25,
        rho_ambient_gcc=1e-24, x_max=50.0)

# Jupiter: M~1.9e30 g, R~7e9 cm, rho~1.33 g/cm^3
add_env('Jupiter', 1.33, 7.0e9, 1.898e30,
        rho_ambient_gcc=1e-24, x_max=50.0)


# ---------------------------------------------------------------------------
# Solve all environments for a given eta
# ---------------------------------------------------------------------------

def solve_all(eta_earth, verbose=False):
    """Solve BVP for all environments with the same potential parameter.

    eta_earth is the eta value for Earth (rho_ref = rho_earth, R = R_earth).
    For other environments, eta scales as (R/R_earth)^2 * (rho_ref/rho_earth_ref).
    But since we use each environment's own rho_ref and R, the eta for each
    environment is eta_earth * (R_env/R_earth)^2 * (rho_env_ref/rho_earth_ref)
    ... actually no. The potential parameter Lambda is FIXED.
    eta = Lambda^4 R^2 / M_Pl^2, so eta_env = Lambda^4 R_env^2 / M_Pl^2.
    Since Lambda is the same: eta_env / eta_earth = (R_env / R_earth)^2.
    """
    R_earth_nat = to_natural_length(6.371e8)

    results = {}
    for env in environments:
        name = env['name']
        R_env = env['R_nat']
        eta_env = eta_earth * (R_env / R_earth_nat)**2

        pot = DimlessCosh(eta_env)
        u_eq = pot.equilibrium(env['lambda'], 1.0)

        # Initial guess
        x_init = np.logspace(-4, np.log10(env['source'].x[-1]), 300)
        u_guess = np.where(x_init <= 1.0, u_eq,
                          (u_eq) * np.exp(-(x_init - 1.0)))

        x, u, u_p, sol = solve_bvp_dimless(pot, env['source'], 0.0,
                                           tol=1e-10, u_init_arr=u_guess,
                                           max_nodes=100000)

        if not sol.success:
            if verbose:
                print(f"  {name}: FAILED - {sol.message}")
            # Try relaxed
            x, u, u_p, sol = solve_bvp_dimless(pot, env['source'], 0.0,
                                               tol=1e-6, u_init_arr=u_guess,
                                               max_nodes=200000)
            if not sol.success:
                results[name] = None
                if verbose:
                    print(f"  {name}: FAILED (relaxed)")
                continue

        obs = extract_observables_dimless(x, u, u_p, env['source'], pot,
                                          0.0, env['lambda'])

        # Additional observables
        u_center = u[np.argmin(np.abs(x - 0.01))]
        u_surface = u[np.argmin(np.abs(x - 1.0))]
        Delta_ln_A = abs(BETA_A) * abs(u_center - 0)  # conformal clock shift
        S_A_shape = Delta_ln_A  # fractional clock-rate shift

        results[name] = {
            'x': x, 'u': u, 'u_p': u_p, 'obs': obs,
            'u_center': u_center, 'u_surface': u_surface,
            'Delta_ln_A': Delta_ln_A,
            'S_Sigma': obs['S_Sigma'],
            'eta_env': eta_env,
            'u_eq': u_eq,
            'converged': True,
        }
        if verbose:
            print(f"  {name}: u_center={u_center:.4e}, S_Sigma={obs['S_Sigma']:.4e}, "
                  f"Delta_ln_A={Delta_ln_A:.4e}")

    return results


# ---------------------------------------------------------------------------
# Run the scan
# ---------------------------------------------------------------------------

print("=" * 80)
print("MULTI-DOMAIN PROFILE SCAN")
print("V = Lambda^4 (cosh(phi/M_Pl) - 1), K=1, beta_A = -1")
print("=" * 80)
print()

print("Environments:")
for env in environments:
    print(f"  {env['name']:15s}: rho={env['rho_mean']:.2e} g/cm^3, "
          f"R={env['R_cm']:.2e} cm, lambda={env['lambda']:.4e}")
print()

# Test several eta values
for log_eta in [-1, 0, 1]:
    eta = 10.0**log_eta
    print(f"\n{'='*80}")
    print(f"eta_earth = 10^{log_eta} = {eta:.4e}")
    print(f"{'='*80}")
    print()

    results = solve_all(eta, verbose=True)

    # Summary table
    print()
    print(f"{'Environment':>15s}  {'lambda':>12s}  {'eta_env':>12s}  "
          f"{'u_center':>12s}  {'S_Sigma':>12s}  {'Delta_ln_A':>12s}  "
          f"{'eps_B_peak':>12s}")
    print("-" * 100)

    for env in environments:
        name = env['name']
        r = results.get(name)
        if r is None:
            print(f"{name:>15s}  {env['lambda']:12.4e}  {'FAILED':>12s}")
        else:
            eps_peak = np.max(r['obs']['eps_B_shape'])
            print(f"{name:>15s}  {env['lambda']:12.4e}  {r['eta_env']:12.4e}  "
                  f"{r['u_center']:12.4e}  {r['S_Sigma']:12.4e}  "
                  f"{r['Delta_ln_A']:12.4e}  {eps_peak:12.4e}")

    # Check constraints
    print()
    sun = results.get('Sun')
    earth = results.get('Earth')
    if sun:
        cassini_ok = abs(sun['S_Sigma']) < 5.75e-6
        print(f"  Cassini (Sun S_Sigma < 5.75e-6): {'PASS' if cassini_ok else 'FAIL'} "
              f"[S_Sigma = {sun['S_Sigma']:.4e}]")
    if earth:
        print(f"  Earth Delta_ln_A (clock shift): {earth['Delta_ln_A']:.4e}")
        print(f"  GNSS observable scale (~10^-16): "
              f"{'in range' if 1e-20 < earth['Delta_ln_A'] < 1e-5 else 'check needed'}")
