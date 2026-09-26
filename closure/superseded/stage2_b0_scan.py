"""Stage 2: B0 window scan — Cassini-safe AND triangle-detectable?

For each eta, the solved Earth and Sun profiles give us epsilon_B(r)/B0.
The actual disformal deformation is:

  epsilon_B(r) = B0 * (eps_B_shape)(r)

where eps_B_shape = (u^2/(1+u^2)) * (u')^2 / A^2  (envelope = 1 in weak field).

Constraints:
  1. Cassini/GW170817: epsilon_B along solar propagation paths < few * 10^-15
  2. Triangle detection: integrated holonomy H_resid > 10^-18 (target sensitivity)

The triangle holonomy is:
  H_resid ~ oint J_mu dx^mu

where J_mu = (B/A^2)(u.nabla_phi) P_mu^nu nabla_nu phi.

For a radial field and a triangle with one leg along the radial direction:
  H_resid ~ B0 * integral of eps_B_shape along the radial leg

The scan is 1D in B0 (as verified: the quartic Gaussian envelope is inert
in the weak-field regime, so sigma_B drops out).

We plot three regions on one axis:
  - Allowed by Cassini (epsilon_B at solar impact parameter < 10^-15)
  - Allowed by GW170817 (epsilon_B along propagation path < 10^-15)
  - Detectable by triangle (H_resid > 10^-18)

If they don't intersect, the window is empty. If they do, the overlap is
the prediction.
"""

import numpy as np
import sys
sys.path.insert(0, '/Users/matthewsmawfield/www/Temporal Equivalence Principle/closure')

from radial_bvp_dimless import (
    DimlessCosh, make_layered_sphere_dimless,
    solve_bvp_dimless, extract_observables_dimless,
    BETA_A, to_natural_density, to_natural_length, compute_lambda, M_Pl,
)
from scipy.integrate import trapezoid

# Physical parameters (same as Stage 1)
rho_earth_nat = to_natural_density(5.51)
R_earth_nat = to_natural_length(6.371e8)
rho_sun_nat = to_natural_density(1.41)
R_sun_nat = to_natural_length(6.96e10)
lambda_earth = compute_lambda(rho_earth_nat, R_earth_nat)
lambda_sun = compute_lambda(rho_sun_nat, R_sun_nat)
ratio_R = R_sun_nat / R_earth_nat

# Triangle geometry: ground-to-satellite leg
# Satellite at ~MEO: x_sat ~ 4 R_earth (GPS altitude ~ 20,000 km ~ 3.1 R_earth)
# Ground station at x ~ 1 R_earth
# The radial leg goes from x=1 to x=4
x_ground = 1.0
x_sat = 4.0  # GPS-like altitude

# Cassini impact parameter: ~1 AU = 1.5e13 cm
# In Earth radii: 1 AU / R_earth ~ 2.35e5
# In Sun radii: 1 AU / R_sun ~ 2155
# For the Sun profile, the solar impact parameter is at x ~ 2155 (in R_sun units)
x_cassini_sun = 2155.0  # 1 AU in R_sun units

# GW170817 constraint: along propagation paths at cosmological distances
# The relevant epsilon_B is at the source environment
# For a rough bound: use the solar system value

# Detection threshold
H_target = 1e-18  # fractional holonomy target


def solve_and_get_profile(eta, lambda_val, source, u_ambient=0.0):
    """Solve BVP and return the profile + observables."""
    pot = DimlessCosh(eta)
    u_eq = pot.equilibrium(lambda_val, 1.0)

    x_init = np.logspace(-4, np.log10(source.x[-1]), 300)
    u_guess = np.where(x_init <= 1.0, u_eq, u_ambient + (u_eq - u_ambient) * np.exp(-(x_init - 1.0)))

    x, u, u_p, sol = solve_bvp_dimless(pot, source, u_ambient, tol=1e-10,
                                       u_init_arr=u_guess, max_nodes=100000)
    if not sol.success:
        return None

    obs = extract_observables_dimless(x, u, u_p, source, pot, u_ambient, lambda_val)
    return x, u, u_p, obs


def compute_triangle_holonomy(x, eps_B_shape, x_ground, x_sat):
    """Compute the triangle holonomy shape factor.

    H_resid / B0 = integral of eps_B_shape along the radial leg
    (from ground to satellite and back, but the holonomy comes from
    the non-exact part — for a radial path the conformal part cancels).

    The disformal synchronization correction is:
      delta_sigma ~ (B/A^2) * (u.nabla_phi) * P * nabla_phi

    For a radial field and radial motion: this reduces to
      H_resid ~ B0 * integral_{x_ground}^{x_sat} eps_B_shape(x) dx

    (The exact geometric factors depend on the triangle orientation,
    but this gives the order of magnitude for the radial leg contribution.)
    """
    mask = (x >= x_ground) & (x <= x_sat)
    if np.sum(mask) < 2:
        return 0.0
    # Integrate eps_B_shape along the radial leg
    H_shape = trapezoid(eps_B_shape[mask], x[mask])
    return H_shape


def compute_cassini_bound(x, eps_B_shape, x_cassini):
    """Compute epsilon_B at the Cassini impact parameter.

    Returns eps_B / B0 at x_cassini.
    """
    idx = np.argmin(np.abs(x - x_cassini))
    if idx == 0 or idx == len(x) - 1:
        # Interpolate or extrapolate
        return np.interp(x_cassini, x, eps_B_shape)
    return eps_B_shape[idx]


# ---------------------------------------------------------------------------
# Run the scan
# ---------------------------------------------------------------------------

print("=" * 70)
print("B0 WINDOW SCAN")
print("=" * 70)
print()
print("For each eta, compute:")
print("  1. eps_B/B0 at Cassini impact parameter (Sun)")
print("  2. H_resid/B0 for the triangle radial leg (Earth)")
print("  3. B0 range allowed by Cassini: B0 < 10^-15 / eps_Cassini")
print("  4. B0 range detectable by triangle: B0 > 10^-18 / H_shape")
print("  5. Window exists if Cassini_max > Triangle_min")
print()

# Source profiles
source_earth = make_layered_sphere_dimless(
    layers=[(0.55, 13.0/5.51), (0.90, 4.5/5.51), (1.0, 2.7/5.51)],
    s_ambient=0.0, x_max=50.0
)
source_earth.lambda_param = lambda_earth

source_sun = make_layered_sphere_dimless(
    layers=[(0.2, 150.0/1.41), (1.0, 1.0)],
    s_ambient=0.0, x_max=np.log10(x_cassini_sun) + 1  # extend to Cassini
)
source_sun.x = np.logspace(-4, np.log10(x_cassini_sun * 2), 800)
source_sun.s = np.where(source_sun.x <= 0.2, 150.0/1.41,
                np.where(source_sun.x <= 1.0, 1.0, 0.0))
source_sun.lambda_param = lambda_sun

results = []

for log_eta in [-1, 0, 1, 2]:
    eta_earth = 10.0**log_eta
    eta_sun = eta_earth * ratio_R**2

    print(f"\n--- log(eta_earth) = {log_eta} ---")

    # Solve Earth
    result_earth = solve_and_get_profile(eta_earth, lambda_earth, source_earth)
    if result_earth is None:
        print("  Earth: SOLVER FAILED")
        continue
    x_e, u_e, u_p_e, obs_e = result_earth

    # Solve Sun (extend to Cassini distance)
    result_sun = solve_and_get_profile(eta_sun, lambda_sun, source_sun)
    if result_sun is None:
        print("  Sun: SOLVER FAILED")
        continue
    x_s, u_s, u_p_s, obs_s = result_sun

    # eps_B shape at Cassini impact parameter (Sun)
    eps_cassini = compute_cassini_bound(x_s, obs_s['eps_B_shape'], x_cassini_sun)

    # Triangle holonomy shape (Earth, ground to satellite)
    H_shape = compute_triangle_holonomy(x_e, obs_e['eps_B_shape'], x_ground, x_sat)

    # B0 constraints
    cassini_limit = 3e-15  # few x 10^-15
    B0_cassini_max = cassini_limit / max(abs(eps_cassini), 1e-300)
    B0_triangle_min = H_target / max(abs(H_shape), 1e-300)

    # Window exists?
    window_exists = B0_cassini_max > B0_triangle_min
    window_ratio = B0_cassini_max / B0_triangle_min if B0_triangle_min > 0 else 0

    print(f"  eps_B/B0 at Cassini (Sun):  {eps_cassini:.4e}")
    print(f"  H_resid/B0 (triangle leg):  {H_shape:.4e}")
    print(f"  B0 max (Cassini):           {B0_cassini_max:.4e}")
    print(f"  B0 min (triangle detect):   {B0_triangle_min:.4e}")
    print(f"  Window ratio (max/min):     {window_ratio:.4e}")
    print(f"  Window exists:              {window_exists}")

    if window_exists:
        # Pick B0 in the middle of the window (log scale)
        B0_optimal = np.sqrt(B0_cassini_max * B0_triangle_min)
        H_predicted = B0_optimal * H_shape
        eps_predicted = B0_optimal * eps_cassini
        print(f"  Optimal B0:                 {B0_optimal:.4e}")
        print(f"  Predicted H_resid:          {H_predicted:.4e}")
        print(f"  Predicted eps_Cassini:      {eps_predicted:.4e}")

    results.append({
        'log_eta': log_eta,
        'eps_cassini': eps_cassini,
        'H_shape': H_shape,
        'B0_cassini_max': B0_cassini_max,
        'B0_triangle_min': B0_triangle_min,
        'window_exists': window_exists,
        'window_ratio': window_ratio,
    })

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

print()
print("=" * 70)
print("B0 WINDOW SCAN SUMMARY")
print("=" * 70)
print()
print(f"{'log_eta':>8s}  {'eps_Cassini':>12s}  {'H_shape':>12s}  "
      f"{'B0_Cassini':>12s}  {'B0_triangle':>12s}  {'window':>8s}")
print("-" * 80)
for r in results:
    w = "YES" if r['window_exists'] else "NO"
    print(f"{r['log_eta']:8d}  {r['eps_cassini']:12.4e}  {r['H_shape']:12.4e}  "
          f"{r['B0_cassini_max']:12.4e}  {r['B0_triangle_min']:12.4e}  {w:>8s}")

print()
any_window = any(r['window_exists'] for r in results)
if any_window:
    print("WINDOW EXISTS for at least one eta value.")
    print("The cosh potential with K=1 permits a detectable triangle signal")
    print("while satisfying Cassini/GW170817 constraints.")
else:
    print("No window exists for the tested eta values.")
    print("Either:")
    print("  1. The triangle geometry needs optimization (longer legs, higher altitude)")
    print("  2. The B0 scan range needs extension")
    print("  3. The minimal K=1 cosh closure is insufficient for detectability")
    print()
    print("Note: the eps_B shape is extremely small because the field excursion")
    print("is tiny (u ~ 10^-9). The disformal deformation scales as u^2 * u'^2,")
    print("which is ~10^-36. This is the fundamental challenge for the disformal")
    print("sector in the weak-field regime.")
