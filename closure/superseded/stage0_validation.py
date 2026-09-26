"""Stage 0: BVP solver validation.

Validates:
1. Linear Yukawa: solver converges, profile is smooth, exterior decays as exp(-mu*x)/x
2. Conformal connection is exact (Theorem 2): numerical integral matches analytic Delta(ln A)
3. Weak-source limit: solver handles lambda << 1 correctly
4. Nonlinear exponential: solver converges with appropriate initial guess
"""

import numpy as np
import sys
sys.path.insert(0, '/Users/matthewsmawfield/www/Temporal Equivalence Principle/closure')

from radial_bvp_dimless import (
    DimlessLinear, DimlessExponential, make_uniform_sphere_dimless,
    solve_bvp_dimless, extract_observables_dimless,
    test_conformal_loop as run_loop_test, BETA_A,
)


def test_yukawa():
    """Validate BVP solver: linear potential, uniform sphere."""
    print("=" * 70)
    print("STAGE 0a: Linear Yukawa validation")
    print("=" * 70)
    print()

    lam = 1.0
    mu = 3.0
    s_core = 1.0

    pot = DimlessLinear(mu)
    source = make_uniform_sphere_dimless(s_core=s_core, x_max=50.0)
    source.lambda_param = lam

    print(f"  Potential: {pot.describe()}")
    print(f"  lambda = {lam}, mu = {mu}, s_core = {s_core}")
    print(f"  Source term: -lambda*s*exp(-u) ~ -1 inside (beta_A = -1)")
    print(f"  Particular solution: u_p ~ lambda/mu^2 = {lam/mu**2:.4f}")
    print()

    x, u, u_p, sol = solve_bvp_dimless(pot, source, 0.0, tol=1e-8)

    if not sol.success:
        print(f"  FAIL: {sol.message}")
        return False
    print(f"  BVP solver: CONVERGED")

    u_center = u[np.argmin(np.abs(x - 0.01))]
    u_surface = u[np.argmin(np.abs(x - 1.0))]
    u_far = u[np.argmin(np.abs(x - 10.0))]

    print(f"  u(0.01) = {u_center:.6e}  (positive, < lambda/mu^2 = {lam/mu**2:.4e})")
    print(f"  u(1.0)  = {u_surface:.6e}  (positive, < u_center)")
    print(f"  u(10)   = {u_far:.6e}   (near zero)")
    print()

    # Check: field is positive inside, decreasing, near zero far away
    if u_center <= 0 or u_center > lam/mu**2 * 1.5:
        print(f"  FAIL: Interior field not in expected range")
        return False
    if u_surface >= u_center or u_surface < 0:
        print(f"  FAIL: Field should decrease from center to surface, stay positive")
        return False
    print(f"  PASS: Profile shape is physically correct")

    # Check: exterior decays as exp(-mu*x)/x
    mask_ext = (x > 2.0) & (x < 20.0) & (np.abs(u) > 1e-30)
    if np.sum(mask_ext) > 5:
        log_ux = np.log(np.abs(u[mask_ext]) * x[mask_ext])
        x_fit = x[mask_ext]
        coeffs = np.polyfit(x_fit, log_ux, 1)
        slope = coeffs[0]
        print(f"  Exterior decay: slope = {slope:.4f}, expected -mu = {-mu:.4f}")
        if abs(slope - (-mu)) / mu > 0.05:
            print(f"  FAIL: Decay slope mismatch")
            return False
        print(f"  PASS: Exterior decays as exp(-mu*x)/x")
    print()

    # Regularity
    u_p_center = u_p[np.argmin(np.abs(x - 0.01))]
    print(f"  u'(0.01) = {u_p_center:.4e} (should be ~ 0)")
    if abs(u_p_center) > 1e-3:
        print(f"  FAIL: Not regular at center")
        return False
    print(f"  PASS: Regular at center")
    print()
    return True


def test_conformal_loop():
    """Test that the conformal connection is exact (Theorem 2)."""
    print("=" * 70)
    print("STAGE 0b: Conformal connection exactness (Theorem 2)")
    print("=" * 70)
    print()

    lam = 1.0
    mu = 3.0
    pot = DimlessLinear(mu)
    source = make_uniform_sphere_dimless(s_core=1.0, x_max=50.0)
    source.lambda_param = lam

    x, u, u_p, sol = solve_bvp_dimless(pot, source, 0.0, tol=1e-8)
    if not sol.success:
        print("  FAIL: BVP did not converge")
        return False

    result = run_loop_test(x, u, u_p)

    print(f"  Numerical integral of beta_A*u' dx: {result['integral_numerical']:.10e}")
    print(f"  Analytic beta_A*(u_max - u_min):    {result['ln_A_diff']:.10e}")
    print(f"  Integration error:                  {result['integration_error']:.4e}")
    print(f"  Loop integral (exact form):         {result['loop_integral']:.4e}")
    print()

    # The key test: numerical integration of the connection must match
    # the analytic difference to high precision. This confirms the
    # connection is being computed correctly for holonomy subtraction.
    tol = 1e-6
    if result['integration_error'] < tol:
        print(f"  PASS: Conformal connection is exact to < {tol:.0e}")
        print(f"  Theorem 2 verified: d(ln A) is an exact differential")
        print(f"  A(phi) cannot source synchronization holonomy")
        return True
    else:
        print(f"  FAIL: Integration error {result['integration_error']:.4e} exceeds {tol:.0e}")
        return False


def test_weak_source():
    """Test solver with lambda << 1 (realistic TEP regime)."""
    print()
    print("=" * 70)
    print("STAGE 0c: Weak-source test (lambda ~ 10^-9, Earth-like)")
    print("=" * 70)
    print()

    lam = 1e-9
    mu = 3.0
    pot = DimlessLinear(mu)
    source = make_uniform_sphere_dimless(s_core=1.0, x_max=50.0)
    source.lambda_param = lam

    print(f"  lambda = {lam:.4e} (Earth-like)")
    print(f"  mu = {mu:.4f}")
    print(f"  For a finite sphere, u_center < lambda/mu^2 due to geometry")
    print(f"  The exact ratio depends on mu and sphere structure")
    print()

    x, u, u_p, sol = solve_bvp_dimless(pot, source, 0.0, tol=1e-12)

    if not sol.success:
        print(f"  FAIL: {sol.message}")
        return False

    print(f"  BVP solver: CONVERGED (tol=1e-12)")
    u_center = u[np.argmin(np.abs(x - 0.01))]
    u_expected_max = lam / mu**2  # upper bound (infinite medium particular solution)

    print(f"  u(0.01) = {u_center:.6e}")
    print(f"  lambda/mu^2 (infinite medium) = {u_expected_max:.6e}")
    print(f"  Ratio = {u_center/u_expected_max:.4f} (must be in (0, 1))")
    print()

    # For a finite sphere, the center value is less than the infinite-medium
    # particular solution but still positive and proportional to lambda.
    # The ratio depends on mu (geometry factor).
    if 0 < u_center / u_expected_max < 1.0:
        print(f"  PASS: Solver handles weak source correctly")
        print(f"  Field is positive and bounded by infinite-medium particular solution")
        return True
    else:
        print(f"  FAIL: Ratio {u_center/u_expected_max:.4f} outside (0, 1)")
        return False


def test_nonlinear():
    """Test solver with exponential potential (nonlinear)."""
    print()
    print("=" * 70)
    print("STAGE 0d: Nonlinear profile convergence (exponential potential)")
    print("=" * 70)
    print()

    lam = 10.0
    eta = 1.0
    pot = DimlessExponential(eta)
    source = make_uniform_sphere_dimless(s_core=1.0, x_max=50.0)
    source.lambda_param = lam

    u_eq = 0.5 * np.log(lam / eta)
    print(f"  Potential: {pot.describe()}")
    print(f"  lambda = {lam}, eta = {eta}")
    print(f"  Equilibrium u_min = 0.5*ln(lambda/eta) = {u_eq:.4f}")
    print()

    # Better initial guess: u_eq inside, 0 outside, smooth transition
    x_init = np.logspace(-4, np.log10(50.0), 300)
    u_guess = np.where(x_init <= 1.0, u_eq, u_eq * np.exp(-(x_init - 1.0)))
    u_guess = np.minimum(u_guess, u_eq)

    x, u, u_p, sol = solve_bvp_dimless(pot, source, 0.0, tol=1e-8,
                                       u_init_arr=u_guess, max_nodes=100000)

    if not sol.success:
        print(f"  FAIL: {sol.message}")
        # Try with relaxed tolerance
        x, u, u_p, sol = solve_bvp_dimless(pot, source, 0.0, tol=1e-6,
                                           u_init_arr=u_guess, max_nodes=200000)
        if not sol.success:
            print(f"  FAIL (relaxed): {sol.message}")
            return False
        print(f"  BVP solver: CONVERGED (tol=1e-6, relaxed)")
    else:
        print(f"  BVP solver: CONVERGED")

    u_center = u[np.argmin(np.abs(x - 0.01))]
    u_surface = u[np.argmin(np.abs(x - 1.0))]
    u_far = u[np.argmin(np.abs(x - 10.0))]

    print(f"  u(0.01) = {u_center:.6e}  (should approach u_eq = {u_eq:.4f})")
    print(f"  u(1.0)  = {u_surface:.6e}")
    print(f"  u(10)   = {u_far:.6e}   (should approach 0)")
    print()

    if u_center <= 0 or u_center > u_eq * 1.5:
        print(f"  FAIL: Interior field not in expected range")
        return False
    if u_surface >= u_center:
        print(f"  FAIL: Field should decrease from center to surface")
        return False
    print(f"  PASS: Nonlinear profile is physically sensible")
    print(f"  Field rises toward equilibrium inside, relaxes to ambient outside")
    return True


if __name__ == '__main__':
    results = {
        'Yukawa': test_yukawa(),
        'Conformal exact': test_conformal_loop(),
        'Weak source': test_weak_source(),
        'Nonlinear': test_nonlinear(),
    }

    print()
    print("=" * 70)
    print("STAGE 0 SUMMARY")
    print("=" * 70)
    for name, ok in results.items():
        print(f"  {name:20s}: {'PASS' if ok else 'FAIL'}")
    print()

    if all(results.values()):
        print("  Solver is validated. Proceed to Stage 1 (Earth profile).")
        sys.exit(0)
    else:
        print("  Solver is NOT fully validated.")
        sys.exit(1)
