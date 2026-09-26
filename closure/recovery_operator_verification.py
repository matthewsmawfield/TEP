"""Recovery-operator Lagrangian verification.

The TEP action is P(X, phi) = K*X + mu^2(phi)*sqrt(2X) - V(phi).
The recovery operator is the non-perturbative solution of the field
equation for a point mass.  The cuscuton term mu^2(phi)*sqrt(2X) has a
non-analytic |grad phi| that produces a non-smooth transition in the
field profile — the steep recovery (k >= 4).

This script solves the field equation numerically for a point mass and
fits the recovery function f(r) = [1 + (R/r)^k]^{-1} to verify k >= 4.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize

# Cosmological parameters (units M_Pl = H_0 = 1)
Om = 0.3153
Or = 9.1e-5
OL = 1 - Om - Or
K = 1.0  # K > 3*Om = 0.946, use K = 1
beta_A = -1.0

# Cosmological branch: V(phi) and mu^2(phi)
def V(phi):
    return -(3*Om*np.exp(-phi) + 3*Or)

def V_prime(phi):
    return 3*Om*np.exp(-phi)

def mu2(phi):
    """mu^2(phi) = -(3*Om*e^{-phi} + 4*Or) / (e^{-phi} * H(z))"""
    A = np.exp(-phi)
    z = np.exp(phi) - 1
    H = np.sqrt(Om*(1+z)**3 + Or*(1+z)**4 + OL)
    return -(3*Om*A + 4*Or) / (A * H)

# Field equation for a point mass M at the origin:
# K * (1/r^2) d/dr(r^2 dphi/dr) + 2*mu^2(phi)/r - V'(phi) = beta_A * M * delta(r) / (4*pi*M_Pl)
#
# Outside the source (r > 0):
# K * (1/r^2) d/dr(r^2 dphi/dr) + 2*mu^2(phi)/r - V'(phi) = 0
#
# Let u = r^2 * dphi/dr, then:
# dphi/dr = u / r^2
# du/dr = (-2*mu^2(phi)*r + V'(phi)*r^2) / K
#
# The Newtonian solution (mu2=0, V'=0) gives u_N = beta_A * M / (4*pi*K)
# The recovery function is f(r) = u(r) / u_N

# Source strength (in units where M_Pl = H_0 = 1)
# For a solar-mass object: M ~ M_sun / M_Pl^2 ~ 10^-38 (too small for numerics)
# Use a dimensionless source strength and scale
M_source = 1.0  # dimensionless
u_Newton = beta_A * M_source / (4 * np.pi * K)  # Newtonian u

def field_eq(r, y):
    """Field equation: y = [phi, u] where u = r^2 * dphi/dr."""
    phi, u = y
    if r < 1e-10:
        return [0, 0]
    dphi = u / r**2
    du = (-2 * mu2(phi) * r + V_prime(phi) * r**2) / K
    return [dphi, du]

# Solve from small r to large r
# At small r, the field is screened: dphi/dr -> 0, so u -> 0
# But we need a non-zero initial condition to avoid trivial solution
# Use u(r_min) = epsilon (small but non-zero)
r_min = 1e-3
r_max = 1e3
phi0 = 0.0  # start at cosmological value
u0 = 1e-6  # small initial gradient (screened)

sol = solve_ivp(field_eq, [r_min, r_max], [phi0, u0],
                method='RK45', max_step=0.1, rtol=1e-8, atol=1e-10,
                dense_output=True)

if sol.success:
    r_arr = np.logspace(np.log10(r_min), np.log10(r_max), 1000)
    phi_arr = sol.sol(r_arr)[0]
    u_arr = sol.sol(r_arr)[1]
    dphi_arr = u_arr / r_arr**2

    # Recovery function: f(r) = u(r) / u_Newton
    f_arr = u_arr / u_Newton

    # Fit f(r) = [1 + (R/r)^k]^{-1}
    # Invert: 1/f - 1 = (R/r)^k => log(1/f - 1) = k*log(R) - k*log(r)
    # Linear fit: log(1/f - 1) = a + b*log(r) where b = -k, a = k*log(R)
    valid = (f_arr > 0.01) & (f_arr < 0.99) & np.isfinite(f_arr)
    if np.sum(valid) > 10:
        log_r = np.log(r_arr[valid])
        log_inv = np.log(1.0 / f_arr[valid] - 1.0)
        # Linear regression
        A_mat = np.column_stack([np.ones(len(log_r)), log_r])
        coef, res, _, _ = np.linalg.lstsq(A_mat, log_inv, rcond=None)
        a_fit, b_fit = coef
        k_fit = -b_fit
        R_fit = np.exp(a_fit / k_fit) if k_fit > 0 else np.nan

        print("=" * 60)
        print("RECOVERY OPERATOR VERIFICATION")
        print("=" * 60)
        print(f"TEP action: P(X, phi) = K*X + mu^2(phi)*sqrt(2X) - V(phi)")
        print(f"  K = {K}, beta_A = {beta_A}")
        print(f"  mu^2(0) = {mu2(0):.4f}, V'(0) = {V_prime(0):.4f}")
        print()
        print(f"Field equation solved for r in [{r_min}, {r_max}]")
        print(f"  Newtonian u = {u_Newton:.6f}")
        print()
        print(f"Recovery function fit: f(r) = [1 + (R/r)^k]^(-1)")
        print(f"  k = {k_fit:.2f}")
        print(f"  R = {R_fit:.4f}")
        print(f"  F4 requirement: k >= 4  =>  {'PASS' if k_fit >= 4 else 'FAIL'}")
        print()

        # Check environmental ordering (F5)
        # Solve with different mu2 (representing different environments)
        # Denser environment -> larger mu2 -> should give larger R
        for env_label, mu2_factor in [("low_density", 0.5), ("midplane", 1.0), ("high_density", 2.0)]:
            def field_eq_env(r, y, mf=mu2_factor):
                phi, u = y
                if r < 1e-10:
                    return [0, 0]
                dphi = u / r**2
                du = (-2 * mu2(phi) * mf * r + V_prime(phi) * r**2) / K
                return [dphi, du]

            sol_env = solve_ivp(field_eq_env, [r_min, r_max], [phi0, u0],
                               method='RK45', max_step=0.1, rtol=1e-8, atol=1e-10,
                               dense_output=True)
            if sol_env.success:
                r_e = np.logspace(np.log10(r_min), np.log10(r_max), 500)
                u_e = sol_env.sol(r_e)[1]
                f_e = u_e / u_Newton
                v = (f_e > 0.01) & (f_e < 0.99) & np.isfinite(f_e)
                if np.sum(v) > 5:
                    lr = np.log(r_e[v])
                    li = np.log(1.0 / f_e[v] - 1.0)
                    A2 = np.column_stack([np.ones(len(lr)), lr])
                    c2 = np.linalg.lstsq(A2, li, rcond=None)[0]
                    k_e = -c2[1]
                    R_e = np.exp(c2[0] / k_e) if k_e > 0 else np.nan
                    print(f"  {env_label:12s} (mu2 x{mu2_factor:.1f}): k={k_e:.2f}, R={R_e:.4f}")

        print()
        print("F5 check: denser environment (larger mu2) should give LARGER R")
        print("  => PASS if R increases with mu2_factor")

        # Save results
        import json
        results = {
            "action": "P(X, phi) = K*X + mu^2(phi)*sqrt(2X) - V(phi)",
            "K": K, "beta_A": beta_A,
            "mu2_at_0": float(mu2(0)), "V_prime_at_0": float(V_prime(0)),
            "k_fit": float(k_fit), "R_fit": float(R_fit),
            "F4_pass": bool(k_fit >= 4),
            "u_Newton": float(u_Newton),
            "note": "Recovery operator is the non-perturbative solution of the TEP field equation. The cuscuton term mu^2(phi)*sqrt(2X) produces the steep recovery through its non-analytic |grad phi|."
        }
        with open("results/recovery_operator_verification.json", "w") as f:
            json.dump(results, f, indent=2)
        print("Results saved to results/recovery_operator_verification.json")
    else:
        print(f"Not enough valid points for fit (only {np.sum(valid)} valid)")
        print(f"f range: [{np.nanmin(f_arr):.4e}, {np.nanmax(f_arr):.4e}]")
else:
    print(f"Integration failed: {sol.message}")
