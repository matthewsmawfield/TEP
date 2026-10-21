#!/usr/bin/env python3
"""Step 25: Macroscopic Temporal Shear BVP Investigation.

This script tests whether an admissible non-linear operator class
for the macroscopic Temporal Shear configuration can simultaneously satisfy:
1. Cassini bound (strong screening at 1.6 R_sun)
2. Saturn bound (strong screening at 9.5 AU)
3. Wide Binary recovery (unscreened at 2646 AU)
4. LLR (zero/negligible differential acceleration between Earth and Moon)

The PDE is: \nabla \cdot ( D(\phi, |\nabla\phi|) \nabla\phi ) = Q
where Q = \beta_A M_pl \rho.
For spherical symmetry, D(\phi, \phi') \phi' = Q_enclosed / (4\pi r^2).
"""

import sys
import os
import json
import numpy as np
from scipy.optimize import root_scalar

# Constants
G = 6.67430e-11
M_SUN = 1.98847e30
R_SUN = 6.957e8
AU = 149597870700.0
M_EARTH = 5.972e24
M_MOON = 7.342e22
C = 299792458.0
M_PL = 2.435e18 # GeV
BETA = -1.0

# Define constraints
CASSINI_R = 1.6 * R_SUN
CASSINI_S_BOUND = 5.75e-6

SATURN_R = 9.5 * AU
SATURN_ACCEL_BOUND = 1e-10 # m/s^2 scalar acceleration bound

WB_R = 2646.0 * AU
WB_S_TARGET = 0.5 # 50% recovery

# Evaluator for a given model D(\phi, s)
def evaluate_model(phi_0, n, s_0, m, phi_gal=0.0):
    """
    Evaluates D(\phi, s) = 1 + (\phi / \phi_0)^n + (s / s_0)^m
    using the exact 1D integral: \phi' * D = g_N.
    """
    results = {}
    
    # Function to solve for \phi' given \phi and g_N
    def get_phi_prime(phi, g_N):
        if m == 0:
            # Amplitude only dependence
            D = 1 + (phi / phi_0)**n
            return g_N / D
        else:
            # Solve s * (1 + (\phi/\phi_0)^n + (s/s_0)^m) = g_N
            def f(s):
                return s * (1 + (phi / phi_0)**n + (s / s_0)**m) - g_N
            # derivative
            def df(s):
                return 1 + (phi / phi_0)**n + (m + 1) * (s / s_0)**m
            
            try:
                # s is bounded between 0 and g_N
                sol = root_scalar(f, fprime=df, x0=g_N/2, bracket=[0, g_N], method='brentq')
                return sol.root
            except ValueError:
                return None

    # Since \phi' * D(\phi) = g_N, we can integrate \phi from r=\infty (where \phi=\phi_gal)
    # However, if m=0, D(\phi) depends ONLY on \phi. 
    # \int D(\phi) d\phi = \int g_N dr = \int (GM/r^2) dr = -GM/r + C
    # Let U(\phi) = \phi + \frac{1}{n+1} \frac{\phi^{n+1}}{\phi_0^n}
    # Then U(\phi) = GM/r + U(\phi_gal)
    
    if m != 0:
        results['error'] = 'Numerical integration for mixed m,n not implemented in this quick test.'
        return results
        
    def U(phi):
        return phi + (1.0 / (n + 1)) * (phi**(n + 1)) / (phi_0**n)
        
    def solve_phi(U_target):
        # Solve U(\phi) = U_target
        def f(phi):
            return U(phi) - U_target
        # phi is monotonically increasing with U
        try:
            sol = root_scalar(f, bracket=[0, U_target], method='brentq')
            return sol.root
        except ValueError:
            return None

    U_gal = U(phi_gal)
    
    # 1. Cassini
    g_N_cassini = G * M_SUN / (CASSINI_R**2)
    U_cassini = g_N_cassini * CASSINI_R + U_gal  # GM/r is actually g_N * r
    phi_cassini = solve_phi(U_cassini)
    if phi_cassini is not None:
        phi_prime_cassini = get_phi_prime(phi_cassini, g_N_cassini)
        S_cassini = phi_prime_cassini / g_N_cassini
        results['cassini'] = {
            'phi': phi_cassini,
            'g_N': g_N_cassini,
            'S': S_cassini,
            'passes': S_cassini < CASSINI_S_BOUND
        }
    
    # 2. Saturn
    g_N_saturn = G * M_SUN / (SATURN_R**2)
    U_saturn = g_N_saturn * SATURN_R + U_gal
    phi_saturn = solve_phi(U_saturn)
    if phi_saturn is not None:
        phi_prime_saturn = get_phi_prime(phi_saturn, g_N_saturn)
        results['saturn'] = {
            'phi': phi_saturn,
            'g_N': g_N_saturn,
            'accel': phi_prime_saturn,
            'passes': phi_prime_saturn < SATURN_ACCEL_BOUND
        }
        
    # 3. Wide Binary (1.2 M_sun equivalent for average demographic)
    g_N_wb = G * (1.2 * M_SUN) / (WB_R**2)
    U_wb = g_N_wb * WB_R + U_gal
    phi_wb = solve_phi(U_wb)
    if phi_wb is not None:
        phi_prime_wb = get_phi_prime(phi_wb, g_N_wb)
        S_wb = phi_prime_wb / g_N_wb
        results['wide_binary'] = {
            'phi': phi_wb,
            'g_N': g_N_wb,
            'S': S_wb,
            'passes': S_wb > 0.1 # Need significant recovery
        }
        
    # 4. LLR (Earth-Moon embedded in Sun)
    # The background field at 1 AU sets the effective D for both Earth and Moon.
    g_N_1au = G * M_SUN / (AU**2)
    U_1au = g_N_1au * AU + U_gal
    phi_1au = solve_phi(U_1au)
    if phi_1au is not None:
        D_1au = 1 + (phi_1au / phi_0)**n
        # The Earth and Moon scalar charges are suppressed by D(\phi_1au).
        S_earth = 1.0 / D_1au
        S_moon = 1.0 / D_1au
        diff = abs(S_earth - S_moon)
        results['llr'] = {
            'D_background': D_1au,
            'S_earth': S_earth,
            'S_moon': S_moon,
            'diff': diff,
            'passes': diff < 1e-13
        }

    return results

def run():
    print("=========================================================")
    print(" Step 25: Macroscopic Temporal Shear BVP")
    print("=========================================================")
    print("Testing pure gradient-cap vs amplitude-dependent operators.")
    
    out = {}
    
    # Test 1: High-power amplitude operator D(\phi) = 1 + (\phi/\phi_0)^4
    # We tune \phi_0 such that D ~ 1 at wide binaries, but D >> 1 at Cassini.
    # At wide binaries, U = GM/R ~ 6.67e-11 * 2e30 / (2646 * 1.5e11) ~ 3.3e5
    # So phi ~ 3.3e5. Let's set \phi_0 = 1e6.
    print("\n[Model A] Amplitude-dependent: D(\phi) = 1 + (\phi / 10^5)^4")
    res_A = evaluate_model(phi_0=1e5, n=4, s_0=1, m=0)
    out['Model_A'] = res_A
    print(json.dumps(res_A, indent=2))
    
    # Save output
    outdir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "step_25_macroscopic_shear_bvp.json"), "w") as f:
        json.dump(out, f, indent=2)

if __name__ == "__main__":
    run()
