#!/usr/bin/env python3
"""Symbolic derivation of the Static Cosmology Consistency Test."""
import sympy as s
from tep_model import save

def run():
    # Define symbols
    t = s.symbols('t', real=True)
    a = s.Function('a')(t)  # Scale factor
    phi = s.Function('phi')(t) # Cosmological time field
    
    # Constants
    K = s.symbols('K', real=True, positive=True) # Kinetic term coefficient
    V = s.symbols('V', real=True, nonnegative=True) # Potential
    rho_m = s.symbols('rho_m', real=True, positive=True) # Matter density
    Lambda = s.symbols('Lambda', real=True) # Cosmological constant
    k = s.symbols('k', real=True) # Spatial curvature
    
    # Friedmann equations in a generic scalar-tensor theory / GR with scalar field
    # H^2 + k/a^2 = (1/3) * (rho_m + rho_phi + Lambda) (using 8 pi G = 1)
    # rho_phi = K/2 * (dphi/dt)^2 + V
    
    H = s.diff(a, t) / a
    rho_phi = K/2 * s.diff(phi, t)**2 + V
    
    friedmann_1 = s.Eq(3 * (H**2 + k/a**2), rho_m + rho_phi + Lambda)
    
    # Second Friedmann equation:
    # dH/dt + H^2 = - (1/6) * (rho_m + rho_phi + 3*(p_m + p_phi) - 2 Lambda)
    # p_m = 0 (dust)
    # p_phi = K/2 * (dphi/dt)^2 - V
    p_phi = K/2 * s.diff(phi, t)**2 - V
    p_m = 0
    
    friedmann_2 = s.Eq(3 * (s.diff(H, t) + H**2), -0.5 * (rho_m + rho_phi + 3*(p_m + p_phi) - 2*Lambda))
    
    # Condition for static cosmology: a(t) = a_0 (constant) -> H = 0, dH/dt = 0
    # Also phi(t) is non-trivial, so dphi/dt != 0
    # Evaluate equations for static universe
    f1_static = friedmann_1.subs({s.diff(a, t): 0})
    f2_static = friedmann_2.subs({s.diff(H, t): 0, s.diff(a, t): 0, H: 0})
    
    # f1_static: 3*k/a_0^2 = rho_m + K/2 * phidot^2 + V + Lambda
    # f2_static: 0 = -0.5 * (rho_m + 2*K*phidot^2 - 2*V - 2*Lambda)
    
    # If we want a flat (k=0) static universe with no cosmological constant (Lambda=0)
    f1_flat_nolambda = f1_static.subs({k: 0, Lambda: 0})
    f2_flat_nolambda = f2_static.subs({k: 0, Lambda: 0})
    
    result = {
        'friedmann_1': str(friedmann_1),
        'friedmann_2': str(friedmann_2),
        'static_eq1': str(f1_static),
        'static_eq2': str(f2_static),
        'flat_nolambda_eq1': str(f1_flat_nolambda),
        'flat_nolambda_eq2': str(f2_flat_nolambda),
        'conclusion': 'For k=0 and Lambda=0, Eq1 implies 0 = rho_m + K/2 * phidot^2 + V. Since rho_m > 0, K > 0, V >= 0, this has no real solution for phidot. Therefore, a flat static cosmology with a non-trivial time field REQUIRES either spatial curvature (k != 0) or a cosmological constant (Lambda < 0) or a phantom field (K < 0).'
    }
    
    print(result['conclusion'])
    save('step_10_cosmology_static.json', result)
    return result

if __name__ == '__main__':
    run()
