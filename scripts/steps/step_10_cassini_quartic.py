#!/usr/bin/env python3
"""Quartic-branch wide-binary test: candidate no-go demonstration.

Tests whether the quartic single-body potential baseline, tuned to the
solar screening regime, can simultaneously recover the wide-binary Temporal
Shear at 2646 AU. Retained as the explicit demonstration that single-body
potential completions cannot reconcile Cassini with wide-binary recovery;
the canonical screening law is the two-body kinetic operator (step_19).
"""
import numpy as np
from tep_model import solve_sphere, diagnostics, M_SUN, R_SUN, AU, save, KG_GEV, M_PL

def run():
    print("Testing wide binary at mu = 10^10...")
    
    mass = M_SUN
    radius = R_SUN
    mass_gev = mass * KG_GEV
    lam = 1e10 / (mass_gev/(4*np.pi*M_PL))**2
    
    # R_WB = 2646 AU
    R_WB_m = 2646 * AU
    x_WB = R_WB_m / R_SUN
    
    # We need to solve out to x > x_WB
    x_max = x_WB * 2
    
    try:
        # Note: running this out to x_max ~ 1e6 might be very slow or stiff.
        solved = solve_sphere(mass_kg=mass, radius_m=radius, lam=lam, rho_bg=1e-24, x_max=x_max, tol=1e-3)
        
        diag_WB = diagnostics(solved, x_WB)
        
        result = {
            'mu': solved['mu'],
            'x_WB': x_WB,
            'S_WB': diag_WB['source_charge_ratio'],
        }
        
        save('step_10_wb_quartic.json', result)
        print(f"S(R_WB) = {result['S_WB']:.2e}")
        
    except Exception as e:
        print(f"Solver failed: {e}")
    
if __name__ == '__main__':
    run()
