#!/usr/bin/env python3
"""Step 24: Compact-body charge closure for TEP binary pulsar constraints.

This script determines the effective scalar charges of a Neutron Star (NS)
and a White Dwarf (WD) using the exact frozen TEP static profile solver
(canonical kinetic term, unified potential with quartic completion).
It then evaluates whether the predicted charge difference is consistent
with the PSR J1738+0333 dipole radiation bound, without refitting.
"""

import sys
import os
import json
import numpy as np

# Add the directory containing tep_model.py to the path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tep_model import (M_SUN, R_SUN, LAMBDA_REFERENCE, solve_sphere,
                       diagnostics, save)

# PSR J1738+0333 parameters (Antoniadis et al. 2012)
M_NS = 1.46 * M_SUN
R_NS = 12000.0  # 12 km approx
M_WD = 0.181 * M_SUN
R_WD = 10000000.0  # 10,000 km approx

# Empirical bound on dipole radiation from J1738+0333
# Freire et al. (2012) gives a bound on the effective scalar charge difference:
# |alpha_NS - alpha_WD| < 2.0e-3 (roughly, depending on the exact ST theory).
# We use this as the order-of-magnitude constraint on the differential charge.
CHARGE_DIFF_BOUND = 2.0e-3

def run():
    print("=========================================================")
    print(" Step 24: Compact-Body Charge Closure (PSR J1738+0333)")
    print("=========================================================")
    print("\n[A] Actual Microscopic Action")
    print("Kinetic term: Canonical -1/2 (\nabla \phi)^2 (frozen)")
    print(f"Potential: Unified Potential with \lambda={LAMBDA_REFERENCE} (frozen)")
    print("Bare coupling: \beta_A = -1 (frozen)")

    results = {
        "action": {
            "kinetic": "canonical",
            "potential": "unified_quartic_plus_plateau",
            "lambda": LAMBDA_REFERENCE,
            "beta_A": -1.0
        },
        "bodies": {}
    }

    print("\n[B] Solving Static Profiles and Extracting Charges")
    
    # 1. Solve Neutron Star
    print(f"Solving Neutron Star (M={M_NS/M_SUN:.2f} M_sun, R={R_NS/1000:.1f} km)...")
    try:
        sol_ns = solve_sphere(M_NS, R_NS, LAMBDA_REFERENCE, x_max=1e6)
        diag_ns = diagnostics(sol_ns, 1e5) # far field charge
        alpha_ns = diag_ns['source_charge_ratio'] * (-1.0) # alpha = beta_A * S = -1 * S
        print(f"  -> Converged! Effective charge \alpha_NS = {alpha_ns:.4e}")
        results["bodies"]["NS"] = {
            "mass_Msun": M_NS/M_SUN,
            "radius_m": R_NS,
            "S_charge_ratio": diag_ns['source_charge_ratio'],
            "alpha_eff": alpha_ns
        }
    except Exception as e:
        print(f"  -> Solver failed for NS: {e}")
        results["bodies"]["NS"] = {"error": str(e)}
        alpha_ns = None

    # 2. Solve White Dwarf
    print(f"Solving White Dwarf (M={M_WD/M_SUN:.3f} M_sun, R={R_WD/1000:.1f} km)...")
    try:
        sol_wd = solve_sphere(M_WD, R_WD, LAMBDA_REFERENCE, x_max=1e6)
        diag_wd = diagnostics(sol_wd, 1e5)
        alpha_wd = diag_wd['source_charge_ratio'] * (-1.0)
        print(f"  -> Converged! Effective charge \alpha_WD = {alpha_wd:.4e}")
        results["bodies"]["WD"] = {
            "mass_Msun": M_WD/M_SUN,
            "radius_m": R_WD,
            "S_charge_ratio": diag_wd['source_charge_ratio'],
            "alpha_eff": alpha_wd
        }
    except Exception as e:
        print(f"  -> Solver failed for WD: {e}")
        results["bodies"]["WD"] = {"error": str(e)}
        alpha_wd = None

    print("\n[C] Pulsar Measurement Requirement")
    if alpha_ns is not None and alpha_wd is not None:
        diff = abs(alpha_ns - alpha_wd)
        print(f"Predicted absolute charge difference: |\alpha_NS - \alpha_WD| = {diff:.4e}")
        print(f"Observational strict bound limit: < {CHARGE_DIFF_BOUND:.1e}")
        
        passes = diff < CHARGE_DIFF_BOUND
        results["constraint"] = {
            "predicted_diff": diff,
            "observational_bound": CHARGE_DIFF_BOUND,
            "passes_without_refitting": passes
        }
        
        if passes:
            print("=> RESULT: The frozen TEP screening architecture NATURALLY SUPPRESSES differential charges below the pulsar constraint.")
        else:
            print("=> RESULT: The predicted differential charge EXCEEDS the pulsar constraint. The screening mechanism is insufficient without modification.")
    else:
        print("=> RESULT: Could not compute differential charge due to solver failure. This indicates the strong-field profiles require a modified solution method or the unified potential breaks down for these compactnesses.")
        results["constraint"] = {"error": "Missing charge data"}

    # Save output
    outdir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "step_24_compact_body_charge_closure.json"), "w") as f:
        json.dump(results, f, indent=2)

if __name__ == "__main__":
    run()
