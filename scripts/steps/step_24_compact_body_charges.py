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
            "beta_A": -1.0,
            "lambda_branch_note": (
                "charges evaluated at the fiducial lambda_ref = 7.526e-71; "
                "at the operative Cassini branch lambda_Cassini = 1e5*lambda_ref "
                "the embedded charge response is ~10-300x smaller (step_67), so "
                "the reference-branch charges quoted here are conservative upper "
                "bounds on the compact-body scalar charges"
            )
        },
        "bodies": {}
    }

    print("\n[B] Solving Static Profiles and Extracting Charges")

    LAMBDAS = {
        "lambda_ref": LAMBDA_REFERENCE,
        "lambda_cassini": LAMBDA_REFERENCE * 1e5,
    }
    branches = {}

    for tag, lam in LAMBDAS.items():
        print(f"\n  -- {tag} (lam = {lam:.3e}) --")
        branch = {}
        for name, M, R in (("NS", M_NS, R_NS), ("WD", M_WD, R_WD)):
            print(f"Solving {name} (M={M/M_SUN:.3f} M_sun, R={R/1000:.1f} km)...")
            try:
                sol = solve_sphere(M, R, lam, x_max=1e6)
                s_ratio = diagnostics(sol, 1e5)['source_charge_ratio']
                alpha = s_ratio * (-1.0) # alpha = beta_A * S = -1 * S
                print(f"  -> Converged! Effective charge \\alpha_{name} = {alpha:.4e}")
                branch[name] = {
                    "mass_Msun": M/M_SUN,
                    "radius_m": R,
                    "S_charge_ratio": s_ratio,
                    "alpha_eff": alpha
                }
            except Exception as e:
                print(f"  -> Solver failed for {name}: {e}")
                branch[name] = {"error": str(e)}
        results["bodies"][tag] = branch
        branches[tag] = branch

    print("\n[C] Pulsar Measurement Requirement")
    results["constraint"] = {}
    for tag, branch in branches.items():
        a_ns = branch.get("NS", {}).get("alpha_eff")
        a_wd = branch.get("WD", {}).get("alpha_eff")
        if a_ns is None or a_wd is None:
            results["constraint"][tag] = {"error": "Missing charge data"}
            continue
        diff = abs(a_ns - a_wd)
        passes = diff < CHARGE_DIFF_BOUND
        print(f"[{tag}] |\\alpha_NS - \\alpha_WD| = {diff:.4e} vs bound {CHARGE_DIFF_BOUND:.1e} -> {'PASS' if passes else 'EXCEEDS'}")
        results["constraint"][tag] = {
            "predicted_diff": diff,
            "observational_bound": CHARGE_DIFF_BOUND,
            "passes_without_refitting": passes,
            "margin_over_bound": CHARGE_DIFF_BOUND / diff if diff > 0 else None,
        }
    results["constraint"]["operative_branch"] = "lambda_cassini"
    op = results["constraint"].get("lambda_cassini", {})
    if op.get("passes_without_refitting"):
        print("=> RESULT: At the operative Cassini coupling the frozen TEP screening architecture naturally suppresses differential charges below the pulsar constraint.")
    else:
        print("=> RESULT: The operative-branch differential charge EXCEEDS the pulsar constraint. The screening mechanism is insufficient without modification.")

    # Save output
    outdir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "step_24_compact_body_charge_closure.json"), "w") as f:
        json.dump(results, f, indent=2)

if __name__ == "__main__":
    run()
