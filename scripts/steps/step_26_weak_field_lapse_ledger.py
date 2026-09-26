#!/usr/bin/env python3
"""Weak-field scalar/lapse normalization ledger for canonical TEP.

The calculation keeps the Einstein-frame Newtonian potential distinct from
the matter-frame potential measured by clocks and freely falling matter.  It
is intentionally algebraic: no observational datum or fitted normalization
enters the result.

Conventions
-----------
    del^2 Phi_N       = rho_*/(2 M_Pl^2),   Phi_N < 0 near positive mass
    del^2 delta_phi   = beta_A rho_*/M_Pl
    A                 = exp(beta_A delta_phi/M_Pl)
    N/N_infinity      = exp(Phi_N) + O(Phi_N^2)
    N_tilde           = A N                         (static conformal branch)

Outputs
-------
    results/step_26_weak_field_lapse_ledger.json
"""

from tep_model import BETA, save


def run():
    beta = float(BETA)

    # Matching the two Poisson equations fixes the scalar response without a
    # fitted coefficient: delta_phi/M_Pl = 2 beta_A Phi_N.
    scalar_field_coefficient = 2.0 * beta
    conformal_lapse_coefficient = beta * scalar_field_coefficient
    einstein_lapse_coefficient = 1.0
    unscreened_matter_lapse_coefficient = (
        einstein_lapse_coefficient + conformal_lapse_coefficient
    )
    unscreened_geff_over_g = unscreened_matter_lapse_coefficient

    if abs(beta + 1.0) > 1e-15:
        raise AssertionError("The corpus-wide canonical coupling beta_A=-1 is required")
    if abs(scalar_field_coefficient + 2.0) > 1e-15:
        raise AssertionError("delta_phi/M_Pl must equal -2 Phi_N for beta_A=-1")
    if abs(conformal_lapse_coefficient - 2.0) > 1e-15:
        raise AssertionError("ln(A/A_inf) must equal 2 Phi_N")
    if abs(unscreened_matter_lapse_coefficient - 3.0) > 1e-15:
        raise AssertionError("ln(N_tilde/N_tilde_inf) must equal 3 Phi_N")

    screening_cases = {}
    for label, source_charge_fraction in {
        "fully_screened": 0.0,
        "cassini_upper_bound": 5.8e-6,
        "unscreened": 1.0,
    }.items():
        lapse_coefficient = (
            einstein_lapse_coefficient
            + conformal_lapse_coefficient * source_charge_fraction
        )
        screening_cases[label] = {
            "S_Sigma": source_charge_fraction,
            "gradient_coefficient_relative_to_Phi_N": lapse_coefficient,
            "G_eff_over_G": lapse_coefficient,
        }

    result = {
        "status": "derived",
        "inputs": {
            "beta_A": beta,
            "poisson_equation": "nabla^2 Phi_N = rho_*/(2 M_Pl^2)",
            "scalar_equation": "nabla^2 delta_phi = beta_A rho_*/M_Pl",
            "einstein_frame_lapse": "ln(N/N_infinity) = Phi_N + O(Phi_N^2)",
        },
        "derived_relations": {
            "delta_phi_over_M_Pl": "2 beta_A Phi_N",
            "delta_phi_over_M_Pl_coefficient": scalar_field_coefficient,
            "ln_A_over_A_infinity": "2 beta_A^2 Phi_N",
            "ln_A_over_A_infinity_coefficient": conformal_lapse_coefficient,
            "unscreened_matter_lapse": "N_tilde/N_tilde_infinity = exp[(1 + 2 beta_A^2) Phi_N] + O(Phi_N^2)",
            "unscreened_matter_lapse_coefficient": unscreened_matter_lapse_coefficient,
            "unscreened_G_eff_over_G": unscreened_geff_over_g,
        },
        "screened_exterior": {
            "relation": "grad ln N_tilde = [1 + 2 beta_A^2 S_Sigma(E)] grad Phi_N",
            "cases": screening_cases,
            "gr_recovery": screening_cases["fully_screened"]["G_eff_over_G"] == 1.0,
        },
        "interpretation": (
            "For beta_A=-1 the bare scalar adds 2 Phi_N to the Einstein-frame "
            "lapse, so the unscreened matter-frame response is 3 Phi_N and "
            "G_eff/G=3. Environmental suppression acts on the exterior source "
            "charge; S_Sigma -> 0 recovers the GR lapse gradient."
        ),
    }
    save("step_26_weak_field_lapse_ledger.json", result)
    print("delta_phi/M_Pl = -2 Phi_N")
    print("ln(A/A_inf) = 2 Phi_N")
    print("ln(N_tilde/N_tilde_inf) = 3 Phi_N (unscreened)")
    print("G_eff/G = 3 (unscreened), -> 1 as S_Sigma -> 0")
    return result


if __name__ == "__main__":
    run()
