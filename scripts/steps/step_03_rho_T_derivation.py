#!/usr/bin/env python3
"""Derive ρ_T from the quartic + geometric identification λ_T = R_T.

The quartic V = λφ⁴/4 at the reference coupling predicts:
  m_eff² = 3λ¹ᐟ³(ρ/M_Pl)²ᐟ³  →  λ_c = 3⁻¹ᐟ² λ⁻¹ᐟ⁶ (ρ/M_Pl)⁻¹ᐟ³
  R_T = (3M/(4πρ))¹ᐟ³  (geometric saturation radius)

The ratio λ_c/R_T = 3⁻¹ᐟ²(4π/3)¹ᐟ³ λ⁻¹ᐟ⁶ (M_Pl/M)¹ᐟ³ is DENSITY-INDEPENDENT.
For Earth: λ_c/R_T ≈ 4.07.

The GNSS correlation length λ_T ≈ 4200 km is identified with R_T (not λ_c).
This identification gives ρ_T = 3M/(4π R_T³) directly from the GNSS observation.
"""
import numpy as np
from tep_model import (M_EARTH, R_EARTH, M_PL, HBAR_C, LAMBDA_REFERENCE,
                       G_CM3_GEV4, equilibrium_varphi, compton_m, save)


def saturation_radius(mass_kg, rho_g_cm3):
    return float((3*mass_kg/(4*np.pi*1000*rho_g_cm3))**(1/3))


def run():
    # GNSS observation (Paper 1/2/6, 25-year multi-centre CODE analysis)
    lambda_T_m = 4200e3

    # === PRIMARY DERIVATION: geometric identification λ_T = R_T ===
    # The GNSS correlation length IS the geometric saturation radius.
    # This gives ρ_T directly:
    rho_T_derived = 3*M_EARTH/(4*np.pi*lambda_T_m**3)/1000

    # Cross-check: R_T from calibrated ρ_T = 20 g/cm³
    rho_calibrated = 20.0
    R_T_from_cal = saturation_radius(M_EARTH, rho_calibrated)
    R_T_match = 100*abs(R_T_from_cal/lambda_T_m - 1)
    rho_T_match = 100*abs(rho_T_derived/rho_calibrated - 1)

    # === QUARTIC UNIVERSAL RATIO (density-independent) ===
    varphi_cal = equilibrium_varphi(rho_calibrated, LAMBDA_REFERENCE)
    lambda_c = compton_m(varphi_cal, rho_calibrated, LAMBDA_REFERENCE)
    ratio = lambda_c / R_T_from_cal

    # Verify density-independence: compute at a different density
    rho_test = 5.0  # g/cm³ (crust density)
    varphi_test = equilibrium_varphi(rho_test, LAMBDA_REFERENCE)
    lambda_c_test = compton_m(varphi_test, rho_test, LAMBDA_REFERENCE)
    R_T_test = saturation_radius(M_EARTH, rho_test)
    ratio_test = lambda_c_test / R_T_test

    # Analytical formula verification
    ratio_analytical = (3**(-0.5) * (4*np.pi/3)**(1/3) *
                        LAMBDA_REFERENCE**(-1/6) * (M_PL / (M_EARTH*5.60958885e26))**(1/3))
    # Note: the analytical formula uses M_Pl/M in GeV units

    # === SECONDARY CHECK: Compton self-consistency (k NOT predicted) ===
    k_inferred = lambda_c / R_EARTH

    # Trial-k self-consistency: solve λ_c(ρ_T) = k R_Earth for ρ_T at an
    # assumed k. With λ_c = 3^(-1/2) λ^(-1/6) (ρ/M_Pl)^(-1/3) (in GeV^-1
    # converted to m), λ_c ∝ ρ^(-1/3), so ρ_T(k) = ρ_ref (λ_c/(k R_Earth))^3.
    k_trial = 2.5
    rho_trial = rho_calibrated * (lambda_c / (k_trial * R_EARTH))**3
    rho_trial_diff_pct = 100 * abs(rho_trial / rho_calibrated - 1)
    # The corresponding radius discrepancy at fixed mass is R ∝ ρ^(-1/3)
    R_trial = saturation_radius(M_EARTH, rho_trial)
    R_trial_diff_pct = 100 * abs(R_trial / R_T_from_cal - 1)
    # Residual of the trial radius against the measured correlation length
    R_trial_vs_lambda_T_pct = 100 * abs(R_trial / lambda_T_m - 1)

    # Corrected (linear-in-S) Cassini-compatible coupling: mu_0 ~ 1e10,
    # i.e. lambda ~ LAMBDA_REFERENCE * 1e5 (step_01 solar scan). The
    # canonical ratio below is evaluated at the reference coupling; at the
    # required coupling it contracts by (1e5)^(-1/6) ~ 0.147.
    lam_required = LAMBDA_REFERENCE * 1e5
    lambda_c_req = compton_m(equilibrium_varphi(rho_calibrated, lam_required),
                             rho_calibrated, lam_required)

    result = {
        'lambda_reference_dimensionless': LAMBDA_REFERENCE,
        'lambda_cassini_compatible_dimensionless': lam_required,
        'primary_derivation': {
            'method': 'Geometric identification: λ_T = R_T (not λ_c)',
            'lambda_T_m': lambda_T_m,
            'identification': 'The GNSS clock covariance correlation length is the geometric saturation radius, not the Compton wavelength. This is a testable prediction: for any source, λ_T(source) = R_T(source) = (3M/(4πρ))^{1/3}.',
            'rho_T_derived_g_cm3': rho_T_derived,
            'rho_T_calibrated_g_cm3': rho_calibrated,
            'rho_T_match_percent': rho_T_match,
            'R_T_from_calibrated_m': R_T_from_cal,
            'R_T_match_percent': R_T_match,
            'classification': 'DERIVED identification, CALIBRATED value — the relation '
                            'λ_T(source) = R_T(source) ∝ M^{1/3} ρ^{-1/3} is the testable '
                            'content; the ρ_T reconstruction returns the calibrated '
                            'density through the same identification and is a '
                            'calibration identity, not an independent confirmation'
        },
        'quartic_universal_ratio': {
            'lambda_c_m': lambda_c,
            'R_T_m': R_T_from_cal,
            'ratio_lambda_c_over_R_T': ratio,
            'ratio_at_rho_5': ratio_test,
            'density_independent': abs(ratio - ratio_test) / ratio < 1e-6,
            'formula': 'λ_c/R_T = 3^(-1/2) (4π/3)^(1/3) λ^(-1/6) (M_Pl/M)^(1/3)',
            'value_for_Earth': ratio,
            'lambda_c_at_cassini_compatible_m': lambda_c_req,
            'ratio_at_cassini_compatible': lambda_c_req / R_T_from_cal,
            'classification': 'DERIVED — universal, density-independent prediction of the quartic'
        },
        'compton_self_consistency': {
            'k_inferred': k_inferred,
            'k_predicted': False,
            'independent_k_derivation_supplied': False,
            'trial_k': k_trial,
            'rho_trial_g_cm3': rho_trial,
            'rho_trial_difference_percent': rho_trial_diff_pct,
            'R_trial_m': R_trial,
            'R_trial_difference_percent': R_trial_diff_pct,
            'R_trial_vs_lambda_T_percent': R_trial_vs_lambda_T_pct,
            'note': 'k = λ_c/R_Earth is inferred from the calibrated density, not predicted ab initio. The geometric identification λ_T = R_T is the primary derivation; the Compton self-consistency is a secondary check only.',
            'classification': 'CALIBRATED — secondary check, k not predicted'
        },
        'derivation_chain': [
            '1. Reference λ = 7.526e-71 (dimensionless, μ₀=10⁵; step_01). FAILS the corrected (linear-in-S) Cassini bound at 1 AU (γ_PPN-1 ≈ -5.8e-3, ~250x over 2.3e-5); the compatible coupling is μ₀ ≳ 10^10 (λ ×10⁵)',
            '2. Quartic predicts m_eff(ρ) ∝ ρ^(1/3) (density-dependent mass)',
            '3. λ_c/R_T = universal ratio ≈ 4.07 (density-independent, derived from quartic)',
            '4. λ_T = R_T (geometric identification — testable prediction of TEP)',
            '5. R_T = 4200 km → ρ_T = 3M⊕/(4π R_T³) ≈ 19.2 g/cm³ (calibration identity)',
            '6. ρ_T ≈ 19.2 vs calibrated 20 → 3.8% (self-consistency of the identification)'
        ]
    }
    save('step_03_rho_t_derivation.json', result)
    print(f"ρ_T derived: {rho_T_derived:.2f} g/cm³ (calibrated: {rho_calibrated}, match: {rho_T_match:.1f}%)")
    print(f"λ_c/R_T = {ratio:.4f} (density-independent: {abs(ratio - ratio_test)/ratio < 1e-6})")
    print(f"R_T match: {R_T_match:.1f}%")
    return result


if __name__ == '__main__':
    run()
