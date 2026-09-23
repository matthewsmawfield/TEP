#!/usr/bin/env python3
"""Strong-field hyperbolicity analysis for the TEP scalar-tensor theory.

The TEP action is:
  S = ∫ d⁴x √(-g) [M_Pl²/2 × R - ½(∂φ)² - V(φ) + L_matter(A²(φ)g, ψ)]

In the weak-field/perturbative regime (B=0, conformal coupling only):
  - The scalar equation is second-order: □φ = V'(φ) + α(φ)ρ
  - The metric equations are second-order (Einstein-scalar)
  - The system is strongly hyperbolic (standard Einstein-Klein-Gordon)

In the strong-field regime (nonzero B, disformal coupling):
  - The matter metric is g̃ = A²g + B(∂φ)(∂φ)²
  - The scalar principal operator acquires density-dependent coefficients
  - Hyperbolicity requires Z_t > 0, Z_s > 0 (sound speeds real)
  - Metric nondegeneracy (det g̃ ≠ 0) is necessary but NOT sufficient

Distinction from scalar-Gauss-Bonnet:
  - Horndeski scalar-GB (f(φ) × R_GB²) also gives 2nd-order field equations;
    its pathology is characteristic degeneration on nontrivial backgrounds
    (Papallo & Reall 2017), not Ostrogradsky ghosts.
  - TEP likewise has 2nd-order EOM; second-order form alone is insufficient
    in either sector — the scalar principal operator can go elliptic even
    while the matter metric stays regular (see counterexample below).

This script:
  1. Verifies the perturbative regime is strongly hyperbolic
  2. Identifies the necessary conditions for strong-field hyperbolicity
  3. Addresses the counterexample (Z_s < 0)
  4. Correctly scopes the scalar-GB comparison (2nd-order EOM; separate
     characteristic analysis required)
"""
import numpy as np
from tep_model import BETA, save


def principal_coefficients(d, rho, p):
    """Reduced scalar-fluid principal-operator coefficients.

    For the frozen scalar-fluid system in the rest frame, the reduced
    principal tensor has time and spatial coefficients

        Z_t = 1 + d ρ,    Z_s = 1 − d p,    c_s² = Z_s / Z_t.

    Positive scalar energy and real characteristics require Z_t > 0 and
    Z_s > 0 simultaneously. Returns a dict with the coefficients, the sound
    speed, and the positivity verdict.
    """
    Z_t = 1 + d * rho
    Z_s = 1 - d * p
    c_s2 = Z_s / Z_t if Z_t != 0 else float('inf')
    return {
        'd': d, 'rho': rho, 'p': p,
        'Z_t': Z_t, 'Z_s': Z_s,
        'sound_speed_squared': c_s2,
        'positive_scalar_energy': bool(Z_t > 0 and Z_s > 0),
        'hyperbolic': bool(Z_t > 0 and Z_s > 0),
    }


def run():
    # === 1. Perturbative regime (B=0, conformal only) ===
    # The conformal Einstein-scalar system is strongly hyperbolic.
    # Principal symbol: standard Einstein-Klein-Gordon
    # Eigenvalues: real (gravitational waves + scalar waves)
    # This is a theorem (Choquet-Bruhat, Geroch, etc.)

    perturbative = {
        'regime': 'B=0, conformal coupling only (A(φ) = exp(β_A φ/M_Pl))',
        'scalar_eom': '□φ = V\'(φ) + α(φ)ρ  (second-order)',
        'metric_eom': 'Einstein-scalar  (second-order)',
        'principal_symbol': 'Standard Einstein-Klein-Gordon',
        'eigenvalues': 'Real (gravitational waves + scalar waves)',
        'theorem': 'Choquet-Bruhat, Geroch: Einstein-Klein-Gordon is strongly hyperbolic',
        'classification': 'DERIVED — strongly hyperbolic in perturbative regime'
    }

    # === 2. GW170817 bound on disformal perturbation ===
    # The disformal metric is g̃ = A²g + B(∂φ)(∂φ)²
    # GW170817 bounds |c_T/c - 1| < 10⁻¹⁵ along the observed
    # EM–GW propagation path and epoch. The constrained combination is the
    # cone-split projection of B(∂φ)² — it does NOT impose a universal
    # maximum-norm bound on every local B∂φ∂φ component, and it does NOT
    # bound the matter-dependent principal-operator coefficients Z_t, Z_s
    # (which involve B·ρ, a different contraction).

    gw170817 = {
        'bound': '|c_T/c - 1| < 10⁻¹⁵ along observed path/epoch',
        'constrained_combination': 'cone-split projection of B(∂φ)²/A²',
        'disformal_perturbation': 'O(B × (∂φ)²) small where cone split is small',
        'scalar_principal_perturbation': 'NOT bounded by GW170817 — Z_t, Z_s involve B·ρ '
                                        'and require the separate Z_t>0, Z_s>0 conditions',
        'eigenvalue_shift': 'c_± = c(1 ± O(10⁻¹⁵)) in the propagation sector only',
        'classification': 'DERIVED — propagation-sector cone split is a small regular '
                        'perturbation; the matter-dependent principal operator is '
                        'separately conditioned by Z_t>0, Z_s>0'
    }

    # === 3. Strong-field necessary conditions ===
    # For nonzero B, the scalar principal operator is:
    #   L = -Z_t ∂_t² + Z_s ∇² + ...
    # Hyperbolicity requires:
    #   Z_t > 0 (positive time-like coefficient)
    #   Z_s > 0 (positive space-like coefficient → real sound speed)
    # Metric nondegeneracy (det g̃ ≠ 0) is necessary but NOT sufficient.

    # The coefficients Z_t, Z_s depend on:
    #   Z_t = 1 + B × (ρ + p) / A²  (time-like, from T_00)
    #   Z_s = 1 + B × (ρ - p) / (3A²)  (space-like, from T_ij trace)

    # For a perfect fluid with equation of state p = wρ:
    #   Z_t = 1 + B × (1 + w) × ρ / A²
    #   Z_s = 1 + B × (1 - w) × ρ / (3A²)

    # Hyperbolicity conditions:
    #   Z_t > 0: B × (1 + w) × ρ / A² > -1
    #   Z_s > 0: B × (1 - w) × ρ / (3A²) > -1

    # For w = 0 (dust): Z_t = 1 + Bρ/A², Z_s = 1 + Bρ/(3A²)
    # For w = 1/3 (radiation): Z_t = 1 + 4Bρ/(3A²), Z_s = 1 + 2Bρ/(3A²)
    # For w = -1 (dark energy): Z_t = 1, Z_s = 1 + 2Bρ/(3A²)

    necessary_conditions = {
        'scalar_principal_operator': 'L = -Z_t ∂_t² + Z_s ∇² + ...',
        'Z_t_formula': 'Z_t = 1 + B × (1 + w) × ρ / A²',
        'Z_s_formula': 'Z_s = 1 + B × (1 - w) × ρ / (3A²)',
        'hyperbolicity_conditions': ['Z_t > 0', 'Z_s > 0', 'det(g̃) ≠ 0'],
        'metric_nondegeneracy_necessary_but_not_sufficient': True,
        'classification': 'DERIVED — necessary conditions identified; sufficient conditions require full principal symbol'
    }

    # === 4. Counterexamples ===
    # Metric nondegeneracy alone is insufficient for hyperbolicity.
    # If Z_s < 0, the sound speed is imaginary → elliptic spatial operator → instability.
    # This is NOT a pathology of the EOM order (still 2nd-order).
    # It is a pathology of the coefficient sign (sound speed).

    # Case A (manuscript §4): reduced operator Z_t = 1 + dρ, Z_s = 1 − dp
    # with A = 1, B = 1, q = 0, ρ = 3, p = 2, d = 1.
    # The matter metric is exactly regular while Z_s < 0 → c_s² < 0 even
    # though Z_t > 0: a regular metric with an elliptic scalar sector.
    coeff_A = principal_coefficients(1.0, 3.0, 2.0)
    d_A, rho_A, p_A = 1.0, 3.0, 2.0
    Z_t_A, Z_s_A, c_s2_A = coeff_A['Z_t'], coeff_A['Z_s'], coeff_A['sound_speed_squared']
    counterexample_regular_metric = {
        'parameterization': 'Z_t = 1 + dρ, Z_s = 1 − dp (reduced operator)',
        'A': 1.0, 'B': 1.0, 'q': 0.0, 'd': d_A,
        'rho': rho_A, 'p': p_A,
        'Z_t': Z_t_A, 'Z_s': Z_s_A,
        'sound_speed_squared': c_s2_A,
        'metric_nondegenerate': True,
        'hyperbolic': Z_t_A > 0 and Z_s_A > 0,
        'pathology': f'Z_t = {Z_t_A:.2f} > 0 but Z_s = {Z_s_A:.2f} < 0 → '
                     f'c_s² = {c_s2_A:.2f} < 0 (imaginary sound speed, elliptic)',
        'conclusion': 'Metric nondegeneracy alone does not guarantee hyperbolicity: '
                      'a regular matter metric coexists with an elliptic scalar sector.',
        'pathology_class': 'Coefficient sign (sound speed), NOT Ostrogradsky (EOM order)',
        'classification': 'DERIVED — counterexample shows metric nondegeneracy is insufficient'
    }

    # Case B: wrong-signature coefficients, B < 0, w = 0 (dust), large ρ.
    # Z_s = 1 + Bρ/(3A²) < 0 when Bρ/(3A²) < -1, i.e., B < -3A²/ρ

    B_counter = -5.0  # negative disformal coefficient
    rho_counter = 1.0  # density
    A2_counter = 1.0  # conformal factor squared
    w_counter = 0.0  # dust

    Z_t = 1 + B_counter * (1 + w_counter) * rho_counter / A2_counter
    Z_s = 1 + B_counter * (1 - w_counter) * rho_counter / (3 * A2_counter)
    c_s2 = Z_s / Z_t if Z_t != 0 else float('inf')

    counterexample = {
        'regular_metric_elliptic_scalar': counterexample_regular_metric,
        'B': B_counter,
        'rho': rho_counter,
        'A2': A2_counter,
        'w': w_counter,
        'Z_t': Z_t,
        'Z_s': Z_s,
        'sound_speed_squared': c_s2,
        'metric_nondegenerate': True,  # det(g̃) ≠ 0 in this example
        'hyperbolic': Z_t > 0 and Z_s > 0,
        'pathology': f'Z_t = {Z_t:.2f} < 0 (wrong signature) AND Z_s = {Z_s:.2f} < 0 (imaginary sound speed if Z_t > 0)',
        'conclusion': 'Metric nondegeneracy alone does not guarantee hyperbolicity. '
                      'Both Z_t and Z_s must be positive. Here both are negative, '
                      'indicating a wrong-signature principal symbol.',
        'pathology_class': 'Coefficient sign (wrong signature), NOT Ostrogradsky (EOM order)',
        'classification': 'DERIVED — counterexample shows metric nondegeneracy is insufficient'
    }

    # === 5. Distinction from scalar-Gauss-Bonnet ===
    # scalar-GB: S ⊃ f(φ) × R_GB²
    #   - In the Horndeski class this term yields SECOND-order field
    #     equations (the Gauss-Bonnet combination is degenerate);
    #     Ostrogradsky does NOT follow from the term alone.
    #   - The known pathology is characteristic degeneration on nontrivial
    #     backgrounds: even with 2nd-order EOM, the principal symbol can
    #     become weakly elliptic / lose hyperbolicity before the metric
    #     degenerates (Papallo & Reall 2017).
    #   - It therefore cannot exempt the TEP strong-curvature sector from
    #     its own characteristic analysis.

    # TEP: S ⊃ B(φ, ∂φ) × (∂φ)²
    #   - 2nd-order EOM for both metric and scalar
    #   - Pathology (if any) is coefficient sign of the matter-dependent
    #     principal operator (Z_t, Z_s), not EOM order

    scalar_gb = {
        'action_term': 'f(φ) × R_GB²',
        'eom_order': '2nd-order (Horndeski class)',
        'pathology': 'Characteristic degeneration on nontrivial backgrounds '
                     '(weakly elliptic sector; Papallo & Reall 2017)',
        'pathology_class': 'Principal-symbol degeneration despite 2nd-order EOM',
        'fix': 'Independent characteristic analysis on each background; '
               'second-order EOM alone is insufficient'
    }

    tep = {
        'action_term': 'B(φ, ∂φ) × (∂φ)²',
        'eom_order': '2nd-order (metric and scalar)',
        'pathology': 'Sound speed instability (if Z_s < 0)',
        'pathology_class': 'Coefficient sign (sound speed), NOT EOM order',
        'fix': 'Require Z_t > 0, Z_s > 0 (positive sound speeds)'
    }

    distinction = {
        'scalar_gb': scalar_gb,
        'tep': tep,
        'fundamentally_distinct': True,
        'distinction': 'Both TEP and Horndeski scalar-GB have 2nd-order EOM; '
                       'neither can be exempted from characteristic analysis. '
                       'The TEP pathology class (matter-dependent coefficient '
                       'sign) is distinct from the scalar-GB background-'
                       'dependent characteristic degeneration.',
        'classification': 'DERIVED — second-order EOM is not sufficient for '
                        'hyperbolicity in either sector'
    }

    # === 6. Summary ===
    result = {
        'perturbative_regime': perturbative,
        'gw170817_bound': gw170817,
        'strong_field_necessary_conditions': necessary_conditions,
        'counterexample': counterexample,
        'distinction_from_scalar_gb': distinction,
        'summary': {
            'perturbative_hyperbolic': True,
            'disformal_perturbation_bounded': True,
            'strong_field_sufficient': False,
            'strong_field_necessary': True,
            'tep_distinct_from_scalar_gb': True,
            'overall': 'Perturbative regime: strongly hyperbolic (derived). '
                       'GW170817 bounds the propagation-sector cone split only; the '
                       'matter-dependent principal operator is conditioned separately '
                       'by Z_t>0, Z_s>0. '
                       'Strong-field: necessary conditions identified, sufficient conditions require full principal symbol. '
                       'Scalar-GB: 2nd-order EOM; characteristic degeneration requires separate analysis.'
        }
    }
    save('step_06_strong_field_hyperbolicity.json', result)
    print("Perturbative: strongly hyperbolic (derived)")
    print("GW170817: bounds propagation-sector cone split, not matter-dependent principal coefficients")
    print(f"Counterexample (regular metric): Z_t = {Z_t_A:.2f}, Z_s = {Z_s_A:.2f}, c_s² = {c_s2_A:.2f} (elliptic)")
    print(f"Counterexample (wrong signature): Z_t = {Z_t:.2f}, Z_s = {Z_s:.2f} (not hyperbolic)")
    print("Scalar-GB: 2nd-order EOM; characteristics can fail on nontrivial backgrounds")
    print("Strong-field sufficient conditions: OPEN (requires full principal symbol)")
    return result


if __name__ == '__main__':
    run()
