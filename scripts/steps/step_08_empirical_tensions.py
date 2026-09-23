#!/usr/bin/env python3
"""Empirical tension analysis: MGEX, JWST, LLR/flyby.

This script analyzes the three empirical tensions using the derived TEP
response coefficients (steps 03-05) and the actual observational data.
No synthetic data is used. All numbers trace to real observations.

Tensions:
  1. MGEX: 1862 ± 112 km vs 4200 km correlation length on a different product type
  2. JWST: mixed evidence structure across comparison spaces
  3. LLR/flyby: post-fit residual circularity
"""
import numpy as np
from tep_model import (M_EARTH, R_EARTH, M_SUN, AU, PC, M_PL, HBAR_C,
                       BETA, LAMBDA_REFERENCE, solve_sphere, diagnostics,
                       equilibrium_varphi, compton_m, save)


def run():
    # === 1. MGEX: 1862 km vs 4200 km ===
    # The canonical GNSS correlation length is λ_T ≈ 4200 km (CODE, 25-year).
    # The MGEX held-out replication (Paper 14, v0.2) returns λ = 1862 ± 112 km
    # on the combined multi-GNSS product, and recovers an anisotropy axis at
    # RA = 180°, Dec = 10°, i.e. within 21.4° of the CMB dipole (exploratory,
    # full-grid permutation p = 0.025). The scale differs by a factor ~2.3;
    # the preferred axis remains CMB-consistent.

    # TEP interpretation: the nested hierarchy
    #   galactic → solar → terrestrial → internal Earth structures
    # The 4200 km scale is R_T (geometric saturation radius, outer structure).
    # The 1862 km scale lies between the inner-core and outer-core hierarchy
    # levels, consistent with a product sensitive to a deeper/denser layer.

    # Inner core radius: ~1220 km (PREM)
    # Outer core radius: ~3480 km (PREM)
    # R_T (from ρ_T = 20): ~4146 km

    # If the MGEX product is sensitive to a different density layer, the
    # correlation length would be R_T for that layer:
    #   R_T(layer) = (3 M_layer / (4π ρ_layer))^(1/3)

    # Inner core: M_inner ~ 1.7e23 kg (1.8% of Earth mass), ρ_inner ~ 13 g/cm³
    M_inner = 0.018 * M_EARTH
    rho_inner = 13.0  # g/cm³
    R_T_inner = float((3 * M_inner / (4 * np.pi * 1000 * rho_inner))**(1/3))

    # Outer core: M_outer ~ 1.8e24 kg (30% of Earth mass), ρ_outer ~ 10 g/cm³
    M_outer = 0.30 * M_EARTH
    rho_outer = 10.0
    R_T_outer = float((3 * M_outer / (4 * np.pi * 1000 * rho_outer))**(1/3))

    # The 1862 km scale vs inner-core and outer-core R_T
    mgex_match_inner = 100 * abs(R_T_inner - 1862e3) / 1862e3
    mgex_match_outer = 100 * abs(R_T_outer - 1862e3) / 1862e3

    # The recovered axis lies 21.4° from the CMB dipole: the same preferred
    # direction as the canonical analysis. The product-dependent quantity is
    # the length scale, not the axis — consistent with a hierarchy-level
    # sensitivity difference rather than a distinct physical component.

    mgex = {
        'canonical_length_km': 4200,
        'mgex_length_km': 1862,
        'mgex_length_err_km': 112,
        'axis_offset_from_cmb_deg': 21.4,
        'tep_interpretation': 'Nested hierarchy: 4200 km = R_T (outer structure), '
                              '1862 km lies between inner-core and outer-core R_T; '
                              'axis remains CMB-consistent',
        'R_T_inner_core_km': R_T_inner / 1000,
        'R_T_outer_core_km': R_T_outer / 1000,
        'inner_core_match_percent': mgex_match_inner,
        'outer_core_match_percent': mgex_match_outer,
        'axis_interpretation': 'Recovered axis within 21.4° of the CMB dipole: '
                               'same preferred direction as the canonical product; '
                               'the discrepancy is in scale, not orientation',
        'classification': 'UNRESOLVED — product-dependent replication difference in scale; '
                        'axis consistent with canonical direction; inner-core '
                        'identification is a hypothesis requiring a prespecified prediction',
        'resolution_test': 'Product-level transfer model: freeze observation operator for '
                          'each product, predict both length and axis from TEP hierarchy'
    }

    # === 2. JWST: mixed evidence structure across comparison spaces ===
    # The JWST analysis (Paper 12, v0.7) reports a covariance-corrected joint
    # comparison favoring TEP at ln BF = +64.5 over the standard mass-plus-z
    # model with four fewer parameters, and an orthogonalized sensitivity set
    # favoring TEP by ln BF = +91 to +138 across eleven alternatives. The
    # unfavorable sub-comparisons are: the conventional raw-mass residual-space
    # comparison (ln BF = -3.9, residual null preferred), the mass-correction
    # transfer onto standard models (mean ln BF = -3.1; AGN -2.8), and the
    # single-observable dust AIC comparison (AGN-threshold model preferred,
    # Akaike weight 0.83).

    # TEP interpretation:
    # At high redshift (z ~ 10), the TEP transfer factor is:
    #   T(z) = (1 + z)^(β_A × Δvarphi)
    # For Δvarphi ~ 3e-7 (galactic field differential):
    #   T(10) = 11^(-3e-7) ≈ 1 - 3e-7 × ln(11) ≈ 1 - 7e-7
    # This is a very small correction (~7e-7), so the decisive evidence comes
    # from the correlation structure, not the SED-level correction.

    delta_ln_Z = -3.9
    bayes_factor = float(np.exp(delta_ln_Z))
    z_example = 10.0
    delta_varphi = 3e-7
    transfer = float((1 + z_example)**(BETA * delta_varphi))
    transfer_correction = abs(transfer - 1)

    jwst = {
        'primary_joint_ln_BF': 64.5,
        'conventional_residual_ln_BF': delta_ln_Z,
        'corrected_residual_ln_BF': 151.4,
        'constrained_agn_residual_ln_BF': 106.4,
        'mass_correction_transfer_mean_ln_BF': -3.1,
        'dust_aic_agn_weight': 0.83,
        'bayes_factor': bayes_factor,
        'classification': 'MIXED — primary covariance-corrected joint evidence favors TEP '
                        '(ln BF = +64.5); conventional residual space favors the null '
                        '(ln BF = -3.9); corrected residual space favors TEP '
                        '(ln BF = +151.4, +106.4 vs constrained AGN)',
        'tep_transfer_factor_z10': transfer,
        'tep_correction_magnitude': transfer_correction,
        'interpretation': 'At z~10, the direct SED correction is ~7e-7 (very small); '
                         'the decisive evidence comes from the correlation structure. '
                         'The unfavorable residual-space result holds only in the '
                         'conventional raw-mass space where signal absorption is expected.',
        'parameter_count_argument': 'AGN models have more free parameters, but parameter '
                                   'count alone does NOT predict evidence ranking. '
                                   'A physically normalized likelihood is required.',
        'resolution_test': 'Lower-redshift analysis where TEP correction is larger; '
                          'prespecify physically normalized SED response'
    }

    # === 3. LLR/flyby: post-fit residual circularity ===
    # The LLR analysis fits the lunar orbit to a model including:
    #   - GM_Earth (absorbs 1/r monopole)
    #   - J₂, J₃, ... (absorbs low-order multipoles)
    #   - Tidal terms (absorbs time-varying components)

    # The TEP perturbation is a Yukawa force with scale R_T ≈ 4200 km.
    # At the lunar distance (r ≈ 384,400 km, r/R_T ≈ 92):
    #   F_Yukawa = exp(-r/R_T) / r² ≈ exp(-92) / r² ≈ 10⁻⁴⁰ / r²
    # This is exponentially suppressed.

    # The GR pipeline absorbs the 1/r monopole. The Yukawa tail is:
    #   F_residual = F_Yukawa - (absorbed 1/r component)
    # For r >> R_T, the Yukawa force is essentially zero, so the residual is zero.

    # BUT: the TEP perturbation is NOT just a Yukawa force. It also includes:
    #   1. Clock-rate effect: δ(ln A) = β_A × δ(varphi)
    #   2. Disformal metric perturbation: O(B × (∂φ)²)
    #   3. Time-varying field: from Earth's motion through the galactic field

    # The clock-rate effect at the lunar distance:
    #   δ(ln A) = β_A × (varphi(r) - varphi_ambient)
    # At r >> R_T, varphi(r) → varphi_ambient, so δ(ln A) → 0.
    # The clock-rate effect is also exponentially suppressed.

    # The disformal metric perturbation is O(10⁻¹⁵) (from GW170817).
    # This is a regular perturbation of the principal symbol, not a new force.

    # The time-varying field from Earth's motion:
    #   d(varphi)/dt ~ varphi_ambient × (v_Earth/c) / r_gal
    # This is a very slow, smooth variation that is absorbed by the tidal terms.

    # Conclusion: at the lunar distance, the TEP perturbation is exponentially
    # suppressed. The LLR residual is essentially zero. The TEP effect is
    # NOT distinguishable at the lunar distance through the standard channels.

    # HOWEVER: the Earth flyby (close approach, r ~ R_Earth) is different:
    #   At r ~ R_Earth, the field is NOT suppressed (exp(-R_Earth/R_T) ≈ exp(-1.5) ≈ 0.22).
    #   The flyby velocity shift is:
    #     δv/v = β_A × varphi_Earth × exp(-r/R_T)
    #   For r = R_Earth: δv/v ≈ (-1) × 1.4e-9 × 0.22 ≈ -3e-10
    #   This is much smaller than the observed flyby anomaly (~10⁻⁶).

    # The flyby anomaly is NOT explained by the conformal coupling alone.
    # The disformal coupling or a different mechanism may be needed.

    r_llr = 384400e3  # lunar distance
    rho_T = 20.0
    R_T = float((3 * M_EARTH / (4 * np.pi * 1000 * rho_T))**(1/3))
    r_over_RT = r_llr / R_T
    yukawa = float(np.exp(-r_over_RT))

    # Flyby at close approach
    r_flyby = R_EARTH
    r_flyby_over_RT = r_flyby / R_T
    yukawa_flyby = float(np.exp(-r_flyby_over_RT))
    varphi_earth = float(M_EARTH * 5.60958885e26 / (4 * np.pi * M_PL**2 * R_EARTH / HBAR_C))
    flyby_delta_v_over_v = abs(BETA) * varphi_earth * yukawa_flyby

    llr = {
        'lunar_distance_km': r_llr / 1000,
        'R_T_km': R_T / 1000,
        'r_over_RT': float(r_over_RT),
        'yukawa_suppression': yukawa,
        'classification': 'NOT CONSTRAINING at lunar distance — Yukawa force exponentially '
                        'suppressed (exp(-92) ≈ 10⁻⁴⁰). LLR cannot distinguish TEP from GR.',
        'flyby_analysis': {
            'r_flyby_km': r_flyby / 1000,
            'r_flyby_over_RT': float(r_flyby_over_RT),
            'yukawa_flyby': yukawa_flyby,
            'varphi_earth': varphi_earth,
            'delta_v_over_v': flyby_delta_v_over_v,
            'observed_anomaly': 1e-6,
            'tep_prediction': flyby_delta_v_over_v,
            'prediction_vs_observed': flyby_delta_v_over_v / 1e-6,
            'classification': 'NOT EXPLAINED — TEP conformal coupling predicts '
                            f'δv/v ≈ {flyby_delta_v_over_v:.1e}, much smaller than '
                            'observed ~10⁻⁶. Disformal coupling or different mechanism needed.'
        },
        'resolution_test': 'Forward-dynamics integration: integrate TEP force model '
                          'through LLR observation times, construct variational/design '
                          'matrices, refit nuisance parameters, evaluate retained signal. '
                          'For flyby: include disformal coupling in the trajectory model.'
    }

    result = {
        'mgex': mgex,
        'jwst': jwst,
        'llr_flyby': llr,
        'summary': {
            'mgex': 'UNRESOLVED — 1862 km lies between inner-core and outer-core R_T; '
                    'axis consistent with CMB direction (21.4°)',
            'jwst': 'MIXED — primary covariance-corrected evidence favors TEP '
                    '(ln BF = +64.5); conventional residual space favors the null',
            'llr': 'NOT CONSTRAINING at lunar distance; flyby NOT EXPLAINED by conformal coupling alone',
            'overall': 'MGEX scale unresolved between hierarchy levels; axis consistent. '
                      'JWST evidence mixed across comparison spaces. '
                      'LLR not constraining (exponential suppression). Flyby requires disformal coupling.'
        }
    }
    save('step_08_empirical_tensions.json', result)
    print(f"MGEX: {mgex['classification']}")
    print(f"  R_T(inner core) = {R_T_inner/1000:.0f} km, match = {mgex_match_inner:.1f}%")
    print(f"JWST: {jwst['classification']}")
    print(f"  TEP correction at z=10: {transfer_correction:.1e}")
    print(f"LLR: {llr['classification']}")
    print(f"  Yukawa suppression: exp(-{r_over_RT:.0f}) = {yukawa:.1e}")
    print(f"  Flyby: δv/v = {flyby_delta_v_over_v:.1e} (observed ~1e-6)")
    return result


if __name__ == '__main__':
    run()
