#!/usr/bin/env python3
"""Fixed-action screening baseline evaluation and physical PPN diagnostics.

Evaluates the static radial scalar equation under the candidate single-body
potential baselines (V=0, quadratic, quartic branches) plus the nested
multi-scale hierarchy decomposition. The quartic branch is retained as a
candidate realization evaluated against the Cassini bound, not as the
canonical screening mechanism — Solar-System screening is performed by the
two-body kinetic operator evaluated in step_19.

The binding Cassini test is the screened-source PPN gamma built from the
far-field (1 AU) screened source charge:
gamma-1 = -4 beta_A^2 S/(1+2 beta_A^2 S), linear in S because the photon
probe is unscreened (Burrage & Sakstein 2018, Eq. 3.31).
Conjunction-point (1.6 R_sun) charge and potential ratios are near-field
diagnostics reported separately; a definitive Cassini-level verdict
requires the full ray-propagation/ephemeris observation operator.
"""
import numpy as np
from tep_model import (M_SUN,R_SUN,M_EARTH,R_EARTH,AU,M_PL,LAMBDA_REFERENCE,BETA,
                       solve_sphere,diagnostics,save,mass_squared,HBAR_C,
                       G_CM3_GEV4,source_parameters)


def run():
    results={'conventions':{'beta_A':-1,'alpha_DEF':-np.sqrt(2),
               'lambda_units':'dimensionless','field':'varphi=phi/M_Pl'},
             'solar_scan':[], 'population_fixed_lambda':[],
             'cassini_fit_performed':False}
    for factor in [1,100,1e4,1e6,1e8]:
        lam=LAMBDA_REFERENCE*factor
        row={'lambda':lam}
        try:
            sol=solve_sphere(M_SUN,R_SUN,lam,x_max=1e6)
            row.update(mu0=sol['mu'],background_varphi=sol['background_varphi'],
                background_compton_m=sol['background_compton_m'],
                max_rms_residual=sol['max_rms_residual'],
                conjunction_impact_diagnostic=diagnostics(sol,1.6),
                at_1AU=diagnostics(sol,AU/R_SUN),
                at_2646AU=diagnostics(sol,2646*AU/R_SUN))
            # Cassini binding test: PPN gamma at 1 AU (far field).
            # The PPN gamma is the far-field parameter that appears in the
            # Shapiro delay formula. The conjunction point is in the near
            # field where the PPN formalism does not apply; the 1 AU source
            # charge gives the correct PPN gamma.
            s_bound=float(2.3e-5/(4*BETA**2))  # ~5.75e-6, screened-source linear bound
            row['cassini_ppn_gamma_at_1AU']=row['at_1AU']['ppn_gamma_minus_one']
            row['cassini_source_charge_at_1AU']=row['at_1AU']['source_charge_ratio']
            row['cassini_bound_S_sigma']=s_bound
            row['cassini_passes']=bool(abs(row['at_1AU']['ppn_gamma_minus_one'])<2.3e-5)
            row['cassini_source_charge_below_bound']=bool(row['at_1AU']['source_charge_ratio']<s_bound)
        except RuntimeError as exc:
            row['solver_failure']=str(exc)
        results['solar_scan'].append(row)
    # Same lambda across all stellar masses; mu0 scales as M^2.
    for name,m,r in [('M dwarf',.15,.2),('Sun',1,1),('F dwarf',1.25,1.15),('Massive',8,5)]:
        row={'name':name,'mass_Msun':m,'radius_Rsun':r,'lambda':LAMBDA_REFERENCE}
        try:
            sol=solve_sphere(m*M_SUN,r*R_SUN,LAMBDA_REFERENCE,x_max=1e6)
            row.update(mu0=sol['mu'],at_1AU=diagnostics(sol,AU/(r*R_SUN)),
                       at_2646AU=diagnostics(sol,2646*AU/(r*R_SUN)))
        except RuntimeError as exc:
            row['solver_failure']=str(exc)
        results['population_fixed_lambda'].append(row)
    bg=2.8e-7
    results['quoted_background_consistency']={
        'varphi':bg,'lambda':LAMBDA_REFERENCE,
        'compton_m_from_Vpp':HBAR_C/np.sqrt(3*LAMBDA_REFERENCE*(M_PL*bg)**2),
        'equilibrium_density_g_cm3':LAMBDA_REFERENCE*M_PL**4*bg**3*np.exp(bg)/G_CM3_GEV4,
        'interpretation':'The quoted zero-potential Galactic amplitude cannot be inserted unchanged into the quartic equilibrium.'}
    results['equations']={
        'quartic_background_increment':'lambda*(3*phi_bg^2*dphi + 3*phi_bg*dphi^2 + dphi^3)',
        'PPN_long_range':'gamma-1 = -4*beta^2*S / (1 + 2*beta^2*S), S = source_charge_ratio; linear in S (unscreened photon probe, Burrage & Sakstein 2018 Eq. 3.31)',
        'cassini_bound':'|gamma-1| = 4*beta^2*S/(1+2*beta^2*S) < 2.3e-5 => S < 5.75e-6',
        'normalization':'mu0=lambda*(M/(4*pi*M_Pl))^2',
        'scope':'PPN gamma is the far-field parameter evaluated at 1 AU (Shapiro delay). The conjunction point (1.6 R_sun) is a near-field diagnostic; the PPN formalism applies at 1 AU, not at the conjunction point.'}
    # Legacy normalized-field sections cited by the manuscript: the
    # three-potential comparison, 3-zone profile, stellar-population scan
    # and nested screening hierarchy.
    import legacy_screening
    results.update(legacy_screening.run())
    save('step_01_radial_ode.json',results)
    for row in results['solar_scan']:
        print(row)
    return results

if __name__=='__main__':
    run()
