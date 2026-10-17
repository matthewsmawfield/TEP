#!/usr/bin/env python3
"""Fail on numerical or scope regressions in the selected consistency repairs.
This validates specified calculations, not the full TEP framework.
"""
import hashlib
import json
from pathlib import Path
import numpy as np
from step_20_landscape_existence import build_landscape, fd_laplacian, solve_landscape, M_PL

ROOT = Path(__file__).resolve().parents[2]


def read(name):
    return json.loads((ROOT/'results'/name).read_text())


def validate():
    checks = {}
    lapse = read('step_21_lapse_solvability.json')
    checks['physical_CMC_initial_slices'] = all(
        x['hamiltonian_converged'] and x['hamiltonian_max_residual'] < 1e-8
        and x['evolving_cmc_initial_slice']['lapse_range'][0] > 0
        and x['evolving_cmc_initial_slice']['lapse_equation_relative_residual'] < 1e-8
        and not x['maximal_slice_diagnostic']['canonical_closure_established']
        for x in lapse['grids'])
    checks['khronon_imposed_stress_lapse_equation_solved'] = all(
        x.get('noncanonical_khronon_static_slice', {}).get('static_maximal_slice_established_on_khronon_branch', False)
        for x in lapse['grids'])
    congruence = read('step_21_congruence_preservation.json')
    checks['negative_kinetic_density_cannot_pass'] = not (
        congruence['joint_solution']['Pi_star2'] < 0 and congruence['verdict'] == 'PASS')
    spectra = read('step_79_inhomogeneous_eigenmodes.json')
    checks['homogeneous_gauge_status_not_overclaimed'] = (
        'UNRESOLVED' in spectra['classification_verdict']['homogeneous_secular_mode'])
    checks['saved_lowest_coupled_eigenpairs_resolved'] = all(
        np.isfinite(value) and value < 1e-4
        for run in (spectra['primary_run_N24'], spectra['convergence_check_N32'])
        for value in run['full_coupled_einstein_scalar_operator']['eigenpair_residual_norms'].values())
    bloch = spectra['supercell_bloch_and_gauge_reduction_audit']
    checks['actual_scalar_Bloch_eigenpairs_obey_bound'] = all(
        x['lowest_scalar_eigenvalue'] >= bloch['all_twist_scalar_lower_bound'] - 1e-8
        and x['eigenpair_absolute_residual'] < 1e-6 and x['hermitian_error'] < 1e-8
        for x in bloch['sampled_bloch_eigenpairs'])
    checks['conditional_Bloch_scope_explicit'] = not (
        bloch['full_coupled_infinite_volume_stability_established']
        or bloch['background_resolved_for_quartic']
        or bloch['gauge_reduction_established'])
    checks['resolved_quartic_coupled_spectrum_positive_unprojected'] = bool(
        spectra.get('resolved_quartic_two_branch_coupled_audit', {}).get('all_grids_strictly_positive_unprojected', False))
    quartic = spectra['resolved_quartic_two_branch_coupled_audit']
    checks['quartic_benchmark_scope_and_residuals_verified'] = (
        not quartic.get('scalar_background_equation_solved', True)
        and not quartic.get('full_action_stability_established', True)
        and all(r['hamiltonian_max_residual'] < 1e-8
                and min(r['unprojected_coupled_lowest_evals']) >= r['all_twist_coupled_lower_bound']
                and all(v < 1e-5 for v in r['eigenpair_residual_norms'].values())
                for r in quartic['grids']))
    joint_q = spectra.get('joint_on_shell_quartic_two_branch_coupled_audit', {})
    checks['joint_on_shell_quartic_coupled_spectrum_positive_unprojected'] = bool(
        joint_q.get('scalar_background_equation_solved', False)
        and joint_q.get('background_resolved_for_quartic', False)
        and joint_q.get('all_grids_strictly_positive_unprojected', False)
        and all(r['scalar_euler_lagrange_max_residual'] < 1e-10
                and r['hamiltonian_max_residual'] < 1e-10
                and r['all_twist_coupled_lower_bound'] > 0
                and min(r['unprojected_coupled_lowest_evals']) >= r['all_twist_coupled_lower_bound']
                and all(v < 1e-7 for v in r['eigenpair_residual_norms'].values())
                for r in joint_q.get('grids', [])))
    comp_c = spectra.get('completed_canonical_reconstructed_coupled_audit', {})
    checks['completed_canonical_reconstructed_matter_second_variation_positive'] = bool(
        comp_c.get('includes_matter_coupling_second_variation', False)
        and comp_c.get('includes_metric_constraint_response', False)
        and not comp_c.get('requires_homogeneous_mode_projection', True)
        and comp_c.get('all_grids_strictly_positive_unprojected', False)
        and all(r['hamiltonian_max_residual'] < 1e-10
                and r['bare_V_rec_pp_mean'] < 0
                and r['complete_V_eff_pp_mean'] > 0
                and r['lowest_unprojected_euclidean_lambda_0'] > 0.4
                and r['lowest_unprojected_coupled_lambda_0'] > 0.4
                and all(v < 1e-7 for v in r['eigenpair_residual_norms'].values())
                for r in comp_c.get('grids', [])))
    bound = read('step_82_stability_lower_bound.json')
    checks['ambient_discrete_bounds_positive'] = all(
        x['branches']['ambient']['operator_lower_bound'] > 0 for x in bound['grids'])
    checks['pure_local_discrete_bounds_positive'] = all(
        x['branches']['pure_local']['operator_lower_bound'] > 0 for x in bound['grids'])
    # Independently check the background sign against the original Hamiltonian
    # equation, rather than accepting a converged flag for the wrong sign.
    background_residuals = {}
    for n in (24, 32):
        data = build_landscape(n, 1., 4, .05, 1., .8, 2., .05, 2., 1., 15.)
        delta = fd_laplacian(n, 1/n)
        sol = solve_landscape(data['S0'], data['g'], delta, verbose=False)
        psi = sol['psi']
        if n == 24:
            saved = spectra['primary_run_N24']
            checks['saved_spectrum_uses_physical_Hamiltonian_background'] = bool(
                abs(sol['lam'] - saved['lambda_offset']) < 1e-9
                and abs(float(psi.min()) - saved['psi_min']) < 1e-9
                and abs(float(psi.max()) - saved['psi_max']) < 1e-9)
        residual = ((delta @ psi.ravel()).reshape(psi.shape)
                    + psi**5*(data['S0']+sol['lam'])/(4*M_PL**2)
                    + psi*data['g']/(8*M_PL**2))
        background_residuals[str(n)] = float(np.max(np.abs(residual)))
    checks['physical_Hamiltonian_sign_verified'] = all(v < 1e-8 for v in background_residuals.values())
    rays = read('step_80_shared_lapse_action_derivation.json')['ppn_and_matter_sector']
    ray = rays['ray_integrated_cassini_shapiro_and_doppler']
    threshold = rays['cassini_bound_abs_gamma_minus_1']
    checks['specified_solar_ray_benchmarks_pass'] = all(
        abs(ray[k]) < threshold for k in
        ('gamma_PPN_minus_1_ray_doppler_observable', 'gamma_PPN_minus_1_ray_shapiro_observable'))
    action_audit = read('step_81_action_contract_audit.json')
    checks['acoustic_shell_flux_and_prescription_values_verified'] = bool(
        action_audit.get('acoustic_shell_reduction', {}).get('matches_bbn_step20_within_1e12', False)
        and action_audit['acoustic_shell_reduction']['flux_law_residual'] < 1e-12)
    shell = action_audit['acoustic_shell_reduction']
    checks['native_shell_hessian_not_conflated_with_effective_prescription'] = (
        shell['native_secant_identity_verified']
        and abs(shell['native_secant_value'] - 1/shell['y_s_sound_horizon']) < 1e-12
        and shell['native_longitudinal_value'] > shell['native_secant_value']
        and not shell['quadratic_perturbation_mapping_established'])
    checks['spherical_shell_hessian_invariants_verified'] = bool(
        shell.get('spherical_rms_ir_limit_equals_sqrt2_over_ys', False)
        and abs(shell['spherical_frobenius_rms_value'] - shell['K_sp_spatial_sound_shell_derived']) / shell['K_sp_spatial_sound_shell_derived'] < 1e-3)
    paths = list((ROOT/'scripts/steps').glob('step_21_*.py')) + [
        ROOT/'scripts/steps'/f for f in (
            'step_20_landscape_existence.py',
            'step_79_inhomogeneous_eigenmodes.py', 'step_80_shared_lapse_action_derivation.py',
            'step_81_action_contract_audit.py', 'step_82_stability_lower_bound.py', Path(__file__).name)]
    paths += [ROOT/'results'/name for name in (
        'step_21_lapse_solvability.json', 'step_21_congruence_preservation.json',
        'step_79_inhomogeneous_eigenmodes.json', 'step_80_shared_lapse_action_derivation.json',
        'step_81_action_contract_audit.json', 'step_82_stability_lower_bound.json')]
    result = {
        'scope': __doc__, 'checks': checks,
        'hamiltonian_max_residuals': background_residuals,
        'status': 'PASS' if all(checks.values()) else 'FAIL',
        'provenance': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'limitations': [
            'Lowest-mode checks do not certify convergence of every higher eigenmode.',
            'Positive finite-grid bounds do not establish a full continuum coupled theory.',
            'CMC initial lapse consistency does not establish eternal zero expansion.',
            'Shared-action observational transfer and absolute loop normalization retain their stated conditional scope.'],
    }
    (ROOT/'results/step_83_consistency_validation.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    if not all(checks.values()):
        raise RuntimeError('Failed checks: '+', '.join(k for k,v in checks.items() if not v))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    validate()
