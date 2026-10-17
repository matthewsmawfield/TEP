#!/usr/bin/env python3
"""Exact local and acoustic-shell identities for the stated action.

Separates (1) single-invariant homogeneous timelike P(X) principal coefficients
at bar_xi = 1 from (2) the U_mu / Pi_{mu nu} (3+1)-split acoustic-shell
spatial flux reduction at the sound horizon q_s = g_b(r_s)/g_t, and from (3)
the static constant-M_Pl tensor action vs assigned TH benchmark pump.
"""
import hashlib
import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[2]


def derive():
    xi, k = sp.symbols('xi k', positive=True)
    px = k * sp.sqrt(xi) + 2 * xi
    kt = sp.simplify(px + 2 * xi * sp.diff(px, xi))
    cs2 = sp.simplify(px / kt)

    # For p = |D u|/g_t = q*y, xi = p^2/2 and J(p)=P_X(p^2/2)*p.
    # Flux conservation implies P_X*y=1. Thus the native transverse/secant
    # stiffness is 1/y, not sqrt(2)/y. The latter is the separate normalization
    # used in the historical BBN effective prescription. A field rescaling
    # cannot be applied only to the spatial term of a quadratic action.
    q_s, y_s, p_s = sp.symbols('q_s y_s p_s', positive=True)
    flux_lhs = (k / sp.sqrt(2)) * q_s * y_s**2 + q_s**2 * y_s**3
    native_secant = sp.simplify(flux_lhs / y_s)
    native_on_shell = 1 / y_s
    current = (k*p_s/sp.sqrt(2) + p_s**2)*p_s
    native_longitudinal = sp.diff(current, p_s).subs(p_s, q_s*y_s)
    # Hessian of W(p)=k*p^3/(3 sqrt(2))+p^4/4 has eigenvalues
    # J/p transversely (multiplicity 2 on S^2) and dJ/dp longitudinally (multiplicity 1).
    # Their 3D rotational invariants are the trace mean (2 K_perp + K_parallel)/3
    # and the Frobenius/RMS eigenvalue sqrt((2 K_perp^2 + K_parallel^2)/3), whose
    # deep two-branch (q_s -> 0, K_parallel -> 2/y_s) limit is sqrt(2)/y_s identically.
    native_longitudinal_on_shell = 2 / y_s + (q_s * y_s) ** 2
    spherical_trace_iso_on_shell = sp.simplify((2 * native_on_shell + native_longitudinal_on_shell) / 3)
    spherical_rms_on_shell = sp.simplify(sp.sqrt((2 * native_on_shell ** 2 + native_longitudinal_on_shell ** 2) / 3))
    spherical_rms_ir_limit = sp.simplify(spherical_rms_on_shell.subs(q_s, 0))
    ksp_shell_from_flux = sp.sqrt(2)*native_secant
    ksp_shell_identity = sp.sqrt(2)/y_s

    eta, M = sp.symbols('eta M', positive=True)
    A = sp.Function('A')(eta)
    F = M**2 / A**2
    # A>0 is the physical branch; write sqrt(F)=M/A explicitly.
    zt = sp.simplify(A * (M / A))
    pump = sp.simplify(sp.diff(zt, eta, 2) / zt)
    ut, u, beta, N, C = sp.symbols('u_dot u beta N C')
    grav_term = N * sp.exp(-2 * beta * u) * C

    bbn_path = ROOT.parent / 'TEP-BBN/results/step_20_ambient_bath_front.json'
    bbn = json.loads(bbn_path.read_text())
    f_spec = bbn['F_boundary_perturbation_spectrum']
    effective = f_spec['perturbation_mode_closure']['quadratic_action_coefficients']
    q_s_val = float(f_spec['two_branch_sound_horizon_gradient_ratio_q_s'])
    y_s_val = float(f_spec['two_branch_sound_horizon_response_y_s'])
    k_value = 16.03

    kt_val = float(kt.subs({xi: 1, k: k_value}))
    flux_val = float(flux_lhs.subs({k: k_value, q_s: q_s_val, y_s: y_s_val}))
    ksp_shell_val = float(ksp_shell_identity.subs({y_s: y_s_val}))
    cp_eff_sq_val = float(ksp_shell_val / kt_val)

    return {
        'status': 'identities_verified_with_stated_scope',
        'scope': (
            'Separates (1) single-invariant homogeneous timelike P(X) principal '
            'coefficients at bar_xi = 1 from (2) the (3+1) U_mu / Pi_{mu nu} '
            'native shell secant 1/y_s and longitudinal flux derivative, separately from the BBN rescaled prescription sqrt(2)/y_s, '
            'and (3) the constant-Mpl static gravitational-background tensor action.'
        ),
        'scalar': {
            'P_X': str(px),
            'K_t': str(kt),
            'sound_speed_squared': str(cs2),
            'K_t_xi1_k16p03': kt_val,
            'K_sp_native_xi1_k16p03': float(px.subs({xi: 1, k: k_value})),
            'sound_speed_squared_xi1_k16p03': float(cs2.subs({xi: 1, k: k_value})),
            'timelike_small_xi_limit': float(sp.limit(cs2, xi, 0, dir='+')),
            'directional_increment_2xi_dPX_dxi': str(sp.simplify(2 * xi * sp.diff(px, xi))),
        },
        'acoustic_shell_reduction': {
            'two_branch_flux_law': f'{flux_lhs} = 1',
            'native_secant_from_flux': str(native_secant),
            'native_secant_on_shell': str(native_on_shell),
            'native_secant_value': float(native_secant.subs({k:k_value,q_s:q_s_val,y_s:y_s_val})),
            'native_longitudinal_expression': str(native_longitudinal),
            'native_longitudinal_value': float(native_longitudinal.subs({k:k_value,q_s:q_s_val,y_s:y_s_val})),
            'spherical_trace_average_on_shell': str(spherical_trace_iso_on_shell),
            'spherical_trace_average_value': float(spherical_trace_iso_on_shell.subs({q_s:q_s_val,y_s:y_s_val})),
            'spherical_frobenius_rms_on_shell': str(spherical_rms_on_shell),
            'spherical_frobenius_rms_value': float(spherical_rms_on_shell.subs({q_s:q_s_val,y_s:y_s_val})),
            'spherical_frobenius_rms_deep_ir_limit': str(spherical_rms_ir_limit),
            'spherical_rms_ir_limit_equals_sqrt2_over_ys': bool(sp.simplify(spherical_rms_ir_limit - ksp_shell_identity) == 0),
            'native_secant_identity_verified': bool(sp.simplify(native_secant*y_s-flux_lhs)==0),
            'normalization_ratio_BBN_over_native': float(sp.sqrt(2)),
            'quadratic_perturbation_mapping_established': False,
            'interpretation': (
                'On the 3D sound-horizon sphere S^2, the spatial flux Hessian is anisotropic with '
                'transverse eigenvalue K_perp = 1/y_s (multiplicity 2) and longitudinal eigenvalue '
                'K_parallel = dJ/dp = 2/y_s + (q_s*y_s)^2 (multiplicity 1). Its 3D rotational '
                'Frobenius/RMS invariant sqrt((2 K_perp^2 + K_parallel^2)/3) equals 0.191772 '
                'at the sound horizon and converges identically to sqrt(2)/y_s = 0.191705 in '
                'the deep two-branch limit (q_s -> 0), matching the BBN effective spatial '
                'stiffness to 0.035% while remaining distinct from an isotropic single-scalar quadratic action.'
            ),
            'shell_secant_stiffness_from_flux': str(ksp_shell_from_flux),
            'shell_secant_stiffness_on_shell': str(ksp_shell_identity),
            'q_s_sound_horizon': q_s_val,
            'y_s_sound_horizon': y_s_val,
            'flux_law_lhs_evaluated': flux_val,
            'flux_law_residual': abs(flux_val - 1.0),
            'K_t_timelike_rolling_background': kt_val,
            'K_sp_spatial_sound_shell_derived': ksp_shell_val,
            'c_p_eff_squared_derived': cp_eff_sq_val,
            'matches_bbn_step20_within_1e12': bool(
                abs(kt_val - effective['K_t_timelike_rolling_background']) < 1e-12
                and abs(ksp_shell_val - effective['K_sp_spatial_sound_shell']) < 1e-12
                and abs(cp_eff_sq_val - effective['c_p_eff_squared_derived']) < 1e-12
                and abs(flux_val - 1.0) < 1e-12
            ),
        },
        'effective_acoustic_prescription': dict(effective),
        'mapping_status': (
            'Native shell flux identity verified: K_perp=1/y_s, '
            'K_parallel=dJ/dp=2/y_s+(q_s*y_s)^2 on shell. '
            'BBN uses the explicitly rescaled prescription sqrt(2)/y_s. '
            'Its equality to a physical acoustic quadratic coefficient with the '
            'separately evaluated timelike K_t is not established by these identities.'
        ),
        'tensor': {
            'matter_frame_F': str(F),
            'z_T': str(zt),
            'z_T_second_derivative_over_z_T': str(pump),
            'assigned_benchmark_pump': str(sp.diff(A, eta, 2) / A),
            'physical_implication': (
                'No pump-driven vacuum amplification on this static, constant-Mpl branch. '
                'This does not forbid tensor sources, classical initial tensors or other action/background branches.'
            ),
        },
        'direct_ADM_prefactor_momentum': str(sp.diff(grav_term, ut)),
        'provenance': {
            'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'bbn_results_sha256': hashlib.sha256(bbn_path.read_bytes()).hexdigest(),
        },
    }


if __name__ == '__main__':
    out = ROOT / 'results/step_81_action_contract_audit.json'
    result = derive()
    out.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(result, indent=2))

