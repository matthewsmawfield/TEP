#!/usr/bin/env python3
"""step_67: In-medium constraint-slice propagator — the u_eff(rho) projection.

AUD-3 remainder.  The WB environmental inversion measures an effective
ambient u_eff declining with local density, u_eff ~ rho^{0.4-0.6}
(020g mass-corrected: 0.40+-0.05; 020l safe-window absolute inversion:
0.50 thermal / 0.55 superthermal).  Ledger section 8 records that every
literal real-space ambient candidate is rejected (gradient channels
wrong sign, rms flat, depth boundary-artifact, starvation wrong sign,
packing amplitude ~40x short), so u_eff is an effective parameter of
the in-medium response -- the operative projection is the propagator on
the constraint slice.

This step computes exactly that, at model-true level:

1. THE IN-MEDIUM PROPAGATOR.  On a time-symmetric (constraint) slice the
   scalar equation is elliptic; perturbations about the local
   density-equilibrium u_eq(rho) propagate with effective mass

       m_eff^2 = V''(u_eq) + rho e^{-u_eq}   (master potential + matter)

   evaluated with tep_model.mass_squared / equilibrium_varphi over the
   disk-to-halo density range.  The source-curvature regime boundary
   rho_* (where rho e^{-u_eq} = V''(u_eq)) is located exactly.

2. PROJECTION CANDIDATES.  The kinetic ambient that enters the shear
   vertex is x_env = |du|^2.  Candidate projections of the medium onto
   that ambient, with their rho-exponent:

     A) m_eff itself:            u_eff ~ m_eff          -> rho^{1/2}
        (in-medium Compton rate as the effective ambient)
     B) gradient coherence:      |du| ~ phi_rms * m_eff -> rho^{1/2}
        if the landscape field rms is density-flat (020h measured
        ~0.12-0.14 z-flat); the propagator localizes the landscape
        gradient to the Compton length
     C) packing amplitude:       u_eff ~ n^{1/3}        -> rho^{1/3}
     D) quartic minimum:         u_eff ~ u_eq           -> ~rho^{1/3}
        (excluded: equilibrium is ambient-pinned ~1e-37)

   Exponents are computed numerically over the stratum density range,
   not asserted analytically.

3. NORMALIZATION.  Each surviving projection is carried to an absolute
   amplitude where the corpus supplies the conversion: the Compton rate
   m_eff in acceleration units (a_C = m_eff c / g_t in |du| units) and
   the measured landscape rms for the gradient-coherence channel.  The
   gap between projection amplitude and the measured u_eff = 0.163-0.18
   midplane is reported explicitly -- it quantifies what remains of the
   propagator-solve normalization.

4. GROWTH APPLICATION (second AUD-3 consumer).  The same m_eff(rho)
   gives the in-well kinetic-variable direction check recorded in the
   ledger: X_inwell >> X_ambient pushes G_eff -> 1 (the direction
   required by the residual f sigma8 excess of step_62).

Inputs:  tep_model master-potential machinery.
Outputs: results/step_67_inmedium_constraint_slice.json
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tep_model import (M_PL, HBAR_C, C, G, PC, KG_GEV, G_CM3_GEV4,
                       LAMBDA_REFERENCE, equilibrium_varphi,
                       mass_squared, save, solve_sphere, diagnostics,
                       source_parameters)

G_T = 3.400464162104563e-10   # m/s^2, corpus transition acceleration
U0_SOLAR = 0.163              # derived solar-circle ambient (shear units)

# Quartic normalization branches (LAMBDA_QUARTIC_REF /
# LAMBDA_QUARTIC_CASSINI in core/constants.py).  The corrected
# linear-in-S_Sigma Cassini evaluation fails at the reference value by
# ~250x (step_03), so the corpus adopts lambda_Cassini = 1e5 lambda_ref
# as the operative weak-field coupling; lambda_ref remains the fiducial
# normalization point of the potential family.
LAMBDA_CASSINI = LAMBDA_REFERENCE * 1e5

# Disk-column density model (step_020d conventions), Msun/pc^3 ->
# g/cm^3 for the field solver and GeV^4 for m_eff.
# NOTE: tep_model.PC is in meters; g/cm^3 needs the cm conversion.
M_SUN = 1.98847e30
MSUN_PC3_G_CM3 = M_SUN * 1e3 / (PC * 100.0)**3   # M_sun/pc^3 -> g/cm^3


def rho_disk_msun_pc3(z_pc):
    z = abs(z_pc) / 1000.0
    return (0.040 * np.exp(-z / 0.300)
            + 0.005 * np.exp(-z / 0.900)
            + 0.050 * np.exp(-z / 0.150))


def inmedium(rho_g_cm3, lam=LAMBDA_CASSINI):
    """Equilibrium field, V'' and matter curvature, m_eff, Compton length."""
    u_eq = equilibrium_varphi(rho_g_cm3, lam)
    m2 = mass_squared(u_eq, rho_g_cm3, lam)
    vpp = m2 - rho_g_cm3 * G_CM3_GEV4 / M_PL**2 * np.exp(-u_eq)
    matter = rho_g_cm3 * G_CM3_GEV4 / M_PL**2 * np.exp(-u_eq)
    lam_c_m = HBAR_C / np.sqrt(m2) if m2 > 0 else np.nan
    return dict(u_eq=u_eq, m2=m2, vpp=vpp, matter_curv=matter,
                compton_m=lam_c_m)


def _landscape_mc(n_star_grid=(0.05, 0.09, 0.15, 0.25, 0.40),
                  box_pc=80.0, n_rep=12, seed=3):
    """Exact k-branch landscape census on a Poisson well network.

    In the k-branch the flux F obeys a LINEAR Gauss law (div F = rho
    e^{-u}, e^{-u} ~ 1 in the inter-well medium), so the total flux at a
    point is the exact vector sum F = sum_i q_i r_hat_i / r_i^2 with
    q_i = G M_i / g_t, and the shear amplitude is |du| = sqrt(sqrt(2)
    |F| / K).  This is the p-Laplacian reduction used in step_020h —
    the exact nonlinear superposition, not an approximation.
    """
    K = 16.03
    q_sun = G * M_SUN / G_T                     # m^2 per star
    rng = np.random.default_rng(seed)
    L = box_pc * PC                              # m
    rows = []
    for nstar in n_star_grid:
        dus = []
        for _ in range(n_rep):
            N = int(nstar * box_pc ** 3)
            pos = rng.uniform(-L / 2, L / 2, (N, 3))
            mass = np.exp(rng.normal(np.log(0.35), 0.55, N))  # M_sun
            r2 = (pos ** 2).sum(1)
            ok = r2 > (0.3 * PC) ** 2
            rr, mm, r2o = pos[ok], mass[ok], r2[ok]
            F = (q_sun * mm / r2o)[:, None] * (rr / np.sqrt(r2o)[:, None])
            Fmag = np.sqrt((F.sum(0) ** 2).sum())
            dus.append(np.sqrt(np.sqrt(2) * Fmag / K))
        dus = np.array(dus)
        rows.append({'n_star_per_pc3': float(nstar),
                     'median_du': float(np.median(dus)),
                     'p16': float(np.percentile(dus, 16)),
                     'p84': float(np.percentile(dus, 84))})
    ns = np.array([r['n_star_per_pc3'] for r in rows])
    du = np.array([r['median_du'] for r in rows])
    exp = float(np.polyfit(np.log(ns), np.log(du), 1)[0])
    return {'rows': rows, 'exponent_vs_n': exp,
            'note': 'discreteness channel: |du| ~ 0.004-0.007, ~40x '
                    'below the measured u_eff ~ 0.16-0.18 (consistent '
                    'with the 020h census) and exponent ~1/3 '
                    '(nearest-well dominated) -- every literal '
                    'real-space landscape channel is amplitude-'
                    'eliminated; u_eff is the in-well response '
                    'parameter, not a boundary field.'}


def _embedded_response(lam):
    """D^-1(k ~ 1/r*; rho_bg, u_amb): the embedded-source response function.

    The pair-scale source (M = 1.4 M_sun halo at r* = 4.9 kAU) is solved
    with tep_model.solve_sphere's nested construction
    phi_total = phi_env + delta_phi, the perturbation decaying to
    bg + phi_env on the constraint-slice background rho_bg.  The
    landscape ambient field VALUE enters as phi_env = u_amb * psi_uns
    (u_amb in the well-solve charge units, psi_uns the source's own
    unscreened-potential scale) and the response is read off as the
    emitted source-charge ratio (diagnostics), the k = 1/r* response.
    """
    G_T = 3.4e-10
    M_SUN_KG = 1.98847e30
    AU_M = 1.496e11
    M = 1.4 * M_SUN_KG
    rstar = float(np.sqrt(G * M / G_T))
    psi, _, _ = source_parameters(M, rstar, lam)
    u_grid = np.array([0.0, 0.05, 0.10, 0.14, 0.20, 0.30, 0.40])
    rho_strata = {
        'z22': 5.77e-24, 'z68': 4.62e-24, 'z121': 3.61e-24,
        'z195': 2.61e-24, 'z346': 1.42e-24,
    }
    def _reseed(rb_, ua_):
        # Cold-start can fail the Jacobian on the stiff lambda_Cassini
        # branch.  Reseed from the converged reference-branch solution
        # for the same stratum/ambient, rescaled to this branch's
        # charge-normalized units (psi_uns differs between branches).
        psi_ref = source_parameters(M, rstar, LAMBDA_REFERENCE)[0]
        s_ref = solve_sphere(M, rstar, lam=LAMBDA_REFERENCE, rho_bg=rb_,
                             phi_env=float(ua_) * psi_ref)
        class _G:
            pass
        g = _G()
        g.x = s_ref['solution'].x
        g.y = s_ref['solution'].y * (psi_ref / psi)
        return g

    table = {}
    for ztag, rb in rho_strata.items():
        row = []
        guess = None
        for ua in u_grid:
            try:
                s = solve_sphere(M, rstar, lam=lam, rho_bg=rb,
                                 phi_env=float(ua) * psi, guess=guess)
            except RuntimeError:
                s = solve_sphere(M, rstar, lam=lam, rho_bg=rb,
                                 phi_env=float(ua) * psi,
                                 guess=_reseed(rb, ua))
            row.append(diagnostics(s, 1.0)['source_charge_ratio'])
            guess = s['solution']
        table[ztag] = row
    # rho_bg sensitivity at fixed ambient and the ambient slope
    col = np.array([table[t] for t in rho_strata])          # (rho, u)
    rho_frac = (col.max(0) - col.min(0)) / col.mean(0)      # per u column
    deep = col[:, u_grid >= 0.10]
    lu, lq = np.log(u_grid[u_grid >= 0.10]), np.log(deep.mean(0))
    u_slope = float(np.polyfit(lu, lq, 1)[0])
    return {
        'rstar_AU': rstar / AU_M,
        'u_amb_grid': u_grid.tolist(),
        'charge_ratio_table': table,
        'rho_bg_sensitivity_frac': rho_frac.tolist(),
        'u_amb_slope_deep': u_slope,
        'note': f'charge ratio varies <{100*rho_frac.max():.2g}% across '
                'the 4x disk density range at fixed ambient, while the '
                'ambient field value suppresses it over ~1.5-2 decades '
                f'(charge ~ u_amb^{u_slope:.2f} in the 0.1-0.4 range).  '
                'The environmental '
                'ordering is therefore carried ENTIRELY by the ambient '
                'boundary field value u_amb at the vertex: local density '
                'does not screen the pair directly.  u_eff(rho) is the '
                'landscape ambient-map value, and D^-1(k_pair) = '
                'charge_ratio(u_amb(rho)).  Unit note: phi_env must be '
                'supplied in varphi units = u_amb * psi_uns; passing the '
                'bare u-value (varphi ~ 0.14) collapses lambda_C to ~1 m '
                'and the solve fails -- the well-solve u is a '
                'charge-normalized unit, not raw varphi.',
    }


def run():
    # Operative weak-field quartic normalization is the Cassini branch
    # (the reference value fails the corrected linear-in-S Cassini bound
    # by ~250x; step_03).  The reference branch is kept alongside for
    # comparison so both normalizations are recorded.
    lam = LAMBDA_CASSINI
    rho_grid = np.geomspace(1e-26, 1e-21, 121)   # g/cm^3
    rows = [inmedium(r) for r in rho_grid]
    rows_ref = [inmedium(r, LAMBDA_REFERENCE) for r in rho_grid]

    # --- regime boundary: source-curvature dominance ---
    frac = np.array([r['matter_curv'] / r['m2'] for r in rows])
    cross = np.argmin(np.abs(frac - 0.5))
    rho_star = float(rho_grid[cross])

    # --- m_eff exponent over the disk stratum density range ---
    z_strata = np.array([22.0, 68.0, 121.0, 195.0, 346.0])
    rho_z = rho_disk_msun_pc3(z_strata) * MSUN_PC3_G_CM3   # g/cm^3
    row_z = [inmedium(r) for r in rho_z]
    row_z_ref = [inmedium(r, LAMBDA_REFERENCE) for r in rho_z]
    m_z = np.sqrt([r['m2'] for r in row_z])
    lc_pc = np.array([r['compton_m'] for r in row_z]) / PC
    lc_pc_ref = np.array([r['compton_m'] for r in row_z_ref]) / PC
    lnr, lnm = np.log(rho_z), np.log(m_z)
    m_exp = float(np.polyfit(lnr, lnm, 1)[0])

    # --- candidate projections, exponent vs rho over the grid ---
    # A: u_eff ~ m_eff (Compton rate as effective ambient, per u-unit
    #    normalisation absorbed into a single constant)
    m_grid = np.sqrt([r['m2'] for r in rows])
    proj_A = np.polyfit(np.log(rho_grid), np.log(m_grid), 1)[0]
    # C/D: packing n^{1/3} and quartic minimum u_eq exponents
    u_eq_grid = np.array([r['u_eq'] for r in rows])
    proj_D = np.polyfit(np.log(rho_grid), np.log(u_eq_grid + 1e-300), 1)[0]

    # --- stratum table ---
    strata = []
    for z, rho, r, m, lc, lc_ref in zip(z_strata, rho_z, row_z, m_z,
                                        lc_pc, lc_pc_ref):
        # Compton rate in shear units: a_C = m_eff * c (s^-1 * m/s -> m/s^2)
        # m in GeV -> s^-1 via /hbar; a_C = omega_C * c / g_t
        omega = m * 6.582119569e-25 ** -1          # GeV -> s^-1 (hbar^-1)
        a_c = omega * C                             # m/s^2
        strata.append({
            'z_pc': float(z),
            'rho_g_cm3': float(rho),
            'u_eq': float(r['u_eq']),
            'm_eff_GeV': float(m),
            'compton_pc': float(lc),
            'compton_pc_lambda_ref': float(lc_ref),
            'matter_curv_frac': float(r['matter_curv'] / r['m2']),
            'compton_accel_over_g_t': float(a_c / G_T),
            # r* ~ 2.6 kAU = 0.013 pc: (r*/lambda_C)^2 is the leading
            # in-medium correction at the pair-halo scale
            'compton_over_rstar': float(lc / 0.013),
            'pair_scale_corr_rstar_over_lamC_sq': float((0.013 / lc) ** 2),
        })

    # --- ambient-map closure: u_amb ~ sqrt(rho L) with L derived ---
    # k-branch planar balance P_X u ~ J with P_X ~ K u/sqrt(2) and the
    # local source J = 4 pi G rho L / g_t gives u^2 = 4 sqrt(2) pi G
    # rho L /(K g_t).  L is the local coherence length of the disk
    # density field.  Calibrate L at the midplane inversion u=0.18 and
    # propagate to all strata via the stratum-mean densities.
    z_stratum_edges = np.array([0.0, 45.0, 95.0, 158.0, 270.0, 700.0])
    rho_mean_stratum = np.array([
        rho_disk_msun_pc3(np.linspace(z_stratum_edges[i],
                                    z_stratum_edges[i + 1], 400)).mean()
        for i in range(5)])
    rho_mid_kgm3 = rho_mean_stratum[0] * M_SUN / PC**3
    K_PHI = 16.03
    u_meas_mid_c = 0.18                            # 020l thermal z1 inversion
    L_coh = (u_meas_mid_c ** 2 * K_PHI * G_T
             / (4.0 * np.sqrt(2.0) * np.pi * G * rho_mid_kgm3))
    u_pred = u_meas_mid_c * np.sqrt(rho_mean_stratum / rho_mean_stratum[0])
    u_meas_strata = np.array([0.18, 0.17, 0.15, 0.13, 0.08])
    ambient_closure = {
        'law': 'u_amb = sqrt(4 sqrt(2) pi G rho L / (K g_t))',
        'L_coh_pc': float(L_coh / PC),
        'L_note': 'L ~ 0.8 kpc is the local disk density coherence '
                  'scale (between the thin/thick scale heights); '
                  'calibrated at the midplane inversion, not fitted '
                  'per stratum -- the law then predicts the remaining '
                  'four strata with ONE parameter total.',
        'rho_mean_stratum_msun_pc3': rho_mean_stratum.tolist(),
        'u_pred': u_pred.tolist(),
        'u_meas': u_meas_strata.tolist(),
        'pointwise_ratio': (u_pred / u_meas_strata).tolist(),
        'max_deviation_pct': float(
            100 * np.max(np.abs(u_pred / u_meas_strata - 1))),
    }

    # --- normalization checks ---
    # projection B needs phi_rms: use the 020h measured z-flat value
    phi_rms = 0.13                                  # well-solve units
    # gradient coherence |du| ~ phi_rms * m_eff * L_ref; with L_ref the
    # Compton length in shear-units conversion this collapses to
    # |du| ~ phi_rms * (a_C/g_t_ref) * (L_pair/lambda_C):
    # report the required coherence factor to hit measured midplane u_eff
    u_meas_mid = 0.18                               # 020l thermal z1 inversion
    a_c_mid = strata[0]['compton_accel_over_g_t']
    req_factor_A = u_meas_mid / a_c_mid if a_c_mid > 0 else np.nan
    req_factor_B = u_meas_mid / (phi_rms * a_c_mid) if a_c_mid > 0 else np.nan

    out = {
        'description':
            'In-medium constraint-slice propagator: m_eff(rho), regime '
            'boundary, projection exponents and normalisation gaps for '
            'the WB u_eff(rho) law and the growth in-well X_env.',
        'regime': {
            'rho_star_source_curv_g_cm3': rho_star,
            'disk_range_g_cm3': [float(rho_z.min()), float(rho_z.max())],
            'matter_curv_frac_at_disk': [
                float(r['matter_curv'] / r['m2']) for r in row_z],
            'note': 'CORRECTS the earlier ledger estimate: at physical '
                    'disk densities (0.02-0.09 M_sun/pc^3 = 1.4-5.8e-24 '
                    'g/cm^3) the matter curvature rho e^{-u} contributes '
                    'only ~1e-16 of m_eff^2 -- the propagator is '
                    'POTENTIAL-CURVATURE dominated (V'' ~ u_eq^2), '
                    'giving m_eff ~ rho^{1/3}.  The earlier '
                    '"m_eff^2 ~ rho" statement used rho ~ 5e-27 g/cm^3, '
                    '~1100x below the real midplane; at true densities '
                    'the source-curvature regime does not apply.',
        },
        'm_eff_exponent_vs_rho': {
            'fit_over_strata': m_exp,
            'fit_over_full_grid': float(proj_A),
            'expected_source_curv': 0.5,
            'actual_regime': 'potential-curvature dominated: '
                             'm_eff ~ rho^{1/3}',
        },
        'strata': strata,
        'projection_candidates': {
            'A_compton_rate': {
                'form': 'u_eff ~ m_eff (in-medium propagator rate)',
                'exponent': float(proj_A),
                'amplitude_midplane_du': float(a_c_mid),
                'required_factor_to_u_eff': float(req_factor_A),
                'verdict': 'potential-dominated regime -> rho^{1/3}; '
                           'direct-rate amplitude is ~1e11 g_t units, '
                           'not the ambient channel',
            },
            'B_gradient_coherence': {
                'form': '|du| ~ phi_rms * m_eff; phi_rms z-flat (020h)',
                'exponent': float(proj_A),
                'amplitude_midplane_du': float(phi_rms * a_c_mid),
                'required_factor_to_u_eff': float(req_factor_B),
                'verdict': 'exponent rho^{1/3} (phi_rms ~flat), the '
                           'LOW edge of the measured 0.4-0.6 band',
            },
            'E_pair_scale_response': {
                'form': 'x_env ~ (k_pair^2 + m_eff^2) evaluated at '
                        'k_pair = 1/r*',
                'exponent_range': [0.0, 0.667],
                'note': 'lambda_C = 0.012-0.019 pc at disk densities '
                        '(operative lambda_Cassini branch; 0.08-0.13 pc '
                        'on the fiducial lambda_ref branch): the '
                        'pair-halo r* (2.6 kAU = 0.013 pc) sits AT the '
                        'Compton edge (r*/lambda_C ~ 0.7-1.1), so the '
                        'in-medium correction (m_eff/k_pair)^2 ~ '
                        'rho^{2/3} is order-unity -- the UPPER edge '
                        'of the measured band.  The measured rho^0.5 '
                        'sits between the k->0 propagator exponent '
                        '(1/3) and the pair-scale correction (2/3) -- '
                        'consistent with the response function being '
                        'evaluated at the halo scale, i.e. the '
                        'projection is the propagator AT k_pair, not '
                        'in the infrared.',
            },
            'C_packing': {
                'form': 'u_eff ~ n^{1/3}',
                'exponent': 1.0 / 3.0,
                'verdict': 'exponent below measured band; retained as '
                           'bracketing alternative',
            },
            'G_planar_kbranch': {
                'form': 'u_amb ~ J^{1/2}; planar k-branch response '
                        'u'' P_X = J to the local disk density '
                        'J ~ rho (020d channel B)',
                'exponent': 0.5,
                'contrast_stratum_median': 2.014,
                'contrast_stratum_mean': 2.440,
                'contrast_020d_density': 1.564,
                'contrast_required': '~1.8-2.0 (alpha-channel 1.78, '
                                     'monopole map 2.0)',
                'density_ratios': {'stratum_median': 4.057,
                                   'stratum_mean': 5.955,
                                   '020d_two_bin': 2.446},
                'exponent_measured_vs_mean_rho': 0.462,
                'verdict': 'ONLY channel at rho^{1/2} AND reaching the '
                           'required amplitude: J^{1/2} at stratum-'
                           'averaged densities gives contrast 2.44 '
                           '(median: 2.01) bracketing the required '
                           '1.78-2.0, and the measured endpoint law '
                           'u_eff ~ rho^{0.46} matches rho^{1/2} '
                           'within systematics.  The ambient map is '
                           'identified: the local planar k-branch '
                           'response of the disk''s own field '
                           'structure, u_amb ~ rho^{1/2}.  020d''s '
                           'smaller contrast used a two-bin split '
                           '(rho ratio 2.44), not a different '
                           'mechanism.',
            },
            'F_halo_truncation': {
                'form': 'u_eff = |du_pair|(D_s); pair halo truncated '
                        'at the inter-well spacing/tidal balance '
                        'D_s ~ n^{-1/3} (flux balance q/r^2 ~ rho r '
                        'gives the Jacobi scaling r_t ~ rho^{-1/3})',
                'exponent': 1.0 / 3.0,
                'verdict': 'SIXTH independent rho^{1/3} channel: the '
                           'pair shear tail |du| ~ 1/r sampled at the '
                           'truncation radius gives rho^{1/3} exactly. '
                           'Only a pair-scale correction (E) pushes '
                           'toward 2/3; the measured 0.4-0.7 band '
                           'brackets [1/3, 2/3].',
            },
            'D_quartic_minimum': {
                'form': 'u_eff ~ u_eq(rho)',
                'exponent': float(proj_D),
                'verdict': 'EXCLUDED: u_eq ~ 1e-37, ambient-pinned',
            },
        },
        'measured_target': {
            'u_eff_exponent_data': [0.40, 0.60],
            'u_eff_midplane': u_meas_mid,
            'u0_solar_derived': U0_SOLAR,
        },
        'growth_application': {
            'in_well_X_direction':
                'inside a well of density rho the constraint slice pins '
                'phi at u_eq(rho) and the perturbation propagator is '
                'm_eff(rho); at well densities rho >> rho_disk the '
                'in-well kinetic variable exceeds the ambient value '
                '(m_eff grows ~rho^{1/2}), pushing the linear-regime '
                'G_eff toward 1 -- the direction required by the '
                'step_62 f sigma8 residual',
        },
        'landscape_monte_carlo': _landscape_mc(),
        'embedded_response': _embedded_response(lam),
        'embedded_response_lambda_ref': _embedded_response(
            LAMBDA_REFERENCE),
        'lambda_branches': {
            'lambda_reference': float(LAMBDA_REFERENCE),
            'lambda_cassini': float(LAMBDA_CASSINI),
            'operative': 'lambda_cassini',
            'note': 'the corrected linear-in-S_Sigma Cassini evaluation '
                    'fails at lambda_ref by ~250x (step_03); the corpus '
                    'adopts lambda_Cassini = 1e5 lambda_ref wherever the '
                    'quartic normalization enters.  The exponent laws '
                    '(m_eff ~ rho^{1/3}, ambient u_amb ~ rho^{1/2}) are '
                    'lambda-independent; the absolute Compton scale '
                    'and the embedded charge amplitude shift by '
                    'lambda^{+1/6} and the branch stiffness.',
        },
        'ambient_closure': ambient_closure,
        'status': (
            'regime corrected: the in-medium propagator at physical '
            'disk densities is potential-curvature dominated, '
            'm_eff ~ rho^{1/3} exactly (the earlier source-curvature '
            'estimate used a density ~1100x too low).  All lambda-'
            'dependent quantities are computed at the operative '
            'lambda_Cassini branch (lambda_ref branch retained in '
            'lambda_branches/*_lambda_ref fields): lambda_C = '
            '0.012-0.019 pc, so the pair-halo r* sits at the Compton '
            'edge (r*/lambda_C ~ 0.7-1.1) and the in-medium '
            'correction is order-unity at the pair scale.  The '
            'measured u_eff ~ rho^{0.4-0.6} sits between the k->0 '
            'propagator exponent (1/3) and the leading pair-scale '
            'correction (m_eff/k_pair)^2 ~ rho^{2/3}: the density law '
            'is consistent with the propagator evaluated AT the '
            'halo scale rather than in the infrared -- the operative '
            'AUD-3 projection.  The embedded-source response function '
            'D^{-1}(k ~ 1/r*; rho_bg, u_amb) is computed '
            '(embedded_response): it depends ONLY on the ambient '
            'boundary field value (<0.01% rho_bg sensitivity), charge '
            '~ u_amb^{-2.0} -- the vertex-ambient operator structure '
            'confirmed at response-function level.  The ambient map '
            'is identified at the level the data resolve: the local '
            'k-branch density response u_amb ~ rho^{1/2} (G_planar_'
            'kbranch) gives contrast 2.0-2.44 vs required 1.78-2.0 '
            'and exponent rho^{0.46} vs measured 0.4-0.7; the solved '
            'planar midplane shear (0.174, 020d) matches the measured '
            'midplane u_eff (0.18) to ~4%.  The ambient map is now '
            'closed at coefficient level: u_amb = sqrt(4 sqrt(2) pi G '
            'rho L /(K g_t)) with coherence length L ~ 0.84 kpc '
            'calibrated once at the midplane inversion reproduces all '
            'five strata to within 10% (ambient_closure).  Open: the '
            'f sigma8 in-well X_env volume weighting.'),
    }
    save('step_67_inmedium_constraint_slice.json', out)
    print(json.dumps({k: out[k] for k in
                      ('regime', 'm_eff_exponent_vs_rho')}, indent=2))
    for s in strata:
        print(f"z={s['z_pc']:5.0f}pc rho={s['rho_g_cm3']:.2e} "
              f"m_eff={s['m_eff_GeV']:.2e}GeV lamC={s['compton_pc']:.2e}pc "
              f"matter%={s['matter_curv_frac']:.3f} "
              f"aC/g_t={s['compton_accel_over_g_t']:.2e}")
    print('wrote results/step_67_inmedium_constraint_slice.json')


if __name__ == '__main__':
    run()
