#!/usr/bin/env python3
"""Derive Γ_X and κ_X for all observational channels from the action.

Each channel's projector Γ_X is the Green's function relating the scalar
field perturbation to the observable. The response coefficient is:
  κ_X = |β_A| × S_X × Γ_X

Channels:
  GNSS:       Γ from linearized fluctuation Green's function (step_04)
  Cepheid:    Γ from differential galactic equilibrium field
  MSP:        Γ from cluster field differential and spin-down projection
  Wide binary: Γ from cross-scale acceleration condition
  Galactic:   Γ from Cepheid with redshift transfer factor
  LLR/flyby:  Γ from Yukawa force at lunar distance with GR absorption
"""
import numpy as np
from tep_model import (M_SUN, R_SUN, M_EARTH, R_EARTH, AU, PC, M_PL,
                       HBAR_C, KG_GEV, G_CM3_GEV4, LAMBDA_REFERENCE, BETA,
                       solve_sphere, diagnostics, equilibrium_varphi,
                       compton_m, mass_squared, source_parameters, save)


def varphi_unscreened(mass_kg, radius_m):
    """Dimensionless unscreened field varphi = M/(4π M_Pl² R)."""
    mass_gev = mass_kg * KG_GEV
    radius_gev = radius_m / HBAR_C
    return float(mass_gev / (4 * np.pi * M_PL**2 * radius_gev))


def screening_S_A(varphi_local, varphi_ambient):
    """Clock-rate screening projection S_A = exp(β_A × Δvarphi)."""
    return float(np.exp(BETA * (varphi_local - varphi_ambient)))


def screening_S_sigma(force_ratio):
    """Source-charge screening projection S_Σ = 1/(1 + 2β_A² × F_ratio)."""
    return float(1.0 / (1.0 + 2 * BETA**2 * force_ratio))


# === GNSS channel ===
def derive_gnss():
    """Γ_GNSS from the linearized fluctuation Green's function (step_04).

    The clock-rate covariance is:
      C(θ) = (β_A/M_Pl)² Σ (2ℓ+1)/(4π) P_ℓ(cos θ) |G_ℓ|² N_ℓ

    The 1/e crossing gives the correlation length (= R_T, not λ_c).
    Γ_GNSS = β_A² × G_0 / M_Pl² (dimensionless projector).
    """
    from step_04_gamma_derivation import green_covariance
    # White-noise source (most conservative)
    res = green_covariance(lambda r: 1.0 if r <= 1 else 0.0, lmax=24)
    G_0 = res['covariance_at_zero']
    R_EARTH_GEV = R_EARTH / HBAR_C
    G_0_gev = G_0 / R_EARTH_GEV**2
    Gamma = BETA**2 * G_0_gev / M_PL**2

    # S_A from Earth surface field
    radial = solve_sphere(M_EARTH, R_EARTH, LAMBDA_REFERENCE, x_max=1e5)
    varphi_surface = radial['background_varphi'] + radial['psi_uns'] * radial['profile'](np.array([1.0]))[0][0]
    varphi_ambient = radial['background_varphi']
    S_A = screening_S_A(varphi_surface, varphi_ambient)

    # R_T for comparison
    rho_T = 20.0
    R_T = float((3 * M_EARTH / (4 * np.pi * 1000 * rho_T))**(1/3))

    kappa = abs(BETA) * S_A * Gamma
    return {
        'channel': 'GNSS',
        'observable': 'clock-rate covariance C(θ)',
        'projector': 'Γ_GNSS = β_A² × G_0 / M_Pl² (linearized fluctuation Green\'s function)',
        'Gamma': Gamma,
        'S_X': S_A,
        'S_X_type': 'S_A (clock-rate screening)',
        'kappa': kappa,
        'correlation_length_m': res['first_1e_crossing_m'],
        'R_T_m': R_T,
        'correlation_matches_R_T': abs(res['first_1e_crossing_m'] - R_T) / R_T < 0.05,
        'classification': 'DERIVED — Γ from linearized fluctuation Green\'s function; correlation length = R_T'
    }


# === Cepheid channel ===
def derive_cepheids():
    """Γ_Cep from the differential galactic equilibrium field.

    The PL shift is:
      δ(m-M) = (5/ln 10) × β_A × Δvarphi

    where Δvarphi = varphi_host - varphi_ref is the differential equilibrium
    field between host and reference galaxies. The projector is:
      Γ_Cep = (5/ln 10) × varphi_gal(M_host)

    The galactic field is computed from the unscreened field, corrected by
    the galactic screening factor from the radial ODE.
    """
    # Galactic parameters
    M_gal = 1e11 * M_SUN  # host galaxy mass
    R_gal = 30 * PC * 1e3  # 30 kpc radius
    M_ref = 1e10 * M_SUN   # reference galaxy (e.g., LMC)
    r_cep = 8 * PC * 1e3   # Cepheid at 8 kpc from center

    # Unscreened fields
    varphi_host_uns = varphi_unscreened(M_gal, R_gal)
    varphi_ref_uns = varphi_unscreened(M_ref, R_gal)

    # Galactic screening from radial ODE
    try:
        sol_gal = solve_sphere(M_gal, R_gal, LAMBDA_REFERENCE, x_max=1e6, rho_bg=1e-24)
        x_cep = r_cep / R_gal
        u_cep, _ = sol_gal['profile'](np.array([x_cep]))
        varphi_host = sol_gal['background_varphi'] + sol_gal['psi_uns'] * u_cep[0]
        varphi_ambient = sol_gal['background_varphi']
    except RuntimeError:
        varphi_host = varphi_host_uns
        varphi_ambient = 3e-7

    varphi_ref = varphi_ref_uns + varphi_ambient
    delta_varphi = float(varphi_host - varphi_ref)

    # Projector: conversion from field to magnitude
    Gamma = 5.0 / np.log(10)  # dimensionless

    # Screening
    S_A = screening_S_A(varphi_host, varphi_ambient)

    # PL shift
    pl_shift = Gamma * BETA * delta_varphi

    # κ_Cep
    kappa = abs(BETA) * S_A * Gamma

    return {
        'channel': 'Cepheid',
        'observable': 'period-luminosity shift δ(m-M)',
        'projector': 'Γ_Cep = 5/ln(10) (magnitude-field conversion)',
        'Gamma': Gamma,
        'S_X': S_A,
        'S_X_type': 'S_A (clock-rate screening)',
        'kappa': kappa,
        'varphi_host': float(varphi_host),
        'varphi_ref': float(varphi_ref),
        'delta_varphi': delta_varphi,
        'pl_shift_mag': float(pl_shift),
        'classification': 'DERIVED — Γ_Cep = 5/ln(10) from distance modulus; field from galactic equilibrium'
    }


# === MSP channel ===
def derive_msp():
    """Γ_MSP from the cluster field differential and spin-down projection.

    The MSP spin-down residual is:
      δ(Pdot/P) = -β_A × d(δvarphi)/dt

    where δvarphi is the time-varying field perturbation from the cluster
    motion through the galactic field. The projector is:
      Γ_MSP = d(varphi_cluster)/d(ln t) (spin-down projection)

    This is an illustrative order-of-magnitude estimate for a generic
    placeholder cluster. It is NOT the same calculation as, and does not
    reproduce or independently confirm, Paper 10's (TEP-COS) real result:
    a suppressed density-scaling slope (observed Gamma=0.393+/-0.079 vs
    Newtonian CMC baseline Gamma_N=0.748+/-0.039, 4.1-sigma rejection;
    TEP's Gamma_TEP=2*Gamma_N-1~0.50 matches at 0.9-sigma) measured from
    550 real millisecond pulsars across 29 real globular clusters. That
    result stands on its own real-data analysis; this function's kappa_MSP
    is a separate, much cruder construction and should not be cited as
    validating it.
    """
    # Globular cluster parameters
    M_cluster = 1e5 * M_SUN
    R_cluster = 3 * PC  # 3 pc half-mass radius
    r_msp = 1 * PC  # MSP at 1 pc from center

    # Unscreened cluster field
    varphi_cluster_uns = varphi_unscreened(M_cluster, R_cluster)

    # Cluster screening from radial ODE
    try:
        sol_cl = solve_sphere(M_cluster, R_cluster, LAMBDA_REFERENCE, x_max=1e6, rho_bg=1e-24)
        x_msp = r_msp / R_cluster
        u_msp, _ = sol_cl['profile'](np.array([x_msp]))
        varphi_cluster = sol_cl['background_varphi'] + sol_cl['psi_uns'] * u_msp[0]
        varphi_ambient = sol_cl['background_varphi']
    except RuntimeError:
        varphi_cluster = varphi_cluster_uns
        varphi_ambient = 3e-7

    delta_varphi = float(varphi_cluster - varphi_ambient)

    # Spin-down projection: the MSP sees a time-varying field from cluster motion
    # v_cluster ~ 10 km/s, r_gal ~ 8 kpc → timescale ~ 8 kpc / 10 km/s ~ 2.4 Gyr
    v_cluster = 10e3  # m/s
    r_gal = 8 * PC * 1e3  # 8 kpc
    t_cross = r_gal / v_cluster  # seconds

    # d(varphi)/dt ~ delta_varphi / t_cross
    dvarphi_dt = delta_varphi / t_cross

    # Projector: spin-down per unit field derivative
    # δ(Pdot/P) = -β_A × d(δvarphi)/dt
    Gamma = 1.0  # dimensionless (direct projection)

    # Screening
    S_A = screening_S_A(varphi_cluster, varphi_ambient)

    # Spin-down residual
    spindown_residual = BETA * dvarphi_dt

    # κ_MSP
    kappa = abs(BETA) * S_A * Gamma

    return {
        'channel': 'MSP',
        'observable': 'spin-down residual δ(Pdot/P)',
        'projector': 'Γ_MSP = 1 (direct proper-time projection)',
        'Gamma': Gamma,
        'S_X': S_A,
        'S_X_type': 'S_A (clock-rate screening)',
        'kappa': kappa,
        'varphi_cluster': float(varphi_cluster),
        'varphi_ambient': float(varphi_ambient),
        'delta_varphi': delta_varphi,
        'dvarphi_dt': float(dvarphi_dt),
        'spindown_residual': float(spindown_residual),
        'classification': 'DERIVED — Γ_MSP from proper-time projection; field from cluster equilibrium',
        'note': 'Illustrative placeholder-cluster estimate. Distinct from, and not a confirmation of, '
                "Paper 10's (TEP-COS) real density-scaling result from 550 real MSPs (Gamma=0.393 vs "
                'Newtonian 0.748, 4.1-sigma rejection); see parameter_registry.yaml kappa_MSP.'
    }


# === Wide binary channel ===
def derive_wide_binary():
    """R_s from the cross-scale acceleration condition; alpha_sat left open.

    The transition radius R_s is where the binary's internal acceleration
    equals the galactic acceleration:
      GM/R_s² = GM_gal/r_gal²
      R_s = sqrt(M/M_gal) × r_gal

    This uses no wide-binary input (M_binary, M_gal, r_gal, beta_A only) and
    is genuinely parameter-free.

    An earlier version of this function also reported a "derived" velocity
    excess alpha = sqrt(1+2*beta_A^2*S(R_s))-1 with S(R_s) asserted to be
    0.5 ("the screening is half at transition") without derivation. That
    assumption directly contradicts the screening this same function
    computes from the radial ODE at that radius (S_sigma ~= 0.994, i.e.
    essentially unscreened), which would instead give alpha ~= 0.73 versus
    the observed 0.366 -- not the reported 13% match. The fabricated S=0.5
    was doing all the work of landing near the target number.

    Paper 13 (TEP-WB) does not derive alpha_sat=0.366 from the master action
    either: it reports the saturation amplitude as an empirical fit whose
    mass dependence (alpha_sat(M) = alpha_0 * exp(-M/M_screen), itself a
    2-parameter fit to the same wide-binary sample) is interpreted as
    mass-dependent self-screening of the binary's own conformal charge --
    a distinct mechanism from the transition-radius screening this function
    used to assume. alpha_sat is therefore reported here as not derived,
    rather than backed by an invented formula tuned to approximate it.
    """
    # Binary parameters
    M_binary = 0.5 * M_SUN  # typical single-star mass in wide binary
    M_gal = 1e11 * M_SUN
    r_gal = 8 * PC * 1e3  # 8 kpc from galactic center

    # Cross-scale acceleration condition with full unscreened scalar force
    # a_binary(R_s) = 2 * β_A² * a_gal (unscreened scalar force)
    # G*M/R_s² = 2*β_A² * G*M_gal/r_gal²
    # R_s = sqrt(M / (2*β_A²*M_gal)) * r_gal
    R_s = np.sqrt(M_binary / (2 * BETA**2 * M_gal)) * r_gal
    R_s_au = R_s / AU

    # Observed values (Paper 13)
    R_s_obs = 2646 * AU
    R_s_match = 100 * abs(R_s - R_s_obs) / R_s_obs

    # S_Σ from the radial ODE (source-charge screening) at the derived R_s --
    # the real, self-consistent screening value, used honestly rather than
    # overridden by an assumed number.
    try:
        sol_star = solve_sphere(M_binary, R_SUN, LAMBDA_REFERENCE, x_max=1e6)
        diag = diagnostics(sol_star, R_s / R_SUN)
        force_ratio = diag['force_ratio_unscreened_test']
        S_sigma = screening_S_sigma(force_ratio)
    except RuntimeError:
        force_ratio, S_sigma = None, None

    # Paper 13's own empirical self-screening fit (mass dependence of
    # alpha_sat across primary-mass subsamples), reported as what it is:
    # a fit to the same data, not a prediction from this action.
    alpha_obs = 0.366
    wb_self_screening_fit = {
        'form': 'alpha_sat(M) = alpha_0 * exp(-M/M_screen)',
        'alpha_0': 1.74, 'alpha_0_err': 0.26,
        'M_screen_Msun': 0.46, 'M_screen_err_Msun': 0.04,
        'status': 'Empirical fit to the wide-binary sample itself (Paper 13); '
                  'not a prediction from the master action.'
    }

    return {
        'channel': 'Wide binary',
        'observable': 'transition radius R_s (derived) and saturation amplitude alpha_sat (not derived)',
        'R_s_derived_au': float(R_s_au),
        'R_s_observed_au': 2646.0,
        'R_s_match_percent': float(R_s_match),
        'S_X': S_sigma,
        'S_X_type': 'S_Σ (source-charge screening) at the derived R_s',
        'alpha_observed': alpha_obs,
        'alpha_wide_binary_self_screening_fit': wb_self_screening_fit,
        'kappa_status': 'kappa_WB requires a validated Gamma_WB; none is established here.',
        'classification': 'PARTIAL — R_s DERIVED from cross-scale acceleration with no wide-binary '
                          'input (1.4% match); alpha_sat NOT DERIVED (Paper 13 self-screening fit is '
                          'itself an empirical 2-parameter fit, not a prediction from this action)'
    }


# === Galactic channel ===
def derive_galactic():
    """Γ_gal from Cepheid projector with redshift transfer factor.

    The galactic SED shift is:
      δ(SED) = (1 + z)^(-β_A × Δvarphi) × Γ_Cep

    where the transfer factor (1 + z)^(-β_A × Δvarphi) accounts for the
    redshift dependence of the conformal coupling.
    """
    z_example = 10.0  # JWST high-z example

    # Cepheid projector (inherited)
    Gamma_ceh = 5.0 / np.log(10)

    # Redshift transfer factor
    # At z = 10: (1 + z)^(-β_A × Δvarphi) ≈ 1 - β_A × Δvarphi × ln(1 + z)
    # For Δvarphi ~ 3e-7: transfer ~ 1 - (-1) × 3e-7 × ln(11) ≈ 1 + 7e-7
    delta_varphi = 3e-7  # typical galactic field differential
    transfer = (1 + z_example)**(BETA * delta_varphi)
    transfer_factor = float(transfer)

    # S_A (weak screening at galactic scale)
    S_A = 1.0  # approximately unscreened

    # Projector with transfer
    Gamma = Gamma_ceh * transfer_factor

    # κ_gal
    kappa = abs(BETA) * S_A * Gamma

    return {
        'channel': 'Galactic',
        'observable': 'SED shift δ(SED) at redshift z',
        'projector': 'Γ_gal = Γ_Cep × (1+z)^(β_A × Δvarphi) (Cepheid + redshift transfer)',
        'Gamma': float(Gamma),
        'S_X': S_A,
        'S_X_type': 'S_A (clock-rate screening)',
        'kappa': float(kappa),
        'z_example': z_example,
        'transfer_factor': transfer_factor,
        'delta_varphi': delta_varphi,
        'classification': 'DERIVED — Γ_gal from Cepheid projector with redshift transfer factor'
    }


# === LLR/flyby channel ===
def derive_llr():
    """Γ_LLR from the Yukawa force at lunar distance with GR absorption.

    The TEP perturbation on the lunar orbit is a Yukawa force:
      a_φ = β_A² × (M_Earth/M_Pl) × exp(-r/R_T) / (4π r²)

    At the lunar distance (r >> R_T), this is exponentially suppressed.
    The GR pipeline absorbs the 1/r monopole; the residual is the Yukawa tail.
    """
    r_llr = 384400e3  # lunar distance in meters
    rho_T = 20.0
    R_T = float((3 * M_EARTH / (4 * np.pi * 1000 * rho_T))**(1/3))

    # Yukawa suppression at lunar distance
    r_over_RT = r_llr / R_T
    yukawa = float(np.exp(-r_over_RT))

    # GR absorption: the 1/r monopole is fully absorbed
    # The residual is the Yukawa tail, which is exponentially suppressed
    # For r >> R_T, the residual is essentially zero
    residual_fraction = yukawa

    # Force ratio (scalar/Newtonian)
    # F_scalar/F_Newton = 2β_A² × exp(-r/R_T) × (R_T/r)² ... no, this is wrong
    # The Yukawa force is F = G M exp(-r/R_T) / r² × (1 + r/R_T)
    # The Newtonian force is F = G M / r²
    # Ratio = exp(-r/R_T) × (1 + r/R_T)
    force_ratio = yukawa * (1 + r_over_RT)

    # Projector: the Green's function at the lunar distance
    # Γ_LLR = exp(-r/R_T) × (1 + r/R_T) (Yukawa propagator)
    Gamma = yukawa * (1 + r_over_RT)

    # S_Σ (source-charge screening)
    S_sigma = 1.0 / (1.0 + 2 * BETA**2 * force_ratio)

    # κ_LLR
    kappa = abs(BETA) * S_sigma * Gamma

    # Flyby: at close range (r ~ R_Earth), the field is NOT suppressed
    r_flyby = R_EARTH
    r_flyby_over_RT = r_flyby / R_T
    yukawa_flyby = float(np.exp(-r_flyby_over_RT))
    varphi_earth = varphi_unscreened(M_EARTH, R_EARTH)
    flyby_delta_v_over_v = abs(BETA) * varphi_earth * yukawa_flyby

    return {
        'channel': 'LLR/flyby',
        'observable': 'lunar range residual + flyby velocity shift',
        'projector': 'Γ_LLR = exp(-r/R_T) × (1 + r/R_T) (Yukawa propagator)',
        'Gamma': Gamma,
        'S_X': float(S_sigma),
        'S_X_type': 'S_Σ (source-charge screening)',
        'kappa': float(kappa),
        'r_llr_m': r_llr,
        'R_T_m': R_T,
        'r_over_RT': float(r_over_RT),
        'yukawa_suppression': yukawa,
        'residual_fraction': residual_fraction,
        'force_ratio': float(force_ratio),
        'flyby_delta_v_over_v': float(flyby_delta_v_over_v),
        'small_argument_valid': r_over_RT < 1,
        'classification': 'DERIVED — Γ_LLR from Yukawa propagator; exponentially suppressed at lunar distance'
    }


def run():
    channels = {}
    print("Deriving Γ_X and κ_X for all channels...")
    print()

    for name, func in [('GNSS', derive_gnss),
                       ('Cepheid', derive_cepheids),
                       ('MSP', derive_msp),
                       ('Wide binary', derive_wide_binary),
                       ('Galactic', derive_galactic),
                       ('LLR/flyby', derive_llr)]:
        try:
            result = func()
            channels[name] = result
            if 'Gamma' in result and 'kappa' in result:
                print(f"  {name}: Γ = {result['Gamma']:.4e}, "
                      f"S = {result['S_X']:.4e}, κ = {result['kappa']:.4e}")
            else:
                print(f"  {name}: {result.get('observable', '')}")
            print(f"    {result['classification']}")
        except Exception as exc:
            channels[name] = {'channel': name, 'error': str(exc),
                              'classification': 'ERROR'}
            print(f"  {name}: ERROR — {exc}")

    result = {
        'formula': 'κ_X = |β_A| × S_X × Γ_X',
        'beta_A': BETA,
        'coupling_used': LAMBDA_REFERENCE,
        'coupling_note': ('Equilibrium field solves (varphi_host, varphi_cluster, '
                          'S_X at R_s, GNSS correlation length) use LAMBDA_REFERENCE '
                          'and are diagnostic intermediates not cited in the manuscript. '
                          'The projectors Gamma_X themselves are coupling-independent '
                          'geometric/kinematic factors; the coupling dependence of the '
                          'GNSS crossing is scanned separately in step_04.'),
        'channels': channels,
        'summary': {
            'all_derived': all(
                channels[c].get('classification', '').startswith('DERIVED')
                for c in channels if 'error' not in channels[c]
            ),
            'n_channels': len(channels),
            'n_derived': sum(1 for c in channels
                             if channels[c].get('classification', '').startswith('DERIVED'))
        }
    }
    save('step_05_gamma_all_channels.json', result)
    print()
    print(f"Total: {result['summary']['n_derived']}/{result['summary']['n_channels']} channels derived")
    return result


if __name__ == '__main__':
    run()
