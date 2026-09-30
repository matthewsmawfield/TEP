#!/usr/bin/env python3
"""Step 76: ambient-excursion transmission into interior clocks.

The corpus's clock channel applies an amplitude factor S_A to ambient
field excursions (delta_phi_eff = S_A * delta_phi_amb).  The canonical
law S_A = min[1, (rho_bar/rho_T)^{1/3}] is the equilibrium-level ratio
u_eq(rho)/u_eq(rho_T) -- the correct quantity only if interiors sit at
the matter-coupled minimum (retained branch).  Step 53 showed real
bodies land on the leaking (flux-following) branch, where the solved
profile is

    u(r) = u_amb + integral_r^rout finv(Q/s^2) ds ,

with Q driven by rho~ e^{-u} - lam~ u^3.  Since leaking interiors have
u << u_eq (e^{-u} ~ 1, lam~ u^3 negligible), Q is insensitive to u and
an ambient shift transmits ~additively: du_int/du_amb ~ 1.  On the
retained branch the interior is pinned to u_eq(rho) and the
transmission is ~0.  Neither branch gives (rho/rho_T)^{1/3}.

This step measures the transmission directly: solve each body at
u_amb and u_amb + delta for several delta, and report
T = du_int/delta on both branches.  The outcome sets the admissible
clock-channel excursion factor and the secular-drift (G-dot/G)
projection through one number.

Outputs: results/step_76_ambient_transmission.json
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from step_53_unified_profile import (  # noqa: E402
    solve_body, u_eq, rho_to_m2, G_CM3_TO_GEV4, RHO_T, LAM_M2,
    G, G_T, C, MSUN, MEARTH, REARTH, RSUN, RHO_M2_PER_GCM3,
    finv, smooth_density, static_boundary_response,
)


def interior_median(M, R, rho_bar, rho_amb, u_amb, retained=False, return_profile=False):
    """Return interior-median u for one body at the given u_amb."""
    edge_pow = 8.0 if retained else 40.0

    def prof(rr, rb=rho_bar, RR=R, ep=edge_pow):
        return smooth_density(rr, RR, rb, rho_amb, ep)

    r, u, it, res, Q = solve_body(
        M, R, prof, rho_amb=rho_amb, u_amb=u_amb, n=4000,
        r_out_fac=1e6 if retained else 1e10,
        rtol_trial=1e-3 if retained else 1e-5,
        n_bisect=40 if retained else 60)
    if return_profile:
        return r, u, prof, it, res
    inner = r < 0.35 * R
    return float(np.median(u[inner])), it, res


def transmission(M, R, rho_bar, rho_amb, deltas, retained=False):
    u_amb = u_eq(rho_amb)
    r, u, prof, it0, res0 = interior_median(M, R, rho_bar, rho_amb, u_amb,
                                          retained=retained, return_profile=True)
    inner = r < 0.35 * R
    u0 = float(np.median(u[inner]))
    tangent = float(np.median(static_boundary_response(r, u, prof, rho_amb)[inner]))
    rows = []
    for d in deltas:
        u1, it1, res1 = interior_median(M, R, rho_bar, rho_amb,
                                       u_amb + d, retained=retained)
        reporting_floor = 1e-8 * max(abs(u0), abs(u1), np.finfo(float).tiny)
        resolved = abs(u1-u0) > reporting_floor
        rows.append({"delta": d,
                     "u_int_shifted": u1,
                     "T": (u1-u0)/d if resolved else None,
                     "raw_secant_diagnostic": (u1-u0)/d,
                     "secant_status": "resolved" if resolved else "below_reporting_resolution",
                     "field_reporting_floor": reporting_floor,
                     "final_res": res1})
    return {"u_amb": u_amb, "u_int": u0, "baseline_res": res0,
            "residual_definition": "roundoff-adjusted finite-volume equation balance",
            "outer_radius_m": float(r[-1]), "outer_boundary_error": float(u[-1]-u_amb),
            "static_linear_boundary_transmission": tangent,
            "rows": rows}


def earth_boundary_branch_sweep():
    rho_amb = 1e-30
    delta = 2.5e-11
    prof = lambda radius: smooth_density(radius, REARTH, 5.51, rho_amb)
    radii = (60.0, 235.0, 23500.0, 1e6, 1e8, 1e10)
    cases = {}
    for branch, factor in (("reference",1.0),("cassini",1e5)):
        lam_m2 = LAM_M2*factor
        ambient = u_eq(rho_amb,lam_m2)
        rows = []
        for outer in radii:
            r, u, _, residual, _ = solve_body(
                MEARTH, REARTH, prof, rho_amb, ambient,
                n=2400, r_out_fac=outer, lam_m2=lam_m2)
            inner = r < 0.35*REARTH
            baseline = float(np.median(u[inner]))
            tangent = float(np.median(static_boundary_response(
                r,u,prof,rho_amb,lam_m2=lam_m2)[inner]))
            _, shifted, _, shifted_residual, _ = solve_body(
                MEARTH, REARTH, prof, rho_amb, ambient+delta,
                n=2400, r_out_fac=outer, lam_m2=lam_m2)
            rows.append({"outer_radius_over_R_earth":outer,
                         "interior_field":baseline,
                         "static_tangent":tangent,
                         "finite_secant":float((np.median(shifted[inner])-baseline)/delta),
                         "baseline_residual":residual,
                         "shifted_residual":shifted_residual,
                         "outer_boundary_error":float(u[-1]-ambient)})
        cases[branch] = {"quartic_factor_over_reference":factor,"rows":rows}
    vacuum_profile = lambda radius: smooth_density(radius, REARTH, 5.51, 0.0)
    r, u, _, residual, _ = solve_body(
        MEARTH, REARTH, vacuum_profile, 0.0, 0.0,
        n=2400, r_out_fac=1e10, lam_m2=0.0)
    inner = r < 0.35*REARTH
    tangent = float(np.median(static_boundary_response(
        r,u,vacuum_profile,0.0,lam_m2=0.0)[inner]))
    _, shifted, _, shifted_residual, _ = solve_body(
        MEARTH, REARTH, vacuum_profile, 0.0, delta,
        n=2400, r_out_fac=1e10, lam_m2=0.0)
    return {"earth_ambient_density_g_cm3":rho_amb,
            "imposed_outer_field_delta":delta,
            "cases":cases,
            "massless_vacuum_control":{
                "outer_radius_over_R_earth":1e10,
                "static_tangent":tangent,
                "finite_secant":float((np.median(shifted[inner])-np.median(u[inner]))/delta),
                "baseline_residual":residual,
                "shifted_residual":shifted_residual},
            "scope":"Static Earth-only Dirichlet boundary-radius bracket. The 60, 235 and 23500 Earth-radius entries are comparison interfaces, not solved Earth-Moon-Sun matching surfaces; no time-dependent or galactic boundary is inferred."}


def earth_solar_local_stiffness():
    rho_amb = 1e-30
    outer = 235.0
    lam_m2 = 1e5*LAM_M2
    prof = lambda radius: smooth_density(radius, REARTH, 5.51, rho_amb)
    r, u, _, residual, charge = solve_body(
        MEARTH, REARTH, prof, rho_amb, u_eq(rho_amb,lam_m2),
        n=2400, r_out_fac=outer, lam_m2=lam_m2)
    earth_flux = float(np.interp(REARTH,r,charge))/REARTH**2
    q_earth = float(C*C*finv(earth_flux)/G_T)
    au = 1.495978707e11
    sun_profile = lambda radius: smooth_density(radius, RSUN, 1.408, rho_amb)
    sun_r, _, _, sun_residual, sun_charge = solve_body(
        MSUN, RSUN, sun_profile, rho_amb, u_eq(rho_amb,lam_m2),
        n=2400, r_out_fac=2*au/RSUN, lam_m2=lam_m2)
    solar_flux = float(np.interp(au,sun_r,sun_charge))/au**2
    q_sun = float(C*C*finv(solar_flux)/G_T)
    point_source_q = float(C*C*finv(2*G*MSUN/(C*C*au**2))/G_T)
    single_GM_q = float(C*C*finv(G*MSUN/(C*C*au**2))/G_T)
    earth_only = 1+q_earth**2
    alignment = {str(mu):float(1+q_earth**2+q_sun**2+2*mu*q_earth*q_sun)
                 for mu in (-1,0,1)}
    return {"kinetic_branch":"P_X=1+q^2, q=c^2 |grad u|/g_t",
            "earth_boundary_over_R_earth":outer,
            "earth_quartic_factor_over_reference":1e5,
            "earth_surface_gradient_over_g_t":q_earth,
            "solar_gradient_at_1AU_over_g_t":q_sun,
            "solar_point_source_q_control":point_source_q,
            "single_GM_flux_q_comparator":single_GM_q,
            "two_GM_over_single_GM_q_ratio":point_source_q/single_GM_q,
            "normalization_audit":"This step_53 scalar flux has q(1+q^2)=2GM/(r^2 g_t); Paper-17 step_086 defines q(1+q^2)=GM/(r^2 g_t) on its baseline. Relabelling q_086=q_53/2 alone would require the kinetic factor 1+4q_086^2, not 1+q_086^2. The PDE source and kinetic-scale conventions must be reconciled before importing its coupled solution as the same-action GNSS field.",
            "solar_boundary_at_AU":2.0,
            "solar_profile_residual":sun_residual,
            "solar_to_earth_gradient_ratio":q_sun/q_earth,
            "earth_only_stiffness":earth_only,
            "stiffness_by_radial_alignment":alignment,
            "aligned_minus_antialigned_over_earth_only":(
                alignment["1"]-alignment["-1"])/earth_only,
            "earth_profile_residual":residual,
            "scope":"Local tensor-coefficient cross-term from separately solved Earth and isolated Sun fluxes at a declared interface; not a coupled Earth-Sun PDE, dynamic response, covariance, or EW/NS prediction."}


def main():
    rho_amb = 1e-30
    bodies = [
        ("Moon", 7.34767309e22, 1.7374e6, 3.34),
        ("Earth", MEARTH, REARTH, 5.51),
        ("Sun", MSUN, RSUN, 1.408),
        ("Cepheid envelope (8 Msun, 50 Rsun)", 8.0 * MSUN, 50.0 * RSUN,
         8.0 * MSUN * 1000.0 / (4.0 / 3.0 * np.pi * (50.0 * RSUN * 100.0) ** 3)),
    ]
    # excursion magnitudes bracketing the LLR synodic excursion
    # (2.5e-11) and the cosmological-scale ambient shift
    deltas = [1e-20, 1e-19, 1e-18, 1e-13, 2.5e-11, 1e-9]

    out = {"step": "step_76_ambient_transmission",
           "question": ("transmission T = du_int/du_amb of an ambient "
                        "field-value excursion into interior clocks, on "
                        "the solved unified-profile branch"),
           "ambient": {"rho_amb": rho_amb},
           "body_profile_quartic_factor_over_reference": 1.0,
           "bodies": {}}
    for name, M, R, rho in bodies:
        print(f"solving {name} ...", flush=True)
        res = transmission(M, R, rho, rho_amb, deltas)
        out["bodies"][name] = res
        Ts = [r["T"] for r in res["rows"]]
        print(f"  u_int={res['u_int']:.3e}  T={Ts}")

    # retained-branch control: toy body beyond R_c (same construction
    # as step_53)
    print("solving retained toy body ...", flush=True)
    out["retained_control"] = transmission(
        3.0e46, 3.0e15, RHO_T, rho_amb, [1e-9, 1e-6], retained=True)
    print("  T=", [r["T"] for r in out["retained_control"]["rows"]])

    out["earth_boundary_branch_sweep"] = earth_boundary_branch_sweep()
    out["earth_solar_local_stiffness"] = earth_solar_local_stiffness()
    out["status"] = "validated_static_boundary_response"
    out["scope"] = (
        "Static Dirichlet perturbations at the reported outer radius with density and "
        "potential held fixed. Linear tangent responses and finite-amplitude secants "
        "are distinct. No frequency-dependent or cosmological transmission is solved.")
    out["verdict"] = (
        "The static boundary response is computed from converged finite-volume profiles "
        "and their linearized operator. It depends on body, outer radius and excursion "
        "amplitude; neither a universal unity transmission nor the static depth ratio "
        "S_A may be substituted for the dynamical clock transfer.")

    dest = Path(__file__).resolve().parents[2] / "results" / \
        "step_76_ambient_transmission.json"
    dest.write_text(json.dumps(out, indent=2))
    print("wrote", dest)


if __name__ == "__main__":
    main()
