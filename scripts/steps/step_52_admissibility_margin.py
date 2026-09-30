#!/usr/bin/env python3
"""Pointwise admissibility-margin ledger for the disformal sector.

The null-cone condition selects B0 >= 0 as the unconditionally
admissible branch; B < 0 is not excluded axiomatically but carries a
pointwise admissibility condition on each realized configuration
(Paper 0, SS4).  This step evaluates that condition numerically on the
corpus's existing realized profiles, so that the branch statement is
supported by a measured margin rather than by assertion.

The three walls, all expressed in the single deformation combination

    D(r) = B(u) (du/dl R_H)^2 = B0 * shape(u) * (du/dl R_H)^2,

the same functional whose peak step_50 caps at 1e-15 via GW170817:

  1. Signature.  For spacelike grad phi the matter metric stays
     Lorentzian iff 1 + D > 0 pointwise; D -> -1 is the degeneracy.
     Wall (negative branch only):  |B0| < 1 / Dbar_max.
  2. Fluid hyperbolicity.  Perfect-fluid principal coefficients
     Z_t = 1 + D (1+w) R,  Z_s = 1 + D (1-w) R / 3,
     with R = rho / (d phi)^2 the ratio of the local matter energy
     density to the scalar gradient energy density (both GeV^4 in
     natural units).  For D < 0 the binding wall is Z_t at stiff
     w -> 1:  |D| < 1 / (2 R).
  3. Causality.  v_par = c / sqrt(1 + D): any D < 0 is superluminal
     with respect to g; the margin is quantitative only.

For B0 > 0 none of the walls binds on spacelike gradients (D > 0
automatically satisfies 1 + D > 0 and Z >= 1); the realized margin at
the GW-saturated B0 = 77.7 is reported so the unconditional branch is
quantified, not merely labelled safe.

The timelike-gradient channel (cosmological drift) is evaluated
separately: there the lapse condition A^2 N^2 > B phi_dot^2 / c^2 caps
the POSITIVE branch (the wall TEP-BBN meets at high redshift), while
the negative branch is lapse-safe but still faces the fluid wall with
rho ~ rho_crit.  The two channels together make the ledger's central
point concrete: the dangerous sign is gradient-character selective,
so no sign choice of B is universally safe.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

# --- constants (same bookkeeping as step_50 / core.constants) ---
C = 2.99792458e8                      # m/s
G = 6.67430e-11                       # m^3 kg^-1 s^-2
MSUN = 1.98847e30                     # kg
PC = 3.085677581e16                   # m
H0 = 70e3 / (1e6 * PC)                # s^-1
R_H = C / H0                          # m
G_T = C * H0 / 2.0                    # m/s^2 (beta_A^2 = 1)
M_PL_GEV = 2.435e18                   # reduced Planck mass
HBAR_C_GEV_M = 1.973269804e-16        # GeV m
HBAR_GEV_S = 6.582119569e-25          # GeV s
G_CM3_TO_GEV4 = 1000.0 * 5.60958885e26 * HBAR_C_GEV_M ** 3
# natural-units conversions used below:
#   (d phi)_nat [GeV^2] = M_PL[GeV] * (du/dl)[m^-1] * hbar c[GeV m]
#   (d_t phi)_nat [GeV^2] = M_PL[GeV] * u_dot[s^-1] * hbar[GeV s]

GW_CEILING = 1.0e-15                  # |delta c / c| path bound
RESONATOR_CEILING = 1.0e-18           # pointwise terrestrial bound
B0_GW = 77.7                          # GW bound transferred (step_50)
B0_BENCH = 1.0                        # canonical benchmark
B0_LENS = 52.0                        # LENS step_62 Refsdal requirement

ENVS = [
    # name, mass kg, radius m, interior mean rho g/cm3, exterior rho
    ("gw170817_host", 1.0e11 * MSUN, 30.0e3 * PC, 1.0e-25, 1.0e-26),
    ("solar",         1.0 * MSUN,    6.96e8,      1.408,    1.0e-12),
    ("earth",         5.972e24,      6.371e6,     5.515,    1.2e-3),
    ("neutron_star",  1.4 * MSUN,    1.2e4,       2.3e14,   1.0e-21),
    ("cluster_halo",  2.0e14 * MSUN, 1.5e6 * PC,  2.0e-28,  1.0e-30),
]


def y_of_g(g):
    g = np.asarray(g, dtype=float)
    qq = (g / G_T) ** 2
    y = np.empty_like(g)
    small = qq < 1e-8
    y[small] = 1.0 - qq[small]
    q = qq[~small]
    disc = (0.5 / q) ** 2 + (1.0 / (3.0 * q)) ** 3
    s = np.sqrt(disc)
    y[~small] = np.cbrt(0.5 / q + s) + np.cbrt(0.5 / q - s)
    return y


def g_uniform(r, mass, radius):
    r = np.asarray(r, dtype=float)
    g = np.empty_like(r)
    inn = r < radius
    g[inn] = G * mass * r[inn] / radius ** 3
    g[~inn] = G * mass / r[~inn] ** 2
    return g


def shape(u):
    u = np.asarray(u, dtype=float)
    return u ** 2 / (1.0 + u ** 2) * np.exp(-0.5 * u ** 4)


def profile(radius_m, mass_kg):
    """u(r), g(r), du/dl(r) on a geomspace grid; same solve as step_50."""
    r = np.geomspace(1e-3 * radius_m, 50.0 * radius_m, 6000)
    g = g_uniform(r, mass_kg, radius_m)
    y = y_of_g(g)
    dudl = y * g / C ** 2
    # u(r) = int_r^inf du/dl dl, integrated inward on the grid
    dr = np.diff(r)
    du = 0.5 * (dudl[1:] + dudl[:-1]) * dr
    u = np.zeros_like(r)
    u[:-1] = np.cumsum(du[::-1])[::-1]
    return r, u, g, dudl


def env_ledger(name, mass, radius, rho_in, rho_out):
    r, u, g, dudl = profile(radius, mass)
    dbar = shape(u) * (dudl * R_H) ** 2          # |D| per unit B0
    i_pk = int(np.argmax(dbar))
    dbar_pk = float(dbar[i_pk])

    rho = np.where(r < radius, rho_in, rho_out) * G_CM3_TO_GEV4

    # --- walls as |B0| limits ---
    b0_sig = 1.0 / dbar_pk
    # fluid wall: |D| < 1/((1+w) R) pointwise; binding at w=1, with
    # R = rho/(d phi)^2 in natural units.  Algebraically
    # R*Dbar = rho * shape(u) * R_H^2 / (M_Pl hbar c)^2:
    # the gradient cancels, so the fluid wall is a pure coupling-density
    # condition |B(u)| rho <~ 1 in natural units.
    with np.errstate(divide="ignore", invalid="ignore"):
        wall_fluid = (M_PL_GEV * HBAR_C_GEV_M) ** 2 / (
            2.0 * rho * shape(u) * R_H ** 2)
    b0_fluid = float(np.min(wall_fluid[wall_fluid > 0]))
    i_f = int(np.argmin(wall_fluid))

    # same-path observational cap if a GW-type bound |D|<1e-15 were
    # evaluated on this environment's own profile; the terrestrial
    # resonator bound |D|<1e-18 applies only at Earth
    b0_local_gw = GW_CEILING / dbar_pk
    dbar_surf = float(np.interp(radius, r, dbar))
    b0_resonator = (RESONATOR_CEILING / dbar_surf
                    if (name == "earth" and dbar_surf > 0) else None)

    b0_neg = min(b0_sig, b0_fluid)
    which = "fluid" if b0_fluid < b0_sig else "signature"

    # margins of the phenomenological positive branch
    phenom = {}
    for tag, b0 in (("benchmark_1", B0_BENCH),
                    ("gw_saturated_77.7", B0_GW),
                    ("lens_required_52", B0_LENS)):
        d_pk = b0 * dbar_pk
        phenom[tag] = {
            "D_peak": float(d_pk),
            "signature_margin_orders":
                float(np.log10(1.0 / d_pk)) if d_pk > 0 else None,
            "v_par_over_c": float(1.0 / np.sqrt(1.0 + d_pk)),
        }

    # negative branch at the same magnitudes: superluminality incurred
    neg = {}
    for tag, b0 in (("neg_1", -B0_BENCH), ("neg_gw", -B0_GW)):
        d_pk = b0 * dbar_pk          # negative
        neg[tag] = {
            "D_peak": float(d_pk),
            "v_par_over_c":
                float(1.0 / np.sqrt(1.0 + d_pk)) if d_pk > -1 else None,
            "signature_ok": bool(d_pk > -1.0),
            "fluid_ok": bool(abs(b0) < b0_fluid),
        }

    return {
        "mass_msun": mass / MSUN, "radius_m": radius,
        "u_surface": float(np.interp(radius, r, u)),
        "u_max": float(u.max()),
        "Dbar_peak": dbar_pk,
        "peak_r_over_R": float(r[i_pk] / radius),
        "walls_negative_B0": {
            "signature": float(b0_sig),
            "fluid_stiff": b0_fluid,
            "fluid_binding_r_over_R": float(r[i_f] / radius),
            "rho_g_cm3_at_binding":
                float(rho[i_f] / G_CM3_TO_GEV4),
            "combined": float(b0_neg),
            "binding_wall": which,
            "note": ("fluid wall is gradient-independent: "
                     "|B0| < (M_Pl hbar c)^2 / (2 rho shape R_H^2)"),
        },
        "observational_caps_B0": {
            "same_path_1e-15": float(b0_local_gw),
            "terrestrial_resonator_1e-18": b0_resonator,
        },
        "walls_over_gw_bound_orders": {
            "signature": float(np.log10(b0_sig / B0_GW)),
            "fluid_stiff": float(np.log10(b0_fluid / B0_GW))
                if np.isfinite(b0_fluid) else None,
        },
        "positive_branch_at_phenom_B0": phenom,
        "negative_branch_spot_checks": neg,
    }


def drift_channel():
    """Timelike-gradient channel: lapse cap on B>0, fluid wall on B<0.

    u drifts at the Hubble rate (u_dot = H0 fiducial); rho = rho_crit.
    """
    u_dot = H0
    gradphi2_t = (M_PL_GEV * u_dot * HBAR_GEV_S) ** 2
    rho_crit_kg_m3 = 3.0 * H0 ** 2 / (8.0 * np.pi * G)
    rho_crit = (rho_crit_kg_m3 / 1000.0) * G_CM3_TO_GEV4  # kg/m3 -> g/cm3
    R_t = rho_crit / gradphi2_t
    out = {"u_dot_over_H0": 1.0,
           "gradphi2_t_GeV4": float(gradphi2_t),
           "rho_crit_GeV4": float(rho_crit),
           "R_t": float(R_t),
           "per_field_value": {}}
    for u_val in (0.3, 0.5, 1.0, 2.0, 5.0, 10.0):
        sh = float(shape(u_val))
        if sh <= 0:
            continue
        out["per_field_value"][f"u_{u_val}"] = {
            "shape": sh,
            "lapse_wall_B0_positive": 1.0 / sh,
            "fluid_wall_B0_magnitude": 1.0 / (2.0 * R_t * sh),
        }

    # --- realized absorber-epoch ambients (Paper 29 gate10 backgrounds):
    # phi_bar/M_Pl = ln(1+z), A = (1+z)^-1, du/dt = H_T(z) ~ 4-5 H0.
    # The covariant lapse ratio is u_bound = B(phi_bar) phi_dot^2/A^2
    # = B0 * shape(u) * (H_T/H0)^2 / A^2; the positive-branch cap is
    # therefore A^2/(shape*(H_T/H0)^2), and the negative branch is
    # bounded by the ambient fluid wall |B|(1+w)rho_m(z) < A^2, with
    # rho_m(z) = Omega_m (1+z)^3 rho_crit. B0 in the corpus
    # dimensionless normalization (B_phys = B0 shape R_H^2/(M_Pl hbar
    # c)^2); identical numbers are computed independently in Paper 29
    # step_10c (results/gate10c_ambient_lapse_margin.json).
    OMEGA_M = 0.3
    realized = []
    for sid, z in (("Q1009+2956", 2.5042), ("PKS1937-101", 3.572),
                   ("J1332+0052", 3.421)):
        u_amb = float(np.log(1.0 + z))
        A2 = float(np.exp(-2.0 * u_amb))
        E2 = OMEGA_M * (1.0 + z) ** 3 + (1.0 - OMEGA_M)
        sh = float(shape(u_amb))
        rho_m = OMEGA_M * (1.0 + z) ** 3 * rho_crit
        realized.append({
            "id": sid, "z_abs": z, "u_ambient": u_amb, "A2": A2,
            "H_T_over_H0": float(np.sqrt(E2)), "shape": sh,
            "lapse_wall_B0_positive": float(A2 / (sh * E2)),
            "fluid_wall_B0_magnitude_w1": float(
                A2 * (M_PL_GEV * HBAR_C_GEV_M) ** 2 / (
                    2.0 * rho_m * sh * R_H ** 2)),
            "rho_m_GeV4": float(rho_m),
        })
    out["realized_absorber_epochs"] = {
        "source": "gate10 backgrounds (phi_bar=ln(1+z), A=e^-phi_bar, "
                  "d0 phi_bar = M_Pl H_T); cross-checked by Paper 29 "
                  "step_10c_ambient_lapse_margin",
        "sightlines": realized,
        "tightest_lapse_cap": float(
            min(r["lapse_wall_B0_positive"] for r in realized)),
        "tightest_fluid_cap": float(
            min(r["fluid_wall_B0_magnitude_w1"] for r in realized)),
        "verdict": (
            "on the realized absorber-epoch ambients the ungated "
            "canonical envelope is ambient-lapse-supercritical for "
            "phenomenological B0 ~ 52-78 by ~10^3 (cap ~0.03), and the "
            "negative branch is capped at |B0| ~ 10^-2 by the ambient "
            "fluid wall: neither sign can carry the phenomenological "
            "amplitude ungated. Ambient closure of the environmental "
            "gate is therefore lapse-mandated — independently of, and "
            "consistent with, the absorber-amplitude requirement of "
            "Paper 29 gate10b (required ambient suppression G_ambient "
            "<~ 4e-4 at B0 = 78 vs the in-well saturation "
            "G ~ 1e-3-3e-3)."),
    }
    return out


def main():
    envs = {name: env_ledger(name, m, rad, ri, ro)
            for name, m, rad, ri, ro in ENVS}

    # cross-reference: step_13's computed negative-branch failure
    s13_path = Path(__file__).resolve().parents[2] / "results" / \
        "step_13_disformal_closure.json"
    s13 = None
    if s13_path.exists():
        s13 = json.loads(s13_path.read_text())
        s13 = {"Q_min_universal_scan":
               s13["scan"]["universal,delta=1e-6"]["Q_min"],
               "physical_status": s13["verdict"]["physical_status"]}

    host = envs["gw170817_host"]
    earth = envs["earth"]
    drift = drift_channel()

    # --- drift-transit analysis: the cosmological field amplitude is
    # not a free choice. Under the adopted clock map
    # A_clock(z) = (1+z)^-1 with beta_A = -1, phi_bar/M_Pl = ln(1+z):
    # the drift transits every amplitude u in [0, ln(1+z_max)] over
    # cosmic history. The REALIZED lapse cap follows the full ratio
    # B0 < A^2(u) / (shape(u) * (H_T/H0)^2) with A^2 = e^{-2u} and
    # (H_T/H0)^2 = E^2(z): the A^2 suppression and the growing drift
    # rate pull the cap minimum to z ~ 3 (u ~ 1.4, cap ~ 0.03), well
    # below both the shape-peak fiducial (u* ~ 0.87, cap ~ 3.1 at
    # u_dot = H0, A = 1) and the phenomenological range. A uniform
    # positive B0 therefore cannot be parked below the wall: the
    # realized configuration is forced through the absorber-epoch
    # minimum. The transit wall is evaded by the absorber-sector
    # structure itself: the B_eff = B(phi)*G(X_local) gate opens at
    # high X inside nonlinear wells and closes on the low-X
    # cosmological ambient — precisely the suppression the drift
    # sector needs (the ambient fluid wall binds the negative branch
    # at the same epochs, |B0| <~ 10^-2, so neither sign escapes).
    uu = np.linspace(0.05, 5.0, 4000)
    ss = shape(uu)
    iu = int(np.argmax(ss))
    u_star = float(uu[iu])
    cap_star = float(1.0 / ss[iu])          # fiducial lapse wall at peak
    # realized cap along the transit: min_z A^2/(shape*E^2)
    zz = np.linspace(0.05, 12.0, 4000)
    u_z = np.log(1.0 + zz)
    E2_z = 0.3 * (1.0 + zz) ** 3 + 0.7
    cap_z = np.exp(-2.0 * u_z) / (shape(u_z) * E2_z)
    iz = int(np.argmin(cap_z))
    z_cap_min = float(zz[iz])
    cap_min = float(cap_z[iz])
    drift["transit"] = {
        "clock_map": "phi_bar/M_Pl = ln(1+z) (A_clock=(1+z)^-1, "
                     "beta_A=-1, ambient phi=0 today)",
        "u_star_shape_peak": u_star,
        "lapse_cap_at_peak_fiducial": cap_star,
        "realized_lapse_cap_min": cap_min,
        "realized_cap_min_redshift": z_cap_min,
        "binding_wall_positive_branch": "lapse (B > 0; the fluid wall "
                        "Z_t = 1 + B(1+w)rho/A^2 binds only B < 0)",
        "transit_unavoidable": True,
        "verdict": (
            "a uniform positive B0 >~ %.3f violates the ambient lapse "
            "during the realized transit (cap minimum at z ~ %.1f, "
            "u ~ %.2f, where A^2 suppression and the H_T ~ 4-5 H0 drift "
            "rate coincide with a still-appreciable envelope shape) — "
            "forced by the clock map, not a choice. The GW-saturated "
            "B0 = 77.7 exceeds it ~%dx and the LENS requirement ~420 "
            "~%dx. The wall is evaded exactly by the absorber-sector "
            "structure: the cosmological ambient sits at low X, where "
            "the B_eff = B(phi)*G(X_local) gate is closed — so the "
            "modulation Paper 29 requires on independent grounds is "
            "also the one this wall selects. Uniform-B0 phenomenology "
            "(LENS, GW transfers) remains excluded at these epochs" %
            (cap_min, z_cap_min, np.log(1.0 + z_cap_min),
             round(B0_GW / cap_min), round(422.0 / cap_min)))}
    out = {
        "step": "step_52_admissibility_margin",
        "deformation": "D = B0 * u^2/(1+u^2) * exp(-u^4/2) * (du/dl R_H)^2",
        "convention_note": ("corpus convention bounds |D| < 1e-15 "
                            "directly; the strict photon mapping is "
                            "|dc/c| = |(1+D)^(-1/2) - 1| ~ D/2, so the "
                            "transferred B0 bound is conservative by a "
                            "factor ~2 (v_par = c/sqrt(1+D), null "
                            "condition on the inverse metric)"),
        "total_field_convention_note": (
            "the environment ledgers evaluate B(u) on the local "
            "halo/profile field with u_amb = 0 (Rule 4 convention, "
            "consistent with step_50's galactic-profile bound "
            "B0 <~ 78). Under the total-field reading "
            "u = u_bar(z) + u_halo the multimessenger path picks up "
            "the cosmic-drift term (du_bar/dx)*R_H = 1, which "
            "tightens the uniform-B0 bound to ~3e-11 (Paper 19 "
            "step_063 ambient-consistent variant); this reinforces "
            "the same conclusion as the realized lapse ledger — a "
            "uniform positive B0 has no phenomenological window, and "
            "the ambient-closed gate is mandatory."),
        "walls": {
            "signature": "1 + D > 0 pointwise (spacelike grad phi)",
            "fluid": "Z_t = 1 + D(1+w)R > 0, Z_s > 0; R = rho/(dphi)^2",
            "causality": "v_par = c/sqrt(1+D) <= c requires D >= 0",
        },
        "environments": envs,
        "drift_channel": drift,
        "step_13_cross_reference": s13,
        "headline": {
            "signature_wall_over_GW_bound_galactic":
                host["walls_over_gw_bound_orders"]["signature"],
            "fluid_wall_over_GW_bound_galactic":
                host["walls_over_gw_bound_orders"]["fluid_stiff"],
            "negative_B0_max_galactic":
                host["walls_negative_B0"]["combined"],
            "negative_B0_max_terrestrial":
                earth["walls_negative_B0"]["combined"],
            "negative_B0_max_neutron_star":
                envs["neutron_star"]["walls_negative_B0"]["combined"],
            "negative_B0_binding_wall_terrestrial":
                earth["walls_negative_B0"]["binding_wall"],
            "lapse_wall_cosmological_u1":
                drift["per_field_value"]["u_1.0"]["lapse_wall_B0_positive"],
            "fluid_wall_cosmological_u1":
                drift["per_field_value"]["u_1.0"]["fluid_wall_B0_magnitude"],
            "lapse_wall_realized_absorber_epochs":
                drift["realized_absorber_epochs"]["tightest_lapse_cap"],
            "fluid_wall_realized_absorber_epochs":
                drift["realized_absorber_epochs"]["tightest_fluid_cap"],
            "note": ("on the realized absorber-epoch ambients "
                     "(phi_bar=ln(1+z)~1.25-1.5, du/dt=H_T~4-5 H0, "
                     "A~0.22-0.29) the positive branch's lapse cap is "
                     "B0<~0.03 and the negative branch's ambient fluid "
                     "cap |B0|<~10^-2: both sit ~10^3 below the "
                     "phenomenological B0~52-78, so ambient closure of "
                     "the environment gate G(X_local) is lapse-mandated "
                     "independently of the absorber-amplitude argument "
                     "(Paper 29 gates 10b/10c). In dense bodies the "
                     "|B|rho fluid wall binds the negative branch "
                     "(|B0|<~0.7 in Earth, <~2e-16 in a neutron star); "
                     "neither wall is a theorem, both are computed "
                     "margins on realized profiles"),
        },
    }
    dest = Path(__file__).resolve().parents[2] / "results" / \
        "step_52_admissibility_margin.json"
    dest.write_text(json.dumps(out, indent=2))
    print(json.dumps(out["headline"], indent=2))


if __name__ == "__main__":
    main()
