#!/usr/bin/env python3
"""Unified-profile derivation: S_Sigma and S_A read off one solved field.

Gate-A item (a1). The corpus handles the two screening channels as
separate operators:

    shear branch:    S_Sigma(g) = P_,X^-1 at the local gradient
                     (two-branch sector P_,X = k sqrt(xi) + 2 xi,
                      k = 4 sqrt(2) g_t/a_0 ~ 16.0 — promoted per
                      MASTER_SECTOR_LEDGER step_54 / WB step_020 /
                      growth step_62; former unit-term baseline
                      P_,X = 1 + 2 xi retained via --operator baseline)
    amplitude branch: S_A(rho) = min[1, (rho_bar/rho_T)^{1/3}]
                     (quartic matter-coupled equilibrium u_min ∝ rho^{1/3},
                      step_27 sec.6 / equilibrium_u)

This step solves the ONE static spherical equation that carries both,

    (1/r^2) d/dr [ r^2 P_,X(u') u' ] = d V_eff / du ,
    V_eff(u; rho) = (lambda/4) M_Pl^4 u^4 + rho e^{-u},
    P_,X(u') = (k/sqrt2)(c|u'|/H0) + (c u'/H0)^2   [two-branch]
             or 1 + (c u'/H0)^2                     [--operator baseline],
    u = phi/M_Pl,

for a stratified body (density profile rho(r)) embedded in an ambient
medium at u_amb = u_eq(rho_amb), and reads both channels off the same
solution u(r):

    density branch:  u(r) tracks u_eq(rho(r)) inside the body
                     => interior amplitude response ∝ rho^{1/3}
                     => S_A(rho_bar) = u_c/u_eq(rho_T) = (rho/rho_T)^{1/3}
    shear branch:    exterior flux response y(r) = a/a_lin, where
                     a = c^2|u'| and a_lin = Q/r^2 is the P_,X=1 gradient,
                     giving y = P_,X(u')^{-1}. The unit-term baseline has
                     y = [1+(a/(c H0))^2]^{-1}; the selected two-branch
                     expression includes the linear-in-|gradient| term.

Solver: integrating once, the equation is algebraic in u':

    r^2 f(u') = -Q(r),   Q(r) = int_0^r [rho~(r') exp(-u) - lam~ u^3] r'^2 dr',
    f(v) = v P_X(v)  [two-branch: (H0/c) a^2 (k/sqrt2 + a),
                      a = c |v| / H0 ; baseline: v(1 + (c v/H0)^2)],

so u'(r) = -finv(Q/r^2) and the profile follows by quadrature,

    u(r) = u_amb + int_r^r_out finv(Q(s)/s^2) ds.

Both the matter coupling and potential contribute to Q; it is not a
conserved bare source charge. The production solver uses conservative
radial finite volumes and a damped Newton solve, with the full source
derivative rho~ exp(-u) + 3 lam~ u^2 in the positive static Hessian.
The shooting implementation remains an optional cross-check. Two regimes:

    leaking (flux-following): interior u << u_eq, Q(r) ~ rho~r^3/3,
        interior excess Delta u ~ sqrt(Q_bare g_t)/c set by matching
        to the exterior tail;
    retained (saturated): interior reaches u_eq, the source dies and
        Q stops growing; possible only for R >~ R_c with
        R_c^3 ~ c^2 u_eq(rho)^2/(rho~ g_t).

Output: results/step_53_unified_profile_{operator}.json
"""
import json
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import brentq

# ---------------------------------------------------------------------------
# constants (SI + natural-unit conversion)
# ---------------------------------------------------------------------------
C = 2.99792458e8
G = 6.67430e-11
MSUN = 1.98847e30
MEARTH = 5.972e24
REARTH = 6.371e6
RSUN = 6.957e8
H0 = 70.0e3 / 3.085677581e22           # s^-1 (corpus convention)
BETA_A = -1.0
G_T = C * H0 / (2.0 * BETA_A ** 2)      # 3.40e-10 m/s^2
M_PL_GEV = 2.435e18
GEV2_TO_M2 = (5.0677307e15) ** 2        # 1 GeV^2 in m^-2
G_CM3_TO_GEV4 = 4.3070e-18              # 1 g/cm^3 in GeV^4
LAMBDA = 7.526e-71                      # quartic coupling (corpus ref)
RHO_T = 20.0                            # g/cm^3, corpus calibration

# Kinetic-sector operator. The promoted corpus realization is the
# two-branch sector P_,X = k sqrt(xi) + 2 xi with xi = |X|/Lam^4
# (MASTER_SECTOR_LEDGER step_54, WB step_020, growth step_62/74); the
# former unit-term baseline 1 + 2 xi is retained for comparison.
# In gradient variable a = c |u'| / H0 the baseline maps xi = a^2/2,
# so the two branches are P_,X = (k/sqrt2) a + a^2.
A0 = 1.2e-10                            # m/s^2, corpus acceleration floor
K2 = 4.0 * np.sqrt(2.0) * G_T / A0      # k ~= 16.0 (zero free parameters)
KP = K2 / np.sqrt(2.0)                  # coefficient of the sqrt branch
OPERATOR = "baseline"                   # set by --operator in main()

LAM_M2 = LAMBDA * M_PL_GEV ** 2 * GEV2_TO_M2   # lambda M_Pl^2 in m^-2


RHO_M2_PER_GCM3 = G_CM3_TO_GEV4 / M_PL_GEV ** 2 * GEV2_TO_M2


def rho_to_m2(rho_g_cm3):
    """rho/(M_Pl^2) converted to m^-2 source units."""
    return rho_g_cm3 * RHO_M2_PER_GCM3


def u_eq(rho_g_cm3, lam_m2=LAM_M2):
    """Local matter-coupled equilibrium: lambda M^4 u^3 = rho e^-u."""
    b = rho_to_m2(rho_g_cm3)
    # solve in log-u to span the huge dynamic range
    g = lambda lu: np.log(lam_m2) + 3.0 * lu - np.log(b) + np.exp(lu)
    lo, hi = np.log(1e-40), np.log(1e3)
    while g(hi) < 0:
        hi += np.log(10.0)
    return float(np.exp(brentq(g, lo, hi)))


def P_X(v):
    """Kinetic stiffness evaluated at gradient v = du/dr [m^-1]."""
    a = C * np.abs(v) / H0
    if OPERATOR == "two_branch":
        return KP * a + a ** 2
    return 1.0 + a ** 2


def J_flux(v):
    """d[v P_X]/dv (radial flux-map Jacobian)."""
    a = C * np.abs(v) / H0
    if OPERATOR == "two_branch":
        return 2.0 * KP * a + 3.0 * a ** 2
    return 1.0 + 3.0 * a ** 2


def finv(q):
    """Root of v P_X(v) = q, vectorised, for the selected operator.

    baseline:    v (1 + (c v/H0)^2) = q
    two-branch:  in a = c|v|/H0, a^2 (KP + a) = |q| c/H0
    Both maps are convex and increasing on the positive half-line, so
    Newton from an upper bound descends monotonically — no divergence
    possible.
    """
    q = np.asarray(q, dtype=float)
    aq = np.abs(q)
    if OPERATOR == "two_branch":
        qq = aq * C / H0
        a = np.maximum(np.sqrt(qq / KP), 1e-300) \
            + np.maximum(qq, 1e-300) ** (1.0 / 3.0)
        for _ in range(80):
            with np.errstate(over="ignore", invalid="ignore"):
                f = a * a * (a + KP) - qq
                fp = 3.0 * a * a + 2.0 * KP * a
                step = np.where(np.isfinite(f / fp), f / fp, 0.0)
            a = np.maximum(a - step, 1e-300)
            if np.max(np.abs(step / a)) < 1e-14:
                break
        return np.sign(q) * a * H0 / C
    alpha = (C / H0) ** 2
    v = np.minimum(np.maximum(aq, 1e-300),
                   np.maximum((aq / alpha) ** (1.0 / 3.0), 1e-300))
    for _ in range(60):
        with np.errstate(over="ignore", invalid="ignore"):
            f = v * (1.0 + alpha * v ** 2) - aq
            fp = 1.0 + 3.0 * alpha * v ** 2
            step = np.where(np.isfinite(f / fp), f / fp, 0.0)
        v = np.maximum(v - step, 1e-300)
        if np.max(np.abs(step / v)) < 1e-14:
            break
    return np.sign(q) * v


def smooth_density(r, radius, body_density, ambient_density, edge_power=40.0):
    x = np.maximum(np.asarray(r, dtype=float) / radius, np.finfo(float).tiny)
    weight = np.exp(-np.logaddexp(0.0, edge_power * np.log(x)))
    return ambient_density + (body_density - ambient_density) * weight


def _radial_balance(r, u, rho, rho_amb, lam_m2):
    faces = (r[1:] + r[:-1]) / 2
    dr = np.diff(r)
    volumes = np.diff(np.r_[0.0, faces**3]) / 3
    slope = np.diff(u) / dr
    flux = -faces**2 * slope * P_X(slope)
    equilibrium = u_eq(rho_amb / RHO_M2_PER_GCM3, lam_m2) if rho_amb > 0 and lam_m2 > 0 else 0.0
    offset = u - equilibrium
    source = ((rho - rho_amb) * np.exp(-u)
              + rho_amb * np.exp(-equilibrium) * np.expm1(-offset)
              - lam_m2 * offset * (u*u + u*equilibrium + equilibrium**2))
    if lam_m2 == 0:
        source = rho * np.exp(-u)
    residual = np.diff(np.r_[0.0, flux]) - source[:-1] * volumes
    conductance = faces**2 * J_flux(slope) / dr
    mass = rho * np.exp(-u) + 3 * lam_m2 * u*u
    diagonal = conductance + np.r_[0.0, conductance[:-1]] + mass[:-1] * volumes
    scale = (np.abs(flux) + np.r_[0.0, np.abs(flux[:-1])]
             + (rho[:-1]*np.exp(-u[:-1]) + lam_m2*np.abs(u[:-1])**3)*volumes)
    roundoff = (8*np.finfo(float).eps
                * (conductance*(np.abs(u[:-1]) + np.abs(u[1:]))
                   + np.r_[0.0, conductance[:-1]*(np.abs(u[:-2]) + np.abs(u[1:-1]))]))
    residual = np.where(np.abs(residual) <= roundoff, 0.0, residual)
    norm = float(np.max(np.abs(residual) / np.maximum(scale, np.finfo(float).tiny)))
    return residual, diagonal, conductance, flux, norm


def _kinetic_density(v):
    """G(v) = int_0^v w P_X(w) dw for the selected operator.

    baseline:   v^2/2 + (c/H0)^2 v^4/4
    two-branch: KP (c/H0) |v|^3/3 + (c/H0)^2 v^4/4
    """
    beta = C / H0
    av = np.abs(v)
    if OPERATOR == "two_branch":
        return KP * beta * av ** 3 / 3.0 + beta ** 2 * v ** 4 / 4.0
    return v ** 2 / 2.0 + beta ** 2 * v ** 4 / 4.0


def _radial_energy_change(r, u, change, rho, lam_m2):
    dr = np.diff(r)
    faces = (r[1:] + r[:-1]) / 2
    volumes = np.diff(np.r_[0.0, faces**3]) / 3
    slope = np.diff(u) / dr
    dslope = np.diff(change) / dr
    kinetic = _kinetic_density(slope + dslope) - _kinetic_density(slope)
    du, field = change[:-1], u[:-1]
    potential = (lam_m2*du*(field**3 + 1.5*field*field*du + field*du*du + du**3/4)
                 + rho[:-1]*np.exp(-field)*np.expm1(-du))
    return float(np.sum(kinetic*faces**2*dr, dtype=np.longdouble)
                 + np.sum(potential*volumes, dtype=np.longdouble))


def solve_body(M_kg, R_m, rho_prof, rho_amb, u_amb, n=20000,
               r_out_fac=1e10, tol=1e-10, rtol_trial=1e-5,
               n_bisect=60, lam_m2=LAM_M2, method="finite_volume"):
    if method == "shooting":
        return _solve_body_shooting(M_kg, R_m, rho_prof, rho_amb, u_amb,
                                   n, r_out_fac, tol, rtol_trial, n_bisect, lam_m2)
    if method != "finite_volume":
        raise ValueError(f"Unknown radial solver: {method}")
    if R_m <= 0 or r_out_fac <= 1 or n < 16 or tol <= 0 or lam_m2 < 0:
        raise ValueError("Invalid radial domain, resolution, tolerance or potential")
    r = np.geomspace(1e-6 * R_m, r_out_fac * R_m, n)
    rho = rho_to_m2(np.array([rho_prof(float(radius)) for radius in r]))
    ambient = rho_to_m2(rho_amb)
    if not np.all(np.isfinite(rho)) or np.any(rho < 0) or not np.isfinite(u_amb):
        raise ValueError("Density and boundary values must be finite; density must be nonnegative")
    faces = (r[1:] + r[:-1]) / 2
    volumes = np.diff(np.r_[0.0, faces**3]) / 3
    bare_charge = np.cumsum((rho[:-1] - ambient) * volumes)
    increments = finv(bare_charge / faces**2) * np.diff(r)
    seed_boundary = u_eq(rho_amb, lam_m2) if rho_amb > 0 and lam_m2 > 0 else 0.0
    if lam_m2 == 0:
        seed_boundary = u_amb
    u = seed_boundary + np.r_[np.cumsum(increments[::-1])[::-1], 0.0]
    if lam_m2 > 0 and np.max(rho) > 0:
        cap = max(u_amb, u_eq(float(np.max(rho) / RHO_M2_PER_GCM3), lam_m2))
        u = np.minimum(u, cap)
    u[-1] = u_amb
    if n > 1600:
        coarse_r, coarse_u, _, _, _ = solve_body(
            M_kg, R_m, rho_prof, rho_amb, u_amb, n=1200,
            r_out_fac=r_out_fac, tol=tol, lam_m2=lam_m2)
        u = np.interp(np.log(r), np.log(coarse_r), coarse_u)
        u[-1] = u_amb
    for iteration in range(1, 501):
        residual, diagonal, conductance, flux, norm = _radial_balance(r, u, rho, ambient, lam_m2)
        if np.isfinite(norm) and norm <= tol:
            charge = np.interp(r, faces, flux, left=0.0, right=flux[-1])
            return r, u, iteration, norm, charge
        if not np.isfinite(norm) or np.any(diagonal <= 0):
            raise RuntimeError("Non-finite or non-elliptic radial boundary-value system")
        jacobian = sp.diags((-conductance[:-1] / diagonal[1:],
                             np.ones(n-1), -conductance[:-1] / diagonal[:-1]),
                            (-1, 0, 1), format="csc")
        correction = spla.spsolve(jacobian, -residual / diagonal)
        directional = float(np.dot(residual, correction))
        merit = float(np.max(np.abs(residual / diagonal)))
        for power in range(30):
            change = np.r_[(0.5**power)*correction, 0.0]
            trial = u + change
            with np.errstate(over="ignore", invalid="ignore"):
                trial_residual, _, _, _, trial_norm = _radial_balance(r, trial, rho, ambient, lam_m2)
                trial_merit = float(np.max(np.abs(trial_residual / diagonal)))
                energy_change = _radial_energy_change(r, u, change, rho, lam_m2)
            if np.isfinite(trial_norm) and (energy_change <= 1e-4*(0.5**power)*directional
                                            or trial_merit < merit or trial_norm <= tol):
                u = trial
                break
        else:
            raise RuntimeError(f"Radial Newton line search stalled: balance residual {norm:.3e}")
    raise RuntimeError(f"Radial boundary-value solve did not converge: residual {norm:.3e}")


def static_boundary_response(r, u, rho_prof, rho_amb, lam_m2=LAM_M2):
    rho = rho_to_m2(np.array([rho_prof(float(radius)) for radius in r]))
    _, diagonal, conductance, _, residual = _radial_balance(
        r, u, rho, rho_to_m2(rho_amb), lam_m2)
    if not np.isfinite(residual) or residual > 1e-8:
        raise RuntimeError("Static response requires a converged background profile")
    jacobian = sp.diags((-conductance[:-1] / diagonal[1:],
                         np.ones(len(r)-1), -conductance[:-1] / diagonal[:-1]),
                        (-1, 0, 1), format="csc")
    rhs = np.zeros(len(r)-1)
    rhs[-1] = conductance[-1] / diagonal[-1]
    response = np.r_[spla.spsolve(jacobian, rhs), 1.0]
    if not np.all(np.isfinite(response)) or np.min(response) < -1e-8 or np.max(response) > 1+1e-8:
        raise RuntimeError("Static boundary response violates the elliptic maximum principle")
    return response


def _solve_body_shooting(M_kg, R_m, rho_prof, rho_amb, u_amb, n=20000,
                         r_out_fac=1e10, tol=1e-10, rtol_trial=1e-5,
                         n_bisect=60, lam_m2=LAM_M2):
    """Shooting solution of the once-integrated algebraic-flux system.

    The integrated equation is a first-order system marched outward
    in t = ln r:

        du/dt  = -r finv(Q/r^2),
        dQ/dt  = r^3 [rho~(r) e^{-u} - lam~ u^3],
        Q(0) = 0,   u(0) = u_c  (unknown),

    with the outer condition u(r_out) = u_amb.  u(r_out) is monotone
    in u_c, so bisection is robust; flux conservation is exact by
    construction since Q is the literal enclosed charge.
    """
    from scipy.integrate import solve_ivp
    r0 = 1e-6 * R_m
    r_out_m = r_out_fac * R_m
    t0, t1 = np.log(r0), np.log(r_out_m)
    t_eval = np.linspace(t0, t1, n)

    rho_hi = float(rho_prof(0.5 * R_m))
    u_cap = 5.0 * u_eq(rho_hi, lam_m2)

    def rhs(t, y):
        r = np.exp(t)
        u, Q = y
        # u_amb is a stable node approached from above; trial
        # trajectories above the separatrix run away in u (lam u^3 eats
        # the charge, then Q<0 lifts u), below it they fall to ~0.
        # Clamping keeps trials finite; only the converged
        # near-separatrix profile is kept.
        ue = min(max(u, 0.0), u_cap)
        v = float(finv(Q / r ** 2))
        dudt = -v * r
        rho_s = rho_to_m2(float(rho_prof(r)))
        dQdt = (rho_s * np.exp(-ue) - lam_m2 * ue ** 3) * r ** 3
        return [dudt, dQdt]

    def hit_cap(t, y):
        return y[0] - u_cap
    hit_cap.terminal = True
    hit_cap.direction = 1

    def fell_through(t, y):
        return y[0] + 1e-3 * u_cap
    fell_through.terminal = True
    fell_through.direction = -1

    def end_u(uc):
        # loose tolerance suffices for the sign decision; terminal
        # events abort runaway branches without marching the full box
        sol = solve_ivp(rhs, (t0, t1), [uc, 0.0], method="RK45",
                        rtol=rtol_trial, atol=1e-40,
                        events=(hit_cap, fell_through))
        return sol.y[0, -1]

    # bisection on sign(u(r_out) - u_amb): below the separatrix the
    # trajectory falls through it, above it stays >= u_amb
    lo, hi = u_amb, u_cap
    it = 0
    for _ in range(n_bisect):
        mid = 0.5 * (lo + hi)
        fm = end_u(mid) - u_amb
        it += 1
        if fm > 0:
            hi = mid
        else:
            lo = mid
    uc = 0.5 * (lo + hi)
    sol = solve_ivp(rhs, (t0, t1), [uc, 0.0], method="RK45",
                    t_eval=t_eval, rtol=1e-9, atol=1e-40,
                    events=(hit_cap, fell_through))
    if not sol.success or sol.status != 0 or not np.isclose(sol.t[-1], t1, rtol=0, atol=1e-12):
        raise RuntimeError("Shooting integration did not reach the outer boundary")
    r = np.exp(sol.t)
    u = sol.y[0]
    Q = sol.y[1]
    boundary_error = abs(u[-1] - u_amb) / max(abs(u_amb), abs(uc), 1e-300)
    rho = rho_to_m2(np.array([rho_prof(float(radius)) for radius in r]))
    res = max(boundary_error, _radial_balance(r, u, rho, rho_to_m2(rho_amb), lam_m2)[-1])
    if not np.all(np.isfinite(u)) or res > max(tol, 1e-5):
        raise RuntimeError(f"Shooting profile failed boundary/flux checks: residual {res:.3e}")
    return r, u, it, res, Q


def verify_residual(r, u, Q, rho_src):
    """Discrete residual: |dQ/dr - src r^2| normalised by the charge
    scale, plus flux-shape check f(u') r^2 + Q == 0."""
    drf = np.diff(r)
    vf = np.diff(u) / drf
    rf = 0.5 * (r[1:] + r[:-1])
    flux_err = np.abs(rf ** 2 * vf * P_X(vf)
                      + 0.5 * (Q[1:] + Q[:-1]))
    qscale = max(float(np.max(np.abs(Q))), 1e-300)
    return float(np.max(flux_err / qscale))


def retention_radius(u_eq_in, rho_in_m2):
    """Radius where R * v_surf(R) = u_eq, v_surf = finv(rho~ R/3).

    Operator-generic form of the leakage -> retention criterion: the
    interior drop u_eq fits inside the body only when the flux-limited
    surface gradient sheds it over < R.
    """
    f = lambda R: R * float(finv(np.array(rho_in_m2 * R / 3.0))) - u_eq_in
    lo, hi = 1e-3, 1e30
    while f(lo) > 0:
        lo /= 10.0
    while f(hi) < 0:
        hi *= 10.0
    return float(brentq(f, lo, hi, rtol=1e-10))


def diagnostics(r, u, R_m, rho_prof, rho_amb, M_kg):
    """Read the two screening branches off the solved profile.

    density branch: interior u vs the matter-coupled equilibrium u_eq.
    shear branch:   the flux leaving the surface, Q_surf = r^2 P_X u',
                    versus the bare source charge Q_bare = rho R^3/3
                    (charge dressing), and the surface scalar-to-
                    Newtonian acceleration ratio versus the canonical
                    cubic response y(g_N/g_t) of the S_Sigma operator.
    """
    rho = rho_prof(r)
    drf = np.diff(r)
    vf = np.diff(u) / drf                      # face du/dr
    f_face = vf * P_X(vf)
    Q = (0.5 * (r[1:] + r[:-1])) ** 2 * f_face  # flux at faces [m]
    inside = r < R_m
    inner = r < 0.35 * R_m
    rho_in = float(np.median(rho[inside]))
    u_int = float(np.median(u[inner]))
    u_eq_in = u_eq(rho_in)
    u_amb = u_eq(rho_amb)

    # interior Compton length at the would-be equilibrium vs R:
    # the CHAMELEON-style retention criterion (not the binding one for
    # this equation — see tail_len below)
    m2_eq = 3.0 * LAM_M2 * u_eq_in ** 2 + rho_to_m2(rho_in)
    lam_c = m2_eq ** -0.5

    # flux leaving the body (last interior face) vs bare charge
    iR = int(np.searchsorted(r, R_m))
    iR = min(max(iR, 1), len(Q) - 1)
    Q_surf = abs(Q[iR - 1])
    Q_bare = rho_to_m2(rho_in) * R_m ** 3 / 3.0
    g_N = G * M_kg / R_m ** 2
    x_surf = g_N / G_T
    # shear response y = a_s/(2 g_N) solves the flux map at x = g_N/g_t:
    #   baseline:   y (1 + y^2 x^2) = 1
    #   two-branch: KP x y^2 + x^2 y^3 = 1
    if OPERATOR == "two_branch":
        y_op = brentq(lambda yy: KP * x_surf * yy ** 2
                      + x_surf ** 2 * yy ** 3 - 1.0,
                      1e-40, 1e6)
    else:
        y_op = brentq(lambda yy: yy * (1.0 + yy ** 2 * x_surf ** 2) - 1.0,
                      1e-40, 1.0)
    a_surf = C ** 2 * abs(vf[iR - 1])
    # tail-matched plateau: in the leaking regime the interior excess
    # is set by the exterior-tail integral Delta u = int_R^inf finv(Q/r^2)dr
    _rt = R_m * np.geomspace(1.0, 1e10, 4000)
    _vt = finv(Q_bare / _rt ** 2)
    plateau_pred = float(np.sum(0.5 * (_vt[1:] + _vt[:-1])
                                * np.diff(_rt)))
    # the ACTUAL retention criterion for this equation: the u_eq -> tail
    # drop must fit inside R at the flux-limited surface gradient
    v_surf = float(finv(np.array(Q_bare / R_m ** 2)))
    tail_len = u_eq_in / v_surf       # radial distance to shed u_eq
    retention_ratio = tail_len / R_m  # >>1 => leaking branch forced
    # saturation radius: solve R finv(rho~ R/3) = u_eq for R
    R_c = retention_radius(u_eq_in, rho_to_m2(rho_in))
    return {
        "rho_in_g_cm3": rho_in,
        "u_eq_analytic": u_eq_in,
        "u_inner_median": u_int,
        "u_int_over_u_eq": u_int / u_eq_in,
        "u_int_over_u_amb": u_int / u_amb,
        "compton_len_at_eq_over_R": lam_c / R_m,
        "retention_ratio_tailLen_over_R": retention_ratio,
        "R_c_over_R": float(R_c / R_m),
        "saturated": bool(retention_ratio < 1.0),
        "S_A_from_solve": u_int / u_eq(RHO_T),
        "S_A_canonical": min(1.0, (rho_in / RHO_T) ** (1.0 / 3.0)),
        "g_N_surface": g_N,
        "g_N_over_g_t": g_N / G_T,
        "y_cubic_at_g_N": float(y_op),
        "a_surf_over_2gN": float(a_surf / (2.0 * g_N)),
        "Q_surf_over_Q_bare": float(Q_surf / Q_bare),
        "plateau_pred_4sqrtQg_t/c": float(plateau_pred),
        "u_int_minus_u_amb": u_int - u_amb,
    }


def main():
    global OPERATOR
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--operator", choices=["baseline", "two_branch"],
                    default="two_branch",
                    help="kinetic-sector operator: promoted two-branch "
                         "(default) or former unit-term baseline")
    OPERATOR = ap.parse_args().operator

    rho_amb = 1e-30                      # nearly empty ambient: clean tail
    u_amb = u_eq(rho_amb)
    u_T = u_eq(RHO_T)
    print("=" * 72)
    print("UNIFIED PROFILE: one solve, two screening branches")
    print(f"  u_amb(rho=1e-30)={u_amb:.3e}  u_eq(rho_T)={u_T:.3e}")
    print("=" * 72)

    bodies = [
        ("Cepheid envelope (8 Msun, 50 Rsun)", 8.0 * MSUN, 50.0 * RSUN,
         8.0 * MSUN * 1000.0 /
         (4.0 / 3.0 * np.pi * (50.0 * RSUN * 100.0) ** 3)),
        ("Sun", MSUN, RSUN, 1.408),
        ("Earth", MEARTH, REARTH, 5.51),
        ("Saturation body rho=rho_T", 1.0e30, 8.0e7, RHO_T),
        ("Dense body rho=1e4", 1.0e30, 3.0e6, 1.0e4),
    ]
    # toy body beyond the retention radius R_c = [3 c^4/(lam~ g_t^2)]^{1/4}
    # ~ 1.2e14 m: the saturated branch exists and is selected by SIZE,
    # not density — its tail is short (retention_ratio ~ 6e-3), so a
    # far smaller box suffices
    toy = ("Toy retained body R=3e15", 3.0e46, 3.0e15, RHO_T)

    eq_px = ("(k/sqrt2)(c|u'|/H0)+(c u'/H0)^2"
             if OPERATOR == "two_branch" else "1+(c u'/H0)^2")
    out = {"step": "step_53_unified_profile",
           "operator": OPERATOR,
           "kinetic_coefficients": {"k": float(K2), "KP": float(KP),
                                    "a0": A0, "g_t": G_T, "H0_s_inverse": H0,
                                    "gradient_unit_m_inverse": H0/C,
                                    "xi_mapping": "xi=(c u_prime/H0)^2/2; Lambda_X^2=M_Pl H0"},
           "equation": (f"(1/r^2)d/dr[r^2 P_,X u'] = dV_eff/du; "
                        f"P_,X={eq_px}; V_eff=lambda u^4/4+rho e^-u"),
           "ambient": {"rho_amb": rho_amb, "u_amb": float(u_amb)},
           "solver": "conservative_radial_finite_volume_damped_newton",
           "potential_lambda": LAMBDA,
           "residual_definition": "finite-volume equation balance, with floating-point flux roundoff accounted for",
           "boundary_conditions": "zero central flux and fixed outer field value",
           "bodies": {}}
    for name, M, R, rho_bar in bodies + [toy]:
        retained = name.startswith("Toy")
        # smooth density edge (logistic): keeps the flux equation
        # consistent at the boundary layer; the retained case uses a
        # softer edge since its shell is resolved anyway
        edge_pow = 8.0 if retained else 40.0

        def prof(rr, rb=rho_bar, RR=R, ep=edge_pow):
            return smooth_density(rr, RR, rb, rho_amb, ep)
        # the leaking tail decays over lam_c(u_plateau) ~ 1e7-1e9 R;
        # r_out must sit beyond it or the Dirichlet cap distorts the
        # tail-matched plateau level.  The retained body's shell/tail
        # sits inside ~R so a 1e6 R box is ample there.
        r_out_fac = 1e6 if retained else 1e10
        r, u, it, res, Q = solve_body(M, R, prof, rho_amb=rho_amb,
                                      u_amb=u_amb, n=20000,
                                      r_out_fac=r_out_fac,
                                      rtol_trial=1e-3 if retained
                                      else 1e-5,
                                      n_bisect=40 if retained
                                      else 60)
        Fnorm = verify_residual(r, u, Q, rho_to_m2(prof(r)))
        d = diagnostics(r, u, R, prof, rho_amb, M)
        d.update({"M_kg": M, "R_m": R, "rho_bar": rho_bar,
                  "iters": it, "final_res": res, "final_Fnorm": Fnorm,
                  "r_out_over_R": float(r[-1] / R),
                  "outer_boundary_error": float(u[-1]-u_amb),
                  "reached_outer_boundary": bool(np.isclose(r[-1], R*r_out_fac, rtol=1e-14))})
        out["bodies"][name] = d
        print(f"\n{name}  rho={rho_bar:.3e}  it={it}  res={res:.1e}  "
              f"|F|={Fnorm:.1e}")
        print(f"  u_int/u_eq={d['u_int_over_u_eq']:.4e}  "
              f"lam_c/R={d['compton_len_at_eq_over_R']:.2e}  "
              f"tail/R={d['retention_ratio_tailLen_over_R']:.2e}  "
              f"sat={d['saturated']}")
        print(f"  u_int/u_amb={d['u_int_over_u_amb']:.3e}  "
              f"S_A_solve={d['S_A_from_solve']:.3e}  "
              f"S_A_can={d['S_A_canonical']:.3f}")
        print(f"  g_N/g_t={d['g_N_over_g_t']:.2e}  "
              f"a/2g_N={d['a_surf_over_2gN']:.3e}  "
              f"y_cubic={d['y_cubic_at_g_N']:.3e}  "
              f"Q_surf/Q_bare={d['Q_surf_over_Q_bare']:.3e}")
        print(f"  plateau: u_int-u_amb={d['u_int_minus_u_amb']:.3e}  "
              f"pred={d['plateau_pred_4sqrtQg_t/c']:.3e}")

    # --- branch sweep: fixed density rho_T, R spanning R_c -------------
    # Separates size selection from density dependence: the leakage ->
    # retention transition should occur at R ~ R_c regardless of the
    # body's other parameters.
    rho_s = RHO_T
    u_eq_s = u_eq(rho_s)
    R_c_s = retention_radius(u_eq_s, rho_to_m2(rho_s))
    print("\n" + "=" * 72)
    print(f"BRANCH SWEEP at rho={rho_s} g/cm3:  R_c={R_c_s:.3e} m")
    out["branch_sweep"] = {"rho_g_cm3": rho_s, "R_c_m": float(R_c_s),
                           "points": []}
    for mult in (0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0):
        R = mult * R_c_s
        M = rho_s * 4.0 / 3.0 * np.pi * (R * 100.0) ** 3 / 1000.0

        def prof(rr, rb=rho_s, RR=R):
            return smooth_density(rr, RR, rb, rho_amb, 8.0)
        big = mult >= 1.0
        r, u, it, res, Q = solve_body(
            M, R, prof, rho_amb=rho_amb, u_amb=u_amb, n=8000,
            r_out_fac=1e6 if big else 1e10,
            rtol_trial=1e-3 if big else 1e-5,
            n_bisect=40 if big else 60)
        d = diagnostics(r, u, R, prof, rho_amb, M)
        out["branch_sweep"]["points"].append({
            "R_over_R_c": mult, "R_m": R,
            "u_int_over_u_eq": d["u_int_over_u_eq"],
            "S_A_from_solve": d["S_A_from_solve"],
            "S_A_canonical": d["S_A_canonical"],
            "Q_surf_over_Q_bare": d["Q_surf_over_Q_bare"],
            "final_res": res})
        print(f"  R={mult:>5.2f} R_c: u_int/u_eq={d['u_int_over_u_eq']:.3e} "
              f" S_A={d['S_A_from_solve']:.3e} "
              f"Q_surf/Q_bare={d['Q_surf_over_Q_bare']:.3e} res={res:.1e}")

    dens = np.array([b["rho_bar"] for b in out["bodies"].values()])
    uint = np.array([b["u_inner_median"] for b in out["bodies"].values()])
    sel = np.isfinite(uint) & (uint > 0)
    if sel.sum() >= 3:
        p = np.polyfit(np.log10(dens[sel]), np.log10(uint[sel]), 1)
        out["density_branch_fit"] = {"exponent": float(p[0]),
                                     "expected": 1.0 / 3.0}
        print(f"\ndensity branch: u_int ~ rho^{p[0]:.3f} (expect 1/3)")

    dest = Path(__file__).resolve().parents[2] / "results" / \
        f"step_53_unified_profile_{OPERATOR}.json"
    dest.write_text(json.dumps(out, indent=2))
    print("wrote", dest)


if __name__ == "__main__":
    main()
