"""Radial BVP solver for the TEP scalar field equation.

Solves the quasistatic scalar equation

    phi'' + (2/r) phi' = V_{,phi}(phi) + rho_*(r) A_{,phi}(phi)

with A(phi) = exp(beta_A phi / M_Pl), beta_A = -1.

This is the TEP scalar equation, NOT a chameleon thin-shell calculation.
Screening emerges from the solved profile phi(r) and the Temporal Topology
operator, not from a density-dependent mass or geometric thin-shell formula.

Observables are extracted directly from the solved profile:
  - Sigma_r(r) = -phi'(r) / M_Pl           (Temporal Shear)
  - Q from exterior: phi(r) ~ phi_env - Q/(4 pi r)  for r >> R
  - S_Sigma = Q / Q_0  where Q_0 = alpha_0 M / M_Pl (unscreened charge)
  - epsilon_B(r) = (B/A^2) |phi'|^2         (disformal deformation)

No thin-shell estimates. No Khoury-Weltman validation. No density-dependent
mass as the screening mechanism. The Temporal Topology is the ontology;
V(phi) is the engine that generates phi(x).
"""

import numpy as np
from scipy.integrate import solve_bvp
from scipy.optimize import brentq
import json

# ---------------------------------------------------------------------------
# Physical constants (natural units, c = hbar = 1)
# ---------------------------------------------------------------------------

M_Pl = 2.435e18       # GeV, reduced Planck mass
GeV_per_g = 5.61e23   # conversion
GeV_inv_per_cm = 5.07e13

def to_natural_density(rho_g_cc):
    """Convert density from g/cm^3 to GeV^4."""
    return rho_g_cc * GeV_per_g / GeV_inv_per_cm**3

def to_natural_length(cm):
    """Convert length from cm to GeV^-1."""
    return cm * GeV_inv_per_cm

def to_meters(GeV_inv):
    """Convert GeV^-1 to meters."""
    return GeV_inv / (GeV_inv_per_cm * 100)

# ---------------------------------------------------------------------------
# Coupling functions (frozen)
# ---------------------------------------------------------------------------

BETA_A = -1.0

def A_factor(phi, M_Pl_val=M_Pl):
    """Conformal factor A(phi) = exp(beta_A phi / M_Pl)."""
    return np.exp(BETA_A * phi / M_Pl_val)

def A_phi(phi, M_Pl_val=M_Pl):
    """dA/dphi = (beta_A / M_Pl) A."""
    return (BETA_A / M_Pl_val) * A_factor(phi, M_Pl_val)

# ---------------------------------------------------------------------------
# Potential families
# ---------------------------------------------------------------------------

class Potential:
    """Base class for scalar self-interaction V(phi).

    Each subclass provides V(phi), V_phi(phi), V_phiphi(phi).
    All in natural units (GeV).
    """

    def V(self, phi):
        raise NotImplementedError

    def V_phi(self, phi):
        raise NotImplementedError

    def V_phiphi(self, phi):
        raise NotImplementedError

    def describe(self):
        raise NotImplementedError


class LinearPotential(Potential):
    """V = m^2 phi^2 / 2.  For Yukawa validation only.

    Gives phi(r) ~ exp(-m r) / r for a point source.
    No chameleon mechanism — just Yukawa suppression.
    """

    def __init__(self, m_GeV):
        self.m = m_GeV
        self.m2 = m_GeV**2

    def V(self, phi):
        return 0.5 * self.m2 * phi**2

    def V_phi(self, phi):
        return self.m2 * phi

    def V_phiphi(self, phi):
        return self.m2

    def describe(self):
        return f"Linear V = (1/2) m^2 phi^2, m = {self.m:.4e} GeV"


class ExponentialPotential(Potential):
    """V = Lambda^4 exp(phi / M_Pl).

    Simplest viable candidate:
      V_{,phi} > 0 for all phi (equilibrium condition satisfied)
      V_{,phi phi} > 0 (stable)
      No singularity
      Connects to phi -> +inf (A -> 0, temporal horizon)
      phi -> -inf: V -> 0 (flat, unscreened)

    Equilibrium: phi_min = (M_Pl/2) ln(rho_* / Lambda^4)
    Effective mass: m_eff^2 = 2 sqrt(Lambda^4 rho_*) / M_Pl^2

    One parameter (Lambda). NOT a chameleon potential — it is a
    TEP scalar self-interaction whose solved profile generates
    the Temporal Topology response.
    """

    def __init__(self, Lambda_GeV):
        self.Lambda = Lambda_GeV
        self.Lambda4 = Lambda_GeV**4

    def V(self, phi):
        return self.Lambda4 * np.exp(phi / M_Pl)

    def V_phi(self, phi):
        return (self.Lambda4 / M_Pl) * np.exp(phi / M_Pl)

    def V_phiphi(self, phi):
        return (self.Lambda4 / M_Pl**2) * np.exp(phi / M_Pl)

    def equilibrium(self, rho_star):
        """Adiabatic equilibrium phi_min(rho_*)."""
        return 0.5 * M_Pl * np.log(rho_star / self.Lambda4)

    def m_eff(self, rho_star):
        """Effective mass at equilibrium."""
        return np.sqrt(2) * self.Lambda**2 * rho_star**0.25 / M_Pl

    def describe(self):
        return f"Exponential V = Lambda^4 exp(phi/M_Pl), Lambda = {self.Lambda:.4e} GeV"


class CosPotential(Potential):
    """V = Lambda^4 (1 - cos(phi / M_Pl)).

    Periodic, bounded below. V_{,phi} = (Lambda^4/M_Pl) sin(phi/M_Pl).
    V_{,phi} > 0 for 0 < phi < pi M_Pl.
    Stable: V_{,phi phi} = (Lambda^4/M_Pl^2) cos(phi/M_Pl) > 0 for |phi| < pi M_Pl/2.
    Connects to phi -> +inf only through periodicity (not monotonic).
    Useful for testing a different local structure.
    """

    def __init__(self, Lambda_GeV):
        self.Lambda = Lambda_GeV
        self.Lambda4 = Lambda_GeV**4

    def V(self, phi):
        return self.Lambda4 * (1 - np.cos(phi / M_Pl))

    def V_phi(self, phi):
        return (self.Lambda4 / M_Pl) * np.sin(phi / M_Pl)

    def V_phiphi(self, phi):
        return (self.Lambda4 / M_Pl**2) * np.cos(phi / M_Pl)

    def describe(self):
        return f"Cosine V = Lambda^4 (1 - cos(phi/M_Pl)), Lambda = {self.Lambda:.4e} GeV"


# ---------------------------------------------------------------------------
# Source profiles
# ---------------------------------------------------------------------------

class SourceProfile:
    """Density profile rho_*(r) for a spherical body.

    All quantities in natural units (GeV).
    """

    def __init__(self, r_array, rho_array, name=""):
        self.r = r_array
        self.rho = rho_array
        self.name = name

    def rho_at(self, r):
        """Interpolate density at radius r. Handles scalar or array input."""
        r = np.asarray(r)
        scalar = r.ndim == 0
        r_flat = np.atleast_1d(r)
        result = np.interp(r_flat, self.r, self.rho)
        result = np.where(r_flat <= self.r[0], self.rho[0], result)
        result = np.where(r_flat >= self.r[-1], self.rho[-1], result)
        if scalar:
            return float(result[0])
        return result


def make_uniform_sphere(rho_core_gev4, R_gev_inv, r_max_factor=100.0,
                        rho_ambient=0.0, name="uniform"):
    """Uniform density sphere with ambient background."""
    r_max = R_gev_inv * r_max_factor
    # Log-spaced grid
    r = np.logspace(np.log10(R_gev_inv * 1e-4), np.log10(r_max), 500)
    rho = np.where(r <= R_gev_inv, rho_core_gev4, rho_ambient)
    return SourceProfile(r, rho, name)


def make_earth_profile():
    """Simplified Earth density profile.

    Core: ~13 g/cm^3, Mantle: ~4.5 g/cm^3, Crust: ~2.7 g/cm^3,
    transitioning to ambient at ~100 AU.
    """
    R_core = to_natural_length(3480e5)   # 3480 km in cm
    R_mantle = to_natural_length(5701e5) # 5701 km
    R_crust = to_natural_length(6371e5)  # 6371 km
    R_atmo = to_natural_length(6471e5)   # +100 km

    rho_core = to_natural_density(13.0)
    rho_mantle = to_natural_density(4.5)
    rho_crust = to_natural_density(2.7)
    rho_ambient = to_natural_density(1e-24)  # interplanetary

    r_max = to_natural_length(1.5e13)  # 1 AU in cm
    r = np.logspace(np.log10(R_core * 1e-3), np.log10(r_max), 800)

    rho = np.zeros_like(r)
    for i, ri in enumerate(r):
        if ri <= R_core:
            rho[i] = rho_core
        elif ri <= R_mantle:
            # Linear interpolation core -> mantle
            frac = (ri - R_core) / (R_mantle - R_core)
            rho[i] = rho_core * (1 - frac) + rho_mantle * frac
        elif ri <= R_crust:
            frac = (ri - R_mantle) / (R_crust - R_mantle)
            rho[i] = rho_mantle * (1 - frac) + rho_crust * frac
        elif ri <= R_atmo:
            frac = (ri - R_crust) / (R_atmo - R_crust)
            rho[i] = rho_crust * (1 - frac) + rho_ambient * frac
        else:
            rho[i] = rho_ambient

    return SourceProfile(r, rho, "Earth")


def make_sun_profile():
    """Simplified Sun density profile.

    Core: ~150 g/cm^3, mean: ~1.4 g/cm^3, transitioning to ambient.
    """
    R_sun = to_natural_length(6.96e10)  # cm
    R_core = to_natural_length(0.2 * 6.96e10)  # 20% of R_sun

    rho_core = to_natural_density(150.0)
    rho_mean = to_natural_density(1.4)
    rho_ambient = to_natural_density(1e-24)

    r_max = to_natural_length(1.5e13)  # 1 AU
    r = np.logspace(np.log10(R_core * 1e-3), np.log10(r_max), 800)

    rho = np.zeros_like(r)
    for i, ri in enumerate(r):
        if ri <= R_core:
            rho[i] = rho_core
        elif ri <= R_sun:
            # Power-law decrease (rough solar model)
            frac = (ri - R_core) / (R_sun - R_core)
            rho[i] = rho_core * (1 - frac)**3 + rho_mean * frac
        else:
            rho[i] = rho_ambient

    return SourceProfile(r, rho, "Sun")


# ---------------------------------------------------------------------------
# BVP solver
# ---------------------------------------------------------------------------

def solve_radial_bvp(potential, source, phi_ambient, r_factor=1.0,
                     verbose=False, tol=1e-6, max_nodes=50000):
    """Solve the radial scalar field BVP.

    Equation: phi'' + (2/r) phi' = V_{,phi}(phi) + rho_*(r) A_{,phi}(phi)

    Boundary conditions:
      phi'(r=0) = 0  (regularity)
      phi(r_max) = phi_ambient  (ambient matching)

    Returns: (r, phi, phi_prime, solution_object)
    """
    r_min = source.r[0]
    r_max = source.r[-1] * r_factor

    # Initial guess: constant phi = phi_ambient everywhere
    # with a small dip near the center where density is high
    r_init = np.logspace(np.log10(r_min), np.log10(r_max), 200)
    phi_init = np.full_like(r_init, phi_ambient)

    # The ODE system: y = [phi, phi'], y' = [phi', phi'']
    def ode(r, y):
        phi = y[0]
        phi_p = y[1]

        # Interpolate density
        rho = source.rho_at(r)

        # V_{,phi} and A_{,phi}
        Vp = potential.V_phi(phi)
        Ap = A_phi(phi)

        # phi'' = V_{,phi} + rho_* A_{,phi} - (2/r) phi'
        # Guard against r -> 0
        r_safe = np.maximum(r, r_min * 0.1)
        phi_pp = Vp + rho * Ap - (2.0 / r_safe) * phi_p

        return np.vstack([phi_p, phi_pp])

    def bc(ya, yb):
        # ya = [phi(0), phi'(0)], yb = [phi(r_max), phi'(r_max)]
        # phi'(0) = 0, phi(r_max) = phi_ambient
        return np.array([ya[1], yb[0] - phi_ambient])

    sol = solve_bvp(ode, bc, r_init, np.vstack([phi_init, np.zeros_like(phi_init)]),
                    tol=tol, max_nodes=max_nodes, verbose=verbose)

    if not sol.success:
        print(f"BVP solver failed: {sol.message}")
        return None, None, None, sol

    # Dense output
    r_dense = np.logspace(np.log10(r_min), np.log10(r_max), 2000)
    y_dense = sol.sol(r_dense)
    phi = y_dense[0]
    phi_p = y_dense[1]

    return r_dense, phi, phi_p, sol


# ---------------------------------------------------------------------------
# Observable extraction (direct from solved profile — no thin-shell)
# ---------------------------------------------------------------------------

def extract_observables(r, phi, phi_p, source, potential, body_mass_gev,
                        r_body_gev_inv, phi_ambient):
    """Extract TEP observables directly from the solved profile.

    No thin-shell estimates. Everything is read off the solution.
    """
    M_Pl_val = M_Pl

    # Temporal Shear: Sigma_r = -phi'/M_Pl (radial component)
    Sigma_r = -phi_p / M_Pl_val

    # Conformal factor profile
    A_profile = A_factor(phi)

    # Effective mass along the profile (local property, not screening mechanism)
    rho_interp = np.array([source.rho_at(ri) for ri in r])
    m_eff_sq = potential.V_phiphi(phi) + rho_interp * A_factor(phi) / M_Pl_val**2
    m_eff = np.sqrt(np.maximum(m_eff_sq, 0))

    # Exterior scalar charge Q from asymptotic fit
    # phi(r) ~ phi_ambient - Q/(4 pi r) for r >> R_body
    mask_ext = r > r_body_gev_inv * 5  # well outside the body
    if np.sum(mask_ext) > 10:
        # Fit: phi_ambient - phi(r) = Q / (4 pi r)
        # => Q = 4 pi r (phi_ambient - phi(r))
        Q_values = 4 * np.pi * r[mask_ext] * (phi_ambient - phi[mask_ext])
        # Use the most stable region (far enough out, not at boundary)
        Q = np.median(Q_values[len(Q_values)//4:len(Q_values)*3//4])
    else:
        Q = 0.0

    # Unscreened charge Q_0 = |d lnA/d phi| * M = |beta_A| * M / M_Pl
    # (dimensionful-phi convention; S_Sigma = Q/Q_0 is convention-free)
    Q_0 = abs(BETA_A) * body_mass_gev / M_Pl_val

    # Source-charge screening factor (directly from solved profile)
    S_Sigma = Q / Q_0 if Q_0 > 0 else 0.0

    # Disformal deformation (B envelope — weak-field limit B ~ B0 (phi/M_Pl)^2)
    # epsilon_B = (B/A^2) |phi'|^2
    # In weak field: B(phi) ~ B0 * (phi/M_Pl)^2 (envelope = 1)
    # So epsilon_B / B0 ~ (phi/M_Pl)^2 * |phi'|^2 / A^2
    # We compute the shape factor; B0 is scanned separately
    u = phi / M_Pl_val  # dimensionless field
    B_shape = u**2 / (1 + u**2)  # rational part (envelope = 1 in weak field)
    epsilon_B_over_B0 = B_shape * phi_p**2 / A_profile**2

    return {
        'r': r,
        'phi': phi,
        'phi_prime': phi_p,
        'Sigma_r': Sigma_r,
        'A': A_profile,
        'm_eff': m_eff,
        'rho': rho_interp,
        'Q': Q,
        'Q_0': Q_0,
        'S_Sigma': S_Sigma,
        'epsilon_B_over_B0': epsilon_B_over_B0,
        'u': u,
    }


# ---------------------------------------------------------------------------
# Conformal loop integral test (Theorem 2 in code)
# ---------------------------------------------------------------------------

def test_conformal_loop_integral(r, phi, phi_p):
    """Verify that oint d(ln A) = 0 around a closed radial loop.

    For a static radial field, ln A = beta_A phi / M_Pl.
    A closed loop out and back along the same radial path gives:
      oint d(ln A) = [ln A(r_out) - ln A(r_in)] + [ln A(r_in) - ln A(r_out)] = 0

    This is trivial for a radial path, but the test checks that the
    numerical integration of the connection reproduces this to
    machine precision. It validates the holonomy subtraction machinery.
    """
    M_Pl_val = M_Pl

    # Connection: d(ln A) = (beta_A / M_Pl) dphi = (beta_A / M_Pl) phi' dr
    connection = (BETA_A / M_Pl_val) * phi_p

    # Integrate outward
    from scipy.integrate import trapezoid
    integral_out = trapezoid(connection, r)

    # Integrate back (reverse direction)
    integral_back = trapezoid(-connection[::-1], r[::-1])

    loop_integral = integral_out + integral_back

    # Also check: ln A(r_max) - ln A(r_min) should equal integral_out
    ln_A_diff = BETA_A * (phi[-1] - phi[0]) / M_Pl_val

    return {
        'loop_integral': loop_integral,
        'integral_out': integral_out,
        'ln_A_diff': ln_A_diff,
        'loop_residual': abs(loop_integral),
        'integration_error': abs(integral_out - ln_A_diff),
    }


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def save_results(obs, filename, metadata=None):
    """Save observable results to JSON."""
    # Convert arrays to lists, handle non-serializable
    result = {}
    for key, val in obs.items():
        if isinstance(val, np.ndarray):
            result[key] = val.tolist()
        elif isinstance(val, (np.floating, np.integer)):
            result[key] = float(val)
        else:
            result[key] = val

    if metadata:
        result['_metadata'] = metadata

    with open(filename, 'w') as f:
        json.dump(result, f, indent=2, default=str)
    print(f"  Saved: {filename}")
