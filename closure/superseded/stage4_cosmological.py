"""Stage 4: Cosmological continuation and diffuse regime.

Tests whether the cosh V permits:
1. A viable cosmological evolution (Friedmann + scalar)
2. The temporal horizon asymptotic (phi -> +inf, A -> 0)
3. A diffuse/void unscreened branch

The cosmological scalar equation in a homogeneous FRW universe:
  phi_ddot + 3H phi_dot + V_{,phi} = -rho_m A_{,phi} = (rho_m A)/M_Pl

In the quasi-static limit (slow roll):
  3H phi_dot + V_{,phi} = (rho_m A)/M_Pl

The Friedmann equation:
  H^2 = (8 pi/3 M_Pl^2) (rho_m + rho_V)
  where rho_V = V(phi) + V_ref

For the temporal horizon: phi -> +inf, A -> 0, rho_m A -> 0
The field rolls up the potential, driven by the matter coupling.
As A -> 0, the matter coupling weakens, and the field asymptotes.

Key question: does the cosh V support a rolling solution toward phi -> +inf?

The cosh V at large phi: V ~ (Lambda^4/2) exp(phi/M_Pl) -> infinity
This is a STEEP WALL. The field cannot roll to phi -> +inf because
the potential gradient pushes it back.

CONCLUSION: The cosh V does NOT naturally produce the temporal horizon.
The TH requires a potential that DECREASES toward large phi (or saturates),
so the matter coupling can push the field forward.

This is a fundamental finding: the cosh V solves the LOCAL closure but
NOT the global/cosmological closure. A different V (or a composite V)
is needed for the TH.

However, the cosh V is still the correct LOCAL EFT. The global V may
differ from the local V if there are additional operators (K != 1, or
higher-order terms) that become important at large phi.
"""

import numpy as np
import sys
sys.path.insert(0, '/Users/matthewsmawfield/www/Temporal Equivalence Principle/closure')

from radial_bvp_dimless import BETA_A, M_Pl, to_natural_density


def analyze_cosmological_continuation():
    """Analyze whether cosh V supports the temporal horizon."""
    print("=" * 80)
    print("COSMOLOGICAL CONTINUATION ANALYSIS")
    print("=" * 80)
    print()

    print("The cosh potential: V = Lambda^4 (cosh(phi/M_Pl) - 1)")
    print()
    print("At large phi/M_Pl:")
    for u in [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
        V = np.cosh(u) - 1
        Vp = np.sinh(u)
        Vpp = np.cosh(u)
        print(f"  u = {u:5.1f}: V/Lambda^4 = {V:12.4e}, "
              f"V_phi*M_Pl/Lambda^4 = {Vp:12.4e}, "
              f"V_phiphi*M_Pl^2/Lambda^4 = {Vpp:12.4e}")

    print()
    print("The potential gradient V_{,phi} = (Lambda^4/M_Pl) sinh(u) -> infinity")
    print("as u -> infinity. This is a STEEP WALL.")
    print()
    print("The matter driving force: rho_m * A / M_Pl = rho_m * exp(-u) / M_Pl")
    print("This DECREASES as u -> infinity (A -> 0).")
    print()
    print("Equilibrium: sinh(u) = (rho_m/Lambda^4) * exp(-u)")
    print("=> u_eq = (1/2) ln(1 + 2*rho_m/Lambda^4)")
    print()
    print("As rho_m -> 0 (cosmic expansion): u_eq -> 0")
    print("The field returns to the vacuum, NOT to the temporal horizon.")
    print()
    print("CONCLUSION: The cosh V does NOT produce the temporal horizon.")
    print("The TH requires phi -> +infinity, but the cosh V pushes phi back")
    print("to 0 as the universe expands. The TH needs a different mechanism.")
    print()

    # What potential WOULD produce the TH?
    print("=" * 80)
    print("WHAT POTENTIAL PRODUCES THE TEMPORAL HORIZON?")
    print("=" * 80)
    print()
    print("Requirements for TH:")
    print("  1. V_{,phi} > 0 for phi > 0 (local screening, equilibrium condition)")
    print("  2. V -> finite or 0 as phi -> +inf (field can roll forward)")
    print("  3. V_{,phi} -> 0 as phi -> +inf (driving force dominates at large phi)")
    print("  4. V_{,phi phi} > 0 at equilibrium (stable)")
    print()
    print("Conditions 1 and 2 together require V to be increasing but bounded.")
    print("This means V_{,phi} > 0 but V_{,phi} -> 0, so V_{,phi phi} < 0")
    print("at large phi (concave). This is UNSTABLE at large phi.")
    print()
    print("The resolution: the TH is NOT a static equilibrium.")
    print("It is a DYNAMICAL ATTRACTOR — the field rolls toward phi -> +inf")
    print("with phi_dot > 0, never reaching equilibrium.")
    print("The stability condition m_eff^2 > 0 applies to perturbations")
    print("around the ROLLING trajectory, not around a static point.")
    print()
    print("Candidate: V = Lambda^4 (1 - exp(-phi/M_Pl))")
    print("  V_{,phi} = (Lambda^4/M_Pl) exp(-phi/M_Pl) > 0  [GOOD]")
    print("  V -> Lambda^4 as phi -> inf  [FINITE — TH possible]")
    print("  V_{,phi} -> 0 as phi -> inf  [driving force dominates]")
    print("  V_{,phi phi} = -(Lambda^4/M_Pl^2) exp(-phi/M_Pl) < 0  [CONCAVE]")
    print()
    print("  At equilibrium: (Lambda^4/M_Pl) exp(-u) = (rho_m/M_Pl) exp(-u)")
    print("  => Lambda^4 = rho_m  [NO density-dependent equilibrium!]")
    print("  This potential has NO density-dependent minimum. It doesn't screen.")
    print()

    # The fundamental tension
    print("FUNDAMENTAL TENSION:")
    print("  Local screening needs: V_{,phi} = rho_m * A / M_Pl (density-dependent)")
    print("  This requires V_{,phi} to track rho_m, which means the equilibrium")
    print("  phi depends on rho_m. This is the density-mass relation.")
    print()
    print("  TH needs: V -> finite as phi -> inf, V_{,phi} -> 0")
    print("  This means V is bounded, so V_{,phi} decreases at large phi.")
    print()
    print("  A bounded V with density-dependent equilibrium:")
    print("  V = Lambda^4 * f(phi/M_Pl) where f is bounded and f' > 0")
    print("  The equilibrium f'(u) = (rho_m/Lambda^4) exp(-u)")
    print("  As rho_m -> 0: f'(u) -> 0, so u -> u_max where f'(u_max) = 0")
    print("  But f' > 0 everywhere means f' never reaches 0 at finite u.")
    print("  So u -> infinity as rho_m -> 0. THIS IS THE TH!")
    print()
    print("  The key: f must be BOUNDED with f' > 0 and f' -> 0 as u -> inf.")
    print("  Then as rho_m -> 0, the equilibrium u -> infinity (TH).")
    print("  And for finite rho_m, the equilibrium is at finite u (screening).")
    print()

    # Test candidate: V = Lambda^4 * tanh(phi/M_Pl)
    print("CANDIDATE: V = Lambda^4 * tanh(phi/M_Pl)")
    print("  f(u) = tanh(u), f'(u) = sech^2(u) > 0, f -> 1 as u -> inf")
    print("  Equilibrium: sech^2(u) = (rho_m/Lambda^4) exp(-u)")
    print()
    for rho_ratio in [1e-5, 1e-3, 1e-1, 1.0, 10.0, 100.0]:
        # Solve sech^2(u) = rho_ratio * exp(-u)
        # 1/cosh^2(u) = rho_ratio * exp(-u)
        # exp(-u) * cosh^2(u) = 1/rho_ratio
        # cosh^2(u) = exp(u) / rho_ratio
        # (exp(u) + exp(-u))^2 / 4 = exp(u) / rho_ratio
        # For large u: exp(2u)/4 ~ exp(u)/rho_ratio => exp(u) ~ 4/rho_ratio
        # u ~ ln(4/rho_ratio)
        from scipy.optimize import brentq
        def eq(u):
            return 1.0/np.cosh(u)**2 - rho_ratio * np.exp(-u)
        try:
            u_eq = brentq(eq, 0.01, 50.0)
            m_eff_sq = 1.0/np.cosh(u_eq)**2 * np.tanh(u_eq)  # -f'' = -(-2 sech^2 tanh) = 2 sech^2 tanh
            # Actually f'' = -2 sech^2(u) tanh(u), so V_{,phi phi} = -2 Lambda^4 sech^2(u) tanh(u) / M_Pl^2
            # This is NEGATIVE for u > 0. UNSTABLE.
            Vpp = -2.0 / np.cosh(u_eq)**2 * np.tanh(u_eq)
            print(f"  rho/Lambda^4 = {rho_ratio:10.4e}: u_eq = {u_eq:8.4f}, "
                  f"V_{,phi phi}*M_Pl^2/Lambda^4 = {Vpp:8.4f} "
                  f"{'UNSTABLE' if Vpp < 0 else 'stable'}")
        except:
            print(f"  rho/Lambda^4 = {rho_ratio:10.4e}: no equilibrium found")

    print()
    print("  PROBLEM: tanh has V_{,phi phi} < 0 for phi > 0. UNSTABLE.")
    print("  This is the generic problem: bounded + increasing = concave = unstable.")
    print()

    # The resolution: composite potential
    print("RESOLUTION: COMPOSITE POTENTIAL")
    print()
    print("  V = Lambda^4 (cosh(u) - 1)  for u < u_*  (local screening, stable)")
    print("  V = Lambda^4 (cosh(u_*) - 1) + Lambda^4 sinh(u_*) (u - u_*)")
    print("       + higher-order terms  for u > u_*  (TH branch)")
    print()
    print("  Or: V = Lambda^4 [cosh(u) - 1] * g(u)")
    print("  where g(u) -> 0 as u -> inf, softening the growth.")
    print()
    print("  Or: the TH is produced by a DIFFERENT sector (noncanonical K,")
    print("  or a second field, or a geometric effect not captured by V alone).")
    print()
    print("  The Jakarta manuscript already states that V is 'part of the")
    print("  action-closure problem.' The cosh V closes the LOCAL problem.")
    print("  The GLOBAL/TH closure requires additional structure.")
    print("  This is NOT a failure — it is the expected scope of a local EFT.")


def analyze_diffuse_regime():
    """Analyze the diffuse/void regime with proper nondimensionalization."""
    print()
    print("=" * 80)
    print("DIFFUSE REGIME ANALYSIS")
    print("=" * 80)
    print()

    # The diffuse regime: rho ~ 10^-29 g/cm^3, R ~ 1 Mpc
    # The issue was that eta_env is enormous and the field is pinned at u=0
    # The physical question: what is the field in the diffuse regime?

    rho_diffuse = to_natural_density(1e-29)  # GeV^4
    R_mpc = 3.086e24 * 5.07e13  # 1 Mpc in GeV^-1
    lam_diffuse = rho_diffuse * R_mpc**2 / M_Pl**2

    print(f"Diffuse void: rho = 1e-29 g/cm^3 = {rho_diffuse:.4e} GeV^4")
    print(f"R = 1 Mpc = {R_mpc:.4e} GeV^-1")
    print(f"lambda = {lam_diffuse:.4e}")
    print()

    # For the cosh V, the equilibrium in the diffuse regime:
    # u_eq = (1/2) ln(1 + 2 rho / Lambda^4)
    # If Lambda ~ 1 GeV (eta_earth = 1 => Lambda^4 = M_Pl^2/R_earth^2):
    # Lambda^4 = M_Pl^2 / R_earth^2 = (2.435e18)^2 / (3.23e22)^2 = 5.93e36 / 1.04e45 = 5.7e-9 GeV^4
    # rho_diffuse = 8.6e-47 GeV^4 (from 1e-29 g/cm^3)
    # rho/Lambda^4 = 8.6e-47 / 5.7e-9 = 1.5e-38
    # u_eq = (1/2) ln(1 + 3e-38) ~ 1.5e-38

    Lambda4 = M_Pl**2 / (to_natural_density(5.51) * to_natural_length(6.371e8)**2 / 1.0)
    # Actually: eta_earth = Lambda^4 R_earth^2 / M_Pl^2 = 1
    # => Lambda^4 = M_Pl^2 / R_earth^2
    R_earth_nat = to_natural_length(6.371e8)
    Lambda4 = M_Pl**2 / R_earth_nat**2
    Lambda = Lambda4**0.25

    print(f"For eta_earth = 1: Lambda^4 = M_Pl^2/R_earth^2 = {Lambda4:.4e} GeV^4")
    print(f"Lambda = {Lambda:.4e} GeV = {Lambda*1e9:.4e} eV")
    print()

    rho_ratio = rho_diffuse / Lambda4
    u_eq_diffuse = 0.5 * np.log(1 + 2 * rho_ratio)
    print(f"rho_diffuse / Lambda^4 = {rho_ratio:.4e}")
    print(f"u_eq (diffuse) = {u_eq_diffuse:.4e}")
    print(f"This is the ambient field value in the diffuse regime.")
    print()

    # The field in the diffuse regime is u ~ 10^-38, essentially zero.
    # This means the diffuse regime is UNSCREENED (field at vacuum).
    # The screening factor S_Sigma ~ 1 (no screening).
    # This is correct: in the void, there's no body to screen.
    print("The diffuse field is at u ~ 10^-38, essentially the vacuum.")
    print("This is correct: in a void, there is no body to screen.")
    print("The ambient field is set by the cosmological solution,")
    print("which for the cosh V is u -> 0 as rho -> 0.")
    print()
    print("For the TH, we need u -> infinity, but the cosh V gives u -> 0.")
    print("This confirms: the cosh V does NOT produce the TH.")
    print("The TH requires a different global structure.")
    print()

    # What about the absorber regime (Paper 29)?
    # Paper 29 uses rho ~ 10^-24 g/cm^3 (ISM density)
    rho_ism = to_natural_density(1e-24)
    rho_ratio_ism = rho_ism / Lambda4
    u_eq_ism = 0.5 * np.log(1 + 2 * rho_ratio_ism)
    print(f"ISM density: rho = 1e-24 g/cm^3 = {rho_ism:.4e} GeV^4")
    print(f"rho_ism / Lambda^4 = {rho_ratio_ism:.4e}")
    print(f"u_eq (ISM) = {u_eq_ism:.4e}")
    print()
    print("In the ISM, the field is u ~ 10^-30, still essentially vacuum.")
    print("The disformal deformation B*(dphi)^2 is suppressed by u^2 ~ 10^-60.")
    print("This is far more suppressed than the constant-B form Paper 29 tested.")
    print("Paper 29 found constant B insufficient; quartic-Gaussian is worse.")
    print("CONSISTENT: both agree perturbative disformal can't close the gap.")


if __name__ == '__main__':
    analyze_cosmological_continuation()
    analyze_diffuse_regime()
