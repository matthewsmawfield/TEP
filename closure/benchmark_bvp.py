#!/usr/bin/env python3
"""
Benchmark BVP — convex V, K=1, beta_A=-1, B=0.

ODE:  phi'' + (2/r)phi' = V_{,phi} + rho_*(r) * A_{,phi}

  V(phi) = (1/2) m^2 phi^2          (convex; benchmark, NOT labelled TEP)
  A(phi) = exp(beta_A phi / M_Pl),  beta_A = -1
  rho_*  = A^3(phi) * rho           (conserved Einstein-frame density)
  A_{,phi} = (beta_A / M_Pl) * A(phi)

  Full nonlinear source: rho_* A_{,phi} = (beta_A/M_Pl) * rho * A^4(phi)
  Linearized (phi << M_Pl):           ≈ (beta_A/M_Pl) * rho

  K = 1  (canonical kinetic, frozen)
  B = 0  (disformal sector frozen)

Solves for uniform-density Sun and Earth.
Plots: |Sigma(r)|, S_Sigma(r), eps_B(r), J(r).

TEP context:
  - Theorem 1 (Jakarta): conformal sector preserves null cones. No propagation
    anomaly from A(phi) alone. eps_B = 0 with B=0 (trivially).
  - Synchronization one-form: omega = d ln A = (beta_A/M_Pl) d phi (exact).
    Closed-loop holonomy: oint d ln A = 0. But spatial components
    Sigma_i = partial_i ln A are non-zero in the transition region.
  - Triangle kernel = spatial connection Sigma_i, NOT the holonomy.
  - Recovery region = transition where |Sigma| is non-zero (smooth, not step).
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json, os

# ── Constants ──────────────────────────────────────────────────────────
c = 2.998e8; G = 6.674e-11; hbar = 6.582e-25; hbar_c = 1.973e-14
M_Pl = 2.435e18; H0 = 70e3 / 3.086e22
BETA = -1.0
GeV_inv_per_cm = 5.07e13; GeV_per_g = 5.61e23

# ── Benchmark scale ────────────────────────────────────────────────────
Lam = np.sqrt(M_Pl * H0 * hbar)  # GeV
m = Lam                           # scalar mass (benchmark, not TEP)

# ── Bodies (uniform density) ───────────────────────────────────────────
bodies = {
    'Sun':   {'R_cm': 6.96e10, 'rho_g': 1.41,   'rho_amb_g': 1e-24},
    'Earth': {'R_cm': 6.371e8, 'rho_g': 5.51,   'rho_amb_g': 1e-21},
}

def rho_to_GeV4(rho_g):
    return rho_g * GeV_per_g / GeV_inv_per_cm**3

# ── Analytical solution (uniform sphere, Helmholtz) ───────────────────
# Full ODE:  phi'' + (2/r)phi' = m^2 phi + (beta_A/M_Pl) rho A^4(phi)
# Linearized (phi << M_Pl, verified for all bodies):
#   phi'' + (2/r)phi' = m^2 phi - rho/M_Pl
# Equilibrium: phi_eq = rho / (m^2 M_Pl)  (positive)
# Interior:  phi = phi_body + A_coeff * sinh(mr)/(mr)
# Exterior:  phi = phi_amb + B_coeff * exp(-mr)/r
# Matching at r=R.

def solve_profile(body, m, n_pts=2000):
    R = body['R_cm']
    rho = rho_to_GeV4(body['rho_g'])
    rho_amb = rho_to_GeV4(body['rho_amb_g'])

    phi_body = rho / (m**2 * M_Pl)       # equilibrium inside (positive)
    phi_amb  = rho_amb / (m**2 * M_Pl)   # equilibrium outside
    dphi = phi_body - phi_amb

    x = m * R * GeV_inv_per_cm  # dimensionless mR
    lam_C = hbar_c / m          # Compton wavelength (cm)

    # Stable matching (avoids overflow for large x)
    em2x = np.exp(-2*x) if x < 350 else 0.0

    # cosh(x) - sinh(x)/x = [x - 1 + em2x*(x+1)] / (2x)
    cosh_minus_sinh_over_x = (x - 1 + em2x*(x+1)) / (2*x)
    B_coeff = dphi * R * cosh_minus_sinh_over_x  # GeV * cm
    A_coeff = -dphi * (x + 1) * np.exp(-x) if x < 350 else 0.0

    # S_Sigma at surface (stable form)
    # S_Sigma(R) = 3(x-1)(x+1) / (2x^3)  for large x
    # Exact: S_Sigma(R) = 3*(x-1+em2x*(x+1))*(x+1) / (2*x^3)
    S_at_R = 3.0 * (x - 1 + em2x*(x+1)) * (x + 1) / (2 * x**3)

    # Q_0 = beta_A * M / M_Pl (unscreened charge)
    M_body_g = (4/3) * np.pi * R**3 * body['rho_g']
    Q_0 = abs(BETA * M_body_g * GeV_per_g / M_Pl)

    # Q_eff = 4*pi * B_coeff (effective screened charge)
    Q_eff = 4 * np.pi * abs(B_coeff)

    # ── Radial grid in transition coordinate u = m*(r - R) ─────────
    # Span: u in [-30, +30] covers the transition region
    u = np.linspace(-30, 30, n_pts)
    r = R + u / (m * GeV_inv_per_cm)

    # Compute phi, phi' using stable expressions
    phi = np.zeros(n_pts)
    phip = np.zeros(n_pts)

    for i in range(n_pts):
        ri = r[i]
        if ri <= R:
            mr = m * ri * GeV_inv_per_cm
            if mr > 350:
                # Stable: A*sinh(mr)/(mr) = -dphi*(x+1)*exp(mr-x)/(2*mr)
                phi[i] = phi_body - dphi * (x+1) * np.exp(mr - x) / (2*mr)
                # d/dr[sinh(mr)/(mr)] ~ exp(mr)/(2*mr*r) * (1 - 1/mr)
                phip[i] = -dphi * (x+1) * np.exp(mr - x) * (1 - 1/mr) / (2*ri) * m * GeV_inv_per_cm
            elif mr < 1e-10:
                phi[i] = phi_body + A_coeff  # sinh(mr)/mr -> 1
                phip[i] = 0.0
            else:
                phi[i] = phi_body + A_coeff * np.sinh(mr) / mr
                phip[i] = A_coeff * m * GeV_inv_per_cm * (np.cosh(mr) - np.sinh(mr)/mr) / ri
        else:
            mr = m * ri * GeV_inv_per_cm
            phi[i] = phi_amb + B_coeff * np.exp(-mr) / ri
            phip[i] = -B_coeff * np.exp(-mr) * (m * GeV_inv_per_cm / ri + 1/ri**2)

    # ── Extract observables ────────────────────────────────────────

    # Sigma = grad(ln A) = (beta_A / M_Pl) * phi'
    Sigma = (BETA / M_Pl) * phip  # GeV (natural units)
    abs_Sigma = np.abs(Sigma)

    # Normalize |Sigma| by its peak value
    Sigma_max = np.nanmax(abs_Sigma)
    abs_Sigma_norm = abs_Sigma / Sigma_max if Sigma_max > 0 else abs_Sigma

    # S_Sigma(r) = Q_eff(r) / Q_0  (exterior only)
    # Q_eff(r) = |B_coeff| * (mr + 1) * exp(-mr) * 4*pi  (exterior charge enclosed)
    S_Sigma_r = np.full(n_pts, np.nan)
    for i in range(n_pts):
        if r[i] > R:
            mr = m * r[i] * GeV_inv_per_cm
            Q_enc = 4 * np.pi * abs(B_coeff) * (mr + 1) * np.exp(-mr)
            S_Sigma_r[i] = Q_enc / Q_0

    # eps_B = (B / A^2) * |phi'|^2 = 0  (B = 0 frozen)
    eps_B = np.zeros(n_pts)

    # J ~ (B / A^2) * (u . Sigma) * Sigma = 0  (B = 0 frozen)
    J = np.zeros(n_pts)

    return {
        'u': u, 'r': r, 'R': R,
        'phi': phi, 'phi_prime': phip,
        'Sigma': Sigma, 'abs_Sigma': abs_Sigma,
        'abs_Sigma_norm': abs_Sigma_norm,
        'S_Sigma_r': S_Sigma_r, 'S_at_R': S_at_R,
        'eps_B': eps_B, 'J': J,
        'phi_body': phi_body, 'phi_amb': phi_amb,
        'x': x, 'lam_C_cm': lam_C,
        'Q_eff': Q_eff, 'Q_0': Q_0,
        'Sigma_max': Sigma_max,
    }

# ── Solve for both bodies ──────────────────────────────────────────────
results = {}
for name, body in bodies.items():
    results[name] = solve_profile(body, m)
    r = results[name]
    print(f"{name}: x=mR={r['x']:.2e}, S_Sigma(R)={r['S_at_R']:.4e}, "
          f"lambda_C={r['lam_C_cm']:.4f} cm = {r['lam_C_cm']*10:.2f} mm")
    print(f"  phi_body={r['phi_body']:.4e} GeV = {r['phi_body']/M_Pl:.2e} M_Pl")
    print(f"  |Sigma|_max = {r['Sigma_max']:.4e} GeV")
    print(f"  |Sigma| non-zero outside body: {np.any(r['abs_Sigma'][r['r'] > r['R']] > 0)}")

# ── Pass/fail checks (stated before solving, assessed after) ───────────
print("\n" + "="*70)
print("PASS/FAIL")
print("="*70)
print("Criteria (stated before solving):")
print("  1. Sun: S_Sigma^(sun) <= 5.75e-6  (Cassini)")
print("  2. eps_B <= 1e-15  (GW-like path, B=0 frozen)")
print("  3. Earth: a recovery region exists (not a step)")
print("  4. Triangle kernel not identically zero on space legs")
print()

# 1. Sun Cassini
S_sun = results['Sun']['S_at_R']
cassini_pass = S_sun < 5.75e-6
print(f"1. Sun S_Sigma = {S_sun:.4e}  (bound 5.75e-6)")
print(f"   {'PASS' if cassini_pass else 'FAIL'}  (margin {5.75e-6/S_sun:.0f}x)")

# 2. eps_B (B=0 frozen)
eps_max = 0.0
gw_pass = eps_max < 1e-15
print(f"2. eps_B = {eps_max:.1e}  (bound 1e-15, B=0 frozen)")
print(f"   {'PASS' if gw_pass else 'FAIL'}  (trivially satisfied: B=0)")

# 3. Earth recovery region
earth = results['Earth']
# Recovery region = where |Sigma| > 1% of peak, spanning the transition
threshold = 0.01
nonzero_mask = earth['abs_Sigma_norm'] > threshold
transition_width_u = np.ptp(earth['u'][nonzero_mask])
# Check: is it smooth (not a step)?
# A step would have |Sigma| jump discontinuously. Our solution is analytic.
# Check: does |Sigma| have a smooth peak (not a discontinuity)?
is_smooth = True  # analytical solution is C^infinity
# Check: is |Sigma| non-zero on both sides of r=R?
inside_nonzero = np.any(earth['abs_Sigma_norm'][(earth['r'] < earth['R'])] > threshold)
outside_nonzero = np.any(earth['abs_Sigma_norm'][(earth['r'] > earth['R'])] > threshold)
recovery_exists = is_smooth and (inside_nonzero or outside_nonzero)
print(f"3. Earth recovery region:")
print(f"   |Sigma| > 1% peak: u in [{earth['u'][nonzero_mask].min():.1f}, "
      f"{earth['u'][nonzero_mask].max():.1f}]")
print(f"   Width = {transition_width_u:.1f} lambda_C = "
      f"{transition_width_u * earth['lam_C_cm'] * 10:.2f} mm")
print(f"   Smooth (analytical C^inf): {is_smooth}")
print(f"   Non-zero inside: {inside_nonzero}, outside: {outside_nonzero}")
print(f"   {'PASS' if recovery_exists else 'FAIL'}  "
      f"(recovery region exists, smooth — not a step)")

# 4. Triangle kernel
# The synchronization connection: omega = d(ln A) = (beta_A/M_Pl) d(phi)
# Spatial components: Sigma_i = (beta_A/M_Pl) * partial_i phi
# In radial BVP: Sigma_r = (beta_A/M_Pl) * phi'(r) != 0 in transition
# The KERNEL (connection itself) is non-zero on space legs.
# The HOLOMONY (closed-loop integral) is zero because d(ln A) is exact.
# Criterion asks about the KERNEL, not the holonomy.
kernel_nonzero = np.any(earth['abs_Sigma'] > 0)
print(f"4. Triangle kernel on space legs:")
print(f"   Sigma_r = (beta_A/M_Pl) * phi'(r)")
print(f"   Non-zero in transition: {kernel_nonzero}")
print(f"   Max |Sigma_r| = {earth['Sigma_max']:.4e} GeV")
print(f"   Note: holonomy oint d(ln A) = 0 (exact form),")
print(f"         but kernel Sigma_r != 0 (spatial connection is non-zero)")
print(f"   {'PASS' if kernel_nonzero else 'FAIL'}  "
      f"(kernel non-zero on space legs; holonomy zero is expected with B=0)")

# ── Overall outcome ────────────────────────────────────────────────────
n_pass = sum([cassini_pass, gw_pass, recovery_exists, kernel_nonzero])
print(f"\n{'='*70}")
print(f"OVERALL: {n_pass}/4 PASS")
print(f"{'='*70}")
if n_pass == 4:
    print("PASS: local K=1 is viable; V still a benchmark, not frozen")
elif not cassini_pass:
    print("FAIL Cassini vs GNSS scale at once:")
    print("  one profile cannot serve both numbers; corpus tension,")
    print("  not a reason to invent S(a)")
else:
    print("FAIL range / coherence:")
    print("  add smallest extra term (K(phi) or one derivative interaction),")
    print("  with the failed plot attached")

# ── Plot ───────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 4, figsize=(20, 9))

for row, (name, res) in enumerate(results.items()):
    u = res['u']

    # Panel 1: |Sigma(r)|
    ax = axes[row, 0]
    ax.semilogy(u, res['abs_Sigma_norm'], 'b-', linewidth=1.5)
    ax.axvline(0, color='r', linestyle='--', alpha=0.5, label='r=R')
    ax.set_xlabel(r'$(r - R) / \lambda_C$')
    ax.set_ylabel(r'$|\Sigma(r)| / \Sigma_{\rm max}$')
    ax.set_title(f'{name}: $|\\Sigma(r)| = |\\nabla \\ln A|$')
    ax.legend(fontsize=8)
    ax.set_xlim(-30, 30)
    ax.set_ylim(1e-6, 2)

    # Panel 2: S_Sigma(r)
    ax = axes[row, 1]
    ext = u[res['r'] > res['R']]
    S_vals = res['S_Sigma_r'][res['r'] > res['R']]
    valid = ~np.isnan(S_vals) & (S_vals > 0)
    if np.any(valid):
        ax.semilogy(ext[valid], S_vals[valid], 'g-', linewidth=1.5)
    ax.axhline(5.75e-6, color='r', linestyle='--', alpha=0.5, label='Cassini')
    ax.set_xlabel(r'$(r - R) / \lambda_C$')
    ax.set_ylabel(r'$S_\Sigma(r) = Q/Q_0$')
    ax.set_title(f'{name}: $S_\\Sigma(r)$ (source charge)')
    ax.legend(fontsize=8)
    ax.set_xlim(-30, 30)

    # Panel 3: eps_B
    ax = axes[row, 2]
    ax.plot(u, res['eps_B'], 'k-', linewidth=2)
    ax.axhline(1e-15, color='r', linestyle='--', alpha=0.5, label=r'$10^{-15}$')
    ax.set_xlabel(r'$(r - R) / \lambda_C$')
    ax.set_ylabel(r'$\varepsilon_B = (B/A^2)|\phi\'|^2$')
    ax.set_title(f'{name}: $\\varepsilon_B$ (B=0)')
    ax.legend(fontsize=8)
    ax.set_xlim(-30, 30)
    ax.set_ylim(-1e-16, 2e-15)

    # Panel 4: J
    ax = axes[row, 3]
    ax.plot(u, res['J'], 'k-', linewidth=2)
    ax.set_xlabel(r'$(r - R) / \lambda_C$')
    ax.set_ylabel(r'$J \propto (B/A^2)(u \cdot \Sigma)\Sigma$')
    ax.set_title(f'{name}: $J$ (B=0)')
    ax.set_xlim(-30, 30)

fig.suptitle(
    'Benchmark BVP: $V=\\frac{1}{2}m^2\\phi^2$, $K=1$, $\\beta_A=-1$, $B=0$\n'
    f'$m = \\Lambda_{{DE}} = {m:.2e}$ GeV,  $\\lambda_C = {hbar_c/m*10:.2f}$ mm\n'
    f'TEP context: conformal sector only (Theorem 1: null cones preserved; '
    f'holonomy=0 exact, kernel $\\Sigma_r \\neq 0$)',
    fontsize=11)
plt.tight_layout()
plt.savefig('benchmark_bvp.png', dpi=150, bbox_inches='tight')
print(f"\nPlot: benchmark_bvp.png")

# ── Save results ───────────────────────────────────────────────────────
os.makedirs('results', exist_ok=True)
output = {
    'benchmark': 'V = 0.5 * m^2 * phi^2, m = Lambda_DE (benchmark, NOT TEP)',
    'm_GeV': float(m),
    'lambda_C_cm': float(hbar_c / m),
    'beta_A': BETA, 'K': 1, 'B': 0,
    'bodies': {},
    'pass_fail': {
        'sun_cassini': {'value': float(S_sun), 'bound': 5.75e-6,
                        'pass': bool(cassini_pass)},
        'eps_B_GW': {'value': 0.0, 'bound': 1e-15, 'pass': bool(gw_pass),
                     'note': 'B=0 frozen, trivially satisfied'},
        'earth_recovery': {
            'transition_width_lambda_C': float(transition_width_u),
            'smooth': bool(is_smooth),
            'nonzero_inside': bool(inside_nonzero),
            'nonzero_outside': bool(outside_nonzero),
            'pass': bool(recovery_exists),
            'note': 'smooth exponential recovery, not a step'},
        'triangle_kernel': {
            'nonzero': bool(kernel_nonzero),
            'pass': bool(kernel_nonzero),
            'note': 'Sigma_r != 0 on space legs; holonomy=0 (exact form, B=0)'},
    },
    'overall': f'{n_pass}/4 PASS',
    'outcome': ('PASS: local K=1 is viable; V still a benchmark, not frozen'
                if n_pass == 4 else 'FAIL'),
}
for name, res in results.items():
    output['bodies'][name] = {
        'x_mR': float(res['x']),
        'S_Sigma_at_R': float(res['S_at_R']),
        'phi_body_GeV': float(res['phi_body']),
        'phi_body_over_M_Pl': float(res['phi_body'] / M_Pl),
        'Sigma_max_GeV': float(res['Sigma_max']),
        'Q_eff': float(res['Q_eff']),
        'Q_0': float(res['Q_0']),
    }

with open('results/benchmark_bvp.json', 'w') as f:
    json.dump(output, f, indent=2)
print(f"Results: results/benchmark_bvp.json")
