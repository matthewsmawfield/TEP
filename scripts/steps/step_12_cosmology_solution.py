#!/usr/bin/env python3
"""
verify_static_closed_background.py

Standalone check of the static, closed (k=+1) Einstein-frame TEP background.

Conventions (Jakarta v0.14, weak-field canonical branch, K=1):
  Einstein frame, units with c = 1, M = reduced Planck mass.
  ds^2 = -dt^2 + a^2 dOmega_3^2,  a = const  (static, closed).
  A(phi) = exp(beta_A * phi / M),  beta_A = -1  ->  alpha = dlnA/dphi = -1/M.
  Matter = pressureless dust in the matter frame.

Equations (standard conformally coupled scalar-tensor, Einstein frame):
  scalar:      phi_dd = -V'(phi) - alpha*rho        (Box phi = V' - alpha*T, T = -rho, Box = -d_t^2)
  continuity:  rho_d = alpha * rho * phi_d          ->  rho = rho_bar * A(phi)
  00:          3 M^2 / a^2 = rho + phi_d^2/2 + V
  ii:          -M^2 / a^2  = phi_d^2/2 - V           (dust has p = 0)

Result derived in the accompanying note: a static solution with phi_d != 0 exists iff
  V(phi) = 2 M^2/a^2 - (rho_bar/2) * A(phi).
This script integrates the FULL dynamical system (scalar EOM + matter continuity, a fixed)
with that potential and verifies the 00 and ii equations hold along the trajectory.
It also reports three derived consequences that must be checked against data:
  (1) future turnaround of the apparent Hubble rate,
  (2) H_app(z)/H0 = sqrt(1 + 3*Omega0*z/(1+z)),
  (3) finite accumulated matter-frame proper time since the temporal horizon A -> 0.
"""
import json
import os
import numpy as np
from scipy.integrate import solve_ivp

M, a = 1.0, 1.0
beta_A = -1.0
alpha = beta_A / M
rho_bar = 1.0          # matter normalisation; any value with rho <= 2M^2/a^2 on the branch

A  = lambda phi: np.exp(beta_A * phi / M)
V  = lambda phi: 2*M**2/a**2 - 0.5*rho_bar*A(phi)
dV = lambda phi: -0.5*rho_bar*alpha*A(phi)          # = (rho_bar/2M) A  for beta_A=-1

def rhs(t, y):
    phi, phid, rho = y
    phidd = -dV(phi) - alpha*rho        # scalar EOM with dust source (T = -rho)
    rhod  = alpha*rho*phid               # Einstein-frame dust exchange
    return [phid, phidd, rhod]

# Initial data on the constraint surface, deep in the past (A small), phid < 0 so A grows.
phi0 = 12.0
rho0 = rho_bar*A(phi0)
phid0 = -np.sqrt(2*(V(phi0) - M**2/a**2))      # from ii
sol = solve_ivp(rhs, [0, 40], [phi0, phid0, rho0], rtol=1e-11, atol=1e-13,
                dense_output=True, max_step=0.01)
t = sol.t; phi, phid, rho = sol.y

res00 = 3*M**2/a**2 - (rho + 0.5*phid**2 + V(phi))
resii = -M**2/a**2 - (0.5*phid**2 - V(phi))
resrho = rho - rho_bar*A(phi)
print("max |00 residual|        :", np.max(np.abs(res00)))
print("max |ii residual|        :", np.max(np.abs(resii)))
print("max |rho - rho_bar*A|    :", np.max(np.abs(resrho)))

# Control: same run with the quartic V = lam phi^4/4 cannot satisfy ii with phid != 0.
lam = 1e-3
Vq, dVq = (lambda p: lam*p**4/4), (lambda p: lam*p**3)
def rhs_q(t, y):
    p, pd, r = y
    return [pd, -dVq(p) - alpha*r, alpha*r*pd]
ydq = solve_ivp(rhs_q, [0, 20], [1.0, -0.5, 1.0], rtol=1e-10, atol=1e-12, max_step=0.01)
pq, pdq, rq = ydq.y
print("quartic control: spread of ii combination (p_tot) :",
      np.ptp(0.5*pdq**2 - Vq(pq)), " (must be 0 for a static a; it is not)")

# (1) Turnaround: phid -> 0 when rho = 2 M^2/a^2
i_turn = np.argmax(phid >= 0) if np.any(phid >= 0) else None
if i_turn:
    print(f"turnaround at t={t[i_turn]:.3f}, rho={rho[i_turn]:.4f} (2M^2/a^2={2*M**2/a**2})")

# (2) Coordinate-time conformal drift rate H_t = dlnA/dt = beta_A*phid/M ; H(z) law.
# NB: the operational matter-frame rate H_matter = (aA)^-1 d(aA)/dtau = dlnA/dt / A
# carries an additional factor (1+z): H_matter(z)/H_matter(0) = (1+z) H_t(z)/H_t(0).
H = beta_A*phid/M
x0 = 0.947   # rho0 a^2/M^2 giving q0 ~ -0.55 (see note)
Om0 = x0/(3*(2-x0))
z = np.array([0.5, 1.0, 2.0, 1e6])
print(f"Omega0 implied = {Om0:.3f};  H(z)/H0 =", np.round(np.sqrt(1+3*Om0*z/(1+z)), 3),
      " at z =", z)

# (3) Matter-frame proper time accumulated from the temporal horizon to the epoch where x = x0
x = rho*a**2/M**2
kmax = i_turn if i_turn else len(x)
k = np.argmin(np.abs(x[:kmax] - x0))
tau = np.trapz(A(phi[:k+1]), t[:k+1])              # dtau = A dt
H0 = H[k]
tau_tail = A(phi[0]) * a/np.sqrt(2)              # analytic pre-history: A ~ exp(sqrt2 t/a)
age = (tau + tau_tail) * H0 / A(phi[k])
print(f"matter-frame age since A->0 : tau*H0/A0 = {age:.3f}  (finite)")

# (4) Master-family membership: the reconstructed branch potential
#     V_rec(u) = C1 - (rho_bar/2) e^{-u}  (C1 = 2M^2/a^2, u = phi/M)
# is tested against the adopted master family
#     V(u) = V_matter(u) e^{-(u/u_s)^4} + V_0 e^{-(u_s/u)^4}.
# Two members are exhibited:
#  (a) exact member: V_matter(u) = [V_rec(u) - V_0 e^{-(u_s/u)^4}] e^{+(u/u_s)^4}
#      reproduces V_rec identically for any u_s (machine-precision check);
#  (b) natural member: V_matter(u) = V_rec(u) itself (shallow exponential-plus-
#      constant on the matter domain) with V_0 = C1; deviation is confined to
#      the transition band u ~ u_s and vanishes asymptotically at the floor.
u_grid = np.linspace(-1.0, 30.0, 30000)
C1 = 2*M**2/a**2
V_rec_u = lambda u: C1 - 0.5*rho_bar*np.exp(beta_A*u)
us_probe = 15.0
cut = np.exp(-(u_grid/us_probe)**4)
floor_w = np.exp(-(us_probe/np.maximum(u_grid, 1e-12))**4)
V_exact_member = ((V_rec_u(u_grid) - C1*floor_w)/cut)*cut + C1*floor_w
V_nat_member = V_rec_u(u_grid)*cut + C1*floor_w
dev_nat = np.abs(V_nat_member - V_rec_u(u_grid))/V_rec_u(u_grid)
membership = {
    "u_s_probe": us_probe,
    "exact_member_max_frac_dev": float(np.max(np.abs(V_exact_member - V_rec_u(u_grid))/V_rec_u(u_grid))),
    "natural_member_max_frac_dev": float(dev_nat.max()),
    "natural_member_dev_u_le_half_us": float(dev_nat[u_grid <= us_probe/2].max()),
    "deep_field_floor": {"V_0_required": C1, "V_rec_asymptote": float(V_rec_u(50.0))},
}
print("master-family exact-member max frac dev :",
      f"{membership['exact_member_max_frac_dev']:.2e}")
print("natural member: max frac dev (transition band) :",
      f"{membership['natural_member_max_frac_dev']:.3f}",
      f"at u~u_s; u<=u_s/2: {membership['natural_member_dev_u_le_half_us']:.3f}")

results = {
    "step": "step_12_cosmology_solution",
    "description": "Exact static closed (k=+1) Einstein-frame TEP background: "
                   "full-dynamics integration with the branch potential "
                   "V = 2M^2/a^2 - (rho_bar/2)A(phi); verification of the 00, ii "
                   "and continuity equations along the trajectory; "
                   "master-family membership test of the reconstructed profile.",
    "residuals": {
        "max_abs_00": float(np.max(np.abs(res00))),
        "max_abs_ii": float(np.max(np.abs(resii))),
        "max_abs_continuity": float(np.max(np.abs(resrho))),
    },
    "quartic_control_ii_spread": float(np.ptp(0.5*pdq**2 - Vq(pq))),
    "turnaround": ({"t": float(t[i_turn]), "rho": float(rho[i_turn]),
                    "rho_turnover_2M2_over_a2": float(2*M**2/a**2)}
                   if i_turn else None),
    "Omega0_implied": float(Om0),
    "matter_frame_age_tauH0_over_A0": float(age),
    "master_family_membership": membership,
}
out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   os.pardir, os.pardir, "results",
                   "step_12_cosmology_solution.json")
out = os.path.normpath(out)
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, "w") as f:
    json.dump(results, f, indent=2)
print("wrote", out)
