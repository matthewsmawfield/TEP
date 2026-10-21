import json
import os
import math
import sympy as sp

# Physical constants
G = 6.67430e-11
c = 299792458.0
M_sun = 1.989e30
R_sun = 6.96e8
M_pl = math.sqrt(c * c / (8 * math.pi * G))  # Reduced Planck mass in kg (roughly)

# J1738+0333 benchmark parameters
# Antoniadis et al. 2012
M_ns = 1.46 * M_sun
M_wd = 0.181 * M_sun
P_b = 8.5 * 3600  # seconds
dot_P_b_obs = -2.59e-14
dot_P_b_GR = -2.77e-14
dot_P_b_err = 0.32e-14

def run():
    print("=========================================================")
    print(" TEP Binary Pulsar Scalar Dipole Derivation (J1738+0333)")
    print("=========================================================")
    print()

    # 1. Action and Perturbation Equations
    print(r"[1] Radiative Modes from the TEP Action")
    print(r"The frozen TEP scalar action is S_phi = \int d^4x \sqrt{-g} [ -1/2 (\nabla \phi)^2 - \lambda \phi^4 / 4 ] + S_m")
    print(r"Equation of motion: \square \phi - \lambda \phi^3 = - Q")
    print(r"where the bare source is Q = (1 / M_pl) T (for universal \beta_A = -1).")
    
    # 2. Linearized propagation
    # In the interstellar medium (ISM), phi -> phi_gal
    # Effective mass m_eff^2 = 3 \lambda \phi_gal^2
    # From Step 01 / 19, the ISM Compton wavelength is ~ 1.33 pc = 4e16 m
    lambda_compton = 1.33 * 3.086e16  # meters
    
    # Orbital wavelength of J1738
    lambda_gw = c * P_b
    
    print(f"\n[2] Scalar Radiation Propagation")
    print(f"ISM Compton wavelength: {lambda_compton:.2e} m")
    print(f"Orbital GW wavelength:  {lambda_gw:.2e} m")
    
    propagates = lambda_gw < lambda_compton
    if propagates:
        print("Result: lambda_gw << lambda_compton. The scalar field is effectively massless on orbital scales.")
        print("Propagating scalar radiation IS kinematically supported by the TEP action.")
    else:
        print("Result: lambda_gw > lambda_compton. Scalar radiation is evanescent and suppressed.")
        
    # 3. Effective Scalar Charges
    print("\n[3] Effective Scalar Charges (Nonlinear Solution)")
    print("Because both the NS and WD are highly compact, their interiors are deep within the Temporal Topology saturation regime (R << R_s).")
    print(r"In the \lambda \phi^4 completion, a saturated body exhibits a thin-shell effect.")
    print(r"The exterior field behaves as \phi(r) ~ \phi_ext + Q_eff / (4 \pi r).")
    print(r"The effective charge is \alpha_eff = Q_eff / M = \alpha_bare * (3 \Delta R / R).")
    print(r"For \lambda \phi^4, the thin-shell factor \Delta R / R \propto \phi_ext / \Phi_N.")
    print(r"Since the Newtonian potential \Phi_N is completely different for a NS (~0.2) and a WD (~10^-4),")
    print(r"their effective charges will be vastly different: \alpha_NS << \alpha_WD << 1.")
    
    # 4. Action-Closure Issue
    print("\n[4] The Action-Closure Bottleneck")
    print(r"To predict the exact \dot{P}_b, we need the exact numerical values of \alpha_NS and \alpha_WD.")
    print(r"This requires knowing the exact asymptotic field \phi_ext at the binary, which depends on the absolute normalization of \lambda and the galactic field.")
    print(r"As noted in the Jakarta manuscript, the absolute normalization of \lambda and the dimensional field \phi is an outstanding action-closure issue.")
    print(r"Without inserting phenomenological screening factors (which we must not do), the fixed action currently dictates:")
    print(r"  a) \alpha_NS != \alpha_WD (because \Phi_N differs by 10^3)")
    print(r"  b) Therefore, scalar dipole radiation will be generated.")
    print(r"  c) But the amplitude cannot be evaluated without closing the absolute dimensional scale of \phi.")
    
    # We report this as a formal theoretical output.
    output = {
        "benchmark": "PSR J1738+0333",
        "orbital_wavelength_m": lambda_gw,
        "compton_wavelength_m": lambda_compton,
        "propagating_modes_supported": propagates,
        "charge_evaluation": {
            "ns_potential_approx": 0.2,
            "wd_potential_approx": 1e-4,
            "charges_equal": False,
            "dipole_radiation_present": True
        },
        "closure_status": "PENDING_ACTION_NORMALIZATION",
        "conclusion": "The TEP action supports propagating scalar waves for this binary. Because the NS and WD have vastly different binding potentials, the nonlinear field profiles yield different effective scalar charges (thin-shell suppression). This sources scalar dipole radiation. However, calculating the precise P_b_dot requires the absolute dimensional normalization of the field \phi to evaluate the thin-shell factor. This is a recognized action-closure issue."
    }
    
    outdir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "step_23_binary_pulsar_dipole.json"), "w") as f:
        json.dump(output, f, indent=2)
    print("\n[SUCCESS] Wrote JSON output to scripts/results/step_23_binary_pulsar_dipole.json")

if __name__ == "__main__":
    run()
