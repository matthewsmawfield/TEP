import numpy as np
from scipy.optimize import root_scalar

G = 6.67430e-11
M_SUN = 1.98847e30
R_SUN = 6.957e8
AU = 149597870700.0

def U_fun(phi, phi_0, n):
    return phi + (phi**(n + 1)) / ((n + 1) * (phi_0**n))

def solve_phi(U_target, phi_0, n):
    def f(phi):
        return U_fun(phi, phi_0, n) - U_target
    try:
        sol = root_scalar(f, bracket=[0, U_target], method='brentq')
        return sol.root
    except:
        return None

def test_model(phi_0, n):
    # Cassini
    g_N_c = G * M_SUN / (1.6 * R_SUN)
    U_c = g_N_c * (1.6 * R_SUN)
    phi_c = solve_phi(U_c, phi_0, n)
    D_c = 1 + (phi_c / phi_0)**n
    S_c = 1.0 / D_c
    pass_c = S_c < 5.75e-6
    
    # Saturn
    g_N_s = G * M_SUN / (9.5 * AU)
    U_s = g_N_s * (9.5 * AU)
    phi_s = solve_phi(U_s, phi_0, n)
    D_s = 1 + (phi_s / phi_0)**n
    a_s = g_N_s / D_s
    pass_s = a_s < 1e-10
    
    # WB
    g_N_wb = G * (1.2 * M_SUN) / (2646 * AU)
    U_wb = g_N_wb * (2646 * AU)
    phi_wb = solve_phi(U_wb, phi_0, n)
    D_wb = 1 + (phi_wb / phi_0)**n
    S_wb = 1.0 / D_wb
    pass_wb = S_wb > 0.1
    
    if pass_c and pass_s and pass_wb:
        print(f"FOUND ONE! n={n}, phi_0={phi_0}")
        print(f"Cassini S: {S_c:.2e}")
        print(f"Saturn accel: {a_s:.2e}")
        print(f"WB S: {S_wb:.2e}")
        return True
    return False

for n in [2, 4, 6, 8, 10, 12, 14, 16, 20]:
    for phi_0 in np.logspace(4, 7, 300):
        if test_model(phi_0, n):
            exit(0)
            
print("None found.")
