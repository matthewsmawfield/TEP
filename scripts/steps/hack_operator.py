import numpy as np

# Constants
G = 6.67430e-11
M_SUN = 1.98847e30
R_SUN = 6.957e8
AU = 149597870700.0

def solve_phase_transition(M, radii, phi_c, D_max):
    # U(\phi) = \phi for \phi < \phi_c
    # U(\phi) = \phi_c + D_max * (\phi - \phi_c) for \phi > \phi_c
    results = []
    for r in radii:
        g_N = G * M / r**2
        U = G * M / r
        if U <= phi_c:
            phi = U
            D = 1.0
        else:
            phi = phi_c + (U - phi_c) / D_max
            D = D_max
        
        s = g_N / D
        S = 1.0 / D
        results.append((r/AU, g_N, s, S))
    return results

print("Phase Transition Model")
res = solve_phase_transition(M_SUN, [1.6*R_SUN, 9.5*AU, 2646*AU], phi_c=1e6, D_max=1e6)
for r, gn, s, S in res:
    print(f"r={r:.2e} AU, g_N={gn:.2e}, s={s:.2e}, S={S:.2e}")

res_wb = solve_phase_transition(1.2*M_SUN, [2646*AU], phi_c=1e6, D_max=1e6)
print("Wide Binary 1.2 M_sun at 2646 AU:")
print(f"S = {res_wb[0][3]:.2e}")
