#!/usr/bin/env python3
import numpy as np

M_Pl = 1.0
lam = 7.5e-66
u_s = 10.0
V_0 = 0.3  # M_Pl^4

def V(u):
    term1 = (lam / 4) * u**4 * np.exp(-(u / u_s)**4)
    term2 = V_0 * np.exp(-(u_s / u)**4)
    return term1 + term2

def dV(u):
    # d/du
    # term1: (lam/4) * [ 4 u^3 e - u^4 * 4(u^3/u_s^4) e ]
    # term2: V_0 * e * [ 4 u_s^4 / u^5 ]
    term1 = (lam / 4) * (4 * u**3 - 4 * u**7 / u_s**4) * np.exp(-(u / u_s)**4)
    term2 = V_0 * (4 * u_s**4 / u**5) * np.exp(-(u_s / u)**4)
    return term1 + term2

print("z", "u=ln(1+z)", "V(u)", "dV(u)")
for z in [0, 1, 10, 100, 1000]:
    u = np.log(1 + z)
    if u == 0:
        val_V, val_dV = 0, 0
    else:
        val_V = V(u)
        val_dV = dV(u)
    print(f"{z:<5} {u:<10.5f} {val_V:<15.5e} {val_dV:<15.5e}")

