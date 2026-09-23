import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar

# Constants
G = 6.67430e-11
M_SUN = 1.98847e30
R_SUN = 6.957e8
AU = 149597870700.0

def get_s(g_N, phi, phi_0, s_0):
    # Solves s + (\phi/\phi_0)^2 * s_0^4 / s^3 = g_N
    # Let x = s. x^4 - g_N * x^3 + (\phi/\phi_0)^2 * s_0^4 = 0
    c = (phi / phi_0)**2 * s_0**4
    if c == 0:
        return g_N
        
    def f(x):
        return x**4 - g_N * x**3 + c
    def df(x):
        return 4*x**3 - 3*g_N*x**2
        
    # The minimum of f(x) is at x = 3/4 g_N. 
    # If f(3/4 g_N) > 0, there is no real root (screening is too strong).
    # But this shouldn't happen if the model is well posed, wait...
    # Actually, if c > (27/256) g_N^4, there's NO real root!
    # This means the operator D = 1 + A/s^4 has a maximum allowed A for a given g_N.
    # If A is too large, the equation has no solution!
    
    # This is a known issue with MOND-like theories with n > 1. 
    # We must use D(s) that doesn't cause a breakdown.
    # What if D = 1 + (\phi/\phi_0)^2 * (s_0/s)? Then m=1.
    # Then p = (1 - 2) / (2 - 2) = -1 / 0 ... wait.
    pass

def check_well_posed(n, m):
    # We need p = 1/3, where p = (m - n) / (2m - n)
    # 1/3 = (m - n) / (2m - n)  => 2m - n = 3m - 3n => m = 2n
    # If m = 2n, then m/n = 2.
    # For m=2, n=1: D = 1 + (\phi/\phi_0) (s_0/s)^2. 
    # Root of s^3 - g_N s^2 + c = 0. Min at x = 2/3 g_N. 
    # Again, max c limit.
    pass
