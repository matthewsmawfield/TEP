import numpy as np

# Constants
G = 6.67430e-11
M_SUN = 1.98847e30
R_SUN = 6.957e8
AU = 149597870700.0

def evaluate_S(r, M, Phi_0, s_0):
    g_bare = G * M / r**2
    Phi_N = G * M / r
    A_Phi = np.sqrt(Phi_N / Phi_0)
    
    g = g_bare / s_0
    a = 1.0
    b = 1.0 + A_Phi - g
    c = -g
    
    discriminant = b**2 - 4*a*c
    if discriminant < 0:
        return None
        
    if b > 0:
        x = (2*c) / (-b - np.sqrt(discriminant))
    else:
        x = (-b + np.sqrt(discriminant)) / (2*a)
    s = x * s_0
    S = s / g_bare
    
    if r > 1000 * AU:
        print(f"DEBUG WB: g_bare={g_bare:.2e}, Phi_N={Phi_N:.2e}, A_Phi={A_Phi:.2f}, g={g:.2f}, b={b:.2f}, c={c:.2f}, x={x:.2f}, s={s:.2e}, S={S}")
        
    return S

def scan():
    phi_0_vals = np.logspace(-15, -5, 100)
    s_0_vals = np.logspace(-12, -2, 100)
    
    found = False
    for Phi_0 in phi_0_vals:
        for s_0 in s_0_vals:
            S_cas = evaluate_S(1.6 * R_SUN, M_SUN, Phi_0, s_0)
            S_sat = evaluate_S(9.5 * AU, M_SUN, Phi_0, s_0)
            S_wb  = evaluate_S(2646 * AU, 1.2 * M_SUN, Phi_0, s_0)
            
            if S_cas is not None and S_sat is not None and S_wb is not None:
                if S_cas < 5.75e-6 and S_sat < 1.5e-6 and S_wb > 0.1:
                    print(f"FOUND: Phi_0={Phi_0:.2e}, s_0={s_0:.2e}")
                    print(f"  Cassini S = {S_cas:.2e}")
                    print(f"  Saturn S  = {S_sat:.2e}")
                    print(f"  WB S      = {S_wb:.2e}")
                    found = True
                    return
                    
    if not found:
        print("No parameters found that satisfy Cassini, Saturn, and WB simultaneously.")
        
        # Print one example of where it fails but S_wb > 0.1
        Phi_0 = 1e-12
        s_0 = 1e-9
        print(f"\\nDiagnostic for Phi_0={Phi_0:.2e}, s_0={s_0:.2e}")
        S_cas = evaluate_S(1.6 * R_SUN, M_SUN, Phi_0, s_0)
        S_sat = evaluate_S(9.5 * AU, M_SUN, Phi_0, s_0)
        S_wb  = evaluate_S(2646 * AU, 1.2 * M_SUN, Phi_0, s_0)
        print(f"  Cassini S = {S_cas:.2e} (Needs < 5.75e-6)")
        print(f"  Saturn S  = {S_sat:.2e} (Needs < 1.5e-6)")
        print(f"  WB S      = {S_wb}")

scan()
