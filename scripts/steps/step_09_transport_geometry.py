#!/usr/bin/env python3
"""Symbolic ADM, clock-transport and cosmological consistency identities."""
import sympy as s
from tep_model import C, save

def run():
    A, B = s.symbols('A B', real=True, positive=True)
    N = s.symbols('N', real=True, positive=True)
    dt_phi = s.symbols('dt_phi', real=True)
    di_phi = s.symbols('di_phi', real=True)
    Ni = s.symbols('Ni', real=True)
    N_up_i = s.symbols('N_up_i', real=True)
    
    h = s.symbols('h', real=True, positive=True)
    
    # 1. Exact ADM expressions
    g00 = -N**2 + h * N_up_i**2
    g0i = h * N_up_i
    gij = h
    
    gtilde_00 = A**2 * g00 + B * dt_phi**2
    gtilde_0i = A**2 * g0i + B * dt_phi * di_phi
    gtilde_ij = A**2 * gij + B * di_phi**2
    
    Ntilde_i = gtilde_0i
    htilde_ij = gtilde_ij
    htilde_up_ij = 1 / htilde_ij
    Ntilde_up_i = htilde_up_ij * Ntilde_i
    
    Ntilde_sq = -gtilde_00 + Ntilde_i * Ntilde_up_i
    sigma_i = gtilde_0i / (-gtilde_00) # standard convention
    
    # 2. Four automatic checks
    
    # Check 1: Pure conformal limit (B -> 0)
    ch1_Ntilde_sq = s.simplify(Ntilde_sq.subs(B, 0) - A**2 * N**2)
    ch1_Ntilde_i = s.simplify(Ntilde_i.subs(B, 0) - A**2 * g0i)
    ch1_sigma_i = s.simplify(sigma_i.subs(B, 0) - g0i/(-g00))
    pure_conformal = bool(ch1_Ntilde_sq == 0 and ch1_Ntilde_i == 0 and ch1_sigma_i == 0)
    
    # Check 2: Vanishing-gradient limit (dt_phi -> 0, di_phi -> 0)
    # The gradient is the full spacetime gradient
    ch2_sigma_i = s.simplify(sigma_i.subs({dt_phi: 0, di_phi: 0}) - g0i/(-g00))
    vanishing_gradient = bool(ch2_sigma_i == 0)
    
    # Check 3: Zero-shift limit (N_up_i -> 0)
    sigma_i_zs = sigma_i.subs(N_up_i, 0)
    # in zero shift, g00 = -N^2, g0i = 0.
    # sigma_i = (B dt_phi di_phi) / (A^2 N^2 - B dt_phi^2)
    target_sigma_i_zs = (B * dt_phi * di_phi) / (A**2 * N**2 - B * dt_phi**2)
    zero_shift = bool(s.simplify(sigma_i_zs - target_sigma_i_zs) == 0)
    
    # Check 4: Consistency with inverse metric from Appendix A4
    # In A4: gtilde^0i = A^-2 * (g^0i - C * q^0 q^i), where C = (B/A^2) / (1 + (B/A^2) q.q)
    # We test the relation: gtilde^0i / gtilde^00 = - Ntilde^i = - htilde^{ij} Ntilde_j
    # We can just compute gtilde^0i and gtilde^00 directly from A4 formula and check the ratio
    
    # g^mu_nu components:
    # g^00 = -1/N^2
    # g^0i = N_up_i / N^2
    # g^ij = 1/h - (N_up_i)^2 / N^2
    
    g_up_00 = -1/N**2
    g_up_0i = N_up_i/N**2
    g_up_ij = 1/h - N_up_i**2/N**2
    
    q0 = dt_phi
    qi = di_phi
    # q^mu = g^mu_nu q_nu
    q_up_0 = g_up_00 * q0 + g_up_0i * qi
    q_up_i = g_up_0i * q0 + g_up_ij * qi
    
    # q.q = q^mu q_mu
    q_dot_q = q_up_0 * q0 + q_up_i * qi
    
    C_factor = (B/A**2) / (1 + (B/A**2) * q_dot_q)
    
    gtilde_up_00 = (g_up_00 - C_factor * q_up_0 * q_up_0) / A**2
    gtilde_up_0i = (g_up_0i - C_factor * q_up_0 * q_up_i) / A**2
    
    # gtilde^0i / gtilde^00 should equal - Ntilde_up_i
    ratio = gtilde_up_0i / gtilde_up_00
    diff_ratio = s.simplify(ratio - (-Ntilde_up_i))
    consistency_A4 = bool(diff_ratio == 0)
    
    result = {
        'ADM': {
            'pure_conformal': pure_conformal,
            'vanishing_gradient': vanishing_gradient,
            'zero_shift': zero_shift,
            'consistency_A4': consistency_A4
        }
    }
    
    print("Checks passed:")
    print("Pure conformal:", pure_conformal)
    print("Vanishing gradient:", vanishing_gradient)
    print("Zero shift:", zero_shift)
    print("Consistency with A4:", consistency_A4)
    
    save('step_09_transport_geometry.json', result)
    return result

if __name__ == '__main__':
    run()
