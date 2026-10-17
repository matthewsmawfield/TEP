#!/usr/bin/env python3
"""Conservative discrete positivity bound for the step79 ambient branch.

For weighted-zero-mean v, write v=v0+mean(v), mean(v0)=0.
Then ||v||² <= (1+Var(w)/mean(w)²)||v0||², w=psi^6.
Use the periodic lattice Poincare constant, minimum face stiffness,
minimum potential, and a norm bound on the attractive Schur correction.
The nonnegative directional stiffness is discarded, making the bound weaker.
This certifies only the specified discretized operators/backgrounds, not the
continuum limit or every action constraint. A nonpositive bound is inconclusive.
"""
import hashlib
import json
from pathlib import Path
import numpy as np
from step_20_landscape_existence import (
    build_landscape, solve_landscape, fd_laplacian, landscape_potential, M_PL)


def calculate(N):
    dx = 1.0/N
    data = build_landscape(N, 1., 4, .05, 1., .8, 2., .05, 2., 1., 15.)
    sol = solve_landscape(data['S0'], data['g'], fd_laplacian(N, dx), verbose=False)
    if not sol['converged']:
        raise RuntimeError('Background did not converge')
    psi, u, rho = sol['psi'], data['u'], data['rho_m']
    eps = 1e-6
    V = lambda x: landscape_potential(x, 2., 1., 15.)
    vp = (V(u+eps)-V(u-eps))/(2*eps)
    vpp = (V(u+eps)-2*V(u)+V(u-eps))/eps**2
    w = psi**6
    poincare = 4*np.sin(np.pi/N)**2/dx**2
    factor = 1+np.var(w)/np.mean(w)**2
    c = -5*psi**4/(4*M_PL**2)*(data['S0']+sol['lam'])-data['g']/(8*M_PL**2)
    lpsi_bound = poincare+float(np.min(c))
    coupling = psi**5/(4*M_PL**2)*(vp-rho*np.exp(-u))
    correction = float(np.max(coupling**2/w)/lpsi_bound) if lpsi_bound > 0 else None
    rhob_half = 1.0
    rho_cell_mean = 2.0 * rhob_half
    branches = {}
    for name, ambient in [('ambient', .32), ('pure_local', 0.)]:
        xi = np.maximum(.5*data['g']+ambient, 1e-30)
        stiffness = psi**2*(16.03*np.sqrt(xi)+2*xi)
        stiff_3d = stiffness.reshape(N, N, N)
        face_stiffness_min = float(min(
            0.5 * np.min(stiff_3d + np.roll(stiff_3d, -1, axis=d))
            for d in range(3)
        ))
        gradient = float(face_stiffness_min*poincare/(np.max(w)*factor))
        gradient_vertex = float(np.min(stiffness)*poincare/(np.max(w)*factor))
        potential = float(np.min(vpp+rho_cell_mean*np.exp(-u)))
        potential_unnorm = float(np.min(vpp+rho*np.exp(-u)))
        bound = gradient+potential-correction if correction is not None else None
        bound_vertex_unnorm = gradient_vertex+potential_unnorm-correction if correction is not None else None
        branches[name] = {'gradient_lower_bound':gradient,
                          'gradient_lower_bound_vertex_min':gradient_vertex,
                          'face_stiffness_min':face_stiffness_min,
                          'potential_lower_bound':potential,
                          'potential_lower_bound_unnormalized_rho':potential_unnorm,
                          'attractive_correction_upper_bound':correction,
                          'operator_lower_bound':bound,
                          'operator_lower_bound_vertex_unnormalized':bound_vertex_unnorm,
                          'status':'positive_discrete_bound' if bound is not None and bound > 0 else 'inconclusive_bound'}
    return {'N':N,'constraint_converged':True,'poincare_constant':float(poincare),
            'lpsi_zero_mean_lower_bound':lpsi_bound,'branches':branches}


if __name__ == '__main__':
    root = Path(__file__).resolve().parents[2]
    result = {'scope':__doc__, 'grids':[], 'provenance':{}}
    for N in [24,32,40]:
        row = calculate(N)
        result['grids'].append(row)
        print(json.dumps(row), flush=True)
    for p in [Path(__file__), Path(__file__).with_name('step_20_landscape_existence.py'),
              Path(__file__).with_name('step_79_inhomogeneous_eigenmodes.py')]:
        result['provenance'][p.name] = hashlib.sha256(p.read_bytes()).hexdigest()
    (root/'results/step_82_stability_lower_bound.json').write_text(json.dumps(result,indent=2)+'\n')
