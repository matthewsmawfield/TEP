#!/usr/bin/env python3
"""Resolved layered-Earth profile: a mean-field diagnostic, not covariance."""
import numpy as np
from scipy.optimize import minimize_scalar
from tep_model import (M_EARTH,R_EARTH,LAMBDA_REFERENCE,solve_sphere,save,
  equilibrium_varphi,mass_squared,compton_m)


def prem_shape(x):
    core=.5*(1-np.tanh((x-3480/6371)/.01))
    surface=.5*(1-np.tanh((x-1)/.01))
    return (4.5+6.5*core)*surface


def run():
    rows=[]
    for factor in [.01,1,100,1e8]:
        solved=solve_sphere(M_EARTH,R_EARTH,LAMBDA_REFERENCE*factor,density_shape=prem_shape)
        x=np.linspace(3480/6371+.03,.98,2001)
        shear=-solved['profile'](x)[1]
        j=int(np.argmax(shear))
        opt=minimize_scalar(lambda z:solved['profile'](z)[1],bounds=(x[max(0,j-1)],x[min(len(x)-1,j+1)]),method='bounded')
        location=float(opt.x)
        rows.append({'lambda':LAMBDA_REFERENCE*factor,'mu':solved['mu'],
          'mantle_compton_km':compton_m(equilibrium_varphi(4.5,
            LAMBDA_REFERENCE*factor),4.5,LAMBDA_REFERENCE*factor)/1000,
          'mantle_interval_maximum_r_km':location*R_EARTH/1000,
          'maximum_at_interval_boundary':j in [0,len(x)-1],
          'normalized_gradient':float(-solved['profile'](location)[1]),
          'surface_varphi':float(solved['background_varphi']+solved['psi_uns']*solved['profile'](1)[0]),
          'max_rms_residual':solved['max_rms_residual']})
    result={'density_model':'Smoothed two-layer shape (11:4.5 core:mantle), normalized to Earth excess mass; not full PREM.',
      'equation':'Same exact quartic and homogeneous equilibrium subtraction as step 01',
      'classification':'Mean-profile geometry. A profile maximum is not a two-point correlation length.', 'scan':rows}
    save('step_02_earth_topology_closure.json',result)
    print(result)
    return result

if __name__=='__main__':
    run()
