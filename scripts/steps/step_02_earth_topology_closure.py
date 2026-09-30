#!/usr/bin/env python3
"""Resolved layered-Earth profile: a mean-field diagnostic, not covariance."""
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import eval_legendre
from tep_model import (G,M_EARTH,R_EARTH,LAMBDA_REFERENCE,solve_sphere,save,
  equilibrium_varphi,mass_squared,compton_m)


def kinetic_multipole_diagnostic():
    g_t=3.4e-10
    r_star=np.sqrt(G*M_EARTH/g_t)
    y=brentq(lambda z:z+(r_star/R_EARTH)**4*z**3-1,0,1)
    modes={str(ell):{'decaying_power':float((1/3-np.sqrt(1/9+4*ell*(ell+1)/3))/2)}
           for ell in (1,2,3,4)}
    return {'branch':'deep kinetic shell, mass term and boundary matching omitted',
      'r_star_km':float(r_star/1000),'surface_gradient_ratio':float(y),
      'surface_transverse_stiffness':float(1/y),
      'surface_radial_to_transverse_stiffness':float(3-2*y),
      'mode_equation':'3p(p-1/3)=ell(ell+1) for exterior perturbations f_ell proportional to r^p',
      'multipoles':modes,
      'classification':'Conditional exterior radial-response asymptote, not a clock covariance or an exponential correlation length.'}


def angular_source_diagnostic():
    lengths={str(ell):float(R_EARTH*brentq(
        lambda angle:eval_legendre(ell,np.cos(angle))-1/np.e,0,np.pi/ell)/1000)
        for ell in (1,2,4,10)}
    return {'single_mode_first_1e_crossing_km':lengths,
      'geometry':'Single statistically isotropic angular multipole: C(theta)/C(0)=P_ell(cos theta)',
      'classification':'Conditional source spectra, not independent predictions. The response operator alone cannot fix the angular covariance length or its amplitude.'}


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
      'classification':'Mean-profile geometry. A profile maximum is not a two-point correlation length.',
      'kinetic_multipole_diagnostic':kinetic_multipole_diagnostic(),
      'angular_source_diagnostic':angular_source_diagnostic(), 'scan':rows}
    save('step_02_earth_topology_closure.json',result)
    print(result)
    return result

if __name__=='__main__':
    run()
