#!/usr/bin/env python3
"""Converged disformal loop diagnostic, not a fixed-action amplitude prediction.

The connection is f(x) S(r_E) grad u, where
f = -B(u) M_Pl^2 u_dot / (A^2 N^2 c), u_dot=-H0 W(u).
W(u) alone cannot generate a curl. In this reduced ansatz the independent
Earth-radius screening S(r_E), in the presence of the Sun, supplies the
non-exactness. Its derivation from a single solved scalar configuration remains
open. B M_Pl^2 = B0 * envelope(u) * (c/H0)^2 is a dimensionless benchmark
normalization, NOT a direct import of the dimensionful strong-field B0.

Repairs: physical solar-radius cutoff instead of one AU; analytic gradient;
negative Newtonian potential in the lapse; visible satellite at the angular
bisector; parametric Gauss-Legendre quadrature and exact-connection controls.
The previous ~6e-13 s headline was an endpoint-quadrature artifact on a loop
whose clipped solar potential was constant. It is not a retained prediction.
"""
from functools import lru_cache
from pathlib import Path
import json
import numpy as np
from scipy.special import roots_legendre

C_SI=2.998e8
G_SI=6.674e-11
H0=2.27e-18
MPC=3.086e22
R_H=C_SI/H0
M_E=5.972e24
R_E=6.371e6
M_S=1.989e30
R_S=6.96e8
AU=1.496e11
KAPPA=1.0  # Preserved Newtonian-tracking diagnostic normalization.
SIGMA_B=1.0

def B_of_u(u,B0):
    return B0*u**2/(1+u**2)*np.exp(-u**4/(2*SIGMA_B**4))

def u_field(x):
    x=np.asarray(x,dtype=float)
    rE=np.linalg.norm(x,axis=-1)
    rS=np.linalg.norm(x-np.array([AU,0.,0.]),axis=-1)
    return KAPPA*G_SI/C_SI**2*(M_E/np.maximum(rE,R_E)+M_S/np.maximum(rS,R_S))

def grad_u_field(x,h=None):
    """Analytic gradient of the declared clipped benchmark; h is legacy API."""
    x=np.asarray(x,dtype=float);g=np.zeros_like(x)
    for center,mass,radius in [(np.zeros(3),M_E,R_E),(np.array([AU,0.,0.]),M_S,R_S)]:
        d=x-center;r=np.linalg.norm(d,axis=-1)
        factor=np.where(r>=radius,-KAPPA*G_SI*mass/(C_SI**2*np.maximum(r,radius)**3),0.)
        g+=factor[...,None]*d
    return g

def roll_weight(x,u_thresh):
    if u_thresh<=0:raise ValueError('u_thresh must be positive')
    return .5*(1+np.tanh((u_field(x)-u_thresh)/(.2*u_thresh)))

def connection_prefactor(x,B0,u_thresh,Pi_bar):
    u=u_field(x);N=1-u/KAPPA
    if np.any(N<=0):raise ValueError('Benchmark requires a positive weak-field lapse')
    return -B_of_u(u,B0)*R_H**2*(-Pi_bar*roll_weight(x,u_thresh))/(np.exp(-2*u)*N**2*C_SI)

def screening(x,s0=1e-6,ell=3.):
    r=np.linalg.norm(x,axis=-1)
    return s0+(1-s0)*(-np.expm1(-np.maximum(r-R_E,0.)/(ell*R_E)))

@lru_cache(maxsize=16)
def quadrature(npts):
    if npts<8:raise ValueError('At least eight quadrature nodes required')
    x,w=roots_legendre(npts)
    return (x+1)/2,w/2

def triangle(alt_m=2e6,lonB_deg=60.):
    angle=np.radians(lonB_deg)
    if alt_m<=0 or not 0<angle<np.pi:raise ValueError('Invalid triangle geometry')
    if angle/2>=np.arccos(R_E/(R_E+alt_m)):
        raise ValueError('Satellite link intersects Earth; increase altitude or reduce separation')
    A=R_E*np.array([1.,0.,0.]);B=R_E*np.array([np.cos(angle),np.sin(angle),0.])
    S=(R_E+alt_m)*np.array([np.cos(angle/2),np.sin(angle/2),0.])
    return A,S,B,angle

def loop_holonomy(B0,u_thresh,Pi_bar,alt_m=2e6,lonB_deg=60.,npts=128,s0=1e-6,reverse=False):
    """Integrate along two visible straight links and a surface transport arc.

    The ground arc is a schematic surface route, not a free-space chord.
    Reverse actually reverses paths and tangents, providing an orientation test.
    """
    A,S,B,angle=triangle(alt_m,lonB_deg);t,w=quadrature(npts)
    if reverse:t=1-t
    paths=[]
    for name,p,q in [('A->S',A,S),('S->B',S,B)]:
        x=p[None,:]+t[:,None]*(q-p)
        tangent=np.broadcast_to(q-p,x.shape).copy()
        paths.append((name,x,tangent))
    theta=angle*(1-t)
    x=R_E*np.stack([np.cos(theta),np.sin(theta),np.zeros_like(theta)],axis=-1)
    tangent=R_E*angle*np.stack([np.sin(theta),-np.cos(theta),np.zeros_like(theta)],axis=-1)
    paths.append(('B->A',x,tangent))
    legs={}
    for name,x,tangent in paths:
        if reverse:tangent=-tangent
        integrand=connection_prefactor(x,B0,u_thresh,Pi_bar)*screening(x,s0)*np.sum(grad_u_field(x)*tangent,axis=-1)
        legs[name]=float(w@integrand)
    return sum(legs.values())/C_SI,legs

def main():
    threshold=.5*float(u_field(np.array([[R_E,0.,0.]]))[0])
    configs={
        'admissible_sign_benchmark':dict(B0=1.,u_thresh=threshold,Pi_bar=H0),
        'negative_sign_diagnostic':dict(B0=-3.2e-3,u_thresh=threshold,Pi_bar=H0),
        'exact_connection_control':dict(B0=1.,u_thresh=threshold,Pi_bar=H0,s0=1.),
        'Pi=0_control':dict(B0=1.,u_thresh=threshold,Pi_bar=0.),
        'B0=0_control':dict(B0=0.,u_thresh=threshold,Pi_bar=H0),
        'reversed':dict(B0=1.,u_thresh=threshold,Pi_bar=H0,reverse=True),
    }
    scan={}
    for name,kw in configs.items():
        values=[loop_holonomy(**kw,npts=n)[0] for n in [32,64,128,256]]
        scan[name]={'H_resid_s':values[-1],'convergence_nodes':[32,64,128,256],'convergence_H_s':values,'last_refinement_difference_s':abs(values[-1]-values[-2])}
    result={'units':'seconds','status':'conditional reduced-connection diagnostic; not a fixed-action prediction',
        'geometry':{'satellite_altitude_m':2e6,'ground_separation_degrees':60.,'satellite_longitude_degrees':30.,'straight_links_clear_Earth':True,'ground_return':'surface arc'},
        'normalization':{'B_Mpl_squared':'B0 * envelope(u) * (c/H0)^2','KAPPA':KAPPA,'s0':1e-6,'screening_length_Earth_radii':3.,'drift_rate_s_inverse':H0},
        'scan':scan,'superseded_claim':'~6e-13 s from the old clipped-Sun endpoint sum; not a converged prediction',
        'mechanism':'S(r_E) is not a function of the composite Earth+Sun u; W(u) alone remains exact. Setting S=1 restores exactness in this benchmark.',
        'remaining_dependencies':['Derive u, its gradient, screening and drift from one action-level solution.','Fix the dimensionful B normalization across sectors.','Compute the experimental synchronization and reference-clock model on the actual link geometry.']}
    out=Path(__file__).resolve().parents[2]/'results/step_15_holonomy_amplitude.json'
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
