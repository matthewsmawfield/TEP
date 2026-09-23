#!/usr/bin/env python3
"""Timestamp-level disformal loop example and proper-time calibration.

Local test-field example: g=-N(x)^2 dt^2+dx^2+dy^2, N=1+a*x,
phi/M_Pl=q*t+k*y, and B*M_Pl^2/A^2=D constant. The massless scalar
satisfies the background wave equation. Backreaction and sources of optical
path guidance are outside this local transport example. It is not a
terrestrial amplitude forecast. All coordinates here are dimensionless.
"""
import numpy as np
from scipy.integrate import quad
from tep_model import save


def circuit(vertices,D,q,k,acceleration):
    total=0.0
    for start,end in zip(vertices,vertices[1:]):
        start=np.asarray(start,float);v=np.asarray(end,float)-start
        def integrand(s):
            x,y=start+s*v
            lapse=1+acceleration*x
            aa=-lapse*lapse+D*q*q
            bb=2*D*q*k*v[1]
            cc=float(v@v)+D*(k*v[1])**2
            if aa>=0:
                raise ValueError('Observer threading is not timelike')
            return (-bb-np.sqrt(bb*bb-4*aa*cc))/(2*aa)
        total+=quad(integrand,0,1,epsabs=1e-12,epsrel=1e-12)[0]
    return total


def loop_example(D=1000.0,q=.001,k=.001,acceleration=.1,length=1.0):
    loop=[(0,0),(length,0),(length,length),(0,length),(0,0)]
    cw=circuit(loop,D,q,k,acceleration)
    ccw=circuit(loop[::-1],D,q,k,acceleration)
    sigma=lambda x:-D*q*k/((1+acceleration*x)**2-D*q*q)
    holonomy=length*(sigma(length)-sigma(0))
    # beta=-1, phi_receiver=q*t: exact receiver-clock primitive.
    proper=lambda t:-np.expm1(-q*t)/q if q else t
    dtau=proper(cw)-proper(ccw)
    # Receiver normalization belongs outside the stationary coordinate loop;
    # here the exact time-dependent receiver lapse is integrated instead.
    return {'input':{'D':D,'q':q,'k':k,'acceleration':acceleration,'length':length},
            'clockwise_arrival':cw,'counterclockwise_arrival':ccw,
            'coordinate_arrival_difference':cw-ccw,
            'connection_loop':holonomy,
            'stationary_identity_error':cw-ccw+2*holonomy,
            'receiver_proper_time_difference':dtau,
            'GR_reference_arrival_difference':0.0,
            'scope':'Dimensionless local test-field transport example; no astrophysical amplitude calibration.'}


def run():
    result={'nonzero_loop':loop_example(),'conformal_null':loop_example(D=0),
            'uniform_lapse_null':loop_example(acceleration=0),
            'static_scalar_null':loop_example(q=0)}
    for case in result.values():
        if abs(case['stationary_identity_error'])>1e-11:
            raise AssertionError('Receiver timestamp and connection integral disagree')
    save('step_11_loop_observable.json',result)
    print(result['nonzero_loop'])
    return result

if __name__=='__main__':
    run()
