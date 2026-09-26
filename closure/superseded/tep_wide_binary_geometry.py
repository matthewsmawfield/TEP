#!/usr/bin/env python3
"""
Wide-binary recovery with projection and eccentricity applied explicitly (CORPUS_PLAN.md v5.2, section 7.3).

Per-system steep recovery f = [1 + (R_i / r)^k]^-1 in TRUE separation r; R_i lognormal(median, sigma_int);
observed projected separation s = r * proj; velocity factor from thermal eccentricity and projection.
Fits median, sigma_int, amplitude for k = 4, 6, 8 against Paper 13 Table 4.2 and reports the fraction of
midplane Sun-like systems whose own transition radius would violate the Saturn ephemeris bound.
"""
import math, random
random.seed(11)
s=[59,83,116,162,227,319,446,625,875,1225,1715,2402,3363,4709,6594,9233,12929,18105,25352]
v=[0.9068,1.0428,0.9948,1.0336,1.0220,1.0270,1.0610,1.0930,1.1253,1.1410,1.1824,1.2153,1.2359,1.2846,1.3291,1.3704,1.3942,1.3860,1.3499]
sg=[0.0795,0.0442,0.0235,0.0158,0.0103,0.0077,0.0063,0.0055,0.0052,0.0050,0.0050,0.0053,0.0058,0.0066,0.0076,0.0088,0.0104,0.0120,0.0147]
N=1500; geo=[]
for _ in range(N):
    cosi=random.uniform(-1,1); sini=math.sqrt(1-cosi*cosi); e=math.sqrt(random.random())
    Mn=random.uniform(0,2*math.pi); E=Mn
    for _ in range(20): E-=(E-e*math.sin(E)-Mn)/(1-e*math.cos(E))
    nu=2*math.atan2(math.sqrt(1+e)*math.sin(E/2),math.sqrt(1-e)*math.cos(E/2))
    proj=math.sqrt(max(1e-6,1-(sini*math.sin(nu+random.uniform(0,2*math.pi)))**2))
    kin=math.sqrt(2-(1-e*math.cos(E)))*math.sqrt(max(0.0,1-random.uniform(-1,1)**2))
    geo.append((proj,kin,random.gauss(0,1)))
def model(med,sig,amp,k):
    raw=[]
    for sb in s:
        vals=sorted(kin*math.sqrt(proj)*math.sqrt(1+amp/(1+((med*math.exp(sig*z))/(sb/proj))**k)) for proj,kin,z in geo)
        raw.append(vals[len(vals)//2])
    b=sum(raw[:5])/5; return [x/b for x in raw]
chi=lambda m: sum(((a-bb)/e)**2 for a,bb,e in zip(v,m,sg))
Phi=lambda q: 0.5*(1+math.erf(q/math.sqrt(2)))
res={}
for k in (4,6,8):
    best=None
    for med in range(1500,4251,250):
        for sig in [x/10 for x in range(6,17,1)]:
            for amp in [x/10 for x in range(7,14,1)]:
                x=chi(model(med,sig,amp,k))
                if best is None or x<best[0]: best=(x,med,sig,amp)
    x,med,sig,amp=best
    need=9.537*((2*amp/1.5e-10))**(1.0/k)
    Rsun=7131*med/2646.0
    res[k]=(x,med,sig,amp,need,Phi(math.log(need/med)/sig),Phi(math.log(need/Rsun)/sig))
    print(f"k={k}: chi2={x:6.1f} median={med} sigma_int={sig} amp={amp}  min safe R_sun={need:6.0f} AU  "
          f"P(unsafe | population)={res[k][5]:.1%}  P(unsafe | midplane-scaled median {Rsun:.0f} AU)={res[k][6]:.1%}")
b=min(r[0] for r in res.values())
print("\nDelta chi2 relative to best k:", {k:round(r[0]-b,1) for k,r in res.items()})
print("chi2_nu of the canonical fit is ~5, so divide Delta chi2 by ~5 for an effective significance.")

import json, os
os.makedirs("results", exist_ok=True)
json.dump({str(k): {"chi2": r[0], "median_R_AU": r[1], "sigma_int": r[2], "amplitude": r[3],
                    "min_safe_R_sun_AU": r[4], "unsafe_fraction_population": r[5],
                    "unsafe_fraction_midplane": r[6]} for k, r in res.items()},
          open("results/tep_wide_binary_geometry.json", "w"), indent=2)
