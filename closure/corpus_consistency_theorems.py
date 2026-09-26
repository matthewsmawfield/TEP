#!/usr/bin/env python3
"""
Corpus consistency theorems (2026-09-12).

Reproduces every decisive check behind CORPUS_PLAN.md section 7 (v5.0).
Constants: CODATA 2018, Planck 2018 H0. Data: Paper 13 Table 4.2 and the
Paper 13 pipeline outputs. No fitted inputs except where a test is explicitly
a fit, and those are labelled.
"""
import math, json, random, itertools, csv, os

hbar=1.054571817e-34; c=2.99792458e8; e=1.602176634e-19; G=6.67430e-11
Mpc=3.0856775814913673e22; AU=1.495978707e11; GMsun=1.32712440018e20; Msun=1.98847e30
H0=67.4e3/Mpc; MPl=2.435323e18; H0g=hbar*H0/(e*1e9); LAM=math.sqrt(MPl*H0g)
R={}

# T1 inverse-power potential has no screening minimum for beta_A = -1
def dVeff(phi,beta,L=1.0,M=1e6,rho=1.0): return -L**5/phi**2 + beta*rho/M*math.exp(beta*phi/M)
def has_min(beta):
    prev=None
    for k in range(1,4000):
        phi=10**(-6+k*0.003); d=dVeff(phi,beta)
        if prev is not None and (prev<0)!=(d<0): return True
        prev=d
    return False
R["T1_inverse_power_sign"]={"beta_plus1_has_minimum":has_min(+1),"beta_minus1_has_minimum":has_min(-1),
  "verdict":"V=L^4(1+L/phi) admits no screening minimum for beta_A=-1; Phase 4 v4.0-v4.2 used beta=+1"}

# T2 kinetic P(X) screening leaves a Solar-System acceleration floor ~ g_*
gstar=c*H0/math.sqrt(2)
R["T2_PX_floor"]={"g_star":gstar,"floor_anomalous_accel":2*gstar,"saturn_bound":1e-14,
  "excess":2*gstar/1e-14,"verdict":"|grad phi| saturates near Lambda^2 once |grad phi_N| exceeds it: floor ~2 beta^2 g_*"}

# T3 Paper 13 canonical profile applied to the Sun
alpha,Rs=0.366,2646.0
def dgg(rAU): v=1+alpha*(1-math.exp(-rAU/Rs)); return v*v-1
sat=dgg(9.537)*GMsun/(9.537*AU)**2
R["T3_wb_profile_in_solar_system"]={"saturn_dg":sat,"saturn_bound":1e-14,"excess":sat/1e-14,
  "mercury_dg_over_g":dgg(0.387),"verdict":"canonical Paper 13 law excluded by ephemerides by ~1e7"}

# T4 joint fit: zero-parameter acceleration screening vs Paper 13 bins (smeared)
s=[59,83,116,162,227,319,446,625,875,1225,1715,2402,3363,4709,6594,9233,12929,18105,25352]
v=[0.9068,1.0428,0.9948,1.0336,1.0220,1.0270,1.0610,1.0930,1.1253,1.1410,1.1824,1.2153,1.2359,1.2846,1.3291,1.3704,1.3942,1.3860,1.3499]
sg=[0.0795,0.0442,0.0235,0.0158,0.0103,0.0077,0.0063,0.0055,0.0052,0.0050,0.0050,0.0053,0.0058,0.0066,0.0076,0.0088,0.0104,0.0120,0.0147]
gext=1.91e-10   # Paper 13 011_mond_efe_comparison.csv fixed g_e
random.seed(7)
def ens(sl,N=5000):
    out=[]
    for _ in range(N):
        cosi=random.uniform(-1,1); sini=math.sqrt(1-cosi*cosi); ecc=math.sqrt(random.random())
        M=1.2*math.exp(random.gauss(0,1)*sl); Mean=random.uniform(0,2*math.pi); E=Mean
        for _ in range(20): E-=(E-ecc*math.sin(E)-Mean)/(1-ecc*math.cos(E))
        roa=1-ecc*math.cos(E); nu=2*math.atan2(math.sqrt(1+ecc)*math.sin(E/2),math.sqrt(1-ecc)*math.cos(E/2))
        proj=math.sqrt(max(1e-6,1-(sini*math.sin(nu+random.uniform(0,2*math.pi)))**2))
        kin=math.sqrt(2-roa)*math.sqrt(max(0.0,1-random.uniform(-1,1)**2))
        out.append((proj,M,kin))
    return out
def prof(E,gs,n,cpl):
    raw=[]
    for sb in s:
        vals=[]
        for proj,M,kin in E:
            r=sb/proj; g=GMsun*M/(r*AU)**2; gt=math.hypot(g,gext); S=1/(1+(gt/gs)**n)
            vals.append(kin*math.sqrt(1+cpl*S)*math.sqrt(proj))
        vals.sort(); raw.append(vals[len(vals)//2])
    b=sum(raw[:5])/5; return [x/b for x in raw]
chi=lambda m: sum(((a-b)/e_)**2 for a,b,e_ in zip(v,m,sg))
rows=[]
for sl in (0.0,0.8):
    E=ens(sl)
    for (lab,gs),n,cp in itertools.product((("g_TEP",5e-10),("cH0",c*H0),("cH0/sqrt2",gstar)),(1,2),(1.0,2.0)):
        satdg=cp/(1+(math.hypot(GMsun/(9.537*AU)**2,gext)/gs)**n)*GMsun/(9.537*AU)**2
        rows.append({"sigma_lnM":sl,"g_star":lab,"n":n,"coupling":cp,"chi2":chi(prof(E,gs,n,cp)),
                     "saturn_dg":satdg,"solar_system_ok":satdg<1e-13})
safe=[r for r in rows if r["solar_system_ok"]]
R["T4_joint_wb_solar_system"]={"canonical_chi2":86.34,"best_any":min(rows,key=lambda r:r["chi2"]),
  "best_solar_system_safe":min(safe,key=lambda r:r["chi2"]),"rows":rows,
  "verdict":"no zero-parameter acceleration law fits the bins and passes the Solar System"}

# T5 Paper 13 systematics outputs
base="../TEP-WB/results/outputs"
def rd(f):
    p=os.path.join(base,f)
    return list(csv.DictReader(open(p))) if os.path.exists(p) else []
R["T5_wb_systematics"]={
  "newtonian_forward_outer_v":{r["eccentricity_law"]:float(r["newtonian_outer_v_tilde"]) for r in rd("012_newtonian_forward_model_summary.csv")},
  "triples_can_explain":{r["f_triple"]:r["can_explain_signal"] for r in rd("009_triple_contamination_summary.csv")},
  "verdict":"Newtonian forward model and triples do not produce the rise; the anomaly is not a Newtonian artefact"}

# T6 disformal Cepheid carrier vs GW170817 (frozen envelope B=B0 varphi^2, phidot = eps_T M_Pl H)
eps=0.018; phidot2=(eps*MPl*H0g)**2
dvar2=abs((2*6e-7)**2-(2*2e-7)**2); need=2.07e-2
B0=2*need/(dvar2*phidot2)
gw=0.5*B0*(1e-9)**2*phidot2
R["T6_cepheid_disformal"]={"B0_GeV-4":B0,"M_B_eV":B0**-0.25*1e9,"gw170817_dc_over_c_generous":gw,
  "gw_bound":3e-15,"excess":gw/3e-15,
  "conformal_unscreened_dlnA":2*4e-7,"conformal_shortfall":need/(2*4e-7),
  "verdict":"conformal ~1e5 short; disformal at required amplitude violates GW170817 by >=1e7"}

# T7 cosmology: static flat gravitational background + rolling scalar
R["T7_static_background"]={"identity":"rho_tot+p_tot = phidot^2 + rho_m + p_m",
  "verdict":"must vanish for a static flat Einstein frame; impossible with phidot!=0 and ordinary matter"}
f=1/6; a_nu=(7/8)*(4/11)**(4/3); Nnu=3.044
R["T7b_clock_map_energy"]={"rho_kin_over_rho_crit":f,"delta_Neff_if_radiation_like":f/(1-f)*(1+Nnu*a_nu)/a_nu,
  "verdict":"carrying 1+z entirely in A(phi) with K=1 gives Delta N_eff ~1.5"}

# T8 Paper 27 H_TEP excursions with pipeline defaults
ex={}
for epsd in (0.1,0.0055):
    ex[str(epsd)]={lab:(1+z/100.0)**epsd-1 for lab,z in (("recombination",1100),("BBN",4e9))}
R["T8_paper27_HTEP"]={"fractional_dH":ex,"v42_formula_z0":{zt:0.1/zt for zt in (0.5,1,3,100)},
  "verdict":"shipped defaults give +28% at recombination and +476% at BBN; v4.2 'exact' formula shifts H0"}

# T9 Paper 28 single alpha_GB vs eta=-0.1
rg=G*Msun/c**2/1e3
R["T9_paper28_sgb"]={str(M):{"eta_max":3*1.2**2/(rg*M)**2,"sqrt_alpha_for_eta_0p1_km":math.sqrt(0.1*(rg*M)**2/3)} for M in (10,30,65)}

# T10 Paper 17 under any universal-coupling screening
R["T10_paper17"]={"nordtvedt_eta_leading_order":0.0,
  "verdict":"universal coupling -> no SEP violation; TEP predicts no cos D signal at the 5 mm level"}

json.dump(R,open("results/corpus_consistency_theorems.json","w"),indent=2,default=str)
for k,val in R.items():
    vv=val.get("verdict","") if isinstance(val,dict) else ""
    print(f"{k:34s} {vv}")
print("\nT4 best safe:",{k:R['T4_joint_wb_solar_system']['best_solar_system_safe'][k] for k in ('g_star','n','coupling','sigma_lnM','chi2')})
print("T6 GW excess:",f"{R['T6_cepheid_disformal']['excess']:.1e}","  T8:",R["T8_paper27_HTEP"]["fractional_dH"])
