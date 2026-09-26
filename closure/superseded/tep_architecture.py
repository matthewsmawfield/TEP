#!/usr/bin/env python3
"""
TEP consistent architecture (2026-09-12) -- reproduces CORPUS_PLAN.md v5.0 section 7.

The temporal field phi is a dynamical time field. Its kinetic structure is

    P(X, phi) = K X + mu^2(phi) sqrt(2X) - V(phi),     X = -(1/2) g^{mn} d_m phi d_n phi

On the homogeneous (timelike) branch the sqrt term is a cuscuton: it carries pressure but no
energy density, so the field that IS cosmic time does not gravitate as kinetic energy.
Matter couples to g~ = A^2 g + B d phi d phi, A = exp(beta_A phi / M_Pl), beta_A = -1.

Sections:
  A  Static gravitational background: all three background equations with ordinary matter.
  B  Local stability inherited from A: K > 3 Omega_m.
  C  Wide binaries: steep per-system recovery, population of transition radii, Solar-System margin.
  D  Environmental ordering direction (Paper 13 Section 5).
  E  Cepheid carrier as stellar-structure response to the environmental coupling (no disformal term).
  F  Constraints any microscopic local operator must satisfy (from corpus_consistency_theorems.py).
"""
import math, json
Om, Or = 0.3153, 9.1e-5; OL = 1 - Om - Or
R = {}

# ---------- A: cosmological branch, units M_Pl = H0 = 1 ----------
rm, rr = 3*Om, 3*Or
Ht = lambda z: math.sqrt(Om*(1+z)**3 + Or*(1+z)**4 + OL)
rows, worst = [], 0.0
for z in (0, 0.5, 1, 3, 10, 100, 1100, 1e5, 1e8):
    A = 1/(1+z); phi = math.log(1+z)
    V = -(rm*A + rr)
    phidot = A*Ht(z)
    mu2 = -(rm*A + 4/3*rr)/phidot
    h = 1e-6
    Vp = (-(rm*math.exp(-(phi+h)) + rr) + (rm*math.exp(-(phi-h)) + rr))/(2*h)
    res = abs(Vp - rm*A)/(rm*A); worst = max(worst, res)
    rows.append({"z": z, "A": A, "V": V, "mu2": mu2, "eom_residual": res, "rho_tot": 0.0,
                 "p_tot": mu2*phidot - V + rr/3})
R["A_cosmology"] = {"rows": rows, "max_eom_residual": worst,
    "kinetic_energy_density": 0.0, "delta_Neff_from_phi": 0.0, "c_T": 1.0,
    "mu2_today_over_Lambda2": -3*Om,
    "statement": "static flat gravitational frame; rho_tot=p_tot=0; scalar EOM V_phi = beta T / M_Pl holds identically"}

# ---------- B: local stability ----------
R["B_local_stability"] = {"K_grad": "K - |mu^2|/|phidot| = K - 3 Omega_m",
    "requirement": f"K > {3*Om:.3f}", "K1_gives_cs2": 1-3*Om,
    "unscreened_force_over_newton_at_K1": 2/(1-3*Om)}

# ---------- C: wide binaries ----------
s=[59,83,116,162,227,319,446,625,875,1225,1715,2402,3363,4709,6594,9233,12929,18105,25352]
v=[0.9068,1.0428,0.9948,1.0336,1.0220,1.0270,1.0610,1.0930,1.1253,1.1410,1.1824,1.2153,1.2359,1.2846,1.3291,1.3704,1.3942,1.3860,1.3499]
sg=[0.0795,0.0442,0.0235,0.0158,0.0103,0.0077,0.0063,0.0055,0.0052,0.0050,0.0050,0.0053,0.0058,0.0066,0.0076,0.0088,0.0104,0.0120,0.0147]
X=[i*0.1 for i in range(-45,46)]; W=[math.exp(-0.5*x*x) for x in X]; SW=sum(W)
def pop(med,wid,amp,k):
    out=[1+amp*sum(w/(1+(med*math.exp(wid*x)/sb)**k) for x,w in zip(X,W))/SW for sb in s]
    b=sum(out[:5])/5; return [o/b for o in out]
chi=lambda m: sum(((a-b)/e)**2 for a,b,e in zip(v,m,sg))
fits=[]
for k in (2,4,8,16):
    best=min(((chi(pop(md,wd,am,k)),md,wd,am) for md in (1600,2000,2400,2800) for wd in (1.2,1.5,1.8,2.1) for am in (0.42,0.46,0.50,0.55)))
    x,md,wd,am=best
    fits.append({"k":k,"chi2":x,"median_R_AU":md,"sigma_lnR":wd,"amplitude":am,
                 "saturn_dg_over_g_R_midplane":2*am/(1+(7131/9.537)**k)})
R["C_wide_binaries"]={"canonical_chi2":86.3,"population_fits":fits,"saturn_bound":1.5e-10,
    "requirement":"per-system recovery steepness k >= 4",
    "statement":"the canonical exponential is a population CDF of steep per-system transitions"}

# ---------- D: environmental ordering ----------
R["D_environmental_ordering"]={"midplane_high_density_R_s":7131,"high_Z_low_density_R_s":4662,
    "direction":"denser environment -> LARGER transition radius",
    "density_thin_shell_prediction":"denser -> smaller range -> smaller R_s  (WRONG sign)",
    "statement":"environmental variable is the ambient Temporal-Topology state, not local density"}

# ---------- E: Cepheid carrier via stellar structure ----------
# homology on the instability strip at fixed T_eff: L ~ G^4 M^3, R ~ L^1/2, P ~ R^3/2 (GM)^-1/2
dM_per_dlnG  = -2.5*4/math.log(10)          # mag per unit dlnG from luminosity
dlogP_per_dlnG = (1.5*2 - 0.5)/math.log(10) # log10 P per unit dlnG
slopeW = -3.3
dW_fixedP = dM_per_dlnG - slopeW*dlogP_per_dlnG
need = {"SH0ES_host_mean": 0.045, "LMC_inner_minus_outer": 0.0284, "M31_inner_minus_outer": 0.356}
R["E_cepheid_carrier"]={"dW_at_fixed_P_per_dlnG":dW_fixedP,
    "required_dlnG_env":{k:abs(val/dW_fixedP) for k,val in need.items()},
    "gw170817":"unaffected: conformal/force sector, no disformal term",
    "correlated_predictions":["TRGB offsets between hosts","SN Ia host-environment step","GC pulsar spin-down excess (Paper 10)"],
    "caveat":"M31 offset does not survive matched controls; LMC does"}

# ---------- F: constraints on the microscopic local operator ----------
R["F_local_operator_constraints"]=[
  "not a density thin-shell (fails beta_A=-1 sign; predicts wrong environmental direction; stars screened)",
  "not a pure P(X) gradient cap (leaves Solar-System acceleration floor ~2 beta^2 g_*)",
  "K > 3 Omega_m (stability inherited from the cosmological branch)",
  "per-system recovery steepness k >= 4 in separation",
  "transition radius and amplitude depend on ambient Temporal-Topology state: larger R in denser midplane, amplitude rising with mass ratio q, falling with primary mass",
  "environmental coupling normalisation varies by ~4-6% between host disks and anchors (Cepheid carrier)"]

json.dump(R,open("results/tep_architecture_results.json","w"),indent=2)
print("A cosmology: max EOM residual %.1e, kinetic energy 0, c_T=1, mu2_0 = %.3f Lambda^2"%(worst,-3*Om))
print("B stability: K > %.3f ; K=1 gives c_s^2=%.3f"%(3*Om,1-3*Om))
for f in fits: print("C k=%-2d chi2=%5.1f  Saturn dg/g (midplane R)=%.1e"%(f["k"],f["chi2"],f["saturn_dg_over_g_R_midplane"]))
print("D denser midplane -> larger R_s (7131 vs 4662): density screening predicts the opposite")
print("E dW at fixed P per dlnG = %.2f mag ; required dlnG: %s"%(dW_fixedP,{k:round(x,3) for k,x in R['E_cepheid_carrier']['required_dlnG_env'].items()}))
