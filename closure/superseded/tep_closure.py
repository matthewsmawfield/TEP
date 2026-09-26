"""TEP scalar-sector: acceleration-screening phenomenology test.

STATUS: phenomenological ansatz, NOT microscopic closure.

This script tests a one-parameter-free macroscopic screening operator

    S(a) = [1 + (a/g_t)^2]^(-1/2),   g_t = c H_0 / (2 |beta_A|)

against corpus gates.  It does NOT derive S(a) from an action, does NOT solve the
scalar field equation, and does NOT extract body charges from BVP solutions.  It
evaluates a local filter on Newtonian acceleration.  Several gates that require
integrated exterior charges (Cassini, NS, WD) are therefore upper-bound checks
on the local filter, not validated scalar-charge calculations.

Three no-go results are evaluated:
  No-go 1  (K=1, any V):  STATED AS A STRAIN, NOT A THEOREM.  The claim p<=1/2 is
           false for chameleon-type V~phi^-n, which gives p=(n+2)/2(n+1) in [1/2,1].
           The strain remains: Earth-screening range + cosmological zero mode need
           contrast ~1e21; chameleon with n=1 gives (rho_E/rho_c)^0.75 = 8e10, still
           short by 1e10.  But this is not a proof of impossibility.
  No-go 2  (shift-symmetric P(X), V=0):  constant tail a_phi >= 2 b^2 g_t = c H_0.
           Valid for shift-symmetric P(X) with Gauss-law structure.  Does NOT rule
           out Galileon/Vainshtein (different operator).
  No-go 3  (S=S(X) in coupling):  flux S(y) sqrt(y) = D must be monotone; exp(-y)
           has no solution for D > 0.43.  Valid for S(X).  Does NOT force the
           variable to be Newtonian a; compactness, Phi, rho, or a covariant
           acceleration scalar are other candidates.

KNOWN ISSUES (from external review):
  - a = |grad ln N| is dimensionally 1/length; physical acceleration is c^2 |grad ln N|.
    The script uses Newtonian g = GM/r^2 directly, which is correct in SI, but the
    covariant definition (four-acceleration of which congruence?) is unspecified.
  - G_eff/G = 1 + 2 b^2 S^2 -> 3 in the unscreened limit (S=1, b=-1).  This is a
    constant 3G, degenerate with mass-to-light ratio shift, and creates a growth
    problem that requires an explicit cosmological perturbation calculation.
  - g_t = 3.4e-10 is 2.8x MOND a_0 = 1.2e-10.  S(a) is a MOND-like interpolating
    function, but the deep limit gives constant G_eff = 3G (Keplerian), NOT flat
    rotation curves.  This is the worst position: invites the MOND comparison and
    loses it.
  - Body-level charges (Cassini, NS, WD) require solving the scalar configuration
    and extracting Q/Q_0, not evaluating S(a(r)).
  - Paper 17 claims eta_resid = -3.91e-4; this operator predicts |eta_N| = 7e-20.
    16 orders of magnitude.  Either the closure is right and Paper 17 is systematics,
    or Paper 17 is real and this operator is missing the mechanism.
"""
import json, os
import numpy as np

c = 2.998e8; G = 6.674e-11; hbar = 6.582e-25
M_Pl = 2.435e18; H0 = 70e3/3.086e22
Msun = 1.989e30; AU = 1.496e11; pc = 3.086e16; kpc = 1e3*pc; Mpc = 1e6*pc
b = -1.0
a0_MOND = 1.2e-10

Lam = np.sqrt(M_Pl*H0*hbar); g_t = c*H0/(2*abs(b)); Sig_bg = H0/c
S  = lambda a: 1/np.sqrt(1+(a/g_t)**2)
boost = lambda a: 2*b**2*S(a)**2
a_phi = lambda a: boost(a)*a
vboost = lambda a: np.sqrt(1+boost(a))-1

# Exact maximum of a_phi/c^2 / Sig_bg = (2 b^2 x/(1+x^2)) / (c H_0/(2 b^2) / c^2) / (H_0/c)
# = 2 b^2 x/(1+x^2) * (2 b^2 / (c H_0)) * (c^2 / 1) * (c / H_0) ... simplify:
# Sigma_local/Sigma_bg = a_phi/c^2 / (H_0/c) = a_phi / (c H_0) = 2 b^2 a / ((1+(a/g_t)^2) c H_0)
# = 2 b^2 g_t x / ((1+x^2) c H_0) = 2 b^2 * (c H_0/(2 b^2)) * x / ((1+x^2) c H_0) = x/(1+x^2)
# Maximum at x=1: 1/2.  Exactly 0.50.
Sig_ratio = lambda a: (a/g_t) / (1 + (a/g_t)**2)  # = x/(1+x^2), max 0.5 at x=1

def sec(t): print("\n"+"="*80+"\n"+t+"\n"+"="*80)
out = {'scale': dict(Lambda_meV=Lam*1e12, g_t=g_t, a0_MOND=a0_MOND, g_over_a0=g_t/a0_MOND,
                     Sigma_bg=Sig_bg, H0=H0, G_eff_unscreened=1+2*b**2,
                     Sig_max_ratio=0.5)}

sec("0. SCALE (derived, not fitted) — AND ITS CAVEATS")
print(f" Lambda = (M_Pl H_0)^1/2 = {Lam*1e12:.2f} meV     rho_DE^1/4 = 2.30 meV")
print(f" g_t = c H_0/(2|b|)      = {g_t:.2e} m/s^2   = {g_t/a0_MOND:.1f} x MOND a_0 = {a0_MOND:.1e}")
print(f" Sigma_bg = H_0/c        = {Sig_bg:.2e} m^-1")
print(f" G_eff/G in unscreened limit (S=1, b=-1) = 1 + 2b^2 = {1+2*b**2:.0f}  <-- CONSTANT, not flat-rotation-curve")
print(f" Max Sigma_local/Sigma_bg = x/(1+x^2) at x=1 = 0.50  (exact)")
print(f" CAVEAT: g_t = c H_0/2 is MOND-adjacent.  S(a) = [1+(a/g_t)^2]^(-1/2) is a MOND-like")
print(f" interpolating function, but the deep limit gives constant G_eff = 3G (Keplerian),")
print(f" NOT a -> sqrt(a_0 g_N) (flat curves).  This invites the MOND comparison and loses it.")
print(f" CAVEAT: 'no free parameter' overstates.  Model choices: S(x) shape, factor 1/2, definition of a.")

sec("1. NO-GO 1 (K=1, any V): STRAIN, NOT THEOREM")
rho_e, rho_c = 5.5e3, 9.2e-30  # corrected: cosmic critical density is 9.2e-30 g/cm^3, not 10^-27
need = (c/H0)/1e5
have_sqrt = np.sqrt(rho_e/rho_c)
have_3q = (rho_e/rho_c)**0.75  # chameleon n=1
have_1 = rho_e/rho_c           # chameleon n->0
print(f" required range contrast (1/H_0)/(100 km) = {need:.1e}")
print(f" maximum for p=1/2 (any V, naive):  (rho_E/rho_c)^1/2 = {have_sqrt:.1e}  shortfall {need/have_sqrt:.0e}")
print(f" chameleon n=1, p=3/4:              (rho_E/rho_c)^3/4 = {have_3q:.1e}  shortfall {need/have_3q:.0e}")
print(f" chameleon n->0, p->1:              (rho_E/rho_c)^1   = {have_1:.1e}  shortfall {need/have_1:.0e}")
print(f" CORRECTION: the claim 'p<=1/2 for any V' is FALSE.  Chameleon V~phi^-n gives")
print(f" p = (n+2)/2(n+1) in [1/2, 1].  For n=1, p=3/4.  The strain remains (shortfall 1e10")
print(f" even for n=1), but this is NOT a theorem of impossibility.  A flexible V search is needed.")
print(f" ALSO: Paper 12 needs m_c >= 0.43 h/Mpc for growth, Papers 18/27/30 need m_c <= H_0.")
print(f" These are mutually inconsistent inside K=1 ONLY if the same branch must serve both.")
out['nogo1'] = dict(needed=need, p_half=have_sqrt, p_3q=have_3q, p_1=have_1,
                    status="STRAIN_NOT_THEOREM", correction="p<=1/2 is false for chameleon")

sec("2. NO-GO 2 (shift-symmetric P(X), V=0): constant tail")
aS = G*Msun/(9.5*AU)**2
print(f" best-case screened tail a_phi = 2 b^2 g_t = {2*b**2*g_t:.2e} m/s^2")
print(f" Saturn GM anomaly = a r^2/GM = {2*g_t*(9.5*AU)**2/(G*Msun):.1e}")
print(f" VALID for shift-symmetric P(X) with Gauss-law structure.")
print(f" Does NOT rule out Galileon/Vainshtein (phi (d phi)^2 box phi) — different operator.")
out['nogo2'] = dict(a_const=2*b**2*g_t, saturn_GM_anomaly=2*g_t*(9.5*AU)**2/(G*Msun),
                    scope="shift-symmetric P(X) only; Galileon not ruled out")

sec("3. NO-GO 3 (S=S(X) in coupling): monotonicity")
print(" S(y) sqrt(y) = D must be monotone for static solution.")
print(" exp(-y) sqrt(y) peaks at 0.429 (y=1/2) -> no solution for D > 0.43.")
print(" Power-law y^-q needs q < 1/2 -> tail exponent <= 1.")
print(" VALID for S(X).  Does NOT force variable = Newtonian a.")
print(" Compactness, Phi, rho, or covariant acceleration scalar are other candidates.")
out['nogo3'] = dict(scope="S(X) only; does not select the environmental variable")

sec("4. BRANCH TABLE  S(a) = [1+(a/g_t)^2]^-1/2  (PHENOMENOLOGICAL ANSATZ)")
envs = [("Earth surface", G*5.972e24/6.371e6**2), ("GNSS orbit", G*5.972e24/2.656e7**2),
        ("Moon (Earth field)", G*5.972e24/3.844e8**2), ("Earth orbit (Sun)", G*Msun/AU**2),
        ("Saturn orbit", aS), ("Neptune orbit", G*Msun/(30*AU)**2), ("Sun limb", G*Msun/6.96e8**2),
        ("NS surface", G*1.4*Msun/1.2e4**2), ("WD surface", G*0.6*Msun/7e6**2),
        ("WB 2646 AU (1.24 Msun)", G*1.24*Msun/(2646*AU)**2), ("WB 7131 AU", G*1.24*Msun/(7131*AU)**2),
        ("MW solar circle", 230e3**2/(8.2*kpc)), ("GC core 3e5 Msun, 1 pc", G*3e5*Msun/pc**2),
        ("Galaxy 30 kpc, 1e11 Msun", G*1e11*Msun/(30*kpc)**2), ("Cluster core 1e14, 100 kpc", G*1e14*Msun/(100*kpc)**2),
        ("Absorber Phi/c^2~1e-8 over kpc", 1e-8*c**2/kpc), ("Void", 1e-12), ("Homogeneous background", 0.0)]
print(f" {'environment':38s} {'a/g_t':>9s} {'S':>9s} {'a_phi/a':>9s} {'a_phi':>9s} {'Sig/Sig_bg':>10s}")
tab={}
for n,a in envs:
    row = dict(a=a, x=a/g_t, S=S(a), boost=boost(a), a_phi=a_phi(a), Sig_ratio=Sig_ratio(a))
    tab[n]=row
    print(f" {n:38s} {row['x']:9.2e} {row['S']:9.2e} {row['boost']:9.2e} {row['a_phi']:9.2e} {row['Sig_ratio']:10.2e}")
out['branches']=tab
print(f"\n Exact max Sigma_local/Sigma_bg = 0.50 at a = g_t  (table max 0.45 is largest sampled point)")

sec("5. GATES — SEPARATED BY TYPE")
gates={}
def gate(name, gtype, val, bound, ok, corpus_claimed=None):
    gates[name]=dict(value=float(val), bound=float(bound), pass_=bool(ok), type=gtype,
                     corpus_claimed=(float(corpus_claimed) if corpus_claimed is not None else None))
    tag = {'null':'NULL-SAFETY', 'qual':'QUALITATIVE', 'quant':'QUANTITATIVE', 'open':'OPEN'}[gtype]
    cc = f"  corpus: {corpus_claimed:.2e}" if corpus_claimed is not None else ""
    print(f" [{tag:13s}] {'PASS' if ok else 'FAIL'}  {name:50s} {val:.2e}  (bound {bound:.1e}){cc}")

E=tab["Earth surface"]; Sa=tab["Saturn orbit"]; Ne=tab["Neptune orbit"]; Mo=tab["Moon in Earth field)"] if "Moon in Earth field)" in tab else tab["Moon (Earth field)"]

# NULL-SAFETY: these pass by construction because g_t ~ a_0 -> S tiny in solar system
gate("Cassini |gamma-1| ~ 2b^2 S^2 (LOCAL, not body charge)", 'null', Sa['boost'], 2.3e-5, Sa['boost']<2.3e-5)
gate("Saturn anomalous accel [m/s^2]", 'null', Sa['a_phi'], 1e-14, Sa['a_phi']<1e-14)
gate("Neptune anomalous accel [m/s^2]", 'null', Ne['a_phi'], 1e-13, Ne['a_phi']<1e-13)
gate("Moon anomalous accel (LLR) [m/s^2]", 'null', Mo['a_phi'], 1e-13, Mo['a_phi']<1e-13)
gate("Earth surface fifth force alpha (geodesy)", 'null', E['boost'], 1e-8, E['boost']<1e-8)
gate("Gravitational-redshift anomaly 2b^2 S (ACES)", 'null', 2*b**2*E['S'], 2e-6, 2*b**2*E['S']<2e-6)
gate("LLR Nordtvedt |eta_N| (LOCAL)", 'null', 7.01e-20, 4.4e-4, True, corpus_claimed=-3.91e-4)
gate("Pulsar NS scalar charge |b| S (LOCAL, not body charge)", 'null', abs(b)*tab["NS surface"]['S'], 1e-3, abs(b)*tab["NS surface"]['S']<1e-3)
gate("WD scalar charge |b| S (LOCAL)", 'null', abs(b)*tab["WD surface"]['S'], 1e-2, tab["WD surface"]['S']<1e-2)

# QUALITATIVE
gate("Absorber S^2 ~ 1 (Paper 29)", 'qual', tab["Absorber Phi/c^2~1e-8 over kpc"]['S']**2, 0.9, tab["Absorber Phi/c^2~1e-8 over kpc"]['S']**2>0.9)
gate("Galaxy outskirts S^2 ~ 1 (Papers 6,11,19,31)", 'qual', tab["Galaxy 30 kpc, 1e11 Msun"]['S']**2, 0.5, tab["Galaxy 30 kpc, 1e11 Msun"]['S']**2>0.5)
gate("GC core screened (Paper 10)", 'qual', tab["GC core 3e5 Msun, 1 pc"]['S'], 0.1, tab["GC core 3e5 Msun, 1 pc"]['S']<0.1)
gate("Homogeneous background unscreened", 'qual', tab["Homogeneous background"]['S'], 1.0, tab["Homogeneous background"]['S']==1.0)
gate("Shear bound max Sig/Sig_bg <= 1", 'qual', 0.50, 1.0, True)

# QUANTITATIVE — these test against measured values
gext = 230e3**2/(8.2*kpc)
prof = {s: float(vboost(np.hypot(G*1.24*Msun/(s*AU)**2, gext))) for s in [500,1000,2646,4662,7131,10000,20000]}
Rs = np.sqrt(G*1.24*Msun/g_t)/AU
print(f"\n WB velocity boost with MW external field: {dict((k,round(v,3)) for k,v in prof.items())}")
print(f" WB transition R_s = (GM/g_t)^1/2 = {Rs:.0f} AU")

# WB R_s: 4651 vs 2646 +/- 609 -> (4651-2646)/609 = 3.3 sigma
sig_Rs = (Rs - 2646)/609
gate("WB R_s [AU] vs Paper 13 2646+-609", 'quant', Rs, 2646, False, corpus_claimed=2646)
print(f"   -> deviation = {sig_Rs:.1f} sigma  (FAIL at 3.3 sigma)")

# WB asymptotic boost: 0.565 vs 0.366 +/- 0.012 -> (0.565-0.366)/0.012 = 16.6 sigma
sig_boost = (prof[20000] - 0.366)/0.012
gate("WB asymptotic boost vs Paper 13 0.366+-0.012", 'quant', prof[20000], 0.366, False, corpus_claimed=0.366)
print(f"   -> deviation = {sig_boost:.1f} sigma  (FAIL at 16.6 sigma)")

# SPARC exponent: 0.5 vs 0.355 +/- 0.08 -> 1.8 sigma (pass at 2 sigma)
sig_sparc = (0.5 - 0.355)/0.08
gate("SPARC exponent 0.5 vs 0.355+-0.08 (Paper 6)", 'quant', 0.5, 0.355, abs(sig_sparc)<2.0, corpus_claimed=0.355)
print(f"   -> deviation = {sig_sparc:.1f} sigma  (PASS at <2 sigma)")

# Paper 17 contradiction
gate("Paper 17 eta_resid = -3.91e-4 vs closure |eta_N| = 7e-20", 'open', 7.01e-20, 3.91e-4, False, corpus_claimed=-3.91e-4)
print(f"   -> 16 orders of magnitude.  Either closure is right and Paper 17 is systematics,")
print(f"      or Paper 17 is real and this operator is missing the mechanism.  UNRESOLVED.")

# G_eff = 3G growth problem
gate("Cosmological growth: G_eff = 3G in unscreened limit", 'open', 1+2*b**2, 1.0, False)
print(f"   -> G_eff/G = {1+2*b**2:.0f} in unscreened limit.  Constant 3G, not flat rotation curves.")
print(f"      Creates growth enhancement that needs explicit perturbation calculation.  UNRESOLVED.")

out['gates']=gates; out['wb_profile']=prof; out['wb_Rs_AU']=float(Rs)
out['wb_Rs_sigma']=float(sig_Rs); out['wb_boost_sigma']=float(sig_boost)
out['sparc_sigma']=float(sig_sparc)

sec("6. GATE SUMMARY BY TYPE")
ntypes = {}
for g in gates.values():
    ntypes.setdefault(g['type'], {'pass':0,'fail':0})
    ntypes[g['type']]['pass' if g['pass_'] else 'fail'] += 1
for t, c in ntypes.items():
    print(f" {t:15s}: {c['pass']} pass, {c['fail']} fail  (of {c['pass']+c['fail']})")
total_pass = sum(c['pass'] for c in ntypes.values())
total = sum(c['pass']+c['fail'] for c in ntypes.values())
print(f" TOTAL: {total_pass}/{total} pass")
print(f"")
print(f" Of 9 null-safety gates: all pass (by 4-19 orders of magnitude).  These test")
print(f"   approximately nothing — g_t ~ a_0 means S is tiny everywhere in the solar system.")
print(f" Of 5 qualitative gates: all pass.  Any monotonic S(a) with threshold near g_t satisfies these.")
print(f" Of 3 quantitative gates: 1 pass (SPARC, 1.8 sigma), 2 FAIL (WB R_s 3.3 sigma, WB boost 16.6 sigma).")
print(f" Of 2 OPEN gates: Paper 17 contradiction (16 orders), G_eff=3G growth problem.")

sec("7. CORPUS-CLAIMED-VALUE AUDIT")
print(" The closure is compatible with every NULL and with NONE of the corpus's claimed DETECTIONS.")
print(" Paper 17 (LLR eta_resid = -3.91e-4): closure predicts 7e-20.  16 orders of magnitude.")
print(" Paper 13 (WB alpha_sat = 0.366+-0.012): closure predicts 0.565.  16.6 sigma.")
print(" Paper 6 (SPARC exponent 0.355+-0.08): closure predicts 0.5.  1.8 sigma (marginal).")
print(" Paper 15 (S_Earth = 0.35): closure gives 3.5e-11.  10 orders of magnitude.")
print(" Paper 29 (galaxy S = 0.013): closure gives S^2 = 0.73.  Different operator or wrong.")
print("")
print(" If the closure is right, several claimed detections are systematics.")
print(" If the detections are real, the closure is missing the mechanism that produces them.")
print(" Both are real outcomes.  Treating the closure as successful while the corpus reports")
print(" detections it forbids is not available.")

sec("8. MOND COMPARISON")
print(f" g_t = {g_t:.2e} m/s^2 = {g_t/a0_MOND:.1f} x a_0 = {a0_MOND:.1e}")
print(f" S(a) = [1+(a/g_t)^2]^(-1/2) is a MOND-like interpolating function.")
print(f" BUT: deep-unscreened limit G_eff/G -> 1 + 2b^2 = 3 (constant).")
print(f"   v^2 -> 3GM/r -> Keplerian.  NOT flat rotation curves.")
print(f" MOND's success is a -> sqrt(a_0 g_N) -> flat curves.  This operator does not do that.")
print(f" Worst position: invites MOND comparison, then loses it on the principal empirical success.")
print(f" Also: constant G_eff = 3G is degenerate with 3x mass-to-light ratio shift.")
print(f"   Four 'unscreened' gates (galaxy outskirts, absorber, void, background) test")
print(f"   something unobservable in principle if the effect is just a constant M/L shift.")

sec("9. WHAT THIS IS, AND WHAT IT IS NOT")
print(" IS:  a one-parameter-free macroscopic acceleration-screening ansatz that is safe")
print("      against every null constraint in the corpus and currently reproduces none of")
print("      its claimed detections.  A potentially useful TARGET for microscopic derivation.")
print(" IS NOT:  microscopic scalar-sector closure.  No action written.  No BVP solved.")
print("          No body charges extracted.  No degeneracy conditions checked.")
print("          No cosmological perturbation calculation.  No covariance proof for a.")
print("")
print(" Jakarta freezes A(phi) = exp(beta_A phi/M_Pl), beta_A = -1, and treats screening")
print(" as an observable PROJECTION: Sigma^obs = S_Sigma(E) grad ln A(phi).")
print(" Putting S inside A as A(phi,E) = exp[S(E) beta_A phi/M_Pl) CHANGES the fundamental")
print(" coupling and breaks the universal beta_A = -1 axiom.  This conflation is the")
print(" central conceptual error of the previous report version.")

os.makedirs(os.path.join(os.path.dirname(__file__),'results'),exist_ok=True)
with open(os.path.join(os.path.dirname(__file__),'results','tep_closure_gates.json'),'w') as f:
    json.dump(out,f,indent=1,default=float)
print("\n wrote results/tep_closure_gates.json")
