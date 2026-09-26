"""Stage 6: Corpus-wide test of the two candidate closure classes.

Class I  (Stages 0-5):  K = 1, V = Lambda^4 (cosh(phi/M_Pl) - 1)   [potential / mass screening]
Class II (this stage):  K(X) noncanonical, V ~ 0, single scale Lambda [gradient / shear screening]

The corpus screening operator (Paper 26, l.96) is
    S_Sigma = [ 1 + (Sigma_mu Sigma^mu / g_t^2)^n + (rho/rho_half)^2 ]^-1
i.e. suppression when the LOCAL temporal shear exceeds a threshold shear g_t.
In the unscreened branch (Paper 31 l.333)  delta phi = 2 beta_A M_Pl Phi_N, so
    Sigma_local = |grad ln A| = 2 beta_A^2 g / c^2 .
The only shear scale in a static eternal universe is the cosmic shear Sigma_bg = H_0/c
(z = H_0 r / c  <=>  ln(1+z) = int Sigma dl).  Setting g_t' = c^2 Sigma_bg/(2 beta_A^2):
    g_t' = c H_0 / (2 beta_A^2).
This fixes the kinetic scale  Lambda^4 = M_Pl^2 H_0^2 (up to O(1)), i.e. Lambda = (M_Pl H_0)^(1/2).

Nothing here edits manuscripts. All numbers are printed for the report.
"""
import numpy as np

# constants (SI unless noted)
c = 2.998e8; G = 6.674e-11; hbar_GeVs = 6.582e-25
M_Pl_GeV = 2.435e18
H0_kms_Mpc = 70.0
H0 = H0_kms_Mpc*1e3/3.086e22            # s^-1
H0_GeV = H0*hbar_GeVs
Msun = 1.989e30; AU = 1.496e11; pc = 3.086e16; kpc = 1e3*pc; Mpc = 1e6*pc
beta_A = -1.0

def sec(t): print("\n"+"="*78+"\n"+t+"\n"+"="*78)

sec("1. THE SINGLE SCALE")
Lam = np.sqrt(M_Pl_GeV*H0_GeV)              # GeV
rhoDE_quarter = 2.3e-12                     # GeV  (rho_Lambda^(1/4) = 2.3 meV)
g_t = c*H0/(2*beta_A**2)
print(f" H_0 = {H0:.3e} s^-1 = {H0_GeV:.3e} GeV")
print(f" Lambda = (M_Pl H_0)^(1/2) = {Lam:.3e} GeV = {Lam*1e12:.2f} meV")
print(f" rho_DE^(1/4)               = {rhoDE_quarter*1e12:.2f} meV   (same scale, ratio {rhoDE_quarter/Lam:.2f})")
print(f" Shear threshold: g_t' = c H_0 / (2 beta_A^2) = {g_t:.2e} m/s^2")
print(f" Corpus fitted g_TEP (Papers 6, 13)         = 5.0e-10 m/s^2   (ratio {5e-10/g_t:.2f})")
print(f" MOND a_0                                   = 1.2e-10 m/s^2   (ratio {1.2e-10/g_t:.2f})")
print(f" Cosmic shear Sigma_bg = H_0/c = {H0/c:.3e} m^-1")

sec("2. SCREENING VARIABLE ACROSS THE CORPUS  (g / g_t')")
rows = [("Earth surface (GNSS ground, LLR, geodesy)", G*5.972e24/6.371e6**2),
        ("GNSS orbit r=26560 km", G*5.972e24/2.656e7**2),
        ("Moon at Earth", G*5.972e24/3.844e8**2),
        ("Sun at Earth orbit (Cassini path)", G*Msun/AU**2),
        ("Sun at Saturn (Cassini conj.)", G*Msun/(9.5*AU)**2),
        ("Sun surface (Shapiro limb)", G*Msun/6.96e8**2),
        ("NS surface (pulsar J0437)", G*1.4*Msun/1.2e4**2),
        ("WD surface", G*0.6*Msun/7e6**2),
        ("Wide binary 2646 AU, M=1.24 Msun", G*1.24*Msun/(2646*AU)**2),
        ("Wide binary 7131 AU", G*1.24*Msun/(7131*AU)**2),
        ("MW at solar circle", 2.0e-10),
        ("Globular cluster core (Paper 10)", G*3e5*Msun/(1*pc)**2),
        ("Galaxy halo edge 30 kpc, 1e11 Msun", G*1e11*Msun/(30*kpc)**2),
        ("Cluster core 100 kpc, 1e14 Msun", G*1e14*Msun/(100*kpc)**2),
        ("Void / cosmic web", 1e-12),
        ("Absorber (Paper 29, Phi/c^2~1e-8, L~kpc)", 1e-8*c**2/kpc)]
print(f" {'environment':45s} {'g [m/s^2]':>10s} {'g/g_t':>10s}  branch")
for n,g in rows:
    r=g/g_t; br = "SCREENED" if r>10 else ("transition" if r>0.1 else "UNSCREENED")
    print(f" {n:45s} {g:10.2e} {r:10.2e}  {br}")

sec("3. SOLAR-SYSTEM GATES FOR A STEEP K(X)")
print(" Suppression of the source charge in the screened branch: S_Sigma = 1/K'(chi).")
print(" For K = chi + K0 chi^m the fifth force/gravity ~ 2 beta^2 (g_t/g)^p with p=(2m-2)/(2m-1) -> 1.")
for p in [2/3, 0.8, 1.0]:
    gsat = G*Msun/(9.5*AU)**2
    S_cass = (g_t/gsat)**p
    gam = 4*beta_A**2*S_cass/(1+2*beta_A**2*S_cass)
    print(f"  p={p:.2f}: Cassini-path S_Sigma(Saturn) = {S_cass:.2e}, gamma-1 = -4 beta^2 S/(1+2 beta^2 S) = {gam:.1e}  "
          f"{'PASS' if gam<2.3e-5 else 'FAIL'} (bound 2.3e-5) ; direct fifth-force ratio 2 beta^2 S = {2*S_cass:.1e}")
print(" NOTE: Cassini measures gamma via light bending/Shapiro; the fifth-force ratio 2 beta^2 S itself is bounded")
print("       by planetary ephemerides at ~1e-5..1e-4 (Mars/Saturn ranging). Under the linear")
print("       screened-source map only the p=1.0 branch satisfies Cassini at the Saturn conjunction.")
gE = G*5.972e24/6.371e6**2
for p in [0.8,1.0]:
    S_E = (g_t/gE)**p
    print(f"  p={p:.2f}: Earth surface S_Sigma = {S_E:.1e}; Yukawa-type alpha = 2 beta^2 S = {2*S_E:.1e}"
          f"  (geodesy bounds |alpha| <~ 1e-3..1e-8 for lambda>100 km: {'PASS' if 2*S_E<1e-8 else 'marginal'})")
    print(f"           clock-rate anomaly relative to GR redshift = 2 beta^2 S = {2*S_E:.1e}  (GP-A 7e-5, ACES 2e-6: PASS)")
    s_E = 6.95e-10*S_E; s_M = 3.13e-11*(g_t/(G*7.342e22/1.737e6**2))**p
    print(f"           Nordtvedt eta_N = 4 beta^2 (s_E - s_M) = {4*(s_E-s_M):.1e}  (LLR 4.4e-4: PASS)")
gNS = G*1.4*Msun/1.2e4**2
for p in [0.8,1.0]:
    print(f"  p={p:.2f}: NS scalar charge alpha_NS = beta_A S = {(g_t/gNS)**p:.1e}  (binary-pulsar dipole bound ~1e-3: PASS)")

sec("4. LOW-ACCELERATION SECTOR (Papers 6, 13): transition radius law")
print(" Gradient screening gives a transition at g(R_s) = g_t'  =>  R_s = (G M / g_t')^(1/2)  ~ M^(1/2)")
print(" Corpus geometric law (rho_T = 20 g/cc)      =>  R_T = (3M/4 pi rho_T)^(1/3)   ~ M^(1/3)")
for M in [1.24, 0.5, 2.0]:
    Rs = np.sqrt(G*M*Msun/g_t)/AU
    print(f"  M={M:4.2f} Msun: R_s(g_t') = {Rs:6.0f} AU  (Paper 13 measured 2646+-609 AU, high-|Z| 4662, midplane 7131)")
print("  External-field effect: MW g_ext ~ 2.0e-10 = 0.6 g_t' keeps the local branch partially screened ->")
print("  the saturation boost is capped below the bare 1+2 beta^2 = 3 (v boost 0.73); Paper 13 finds 0.366.")
print("  SPARC exponent: M^(1/2) predicts 0.50; Paper 6 measures 0.355 +- 0.043 (stat) +- 0.07 (def) -> 1.8 sigma.")
print("  => the low-g sector is compatible with a shear threshold; M^(1/2) vs M^(1/3) is the decisive test.")

sec("5. WHY CLASS I (K=1, cosh V) CANNOT WORK CORPUS-WIDE")
rho_e=5.5; rho_cos=9.2e-30
print(f" Canonical scalar: m_eff ~ rho^p, p <= 1/2 for any V.  rho_Earth/rho_cosmic = {rho_e/rho_cos:.1e}")
print(f" Max range ratio = {np.sqrt(rho_e/rho_cos):.1e}: an Earth-scale 4000 km range -> {4e6*np.sqrt(rho_e/rho_cos)/Mpc:.1e} Mpc in voids")
print(" Corpus needs: Gpc coherence (Papers 18,26,27,30), kpc gradients delta phi = 2 beta M_Pl Phi_N (Papers 11,19,29,31),")
print("               lambda_C <~ 14.6 Mpc on sigma_8 scales (Paper 12).  A 6000/sqrt(eta) km range fails all of these.")
print(" Local failures of the Stage-5 window (log eta in [-0.5,1.4]):")
for S in [0.285,0.145,0.05,8e-3]:
    print(f"   S_Earth={S:.3f}: Yukawa alpha=2S={2*S:.2f} at range ~R_Earth  -> EXCLUDED by geodesy (|alpha|<<1e-3)")
print("   NS unscreened S~0.65-1 -> alpha_NS ~ 1  -> EXCLUDED by binary-pulsar dipole radiation (~1e-3)")
print("   The Stage-5 GNSS observable uses the linear screened-source PPN map 2 beta^2 S/(1+2 beta^2 S);")
print("   the clock-rate observable is likewise linear: 2 beta^2 S.")

sec("6. COSMOLOGICAL SECTOR")
print(f" Lambda^4 = M_Pl^2 H_0^2 = {(M_Pl_GeV*H0_GeV)**2:.2e} GeV^4 ;  rho_DE = {rhoDE_quarter**4:.2e} GeV^4  (ratio {(M_Pl_GeV*H0_GeV)**2/rhoDE_quarter**4:.2f})")
print(" Background: the unscreened branch carries Sigma_bg = H_0/c, giving ln(1+z) = int Sigma dl (Paper 26) with a_m = 1.")
print(" Early-universe suppression S(z)=exp(-(z/z_T)^n) (Papers 10, 18) = high-gradient/high-density branch of the same operator.")
print(" Paper 27's Z(chi), U(chi) is the homogeneous reconstruction of this kinetic sector (massless, c_A^2 = 1).")

sec("7. DISFORMAL SECTOR WITH M_* = Lambda  (Paper 31: B = beta_B/M_*^4 F(X/M_*^4))")
print(" GW170817 along unscreened paths where |X| ~ Lambda^4:  eps_B = B X / A^2 ~ beta_B F(1) <~ 1e-15  => beta_B <~ 1e-15.")
print(" Screened terrestrial paths have |grad phi| suppressed by 1/K' -> eps_B smaller still.")
print(" Paper 28 quartic-Gaussian B(varphi) and Paper 29 constant/inverse-field B are then the |varphi|>>1 and |X|<<Lambda^4")
print(" limits of ONE function B(varphi, X); the weak-field Stage-2 window (B_0 ~ 1e18-1e21 dimensionless) must be re-derived")
print(" in this normalization before any holonomy forecast is quoted.")

sec("8. WHAT THIS CLOSURE DOES NOT EXPLAIN")
print(" - The GNSS/SLR covariance length lambda_T ~ 1400-4500 km is NOT a Compton or Vainshtein length in Class II.")
print("   (Vainshtein normalized to r_V(Earth)=4146 km puts r_V(Sun)=0.41 R_sun: Cassini fails by 1e4.)")
print("   Consistent with Jakarta l.514 / Paper 1 l.319-328 / Paper 6 l.193: lambda_T is a covariance length, not a mass.")
print(" - The rho_T = 20 g/cc geometric law (Papers 6,13,15) and the shear law g_t' both give R ~ M^q; q=1/3 vs 1/2 is testable.")
print(" - Paper 15 internal: S_Earth = 0.35 (l.209) with beta_eff = beta_A S_Earth is inconsistent with its own gamma bound;")
print("   its fitted b_flyby = 2.6e-3 is the number consistent with Cassini.")
print(" - Paper 28 (sGB hair) and Paper 29 (10^4 amplification) still require a coupled K(X)+sGB / K(X)+transport calculation.")
