# Master-Sector Admissibility and Cross-Domain Constraint Analysis

This document presents the systematic admissibility analysis of the Temporal Equivalence Principle (TEP) master action across the kinetic function $P(X,\phi)$, the scalar potential $V(\phi)$, the disformal coupling $B(\phi)$, and the cosmological ambient drift $A_{\rm sp}(\eta)$. Rather than postulating separate phenomenological screening laws for individual astronomical regimes, the functional forms of $\{P(X,\phi), V(\phi), B(\phi)\}$ are determined uniquely from mathematical well-posedness, elliptic/hyperbolic stability, and cross-domain observational constraints (`scripts/steps/step_54_master_sector_admissibility.py` $\rightarrow$ `results/step_54_master_sector_admissibility.json`).

---

## 1. Foundational Inputs and Invariant Conventions

The master-sector analysis holds the foundational geometric postulates of Paper 0 fixed across all domains:

| Quantity | Mathematical Form / Value | Reference |
| :--- | :--- | :--- |
| Causal matter metric | $\tilde g_{\mu\nu} = A^2(\phi)\,g_{\mu\nu} + B(\phi)\,\nabla_\mu\phi\nabla_\nu\phi$ | Paper 0, Section 1 |
| Universal conformal coupling | $A(\phi) = \exp(\beta_A\phi/M_{\rm Pl}),\quad \beta_A = -1$ | Paper 0, Section 2 |
| Dimensionless field convention | $u \equiv \phi/M_{\rm Pl} > 0$ in gravitational wells; $u = 0$ in the present ambient | Paper 0, Section 2 |
| Kinetic gradient invariant | $\xi \equiv \lvert X\rvert/\Lambda_X^4,\quad \Lambda_X^4 = M_{\rm Pl}^2 H_0^2\;\;(\Lambda_X \simeq 1.91\text{ meV})$ | Paper 0, Section 3 |
| Cosmic temporal-gradient scale | $g_t = \frac{c H_0}{2\lvert\beta_A\rvert} \simeq 3.40 \times 10^{-10}\text{ m s}^{-2}\;\;(H_0 = 70\text{ km s}^{-1}\text{ Mpc}^{-1})$ | Paper 0, Section 3 |
| Linearized lapse response | $\delta\phi/M_{\rm Pl} = 2\beta_A\Phi_N/c^2;\quad \nabla\ln\tilde N = (1 + 2\beta_A^2\mathcal{S}_\Sigma)\nabla\Phi_N/c^2$ | `step_26_weak_field_lapse_ledger.py` |

---

## 2. Kinetic Sector $P(X,\phi)$: Uniqueness of the Two-Branch Completion

### 2.1 Admissibility Criteria

The kinetic derivative $P_{,X}(\xi)$ governs both the static screening factor $\mathcal{S}_\Sigma(\xi) = P_{,X}(\xi)^{-1}$ on spacelike gradients and the cosmological drift stiffness on timelike gradients. Any admissible candidate must simultaneously satisfy twelve structural and observational conditions:

| Criterion | Physical Requirement |
| :--- | :--- |
| P-1 | Strong hyperbolicity and elliptic well-posedness across both signs of $X$ |
| P-2 | Absence of ghost instabilities: positive perturbation kinetic weight $P_{,X} > 0$ for $\xi > 0$ |
| P-3 | Legendre convexity and stability: $Z_\parallel \equiv P_{,X} + 2\xi P_{,XX} > 0$ |
| P-4 | High-gradient Solar-System screening: $\mathcal{S}_\Sigma = P_{,X}^{-1} \to (2\xi)^{-1} \to 0$ for $\xi \gg \xi_*$ |
| P-5 | Low-acceleration gravitational recovery: $\mathcal{S}_\Sigma$ enhances the effective force at $g \lesssim g_t$ |
| P-6 | Preservation of Solar-System PPN bounds ($\lvert\gamma_{\rm PPN}-1\rvert < 2.3 \times 10^{-5}$ at Cassini conjunction) |
| P-7 | Asymptotically flat galactic rotation curves: $a_\phi \propto 1/r$ ($g_{\rm obs} \simeq \sqrt{a_{\rm eff} g_N}$) for $\xi \ll \xi_*$ |
| P-8 | Parameter-free acceleration scale: $a_{\rm eff} = a_0 \simeq 1.20 \times 10^{-10}\text{ m s}^{-2}$ without free fitting |
| P-9 | Consistent cosmological timelike branch ($\xi_{\rm cosmo} = H^2(z)/(2H_0^2)$) |
| P-10 | Single universal kinetic energy scale $\Lambda_X^4 = M_{\rm Pl}^2 H_0^2$ |
| P-11 | Mutual two-body transition profile $y(s) \propto (r/r_*)^{4/3}$ in the screened interior |
| P-12 | Regularity of deep compact-object temporal wells as $A \to 0$ |

### 2.2 Tail-Exponent Theorem and Zero-Parameter Coefficient Derivation

For a general low-gradient power law $P_{,X}(\xi) \sim k\xi^p$, spherical flux conservation $\nabla\cdot(P_{,X}\nabla\phi) = 0$ outside an isolated baryonic mass $M_b$ implies $r^2 P_{,X}(\xi)\,\phi'(r) = \text{const}$, which yields the radial scaling:

$$\phi'(r) \propto r^{-\frac{2}{2p+1}}.$$

An asymptotically flat rotation curve ($v_{\rm circ}^2 = r\,a_\phi = \text{const}$, i.e., $a_\phi \propto r^{-1}$) requires $2/(2p+1) = 1$, which uniquely selects the exponent $p = 1/2$ within the power-law family (numerical verification across $p \in \{1/3, 1/2, 2/3, 1\}$ gives logarithmic slopes $-1.200, -1.000, -0.857, -0.667$).

Combining the $p = 1/2$ low-gradient branch with the quadratic $2\xi$ high-gradient screening branch yields the two-branch kinetic function:

$$P(X,\phi) = -V(\phi) + \frac{2k}{3\Lambda_X^2}\lvert X\rvert^{3/2} + \frac{1}{\Lambda_X^4}X^2
\quad\Longrightarrow\quad
P_{,X}(\xi) = k\sqrt{\xi} + 2\xi.$$

Retaining any canonical unit term ($P_{,X} = 1 + 2\xi$ or $1 + k\sqrt{\xi} + 2\xi$) restores a linear Laplace operator as $\xi \to 0$, forcing a Keplerian $1/r^2$ scalar tail and violating Criterion P-7. In the low-gradient limit ($\xi \ll \xi_*$) of the unit-free two-branch operator, the spherical flux equation gives:

$$a_\phi^2 = a_{\rm eff}\,g_N,
\qquad
a_{\rm eff} = \frac{2\sqrt{2}\,\lvert\beta_A\rvert^3\Lambda_X^2}{k M_{\rm Pl}} = \frac{4\sqrt{2}\,\lvert\beta_A\rvert g_t}{k}.$$

Matching $a_{\rm eff}$ to the universal acceleration scale $a_0 = 1.20 \times 10^{-10}\text{ m s}^{-2}$ (or equivalently setting the natural branch crossover $k = 16$, $\xi_* = (k/2)^2 = 64$, which gives $a_{\rm eff} = cH_0/(4\sqrt{2}) = 1.202 \times 10^{-10}\text{ m s}^{-2}$) determines the dimensionless coefficient:

$$k = \frac{4\sqrt{2}\,g_t}{a_0} \simeq 16.03,$$

deriving Milgrom's $a_0 \sim cH_0$ relation directly from the cosmological time field with zero free parameters.

### 2.3 Candidate Comparison Across Surveyed Kinetic Forms

| Kinetic Candidate $P_{,X}(\xi)$ | Asymptotic Slope | $P_{,X}>0$ | $Z_\parallel>0$ | Sound Speed $c_s^2$ | Admissibility Verdict |
| :--- | :---: | :---: | :---: | :---: | :--- |
| Baseline $1 + 2\xi$ | $-2.000$ | Yes | Yes | $[1/3, 1]$ | Rejected — Keplerian tail fails RAR (Criterion P-7) |
| Two-branch $k\sqrt{\xi} + 2\xi$ | $-1.000$ | Yes | Yes | $[1/3, 1/2]$ | Selected — unique minimal power-law completion |
| Unit-plus-two-branch $1 + k\sqrt{\xi} + 2\xi$ | $-2.000$ | Yes | Yes | $[1/3, 1]$ | Rejected — unit term restores $1/r^2$ tail at small $\xi$ |
| Exponential interpolation $1 - e^{-k\sqrt{\xi}} + 2\xi$ | $-1.000$ | Yes | Yes | $[1/3, 0.78]$ | Excluded by data — fails wide-binary amplitude ($30\sigma$) and $f\sigma_8$ growth ($\chi^2 = 44.2$) |
| Pure AQUAL $k\sqrt{\xi}$ | $-1.000$ | Yes | Yes | $1/2$ | Rejected — lacks high-$\xi$ quadratic screening (Criterion P-4) |

### 2.4 Regime Crossover and Cosmological Growth Stiffness

The crossover between the small-$\xi$ ($k\sqrt{\xi}$) and large-$\xi$ ($2\xi$) branches occurs at $\xi_* = (k/2)^2 \simeq 64.2$:

| Physical Regime | Typical $\xi = \lvert X\rvert/\Lambda_X^4$ | Active Branch | Screening Factor $\mathcal{S}_\Sigma(\xi)$ |
| :--- | :---: | :---: | :--- |
| Solar System (Cassini $1.6\,R_\odot$ to $1\text{ AU}$, Earth, LLR) | $4 \times 10^{10} \text{ – } 10^{23}$ | Large-$\xi$ ($2\xi$) | $10^{-23} \text{ – } 10^{-11}$ (deeply screened) |
| Galactic Center core | $4.2 \times 10^4$ | Large-$\xi$ ($2\xi$) | $1.2 \times 10^{-5}$ |
| Wide binaries ($s \simeq 2{,}646\text{ AU}$) | $9.5$ | Small-$\xi$ ($k\sqrt{\xi}$) |Transition regime ($\xi < \xi_*$) |
| Milky Way solar circle ambient | $0.35$ | Small-$\xi$ ($k\sqrt{\xi}$) | $0.098$ (external-field suppression) |
| Cosmic voids | $8.6 \times 10^{-6}$ | Small-$\xi$ ($k\sqrt{\xi}$) | Enhanced response (clock acceleration in voids) |

On the cosmological timelike branch ($\xi_{\rm cosmo}(z) = H^2(z)/(2H_0^2)$), the two-branch stiffness suppresses late-time unscreened structure growth from $G_{\rm eff}(0)/G = 2.00$ (baseline) to $G_{\rm eff}(0)/G = 1.16$:

| Redshift $z$ | $\mathcal{S}_\Sigma$ (Baseline) | $G_{\rm eff}/G$ (Baseline) | $\mathcal{S}_\Sigma$ (Two-Branch) | $G_{\rm eff}/G$ (Two-Branch) |
| :---: | :---: | :---: | :---: | :---: |
| $0$ | $0.500$ | $2.00$ | $0.081$ | $1.16$ |
| $1$ | $0.244$ | $1.49$ | $0.043$ | $1.09$ |
| $3$ | $0.048$ | $1.10$ | $0.014$ | $1.03$ |
| $10$ | $2.5 \times 10^{-3}$ | $1.005$ | $1.6 \times 10^{-3}$ | $1.003$ |
| $100$ | $3.2 \times 10^{-6}$ | $1.000$ | $3.2 \times 10^{-6}$ | $1.000$ |

---

## 3. Potential Sector $V(\phi)$: Master Potential Family

### 3.1 Structural Conditions on $V(\phi)$

| Criterion | Physical Requirement |
| :--- | :--- |
| V-1 | Existence of a stable matter-coupled minimum of $V_{\rm eff}(u;\rho) = V(u) + \rho e^{-u}$ on the $\beta_A = -1$ branch |
| V-2 | Environmental density-dependent effective mass $m_{\rm eff}(\rho) \propto \rho^p$ ($p > 0$) |
| V-3 | Locally flat interior plateau at large field depth $u \gg 1$ ($V_{,u} \to 0$) |
| V-4 | Asymptotic potential floor $V(u \to \infty) = V_0$ supporting regular de Sitter-like core geometry |
| V-5 | Compatibility with cosmological clock-rate evolution $A_{\rm sp}(t)$ across cosmic history |
| V-6 | Regular temporal-well interior with strictly positive lapse $N > 0$ at all finite radii |
| V-7 | Positive curvature $V''(u_{\min}) > 0$ at matter-domain equilibria |

### 3.2 Evaluation of Candidate Potentials

| Candidate Form $V(u)$ | Matter-Domain Minima | Mass Scaling $m_{\rm eff} \propto \rho^p$ | Admissibility Verdict |
| :--- | :---: | :---: | :--- |
| Quadratic $\frac{1}{2}m^2 u^2$ | $4/4$ | $p \approx 0.015 \simeq 0$ | Rejected — density-independent correlation length $\lambda_C$ |
| Inverse-power $\Lambda^4(1 + u^{-1})$ | $0/4$ (runaway) | — | Rejected — no minimum for $\beta_A = -1$ ($V_{,u}$ and $\rho A_{,u}$ have same sign) |
| Quartic $\frac{1}{4}\lambda M_{\rm Pl}^4 u^4$ | $4/4$ | $p = 1/3$ | Viable weak/intermediate limit ($u_{\min} \propto \rho^{1/3}$, logarithmic saturation at extreme $\rho$) |
| Exponential $\Lambda_V^4 e^u$ | $4/4$ | $p \approx 0.4$ | Viable auxiliary candidate |
| Symmetron $-\frac{1}{2}\mu^2 u^2 + \frac{1}{4}\lambda u^4$ | $4/4$ | $p \approx 0.15$ | Rejected — nonzero vacuum expectation value violates $u_{\rm amb} = 0$ convention |
| Master family $\frac{1}{4}\lambda M_{\rm Pl}^4 u^4 e^{-(u/u_s)^4} + V_0 e^{-(u_s/u)^4}$ | $4/4$ | $p \in [1/3, 1/2]$ | Selected — satisfies Criteria V-1 through V-7 simultaneously |

In the master potential family, the quartic term $\frac{1}{4}\lambda M_{\rm Pl}^4 u^4$ governs weak and moderate field depths ($u < u_c \approx 2.8\text{–}3.3$), while the flat floor $V_0 e^{-(u_s/u)^4}$ governs deep temporal wells ($u > u_c$, $u_s \approx 9\text{–}15$). The canonical quartic couplings in the shared corpus registry are $\lambda_{\rm ref} = 7.526 \times 10^{-71}$ and the Cassini-compatible weak-field branch $\lambda_{\rm Cassini} = 10^5\,\lambda_{\rm ref} = 7.526 \times 10^{-66}$.

---

## 4. Disformal Sector $B(\phi)$ and Ambient Drift Junction $A_{\rm sp}(\eta)$

### 4.1 Disformal Envelope Constraints

The disformal coupling $B(\phi)$ is governed by five physical requirements:
1. Lorentzian signature preservation: $B(\nabla\phi)^2 > -A^2$ pointwise (`step_52_admissibility_margin.py`).
2. Hyperbolicity of matter-fluid perturbations: $Z_t > 0$ and $Z_s > 0$, which hold unconditionally on the $B \ge 0$ branch and exclude negative-branch reconstructions where $Q \to -1$ (`step_13_disformal_closure.py`, `step_52_admissibility_margin.py`).
3. Multimessenger cone tightness: $B(\phi \simeq 0) = 0$ in the intergalactic medium so that $\lvert c_{\rm GW} - c_\gamma\rvert/c \lesssim 10^{-15}$ along the GW170817 path (`step_68_gw_ambient_baseline.py`).
4. Nonvanishing closed-loop synchronization holonomy: $d(\delta\tilde\sigma) \neq 0$ in spatially inhomogeneous gravitational wells.
5. Pure transverse-traceless tensor polarization in the gravitational frame.

Monomial couplings ($B \propto \phi^2$ or $B \propto \phi^5$) fail Solar-System cone bounds or high-redshift signature margins. The admissible class is the localized bump or envelope family,

$$B(u) = B_0\,\frac{u^2}{1+u^2}\,\exp\!\left(-\frac{u^4}{2\sigma_B^4}\right) \quad \text{or} \quad B(u) = B_0\,u^2\exp\!\left(-\frac{u^2}{u_c^2}\right),$$

gated by the local environmental gradient invariant $G(X_{\rm local})$ (`step_68_gw_ambient_baseline.py`).

### 4.2 Piecewise $C^1$ Ambient Drift Junction

In `step_64_ambient_junction.py`, matching the late-time conformal clock history $A_{\rm clock}(z) = (1+z)^{-1}$ with $C^1$ continuity of the field and its derivative onto an asymptotically regular past temporal-horizon power law $A_{\rm clock}(\eta) = C_{\rm tail}\,\eta^{-p}$ determines the effective index at the junction redshift $z_j$:

$$p_{\rm eff}(z_j) = \frac{E(z_j)}{1+z_j}\,\chi(z_j),
\qquad
\chi(z_j) = \int_0^{z_j}\frac{dz'}{E(z')}.$$

Temporal-horizon curvature regularity ($0 < p \le 1/2$) places a strict upper bound on the junction epoch:

$$z_j \le 0.687 \quad (p = 0.447,\;\; C_{\rm tail} = 0.4646\;\;\text{at the reference junction } z_j = 0.60).$$

Furthermore, the deep-ambient well-nesting analysis (`step_65_acoustic_junction.py`) demonstrates that local potential wells embedded in a deep ambient background ($u_{\rm amb} > 0$) experience exponential source starvation $u_{\rm loc} \simeq u_{\rm loc}(0)\,e^{-0.45\,u_{\rm amb}}$, so that emitters at $z \simeq 1100$ ($u_{\rm emit} = \ln 1101 \simeq 7.00$) are $99.3\%$ ambient-dominated. Because the acoustic sound horizon $r_s$ in matter units is relational within the emitter environment and cosmological redshift depends only on the endpoint ratio $1+z = A_{\rm sp}(t_o)/A_{\rm sp}(t_e)$, the acoustic sector imposes no lower bound on $z_j$; regularity ($z_j \le 0.687$) is the sole binding constraint.

---

## 5. Cross-Domain Empirical Verification of the Selected Sector

### 5.1 Galactic Rotation Curves and the Baryonic Tully–Fisher Relation (`TEP-UCD`)

1. Radial Acceleration Relation (RAR, `TEP-UCD` `step_10_two_branch_sparc_rar.py`): Evaluated across 3,389 points in 175 SPARC galaxies using the regulated two-branch operator ($k = 16.03$ fixed externally by $a_0$, zero per-galaxy free parameters, $\varepsilon_{\rm reg} = 10^{-3}$), the two-branch sector reduces the log-residual RMS from $0.2812\text{ dex}$ (unit-term baseline) and $0.2749\text{ dex}$ (exponential interpolation) to $0.2185\text{ dex}$ at $\Upsilon_d = 0.5$ and $0.2004\text{ dex}$ at $\Upsilon_d = 0.35$ under action-derived scalar addition ($g_{\rm obs} = g_{\rm bar} + a_\phi$), and $0.1974\text{ dex}$ under the empirical quadrature comparator.
2. Outer-Edge Baryonic Tully–Fisher Relation (BTFR, `TEP-UCD` `step_13_two_branch_btfr.py`): Predicting outer rotation velocities across 135 SPARC galaxies directly from the baryonic mass profile via the exact two-branch flux law yields a median ratio $v_{\rm obs}/v_{\rm pred} = 0.910$ (scalar addition at $\Upsilon_d = 0.5$, RMS $0.0925\text{ dex}$) and $1.0109$ ($+0.0047\text{ dex}$ bias, RMS $0.0812\text{ dex}$) under quadrature.
3. Acceleration Onset Radius Law (`TEP-UCD` `step_12_onset_radius_law.py`): Across 54 galaxies spanning the transition $g_{\rm bar} = a_{\rm eff}$, the zero-parameter onset-radius law $R_{\rm onset} = \sqrt{GM_{\rm bar}/a_{\rm eff}}$ reproduces the measured intercept to $0.022\text{ dex}$ ($-4.445$ observed vs. $-4.467$ predicted, RMS $0.276\text{ dex}$).

### 5.2 Wide-Binary Orbital Dynamics and Environmental Ordering (`TEP-WB`)

1. Embedded Propagator and Solar-Circle Ambient (`step_56_embedded_two_branch.py`, `TEP-WB` `step_020`): Solving the nonlinear elliptic PDE $\nabla\cdot[P_{,X}(\lvert\nabla\phi\rvert)\nabla\phi] = \rho$ for a binary embedded in the self-consistent solar-circle scalar shear $u_0 = 0.163\,g_t$ (where $u_0[P_{,X}(u_0^2/2) + 2\beta_A^2] = g_{{\rm obs},\odot}/g_t = 0.631$, implying a Milky Way baryonic circular-velocity share $v_b \simeq 162\text{ km s}^{-1}$) and convolving with isotropic 3D deprojection ($\langle d\rangle \simeq 1.57\,s$) yields:
   - Two-branch prediction: $\alpha_{\rm pred} = 0.371$ ($0.44\sigma$ vs. observed $0.366 \pm 0.012$) and $R_s = 3{,}676\text{ AU}$ ($1.69\sigma$, joint $\chi^2 = 3.06$), improving to $\alpha_{\rm pred} = 0.347$ ($1.57\sigma$) and $R_s = 2{,}940\text{ AU}$ ($0.48\sigma$) when integrating the isotropic two-centre mutual force $\langle y\rangle = \frac{2}{3}y_\perp + \frac{1}{3}y_\parallel$ over per-star screening radii $r_{*,i} = \sqrt{Gm_i/g_t}$ (`step_59_two_center.py`, `step_75_offaxis_pair_3d.py`).
   - Baseline and exponential-interpolation exclusion: Both unit-term candidates predict $\alpha_{\rm pred} = 0.719\text{–}0.731$, excluded by the wide-binary plateau at $29\text{–}30\sigma$ ($\chi^2 > 860$).
2. Vertical Disk Stratification and In-Medium Planar Law (`step_67_inmedium_constraint_slice.py`, `TEP-WB` `step_020g`–`step_020o`):
   - Across five vertical height strata ($\lvert Z\rvert_{\rm med} = 22\text{ to }347\text{ pc}$), denser midplane environments exhibit stronger screening ($\alpha = 0.231 \to 0.446$ from midplane to halo).
   - In the selection-safe separation window ($s < 9\text{ kAU}$, unaffected by the `R_chance < 0.01` outer-tail truncation identified in `step_020k` and verified via continuous purity weighting $w = 1 - R_{\rm chance}$ across $N = 358{,}359$ pairs in `step_020n`), the covariate-hardened vertical gradient is positive across all $30/30$ distance, noise-floor, and mass decile slices ($t = +12.7$, bootstrap contrast $+6.9\sigma$).
   - The in-medium planar $k$-branch balance $P_{,X}u_{\rm amb} = J$ with local disk column source $J = 4\pi G\rho L/g_t$ derives the environmental scaling law:

$$u_{\rm amb}(\rho) = \sqrt{\frac{4\sqrt{2}\pi G\rho L}{k g_t}} \propto \rho^{1/2}.$$

   - Calibrating the disk coherence length once at the midplane ($L = 837\text{ pc}$) predicts the effective ambient across all five vertical strata to within $10.1\%$ ($u_{\rm pred} = 0.180, 0.161, 0.141, 0.117, 0.074$ vs. measured $0.18, 0.17, 0.15, 0.13, 0.08$). Jacobi tidal-radius analysis (`step_020o`, median $r_J = 279\text{ kAU}$, median $s/r_J = 0.023$ at $s = 5\text{–}9\text{ kAU}$) rules out differential Galactic tidal truncation by a factor of forty.

### 5.3 Solar-System PPN, Binary Pulsars, Standard Sirens, and Structure Growth

1. Solar-System PPN (`step_57_ppn_two_branch.py`, `step_80_shared_lapse_action_derivation.py`): The Sun's branch crossover sits at $r_c = (\sqrt{2}/k^3)^{1/2} r_* \simeq 78\text{ AU}$, placing all inner Solar-System probes deep in the quadratic $2\xi$ regime. At Cassini conjunction ($1.6\,R_\odot$), the local response is $y = 2.16 \times 10^{-8}$, yielding point-evaluated $\lvert\gamma_{\rm PPN}-1\rvert = 4.32 \times 10^{-8}$ (unstretched space) / $8.63 \times 10^{-8}$ (isotropic benchmark) and ray-integrated Shapiro/Doppler observables of $3.74 \times 10^{-7}$ and $9.01 \times 10^{-8}$, passing the Cassini bound ($2.3 \times 10^{-5}$) by margins of $60\times\text{ to }530\times$.
2. Binary-Pulsar Scalar Dipole Transmission and Shell-Impedance Screening (`step_69_dipole_transmission.py`): Separates (i) the finite-boundary $\ell = 1$ radial flux-state control for $P = X + X\lvert X\rvert/\Lambda_X^4$ with unscreened compact-body sensitivity charges prescribed at the orbital radius $r = a$ ($T_{\rm amp} = 0.0121$ and tightest sensitivity-bracket margin $0.529$ on PSR J1738+0333) from (ii) the selected two-branch sector $P_{,X} = k\sqrt{\xi} + 2\xi$ ($k = 16.03$) with shell-impedance screening across the compact-body-to-wave-zone profile ($T_{\rm amp}^{\rm imp} \simeq a/r_* = 9.2 \times 10^{-7}\text{ to }1.4 \times 10^{-5}$, $2.26 \times 10^{-6}$ on J1738+0333), under which PSR J1738+0333, J0348+0432, J0737$-$3039, and J0437$-$4715 clear observational $\dot P_b$ bounds across the full EOS sensitivity bracket by amplitude margins of $2.8 \times 10^3\text{ to }3.9 \times 10^4$ at $K_D = 0.5$ (tightest margin $2{,}817\times$ on J1738+0333 $s_{\rm hi}$; $996\times\text{ to }1.39 \times 10^4$ at $K_D = 4.0$), and the PSR J0337+1715 hierarchical $\mathcal R = y^3$ strong-equivalence-principle bound is satisfied at $\eta_{\rm pred} = 2.54 \times 10^{-14} \ll 2.6 \times 10^{-6}$ (margin $\sim 10^8\times$).
3. Standard-Siren Invariant and Supernova Luminosity Bound (`step_70_siren_invariant.py`, `step_71_intrinsic_mass_drift.py`, `step_73_sn_luminosity_bound.py`): Tracking the dimensionless inspiral invariant $I = G_{\rm loc}\mathcal{M}_c f / c^3$ confirms that a screened binary is locally standard ($I_{\rm det}/I_{\rm src} \to 1$) while preserving the transport distance ratio $\Xi(z) = d_L^{\rm GW}/d_L^{\rm EM} = 1/(1+z)$ and detector chirp-mass scaling $\mathcal{M}_{\rm det} \propto (1+z)^{3\epsilon - 1}$. Confronting the Chandrasekhar peak-luminosity scaling $\Delta m(z) = -7.5\epsilon\log_{10}(1+z)$ with the full 1,701-supernova Pantheon+ covariance matrix (`TEP-C0` `step_03_11` and `step_03_13`) excludes full interior tracking ($\epsilon = 1$) at $\Delta\chi^2 = 805.8$ and bounds $\epsilon < 0.150$ at $3\sigma$, preserving the inverse chirp-mass discriminant tested against the GWTC binary-black-hole catalog (`TEP-LVK`).
4. In-Well Convolved Cosmological Structure Growth (`step_62_growth_fsigma8.py`, `step_74_fs8_inwell_convolution.py`): Convolving the cosmic-web growth rate with the Press–Schechter collapsed halo fraction (within which $G_{\rm eff} \to G$) reduces the 7-point BOSS/eBOSS $f\sigma_8$ + DES/KiDS $S_8$ $\chi^2$ from $44.2$ (baseline) to $26.8\text{–}27.0$ ($f\sigma_8(0.61) = 0.483\text{–}0.485$), matching the $\Lambda$CDM floor ($\chi^2 = 26.0$).

### 5.4 Compact-Body Temporal Wells and Boundary Transmission (`step_53`, `step_58`, `step_61`, `step_76`)

1. Asymptotic Temporal-Well Regularity (`step_58_well_family.py`): Sweeping uniform-sphere compaction over eight decades at fixed charge shows that the central field depth grows double-logarithmically, $u_c(R) \propto \ln\ln(1/R)$, so the central lapse $N_{\rm min} = e^{-u_c}$ approaches zero asymptotically as an inverse power of $\ln(1/R)$ without reaching a singular halt at any finite radius, while $V_0$ supplies a regular, scale-invariant interior core in $g_{\mu\nu}$.
2. Elliptic and Hyperbolic Stability Certificate (`step_61_stability_certificate.py`): Across $\xi \in [10^{-12}, 10^6]$ on the static branch and out to $z = 100$ on the cosmological timelike branch, the two-branch operator satisfies $P_{,X} > 0$, $Z_\parallel > 0$, $Z_\perp > 0$, and subluminal propagation speed $c_s^2 = P_{,X}/Z_\parallel \in [1/3, 1/2]$.
3. Unified Radial Profiles and Static Boundary Response (`step_53_unified_profile.py`, `step_76_ambient_transmission.py`): Conservative radial finite-volume solutions confirm that Solar-System and stellar bodies reside on the flux-limited (leaking) branch ($R \ll R_c \simeq 6.5 \times 10^{13}\text{ m}$). On the Cassini-compatible screened quartic branch, static interior tangent responses are suppressed to $T_{\rm tan} \simeq 1.2 \times 10^{-8}$ at $10^8\,R_\oplus$, bounding the operative Solar-System gravitational drift to $\lvert\dot G_{\rm loc}/G_{\rm loc}\rvert_{\rm LLR} \lesssim 1.7 \times 10^{-18}\text{ yr}^{-1}$, more than four orders of magnitude below Lunar Laser Ranging limits.
