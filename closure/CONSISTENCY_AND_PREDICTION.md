# Mathematical Consistency and Observational Prediction Specification

This document specifies the mathematical scope and numerical verification of the core field-theoretic and observational calculations in the Temporal Equivalence Principle (TEP) pipeline. It establishes how the structural properties of the canonical scalar sector motivate the noncanonical two-branch $P(X)$ kinetic sector and the spatial-ruler-drift architecture.

## 1. Unified Metric Formulation and Action Identities

Paper 0 formulates the unified causal matter metric

$$\tilde g_{\mu\nu} = -A_{\rm loc}^2(u_{\rm grav}(\mathbf x))\,U_\mu U_\nu + A_{\rm sp}^2(u_{\rm drift}(t,\mathbf x))\,\Pi_{\mu\nu} + B(\phi)\,\Pi_\mu^\alpha \Pi_\nu^\beta \nabla_\alpha\phi \nabla_\beta\phi$$

together with the spatially covariant preferred-foliation (khronon/ADM) gravitational action on static space $\gamma_{ij}(\mathbf x)$. Symbolic verification (`step_80_shared_lapse_action_derivation.py`) confirms the weak-field and quadratic action identities, the exact equivalence at quadratic order between the local khronon acceleration invariant $a_i a^i = \gamma^{ij} D_i\ln N_{\rm GR}\,D_j\ln N_{\rm GR}$ and the longitudinal scalar 3-curvature, and the shared electromagnetic and gravitational-wave principal operator $\Box_{\rm TEP}$, which enforces $c_{\rm GW} = c_\gamma$ identically and eliminates WKB Hubble friction on static space ($\Xi(z) = d_L^{\rm GW}/d_L^{\rm EM} = 1/(1+z)$).

The companion action-contract audit (`step_81_action_contract_audit.py`) separates the homogeneous timelike principal coefficients from the spatial shell flux. At $\bar\xi=1$, $K_t=2k+6=38.06$ and the native homogeneous spatial coefficient is $P_{,X}(1)=18.03$. On the shell, put $p=q_s y_s$ and $\xi=p^2/2$. The flux law implies

$$P_{,X}(p^2/2)y_s=1,\qquad K_\perp=P_{,X}=\frac1{y_s},\qquad K_\parallel=\frac{d}{dp}[P_{,X}(p^2/2)p]=\frac2{y_s}+(q_s y_s)^2.$$

At the saved shell, $K_\perp=0.135556$ (multiplicity 2 on $S^2$) and $K_\parallel\simeq0.271254$ (multiplicity 1 along $\hat r$). These are the transverse and longitudinal eigenvalues of the spatial flux Hessian $\mathcal H_{ij}$. Their 3-D rotational trace average is $K_{\rm iso}=(2K_\perp+K_\parallel)/3=4/(3y_s)+(q_s y_s)^2/3=0.180789$ ($c_{p,\rm iso}^2=0.004750$), and their 3-D Frobenius/RMS principal stiffness is $K_{\rm rms}=\sqrt{(2K_\perp^2+K_\parallel^2)/3}=0.191772$ ($c_{p,\rm rms}^2=0.005039$), whose exact deep-IR ($q_s\to0$) limit is $\sqrt2/y_s=0.191705$. The BBN reduced calculation uses the explicitly rescaled effective prescription $K_{\rm sp}^{\rm eff}=\sqrt2 K_\perp=0.191705$, with $c_{p,\rm eff}^2=K_{\rm sp}^{\rm eff}/38.06=0.005037$, and evaluates the acoustic peak sensitivity across the full native Hessian interval $[K_\perp/K_t,K_\parallel/K_t]$ (`TEP-BBN` Step 20, maximum peak shift $\Delta\ell=0.12$ multipole). Its numerical results remain conditional on that prescription: the flux identity alone does not derive its normalization or its pairing with the separately evaluated timelike coefficient in a common quadratic action. No homogeneous negative direction is removed by this calculation.

In the matter frame, spatial curvature transforms through derivatives of $A_{\rm sp}$ rather than $A_{\rm loc}$; for spatially uniform $A_{\rm sp}$, $R_{\rm matter}=A_{\rm sp}^{-2}R_{\rm gravity}$. Step 20's $A=e^{-u}$ curvature calculation is a separate conformal benchmark.

## 2. Galactic Radial Acceleration Law

In the weak-field static regime, the radial acceleration is additive: $g_{\rm obs} = g_{\rm bar} + a_\phi$. In `TEP-UCD` (`step_10_two_branch_sparc_rar.py`), this relation is evaluated across 3,389 rotation-curve data points in 175 SPARC galaxies using the regulated two-branch kinetic operator with $k = 16.03$ (fixed by $4\sqrt{2}\,g_t/a_0$ using the adopted empirical acceleration scale $a_0$, with no per-galaxy adjustment of $k$), regulator $\varepsilon_{\rm reg} = 10^{-3}$, and standard disk mass-to-light ratio $\Upsilon_* = 0.5\,M_\odot/L_\odot$, yielding an RMS residual of $0.2185\text{ dex}$ (and $0.2004\text{ dex}$ under the $\Upsilon_* = 0.35\,M_\odot/L_\odot$ population-sensitivity check). Quadrature composition is evaluated separately as an empirical reference comparator rather than an orthogonal metric identity. Both evaluations have lower RMS residuals than the unit-term baseline ($0.281\text{ dex}$); these RMS comparisons alone are not a likelihood-based model-selection significance.

## 3. Initial-Slice Lapse Construction and Canonical-Scalar Obstruction

On the time-symmetric inhomogeneous Hamiltonian slice ($K_{ij} = 0$, $\Pi = 0$, pressureless matter at rest), the momentum constraint is identically satisfied and the Lichnerowicz equation determines the positive conformal factor $\psi > 0$. The lapse analysis (`step_21_lapse_solvability.py`) separates two distinct mathematical problems on this background:

1. Canonical static maximal-slicing obstruction: Imposing $\partial_t K = 0$ via a uniform spectral shift of the canonical scalar stress tensor requires $\Pi_*^2 < 0$ (negative kinetic energy density), which violates physical admissibility for a real canonical scalar field. This result proves that an inhomogeneous periodic slice cannot be held eternally static by a uniform canonical kinetic offset, providing the direct theoretical motivation for the noncanonical $P(X)$ kinetic sector and secular ruler drift $A_{\rm sp}(t)$.
2. Constructive evolving CMC initial slice: Retaining $\Pi = 0$ and the unmodified positive source $q = (\rho_m - 2V_{\rm eff})/(2M_{\rm Pl}^2) > 0$ on the initial $K_{ij} = 0$ slice yields the constant-mean-curvature (CMC) lapse equation

$$(-D_h^2 + q)\,N = C, \qquad \langle N \rangle_{\rm vol} = 1,$$

where $D_h^2 N = \psi^{-6}\nabla\cdot(\psi^2\nabla N)$ is the covariant conformal Laplacian. Because $q > 0$, the discrete face-flux operator satisfies the maximum principle, guaranteeing a unique, strictly positive lapse $N > 0$ with spatially uniform initial mean-curvature rate $\partial_t K = C > 0$. Both $16^3$ and $48^3$ lattice realizations converge with $N > 0$ and relative lapse residuals below $10^{-10}$.
3. Preferred-foliation khronon/cuscuton static maximal-slicing closure: On the preferred cosmic foliation with noncanonical drift term $P_{\rm drift}(X) = -\mu_{\rm kh}^2\sqrt{2X}$, the homogeneous drift contributes zero Hamiltonian energy density ($2X P_{{\rm drift},X} - P_{\rm drift} \equiv 0$, leaving the Lichnerowicz solution $\psi$ and $\lambda_*$ unchanged) while supplying a negative spatial stress trace $S_{\rm drift} = -3\mu_{\rm kh}^2\Pi < 0$ for real timelike drift $\Pi > 0$. Setting $\mu_{\rm kh}^2\Pi = \tfrac23 M_{\rm Pl}^2\mu_0 > 0$ solves the exact unshifted static maximal-slicing equation $(-D_h^2 + q_{\rm kh})N = 0$ ($\partial_t K = 0$) with strictly positive lapse $N(\mathbf x) \in [0.9990, 1.0005] > 0$ ($N=16$; $[0.9992, 1.0004]$ at $N=48$) and relative residual below $10^{-10}$.

The noncanonical stress-shift construction checks the lapse equation on one slice. Its scalar drift equation and preservation under time evolution are not solved by this spectral shift; an eternal static cosmology is not established by that calculation.

The companion maximal-slicing diagnostic (`step_21_congruence_preservation.py`) enforces $\Pi^2 \ge 0$ and $N > 0$ as strict admissibility criteria.

## 4. Linearized Fluctuation Spectra and Stability Bounds

Fluctuation stability around the solved inhomogeneous landscape is evaluated across three operator levels in `step_79_inhomogeneous_eigenmodes.py` and verified independently in `step_82_stability_lower_bound.py`:

1. Weighted zero-mean coupled spectrum: The background Lichnerowicz solver uses the negative-semidefinite Laplacian $\Delta$, whereas the second-variation fluctuation operator uses the positive-semidefinite operator $-\Delta$. Restricted to the metric-weighted zero-mean subspace $\int \sqrt{h}\,\delta u\,d^3x = \int \psi_0^6\,\delta u\,d^3x = 0$, the reduced Euclidean scalar operator ($\lambda_1 = +38.81$), the canonical Lichnerowicz-coupled Einstein–scalar operator ($\lambda_1 = +38.96$ at $N=24$, $+39.06$ at $N=32$), and the two-branch Lichnerowicz-coupled operator ($\lambda_1 = +555.35$ at $N=24$, $+551.66$ at $N=32$) are all strictly positive.
2. Independent discrete lower bounds: Using the periodic lattice Poincaré constant, minimum face stiffness, and a norm bound on the attractive gravitational Schur complement, `step_82_stability_lower_bound.py` establishes strictly positive conservative lower bounds on the ambient two-branch operator at $N = 24, 32,$ and $40$ ($+377.40, +378.34,$ and $+378.81$).
3. Homogeneous mode, complete matter-coupling second variation, twisted-boundary Bloch operator, and jointly on-shell quartic coupled stability: In the bare reduced canonical spectrum (`solve_spectrum`), the single unprojected negative eigenvalue ($\lambda_0 \approx -0.4389$, spatial variation $2.19 \times 10^{-5} \ll 1$) arises from differentiating only the bare reconstructed potential $V_{\rm rec}''(u) = -\tfrac12\bar\rho_* e^{-u} < 0$ (`rhob_half` $= \bar\rho_*/2 = 1.0$ from `step_12_cosmology_solution.py`) while omitting the second derivative of the conformal matter coupling $\rho_*(x)A(u) = \rho_*(x)e^{-u}$ in the total effective potential $\mathcal V_{\rm eff}(u,x) = V_{\rm rec}(u) + \rho_*(x)e^{-u}$. When the complete second variation $\mathcal V_{\rm eff}''(u,x) = V_{\rm rec}''(u) + \rho_*(x)e^{-u} = (\rho_*(x) - \tfrac12\bar\rho_*)e^{-u}$ is retained with the consistent cell-mean conserved density $\langle\rho_*(x)\rangle = \bar\rho_* = 2\,\text{rhob\_half} = 2.0$ (`completed_canonical_reconstructed_coupled_audit` in Step 79), the cell-mean effective mass squared flips from $-0.4389$ to $+0.4102 > 0$, and both the unprojected Euclidean operator ($\lambda_0 = +0.4087$ at $N=16$, $+0.4093$ at $N=24$) and the full unprojected Lichnerowicz-coupled two-branch operator ($\lambda_0 = +0.4123$ at $N=16$, $+0.4127$ at $N=24$; $\lambda_1 = +367.02$ and $+350.97$) are strictly positive without projecting out the homogeneous mode. Separately, for the variable-coefficient Hermitian Bloch operator

$$A_\theta = \sum_{i=1}^3 D_{\theta,i}^\dagger W_{{\rm face},i} D_{\theta,i} + \mathrm{diag}(m_{\rm eff}^2(x))$$

evaluated on the quartic benchmark potential, positive face stiffness $W_{\rm face} > 0$ and positive effective mass squared imply $A_\theta \ge \min(m_{\rm eff}^2) > 0$ for every Bloch wavevector $\theta \in [-\pi,\pi]^3$, confirmed numerically across sampled supercell twists. The additional quartic benchmark (`resolved_quartic_two_branch_coupled_audit`) re-solves the Hamiltonian geometry while retaining a prescribed scalar profile and uniform integrability offset; its scalar perturbations are unprojected and its metric response is constrained by $q_{\rm vol}^T\delta\psi=0$ using an exact bordered inverse. Furthermore, the companion jointly on-shell audit (`joint_on_shell_quartic_two_branch_coupled_audit` in Step 79) solves both the nonlinear two-branch scalar Euler–Lagrange equation $R_u[u_*,\psi_*]=-\nabla\cdot[\psi_*^2 P_{,X}(\xi)\nabla u_*]+\psi_*^6(\tilde\lambda u_*^3-\rho_m e^{-u_*})=0$ and the nonlinear Lichnerowicz Hamiltonian constraint $R_\psi[u_*,\psi_*,\lambda_*]=0$ simultaneously to $<3.1\times10^{-12}$ and evaluates the exact on-shell Hessian $J_u[u_*,\psi_*]$ coupled to the attractive metric Schur complement $-B_q L_\psi^{-1} B_q$ ($q_{\rm vol}^T\delta\psi=0$) across the full unprojected scalar Hilbert space, yielding strictly positive eigenvalues ($\lambda_0=+0.5624$ at $N=16$, $+0.5633$ at $N=24$; $\lambda_1=+378.95$ and $+381.67$) and an all-twist coupled lower bound $+0.4662$ ($N=16$) and $+0.4668>0$ ($N=24$) using the half-wave lattice gap $4\sin^2(\pi/(2N))/\Delta x^2$ adjusted for the angle between $q_{\rm vol}$ and the constant vector in periodic Bloch gauge. These constrained finite-grid operators do not constitute a full-action continuum stability theorem.

## 5. Solar-System Ray Propagation

In `step_80_shared_lapse_action_derivation.py`, the two-branch solar scalar profile is integrated from each impact radius to the outer matching boundary via reverse cumulative trapezoidal quadrature (`radial_tail_integral`), validated against exact analytic integrals in `tests/test_consistency_repairs.py`. The resulting ray-integrated Cassini observables are $|\gamma_{\rm PPN} - 1|_{\rm Shapiro} = 3.74 \times 10^{-7}$ (between $1.6\,R_\odot$ and $10\,R_\odot$) and $|\gamma_{\rm PPN} - 1|_{\rm Doppler} = 9.01 \times 10^{-8}$, satisfying the Cassini bound $2.3 \times 10^{-5}$ by margins exceeding sixty-fold.

## 6. Verification Commands

From the corpus root, the complete verification suite is executed via:

```bash
OPENBLAS_NUM_THREADS=1 python3 -m pytest TEP/tests -q
OPENBLAS_NUM_THREADS=1 python3 TEP/scripts/steps/step_21_lapse_solvability.py
OPENBLAS_NUM_THREADS=1 python3 TEP/scripts/steps/step_79_inhomogeneous_eigenmodes.py
OPENBLAS_NUM_THREADS=1 python3 TEP/scripts/steps/step_82_stability_lower_bound.py
python3 TEP/scripts/steps/step_80_shared_lapse_action_derivation.py
python3 TEP/scripts/steps/step_83_consistency_validation.py
```

## 7. Controlled pressure sensitivity

HC Step 26 replaces the unsupported exact $w=0$ dark-energy-fluid substitution with a supported sequence $w=-10^{-2},\ldots,-10^{-7}$. It tests both a zero-sound-speed control against CDM and convergence of the prescribed finite-pressure sequence. The result is a phenomenological near-dust sensitivity calculation, not a derived replacement for the TEP perturbation sector. Current numerical comparisons and tolerances are in `TEP-HC/results/step_26_near_dust_pressure_limit.json`.
