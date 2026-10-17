# Covariant Mathematical Foundations and Field-Equation Derivations

This document presents the formal differential-geometric and variational derivations underlying the Temporal Equivalence Principle (TEP). Each theorem and identity follows directly from the two-metric postulate and the covariant action principle.

---

## 1. Geometric Objects and Two-Metric Structure

The spacetime manifold $M$ is equipped with a gravitational metric $g_{\mu\nu}$, a scalar time field $\phi$, and a universal causal matter metric $\tilde g_{\mu\nu}$ related by the conformal–disformal map:

$$
\tilde g_{\mu\nu}=A^2(\phi)\,g_{\mu\nu}+B(\phi)\,\nabla_\mu\phi\nabla_\nu\phi,
\qquad
A(\phi)=\exp\!\left(\frac{\beta_A\phi}{M_{\rm Pl}}\right),
$$

with universal dimensionless conformal coupling $\beta_A = -1$. Under the unified cosmological and local decomposition (`step_80_shared_lapse_action_derivation.py`), the matter metric separates into the stationary baryon-sourced lapse factor $A_{\rm loc}(u_{\rm grav}(\mathbf x))$ along the unit timelike foliation normal $U_\mu = -\nabla_\mu\bar\phi/\sqrt{2\bar X}$ ($g^{\mu\nu}U_\mu U_\nu = -1$) and the secular spatial-ruler factor $A_{\rm sp}(u_{\rm drift}(t,\mathbf x))$ on the orthogonal spatial projector $\Pi_{\mu\nu} = g_{\mu\nu} + U_\mu U_\nu$:

$$
\tilde g_{\mu\nu} = -A_{\rm loc}^2(u_{\rm grav}(\mathbf x))\,U_\mu U_\nu + A_{\rm sp}^2(u_{\rm drift}(t,\mathbf x))\,\Pi_{\mu\nu} + B(\phi)\,\Pi_\mu^\alpha \Pi_\nu^\beta \nabla_\alpha\phi \nabla_\beta\phi.
$$

The fundamental scalar gradient invariants and coupling derivatives are defined as:

$$
\alpha(\phi)\equiv\frac{d\ln A}{d\phi}=\frac{\beta_A}{M_{\rm Pl}},
\qquad
\Sigma_\mu\equiv\nabla_\mu\ln A=\alpha(\phi)\nabla_\mu\phi,
\qquad
X\equiv g^{\mu\nu}\nabla_\mu\phi\nabla_\nu\phi.
$$

Metric invertibility and preservation of Lorentzian signature require:

$$
A>0,
\qquad
1+\frac{B}{A^2}X\neq 0,
\qquad
B\,X>-A^2.
$$

By the Sherman–Morrison identity, the exact contravariant inverse matter metric is:

$$
\tilde g^{\mu\nu}=A^{-2}\left[
g^{\mu\nu}
-\frac{(B/A^2)\,\nabla^\mu\phi\,\nabla^\nu\phi}{1+(B/A^2)X}
\right].
$$

On a homogeneous cosmological background ($\phi=\phi(t)$ with timelike gradient), the effective lapse satisfies:

$$
N^2(t)=A^2-B\dot\phi^2>0.
$$

---

## 2. Action Principle and Coupled Field Equations

In the gravitational frame, the general noncanonical scalar–tensor action takes the form:

$$
S=\int d^4x\sqrt{-g}\left[
\frac{M_{\rm Pl}^2}{2}R
+P(X,\phi)
\right]
+S_m[\psi_i,\tilde g_{\mu\nu}],
$$

where $P(X,\phi) = -V(\phi) + K(X,\phi)$ reduces in the canonical benchmark limit to $-\frac{1}{2}g^{\mu\nu}\nabla_\mu\phi\nabla_\nu\phi - V(\phi)$, and in the promoted two-branch kinetic completion to $P(X,\phi) = -V(\phi) + \frac{2k}{3\Lambda_X^2}|X|^{3/2} + \frac{1}{\Lambda_X^4}X^2$.

### 2.1 Gravitational Field Equations

Variation with respect to $g^{\mu\nu}$ yields the Einstein equations:

$$
G_{\mu\nu}=\frac{1}{M_{\rm Pl}^2}\left(T_{\mu\nu}^{(\phi)}+T_{\mu\nu}^{(m)}\right).
$$

For the canonical benchmark kinetic term, the scalar stress–energy tensor is:

$$
T_{\mu\nu}^{(\phi)}
=\nabla_\mu\phi\nabla_\nu\phi
-g_{\mu\nu}\left[\tfrac12(\nabla\phi)^2+V(\phi)\right],
$$

while for a general kinetic function $P(X,\phi)$, $T_{\mu\nu}^{(\phi)} = 2P_{,X}\nabla_\mu\phi\nabla_\nu\phi + P(X,\phi)g_{\mu\nu}$.

The gravitational-frame matter stress–energy tensor is the pullback of the physical matter-frame tensor $\tilde T_{\alpha\beta}^{(m)}$:

$$
\tilde T_{\alpha\beta}^{(m)}
=-\frac{2}{\sqrt{-\tilde g}}\frac{\delta S_m}{\delta\tilde g^{\alpha\beta}},
\qquad
T_{\mu\nu}^{(m)}
=\frac{\sqrt{-\tilde g}}{\sqrt{-g}}\,
\tilde T_{\alpha\beta}^{(m)}\,
\frac{\partial\tilde g^{\alpha\beta}}{\partial g^{\mu\nu}}.
$$

In the small-disformal regime ($|BX/A^2| \ll 1$), this pullback reduces to:

$$
T_{\mu\nu}^{(m)}=A^2(\phi)\,\tilde T_{\mu\nu}^{(m)}+O(B).
$$

### 2.2 Scalar Field Equation

Defining the contravariant matter stress densities:

$$
\tilde T^{\mu\nu}
=\frac{2}{\sqrt{-\tilde g}}\frac{\delta S_m}{\delta\tilde g_{\mu\nu}},
\qquad
\mathcal T^{\mu\nu}
=\frac{\sqrt{-\tilde g}}{\sqrt{-g}}\,\tilde T^{\mu\nu},
$$

and varying the metric map with respect to $\phi$:

$$
\delta_\phi\tilde g_{\mu\nu}
=2AA_{,\phi}g_{\mu\nu}\,\delta\phi
+B_{,\phi}\nabla_\mu\phi\nabla_\nu\phi\,\delta\phi
+B\bigl(\nabla_\mu\delta\phi\,\nabla_\nu\phi+\nabla_\mu\phi\,\nabla_\nu\delta\phi\bigr),
$$

integration by parts yields the covariant scalar field equation:

$$
\nabla_\mu\bigl(P_{,X}\nabla^\mu\phi\bigr)-V_{,\phi}=-\mathcal Q,
$$

where the effective matter source $\mathcal Q$ is given by:

$$
\mathcal Q
=AA_{,\phi}\,g_{\mu\nu}\mathcal T^{\mu\nu}
+\frac12 B_{,\phi}\,\mathcal T^{\mu\nu}\nabla_\mu\phi\nabla_\nu\phi
-\nabla_\mu\bigl(B\,\mathcal T^{\mu\nu}\nabla_\nu\phi\bigr).
$$

In the conformal limit ($B\to 0$), the source simplifies to:

$$
\nabla_\mu\bigl(P_{,X}\nabla^\mu\phi\bigr)-V_{,\phi}=-AA_{,\phi}\,g_{\mu\nu}\mathcal T^{\mu\nu}
=-\alpha(\phi)\,T_g,
$$

where $T_g = g_{\mu\nu}T_g^{\mu\nu}$. Nonrelativistic baryonic matter ($T_g \simeq -\rho$) acts as the primary source of $\phi$, whereas relativistic radiation ($T_g \simeq 0$) decouples from the conformal trace.

### 2.3 Covariant Conservation and Energy–Momentum Exchange

Diffeomorphism invariance of the matter action $S_m[\psi_i,\tilde g_{\mu\nu}]$ guarantees exact geodesic conservation in the physical matter frame:

$$
\tilde\nabla_\mu\tilde T^{\mu\nu}_{(m)}=0.
$$

In the gravitational frame, energy and momentum are exchanged between the matter and scalar sectors according to:

$$
\nabla_\mu T^{\mu\nu}_{(m)}=\mathcal Q\nabla^\nu\phi,
\qquad
\nabla_\mu T^{\mu\nu}_{(\phi)}=-\mathcal Q\nabla^\nu\phi,
$$

preserving exact conservation of the total stress–energy tensor $\nabla_\mu(T^{\mu\nu}_{(m)} + T^{\mu\nu}_{(\phi)}) = 0$. The resulting gravitational-frame fifth force $\mathcal Q\nabla_\mu\phi$ is environmentally suppressed in high-gradient regimes by the nonlinear kinetic stiffness $P_{,X} \gg 1$.

---

## 3. Local Lorentz Invariance and the Invariant Speed of Light

Physical clocks measure proper time along timelike worldlines of the matter metric:

$$
d\tau^2=-\frac{1}{c^2}\tilde g_{\mu\nu}dx^\mu dx^\nu.
$$

Within a local freely falling laboratory where $\partial_0\phi\ll|\nabla\phi|$ and $B$ is small:

$$
\frac{d\tilde\tau}{d\tau_g}=A(\phi)+O(B).
$$

Because a conformal transformation rescales local temporal intervals and spatial rods by the identical factor $A(\phi)$, the null cones of $g_{\mu\nu}$ and $\tilde g_{\mu\nu}$ coincide whenever $B=0$. Consequently, the invariance of the local speed of light $c$ in every freely falling laboratory is a mathematical theorem of universal coupling to $\tilde g_{\mu\nu}$.

In the quantum mechanical tangent limit, the Schrödinger evolution with respect to local proper time $\tau$ transforms into coordinate time $\tilde t$ as:

$$
i\hbar\frac{d|\psi\rangle}{d\tau}=\hat H|\psi\rangle
\quad\Longrightarrow\quad
i\hbar\frac{d|\psi\rangle}{d\tilde t}=A(\phi)\,\hat H|\psi\rangle,
$$

producing a universal clock-frequency shift across atomic transitions:

$$
\delta\ln\nu=\delta\ln A+K_\alpha\delta\ln\alpha+K_\mu\delta\ln\mu+K_q\delta\ln X_q+\cdots,
$$

with vanishing differential composition-dependent residuals ($K_i \delta\ln X_i = 0$) in the pure metric coupling limit.

---

## 4. Propagation Theorems

### 4.1 Theorem 1 — Conformal and Shared-Foliation Null-Cone Invariance

*Statement:* If $\tilde g_{\mu\nu}=A^2(\phi)g_{\mu\nu}$ (or more generally if the gravitational and electromagnetic sectors share the principal wave operator $\Box_{\rm TEP} = -(N_{\rm TEP}c)^{-2}\partial_t^2 + A_{\rm sp}^{-2}\gamma^{ij}\partial_i\partial_j$), then $g_{\mu\nu}k^\mu k^\nu=0$ if and only if $\tilde g_{\mu\nu}k^\mu k^\nu=0$.

*Proof:* For $A>0$, $\tilde g_{\mu\nu}k^\mu k^\nu=A^2 g_{\mu\nu}k^\mu k^\nu$. Because classical four-dimensional Maxwell electrodynamics is conformally invariant, photon rays follow $\tilde g$-null geodesics and coincide identically with $g$-null graviton rays ($c_{\rm GW}=c_\gamma$). Under the unified ADM/DHOST gravitational action on static space, the reciprocal weights $A_{\rm sp}^{+1}$ on the extrinsic-curvature kinetic term and $A_{\rm sp}^{-1}$ on the transverse-traceless 3-curvature preserve $c_{\rm GW}(t,\mathbf x) = c_\gamma(t,\mathbf x)$ across both local potential wells and cosmological epochs.

Consequently, the multimessenger arrival-time constraint from GW170817 ($|c_{\rm GW}-c_\gamma|/c \lesssim 10^{-15}$) bounds only the unshared disformal cone split $B(\phi)X/A^2$ along the propagation path, leaving the conformal and shared-cone sectors completely unconstrained.

### 4.2 Theorem 2 — Cancellation of First-Order One-Way Anisotropy in Static Gradients

*Statement:* For a stationary scalar field ($\partial_0\phi=0$), the linear-in-$\nabla\phi$ contribution to the forward and backward one-way light travel time along a common spatial path cancels identically.

*Proof:* Parameterizing the spatial path by arc length $s\in[0,L]$, the forward and reverse coordinate transit times to first order are:

$$
t_\to=\int_0^L\frac{ds}{c}\,A\bigl(\phi_0+s\,\partial_\parallel\phi\bigr),
\qquad
t_\leftarrow=\int_0^L\frac{ds}{c}\,A\bigl(\phi_0+(L-s)\,\partial_\parallel\phi\bigr).
$$

Expanding $A(\phi_0+s\partial_\parallel\phi)=A(\phi_0)+A'(\phi_0)\,s\,\partial_\parallel\phi+O((\nabla\phi)^2)$, the first-order integrals $\int_0^L s\,ds$ and $\int_0^L (L-s)\,ds$ are equal to $\frac{1}{2}L^2$. Any forward–backward asymmetry along the same link begins at $O((\nabla\phi)^2)$ or requires explicit time dependence ($\partial_t\phi \neq 0$) or non-exact closed-loop transport.

---

## 5. Synchronization Holonomy and the Disformal Connection

### 5.1 Threading Decomposition and Gauge Invariance

To define clock synchronization along a physical network with timelike four-velocity congruence $u^\mu$ ($\tilde g_{\mu\nu}u^\mu u^\nu = -1$), the matter metric is threaded in $(1+3)$ form:

$$
d\tilde s^2
=\tilde g_{00}\bigl(dt+\tilde\sigma_i\,dx^i\bigr)^2
+\tilde h_{ij}\,dx^i dx^j,
$$

where the coordinate synchronization connection one-form is given by the mixed threading component:

$$
\tilde\sigma_i=\frac{\tilde g_{0i}}{\tilde g_{00}}.
$$

Under a general coordinate re-synchronization of the clock worldlines,

$$
t'=t+\chi(x^i),
\qquad
dt=dt'-\partial_i\chi\,dx^i,
$$

the connection transforms by an exact differential:

$$
\tilde\sigma \longrightarrow \tilde\sigma - d\chi.
$$

Because the General Relativistic reference connection $\sigma_{\rm GR}$ (incorporating Sagnac, Lense–Thirring, Shapiro, and gravitational redshift terms on the same congruence) transforms by the identical exact differential $-d\chi$, the residual synchronization one-form

$$
\Delta\sigma\equiv\tilde\sigma-\sigma_{\rm GR}
$$

is strictly invariant under coordinate gauge transformations, and its closed-loop integral defines the gauge-invariant synchronization holonomy:

$$
H_{\rm resid}(C)=\oint_C\Delta\sigma=\iint_\Sigma d(\Delta\sigma),
\qquad
C=\partial\Sigma.
$$

### 5.2 Exactness of the Conformal Sector

In the purely conformal limit ($B=0$), the clock-rate contribution to the transport one-form is exact:

$$
\omega^{(A)}=d\ln A(\phi)
\quad\Longrightarrow\quad
\oint_C\omega^{(A)}=0
$$

around any contractible closed loop $C$ on which $\phi$ is smooth and single-valued. Thus, while the conformal factor $A(\phi)$ governs open-path gravitational and cosmological redshift as well as two-point clock-rate covariance $C_A(x,x')$, a nonvanishing closed-loop holonomy $H_{\rm resid}(C) \neq 0$ requires $d(\Delta\sigma) \neq 0$, generated by the disformal sector $B(\phi)\nabla_\mu\phi\nabla_\nu\phi$ or non-stationary kinematics.

### 5.3 Covariant and ADM Derivation of the Disformal Connection

Decomposing thematter metric into conformal and disformal parts,

$$
\tilde g_{\mu\nu}
=\underbrace{A^2(\phi) g_{\mu\nu}}_{\tilde g_{\mu\nu}^{(A)}}
+\underbrace{B(\phi)\nabla_\mu\phi\nabla_\nu\phi}_{\delta_B\tilde g_{\mu\nu}},
$$

and introducing the spatial projector $P_\mu{}^\nu=\delta_\mu{}^\nu+u_\mu u^\nu$ orthogonal to the clock congruence $u^\mu$, the scalar gradient decomposes as:

$$
\nabla_\mu\phi
=D_\mu\phi-(u\cdot\nabla\phi)\,u_\mu,
\qquad
D_\mu\phi\equiv P_\mu{}^\nu\nabla_\nu\phi.
$$

Projecting the disformal perturbation onto the mixed time–space directions of the congruence yields:

$$
P_\mu{}^\alpha u^\beta\,\delta_B\tilde g_{\alpha\beta}
=B\,(P_\mu{}^\alpha\nabla_\alpha\phi)\,(u^\beta\nabla_\beta\phi)
=B\,(u\cdot\nabla\phi)\,D_\mu\phi.
$$

In local rest-frame coordinates of the congruence ($u^\mu=(1,0,0,0)$, $\tilde g_{0i}^{(A)}=0$, $\tilde g_{00}^{(A)}=-A^2$), expanding $\tilde\sigma_i = (\tilde g_{0i}^{(A)} + \delta_B\tilde g_{0i})/(\tilde g_{00}^{(A)} + \delta_B\tilde g_{00})$ to leading order in $B$ gives:

$$
\delta\tilde\sigma_\mu
\simeq
-\frac{B(\phi)}{A^2(\phi)}(u\cdot\nabla\phi)\,P_\mu{}^\nu\nabla_\nu\phi
=
-\frac{B(\phi)}{A^2(\phi)}(u\cdot\nabla\phi)\,D_\mu\phi.
$$

In $(3+1)$ ADM variables of the gravitational metric,

$$
ds^2=-N^2 dt^2+h_{ij}(dx^i+N^i dt)(dx^j+N^j dt),
$$

with unit normal $n^\mu=N^{-1}(1,-N^i)$ and $\partial_0\phi = N(n\cdot\partial\phi) + N^i\partial_i\phi$, the matter-frame ADM components to first order in $B$ are:

$$
\tilde N\approx AN\left[1-\frac{B}{2A^2N^2}(n\cdot\partial\phi)^2\right],
\qquad
\tilde N_i\approx A^2 N_i + B(\partial_i\phi)\partial_0\phi.
$$

For hypersurface-orthogonal clocks ($u^\mu=n^\mu$) and negligible shift ($N_i\to 0$), $\tilde g_{00}^{(A)}\simeq -A^2 N^2$ and $\delta_B\tilde g_{0i}\simeq B N(n\cdot\partial\phi)\partial_i\phi$, yielding the coordinate-time connection:

$$
\delta\tilde\sigma_i
\simeq
-\frac{B(\phi)}{A^2(\phi) N}(\partial_i\phi)\,(n\cdot\partial\phi).
$$

Writing $\delta\tilde\sigma = f\,D\phi$ with $f \equiv -(B/A^2)(u\cdot\nabla\phi)$, the exterior derivative on a spatial slice is:

$$
d(\delta\tilde\sigma) = df\wedge D\phi + f\,d(D\phi).
$$

Whenever a satellite or terrestrial clock congruence moves with velocity $v^i$ across a spatial scalar gradient $\partial_i\phi$, the contraction $u\cdot\nabla\phi = u^j\partial_j\phi \neq 0$ even for a static field ($\partial_t\phi = 0$). Because $u^j\partial_j\phi$ varies along the orbit independently of $\phi$, $df\wedge D\phi \neq 0$, generating a nonvanishing closed-loop synchronization holonomy $H_{\rm resid}(C) = \iint_\Sigma d(\delta\tilde\sigma)$ across Earth–satellite triangles.

---

## 6. Kinetic Screening and Parametrized Post-Newtonian Limit

### 6.1 Environmental Screening Operator

The observable Temporal Shear $\Sigma_\mu^{\rm obs}$ is governed by the inverse kinetic stiffness of the scalar field:

$$
\Sigma_\mu^{\rm obs}=\mathcal{S}_\Sigma(\xi)\,\Sigma_\mu,
\qquad
\mathcal{S}_\Sigma(\xi) = \frac{1}{P_{,X}(\xi)},
\qquad
\xi \equiv \frac{|X|}{\Lambda_X^4}.
$$

In the parameter-free two-branch kinetic sector,

$$
P_{,X}(\xi) = k\sqrt{\xi} + 2\xi,
\qquad
k = \frac{4\sqrt{2}\,g_t}{a_0} \simeq 16.03,
\qquad
\Lambda_X^4 = M_{\rm Pl}^2 H_0^2,
$$

so $\mathcal{S}_\Sigma(\xi) \to (2\xi)^{-1} = (g_t/g)^2 \ll 1$ in high-gradient Solar-System environments ($\xi \gg \xi_* = (k/2)^2 \simeq 64$), and $\mathcal{S}_\Sigma(\xi) \to (k\sqrt{\xi})^{-1}$ in the low-acceleration galactic regime ($\xi \ll \xi_*$), recovering flat galactic rotation curves ($g_{\rm obs}^2 \simeq a_0 g_{\rm bar}$) with zero free parameters.

### 6.2 Parametrized Post-Newtonian (PPN) Slip

In the Damour–Esposito-Farèse (DEF) normalization $\alpha_0 = \sqrt{2}\,\beta_A$, a body with screened exterior scalar charge $\alpha_{\rm eff} = \mathcal{S}_\Sigma\,\alpha_0$ induces the weak-field Jordan-frame metric components:

$$
g_{00}^{\rm J}=-1+\frac{2GM}{r c^2}\left(1+\alpha_0\,\alpha_{\rm eff}\right),
\qquad
g_{rr}^{\rm J}=1+\frac{2GM}{r c^2}\left(1-\alpha_0\,\alpha_{\rm eff}\right)
\quad (\text{isotropic conformal benchmark}),
$$

or $g_{rr}^{\rm J} = 1 + 2GM/(rc^2)$ on the unstretched spatial metric ($A_{\rm sp} = 1$). Because the photon probe is unscreened while the solar source charge is screened by $\mathcal{S}_\Sigma$, the post-Newtonian deflection and Shapiro parameter $\gamma_{\rm PPN}$ is linear in $\mathcal{S}_\Sigma$:

$$
\gamma_{\rm PPN}-1 = -\frac{2\beta_A^2\,\mathcal{S}_\Sigma}{1+2\beta_A^2\,\mathcal{S}_\Sigma}
\quad (\text{unstretched spatial metric}),
\qquad
\gamma_{\rm PPN}-1 = -\frac{4\beta_A^2\,\mathcal{S}_\Sigma}{1+2\beta_A^2\,\mathcal{S}_\Sigma}
\quad (\text{isotropic conformal benchmark}).
$$

As $\mathcal{S}_\Sigma \to 0$ in the Solar System ($\mathcal{S}_\Sigma(1.6\,R_\odot) = 2.16 \times 10^{-8}$), $\gamma_{\rm PPN} \to 1$ and $\beta_{\rm PPN} \to 1$ to well within two orders of magnitude of the Cassini constraint $|\gamma_{\rm PPN}-1| < 2.3 \times 10^{-5}$, without requiring the fundamental coupling $\beta_A = -1$ to vanish.

---

## 7. Mount Wilson Cosmological Redshift Theorem and Frame Obstruction

### 7.1 Endpoint Redshift Identity

Consider a spatially homogeneous, static gravitational geometry $g_{ij}(\mathbf x)$ ($a_{\rm grav} \equiv 1$) on which the matter-frame spatial ruler factor $A_{\rm sp}(t)$ evolves secularly with cosmic time. Spatial homogeneity provides three translational Killing vectors, conserving the covariant spatial wavevector $q_i = A_{\rm sp}^2 \delta_{ij} k^j$ along any null ray. From the null condition $\tilde g_{\mu\nu}k^\mu k^\nu = 0$, the frequency measured by a comoving matter observer with four-velocity $\tilde u^\mu$ is:

$$
\omega = -\tilde g_{\mu\nu}\tilde u^\mu k^\nu = \frac{c\,q}{A_{\rm sp}(t)}.
$$

Conservation of $q$ between emission at $t_e$ and observation at $t_o$ yields the exact cosmological redshift identity:

$$
1+z = \frac{\omega_e}{\omega_o} = \frac{A_{\rm sp}(t_o)}{A_{\rm sp}(t_e)}.
$$

Thus, cosmological redshift arises on a static gravitational manifold as the ratio of matter-frame clock and ruler standards between emission and observation.

### 7.2 Two-Sector Frame Obstruction

Under an arbitrary conformal redefinition of the gravitational metric, $g'_{\mu\nu} = \Omega^2(\phi)g_{\mu\nu}$, the matter metric reads $\tilde g_{\mu\nu} = (A^2/\Omega^2)g'_{\mu\nu} + B\nabla_\mu\phi\nabla_\nu\phi$:

- Choosing $\Omega = A$ absorbs the conformal factor into $g'_{\mu\nu}$, rendering the gravitational sector time-dependent.
- Choosing $\Omega = 1$ preserves the static gravitational geometry $g_{ij}(\mathbf x)$ while retaining the evolving matter-frame calibration $A_{\rm sp}(t)$.

Because the spatially covariant gravitational action $S_{\rm grav}$ on static space carries constant Planck mass $M_{\rm Pl}^2$ and unit WKB tensor weight $\sqrt{A_{\rm sp}^{+1}\cdot A_{\rm sp}^{-1}} \equiv 1$ in conformal time $d\eta = dt/A_{\rm sp}(t)$, no single coordinate or conformal transformation can simultaneously map both the gravitational wave sector and the electromagnetic matter sector onto a minimally coupled expanding FLRW spacetime. The two sectors are physically distinguishable through the diffeomorphism-invariant standard-siren distance ratio:

$$
\Xi(z) \equiv \frac{d_L^{\rm GW}(z)}{d_L^{\rm EM}(z)} = \frac{1}{1+z} \neq 1.
$$

---

## 8. Temporal Horizon Regularity

Identifying the asymptotic past clock-rate scaling in conformal coordinate $\eta \to \infty$ with the power-law profile

$$
A_{\rm clock}(\eta) = C\,\eta^{-p},
$$

geodesic and curvature analysis (`TEP-TH`) establishes that for $0 < p \le 1/2$:

1. All polynomial curvature invariants of $\tilde g_{\mu\nu}$ remain regular and vanish as $A_{\rm clock} \to 0$;
2. Matter-frame timelike proper time $\int A_{\rm clock}\,d\eta$ diverges, ensuring timelike geodesic completeness;
3. Matter-frame null affine parameter diverges, ensuring null geodesic completeness.

Consequently, the apparent initial singularity ($a \to 0$) of standard FLRW cosmology is replaced in the canonical TEP architecture by an asymptotically regular past temporal horizon $\mathscr{T}^-$.

---

## 9. Effective-Field-Theory (Bellini–Sawicki) Mapping

On the conformal cosmological branch ($B=0$, $\alpha_T=0$), the linear perturbation equations in the effective field theory of dark energy (Bellini & Sawicki parametrization, `TEP-HC`) map onto the conformal rate $\alpha_A \equiv d\ln A / d\ln\tilde a = -d\ln A / d\ln(1+z)$ and the Einstein-frame kinetic stiffness $Z_E = P_{,X} + 2XP_{,XX}$ via $w_c \equiv Z_E - 6\beta_A^2$:

$$
\alpha_M = -2\alpha_A,
\qquad
\alpha_B = 2\alpha_A,
\qquad
\alpha_K = w_c\,\alpha_A^2
\quad (w_c = -5 \text{ for the canonical } Z_E = 1 \text{ image}),
\qquad
\alpha_T = 0.
$$

The resulting kinetic no-ghost discriminant is strictly nonnegative across all admissible branches ($Z_E > 0 \iff w_c > -6$):

$$
D = \alpha_K + \frac{3}{2}\alpha_B^2 = (w_c + 6)\,\alpha_A^2 = Z_E\beta_A^{-2}\alpha_A^2 \ge 0,
$$

with $D = \alpha_A^2 \ge 0$ on the canonical $Z_E = 1$ ($w_c = -5$) image and linear perturbation stability achieved on the kinetically screened noncanonical branch $w_c \gtrsim 7 \times 10^2$ (`TEP-HC`).

---

## 10. Summary of the Deductive Chain

The foundational mathematical structure of the corpus is summarized by the implications:

$$
\begin{aligned}
&\tilde g_{\mu\nu} = -A_{\rm loc}^2 U_\mu U_\nu + A_{\rm sp}^2 \Pi_{\mu\nu} + B\,\Pi_\mu^\alpha \Pi_\nu^\beta \nabla_\alpha\phi\nabla_\beta\phi
\\
&\Longrightarrow\quad
\tilde\nabla_\mu\tilde T^{\mu\nu}_{(m)}=0,\quad
\nabla_\mu(P_{,X}\nabla^\mu\phi)-V_{,\phi}=-\mathcal Q,\quad
c_{\rm loc} = c\;\text{invariant},\quad
c_{\rm GW} = c_\gamma
\\
&\Longrightarrow\quad
\oint_C d\ln A = 0,\quad
H_{\rm resid}(C) = \oint_C (\tilde\sigma - \sigma_{\rm GR}) = \iint_\Sigma d(\delta\tilde\sigma)
\\
&\Longrightarrow\quad
1+z = \frac{A_{\rm sp}(t_o)}{A_{\rm sp}(t_e)},\quad
\Xi(z) = \frac{d_L^{\rm GW}}{d_L^{\rm EM}} = \frac{1}{1+z}\;\text{on static space }\gamma_{ij}(\mathbf x)
\\
&\Longrightarrow\quad
\gamma_{\rm PPN} - 1 \propto -\beta_A^2\,\mathcal{S}_\Sigma(\xi) \to 0\;\text{for }\xi \gg \xi_*,\quad
g_{\rm obs}^2 \to a_0 g_{\rm bar}\;\text{for }\xi \ll \xi_*.
\end{aligned}
$$
