These are the **closed derivations** in the corpus: statements that follow from the axioms plus differential geometry, not from a later \(\kappa_X\) fit. Phenomenological identifications that are *not* derived from the action are marked at the end.

---

## 0. Objects

One manifold, two metrics, one scalar:

\[
\tilde g_{\mu\nu}=A^2(\phi)\,g_{\mu\nu}+B(\phi)\,\nabla_\mu\phi\nabla_\nu\phi,
\qquad
A(\phi)=\exp\!\left(\frac{\beta_A\phi}{M_{\rm Pl}}\right).
\]

Definitions used everywhere:

\[
\alpha(\phi)\equiv\frac{d\ln A}{d\phi}=\frac{\beta_A}{M_{\rm Pl}},
\qquad
\Sigma_\mu\equiv\nabla_\mu\ln A=\alpha(\phi)\nabla_\mu\phi,
\qquad
X\equiv g^{\mu\nu}\nabla_\mu\phi\nabla_\nu\phi.
\]

Invertibility and Lorentzian signature:

\[
A>0,
\qquad
1+\frac{B}{A^2}X\neq 0,
\qquad
B\,X>-A^2.
\]

Inverse (exact):

\[
\tilde g^{\mu\nu}=A^{-2}\left[
g^{\mu\nu}
-\frac{(B/A^2)\,\partial^\mu\phi\,\partial^\nu\phi}{1+(B/A^2)X}
\right].
\]

On the homogeneous cosmological branch (\(\phi=\phi(t)\), mostly timelike \(\nabla\phi\)):

\[
N^2(t)=A^2-B\dot\phi^2>0.
\]

---

## 1. Action and field equations

Einstein-frame action:

\[
S=\int d^4x\sqrt{-g}\left[
\frac{M_{\rm Pl}^2}{2}R
-\frac12(\nabla\phi)^2
-V(\phi)
\right]
+S_m[\psi_i,\tilde g_{\mu\nu}].
\]

### Einstein equation

\[
G_{\mu\nu}=\frac{1}{M_{\rm Pl}^2}\left(T_{\mu\nu}^{(\phi)}+T_{\mu\nu}^{(m)}\right),
\]

\[
T_{\mu\nu}^{(\phi)}
=\nabla_\mu\phi\nabla_\nu\phi
-g_{\mu\nu}\left[\tfrac12(\nabla\phi)^2+V(\phi)\right].
\]

Matter stress in the Einstein frame is the pullback of the matter-frame tensor

\[
\tilde T_{\alpha\beta}^{(m)}
=-\frac{2}{\sqrt{-\tilde g}}\frac{\delta S_m}{\delta\tilde g^{\alpha\beta}},
\qquad
T_{\mu\nu}^{(m)}
=\frac{\sqrt{-\tilde g}}{\sqrt{-g}}\,
\tilde T_{\alpha\beta}^{(m)}\,
\frac{\partial\tilde g^{\alpha\beta}}{\partial g^{\mu\nu}}.
\]

Small-\(B\) reduction:

\[
T_{\mu\nu}^{(m)}=A^2(\phi)\,\tilde T_{\mu\nu}^{(m)}+O(B).
\]

### Scalar equation

Define

\[
\tilde T^{\mu\nu}
=\frac{2}{\sqrt{-\tilde g}}\frac{\delta S_m}{\delta\tilde g_{\mu\nu}},
\qquad
\mathcal T^{\mu\nu}
=\frac{\sqrt{-\tilde g}}{\sqrt{-g}}\,\tilde T^{\mu\nu}.
\]

Variation of the map:

\[
\delta_\phi\tilde g_{\mu\nu}
=2AA_{,\phi}g_{\mu\nu}\,\delta\phi
+B_{,\phi}\nabla_\mu\phi\nabla_\nu\phi\,\delta\phi
+B\bigl(\nabla_\mu\delta\phi\,\nabla_\nu\phi+\nabla_\mu\phi\,\nabla_\nu\delta\phi\bigr).
\]

After integration by parts:

\[
\Box\phi-V_{,\phi}=-\mathcal Q,
\]

\[
\mathcal Q
=AA_{,\phi}\,g_{\mu\nu}\mathcal T^{\mu\nu}
+\frac12 B_{,\phi}\,\mathcal T^{\mu\nu}\nabla_\mu\phi\nabla_\nu\phi
-\nabla_\mu\bigl(B\,\mathcal T^{\mu\nu}\nabla_\nu\phi\bigr).
\]

Conformal limit \(B\to 0\):

\[
\Box\phi-V_{,\phi}=-AA_{,\phi}\,g_{\mu\nu}\mathcal T^{\mu\nu}
=-\alpha(\phi)\,g_{\mu\nu}\mathcal T^{\mu\nu}.
\]

Nonrelativistic matter: \(T\simeq-\rho\) sources \(\phi\). Radiation: \(T\simeq 0\) sources the conformal sector only weakly.

### Conservation

Diffeomorphism invariance of \(S_m[\tilde g]\) gives matter-frame geodesic motion:

\[
\tilde\nabla_\mu\tilde T^{\mu\nu}_{(m)}=0.
\]

Einstein-frame exchange:

\[
\nabla_\mu T^{\mu\nu}_{(m)}=\mathcal Q\nabla^\nu\phi,
\qquad
\nabla_\mu T^{\mu\nu}_{(\phi)}=-\mathcal Q\nabla^\nu\phi,
\]

so the total Einstein-frame tensor is conserved. The Einstein-frame fifth force on matter is \(\mathcal Q\nabla\phi\). That is not optional. Screening must hide it locally.

---

## 2. Local Lorentz invariance and emergent \(c\)

Clocks measure

\[
d\tau^2=-\tilde g_{\mu\nu}dx^\mu dx^\nu/c^2.
\]

In a local lab with \(\partial_0\phi\ll|\nabla\phi|\) and small \(B\):

\[
\frac{d\tilde\tau}{d\tau_g}=A(\phi)+O(B).
\]

Conformal rescaling multiplies clocks and rods together. Null cones of \(g\) and \(\tilde g\) coincide when \(B=0\). Local \(c\) is therefore a **theorem** of A2, not an extra postulate.

Quantum clocks in the same limit:

\[
i\hbar\frac{d|\psi\rangle}{d\tau}=\hat H|\psi\rangle
\quad\Longrightarrow\quad
i\hbar\frac{d|\psi\rangle}{d\tilde t}=A(\phi)\,\hat H|\psi\rangle,
\]

\[
\delta\ln\nu=\delta\ln A+K_\alpha\delta\ln\alpha+K_\mu\delta\ln\mu+K_q\delta\ln X_q+\cdots.
\]

---

## 3. Theorem 1 — conformal null-cone invariance

**Statement.** If \(\tilde g_{\mu\nu}=A^2(\phi)g_{\mu\nu}\), then \(g_{\mu\nu}k^\mu k^\nu=0\) iff \(\tilde g_{\mu\nu}k^\mu k^\nu=0\).

**Proof.** \(\tilde g_{\mu\nu}k^\mu k^\nu=A^2 g_{\mu\nu}k^\mu k^\nu\) and \(A>0\). Maxwell in 4D is conformally invariant, so photon trajectories stay \(\tilde g\)-null and therefore \(g\)-null. Gravitons propagate on \(g\). Hence \(c_g=c_\gamma\) identically in the conformal subclass.

GW170817 then bounds the *disformal* combination \(B(\phi)X\) along the observed path. It does not force \(B\equiv 0\) in every regime, and it does not constrain common-mode conformal clock structure.

---

## 4. Theorem 2 — no first-order one-way anisotropy for static \(\phi\)

**Statement.** If \(\partial_0\phi=0\), the linear-in-\(\nabla\phi\) piece of the forward/back one-way light time along the same path cancels.

**Proof.** Parameterize the path by \(s\in[0,L]\). To first order,

\[
t_\to=\int_0^L\frac{ds}{c}\,A\bigl(\phi_0+s\,\partial_\parallel\phi\bigr),
\qquad
t_\leftarrow=\int_0^L\frac{ds}{c}\,A\bigl(\phi_0+(L-s)\,\partial_\parallel\phi\bigr).
\]

Expand \(A(\phi_0+s\partial_\parallel\phi)=A(\phi_0)+A'(\phi_0)\,s\,\partial_\parallel\phi+O((\nabla\phi)^2)\). The linear integrals are equal. The difference starts at \(O((\nabla\phi)^2)\) or requires \(\partial_t\phi\) or kinematics.

Consequence: two-way and reciprocity-even tests can be null while a *loop* observable remains open. Large first-order “variable-\(c\)” claims are excluded.

---

## 5. Synchronization holonomy (the invariant)

A one-way time is not \(c=dx/d\tau_\gamma\) (\(d\tau_\gamma=0\)). It is an elapsed clock interval after a synchronization convention.

Let \(\tilde\sigma\) be the matter-frame synchronization one-form of a specified clock congruence \(u^\mu\), and \(\sigma_{\rm GR}\) the same object computed in GR (Sagnac, Lense–Thirring, Shapiro, redshift, station motion, clock-scale, frame).

Residual connection and holonomy:

\[
\Delta\sigma\equiv\tilde\sigma-\sigma_{\rm GR},
\qquad
H_{\rm resid}(C)=\oint_C\Delta\sigma=\iint_\Sigma d(\Delta\sigma),\quad C=\partial\Sigma.
\]

**Gauge invariance.** A re-gauging \(t\to t+\chi(x^i)\) adds the same exact form to \(\tilde\sigma\) and \(\sigma_{\rm GR}\). Therefore \(\Delta\sigma\) is invariant, and \(\oint d\chi=0\) for single-valued \(\chi\).

**Conformal exactness.**

\[
\omega^{(A)}=d\ln A
\quad\Longrightarrow\quad
\oint_C\omega^{(A)}=0
\]

on a contractible region where \(A\) is smooth and single-valued. Pure \(A(\phi)\) changes rates, open-path redshift, and covariance \(C_A\). It cannot source leading-order \(H_{\rm resid}\).

Nonzero leading holonomy requires \(d(\Delta\sigma)\neq 0\): disformal \(B\neq 0\), non-metricity, or another non-exact term.

### Disformal connection

Clock congruence: \(u^\mu u_\mu=-1\), projector \(P_\mu{}^\nu=\delta_\mu{}^\nu+u_\mu u^\nu\). Leading disformal piece (Paper 0, App. A3):

\[
\delta\tilde\sigma_\mu\simeq-\frac{B}{A^2}(u\cdot\nabla\phi)\,P_\mu{}^\nu\nabla_\nu\phi.
\]

Hypersurface-orthogonal case \(u=n\):

\[
\delta\tilde\sigma_i\simeq-\frac{B}{A^2 N}(\partial_i\phi)\,(n\cdot\partial\phi).
\]

ADM expansion used to get that:

\[
\tilde N\approx AN\left[1-\frac{B}{2A^2N^2}(n\cdot\partial\phi)^2\right],
\qquad
\tilde N_i\approx A\left[N_i+\frac{B}{A^2}(\partial_i\phi)(n\cdot\partial\phi)\right],
\]

\[
\tilde\sigma_i=\frac{\tilde g_{0i}}{\tilde g_{00}}\simeq-\frac{\tilde N_i}{\tilde N^2}.
\]

\(d(\delta\tilde\sigma)\neq 0\) when the prefactor varies independently around the loop (time dependence, anisotropic boundaries, inhomogeneous \(\phi\), lapse/shift, spatially varying \(B\)-response).

If screening drives \(\nabla\phi\to 0\) everywhere on the loop, holonomy vanishes with the shear. Conformal screening and disformal screening are *not* forced to be the same function of \(\mathcal{E}\).

---

## 6. Screening, PPN, response

Canonical dictionary:

\[
\Sigma_\mu=\nabla_\mu\ln A,
\qquad
C_A(x,x')=\langle\delta\ln A(x)\,\delta\ln A(x')\rangle,
\]

\(\lambda_T\) = measured decay length of \(C_A\), not an algebraic combination of local fields.

Observable shear:

\[
\Sigma_\mu^{\rm obs}=\mathcal{S}_\Sigma(\mathcal{E})\,\Sigma_\mu,
\qquad
\mathcal{S}_\Sigma\to 0
\text{ in screened regimes.}
\]

\(\mathcal{E}\) includes \(\rho\), \(\Phi/c^2\), \(\nabla\rho\), \(\nabla\Phi\), proximity, coherence volume, boundary geometry. Domain quantities \(\rho_T\), \(R_T(M)\), \(S_\oplus(r)\), compactness, \(\lambda_T\) are **projections of \(\mathcal{E}\)**, not extra mechanisms.

Channel response (definition, not a solution of \(\Box\phi\)):

\[
\Delta O_X
=\kappa_X\cdot\mathcal{S}_X(\mathcal{E})\cdot
\mathcal{F}_X[\Delta\ln A,\Sigma_\mu,C_A;\Phi,\rho,z].
\]

\(\kappa_X\) is not \(\beta_A\) and not a PPN parameter.

### PPN (DEF)

Unscreened:

\[
\gamma_{\rm PPN}-1\approx-\frac{2\alpha_0^2}{1+\alpha_0^2}\approx-2\alpha_0^2,
\qquad
\alpha_0=\alpha(\phi_\infty).
\]

Cassini \(|\gamma-1|<2.3\times 10^{-5}\) \(\Rightarrow\) \(\alpha_0\lesssim 3.4\times 10^{-3}\) *if unscreened*.

Screened exterior charge:

\[
\alpha_{\rm eff}=\mathcal{S}_\Sigma(\mathcal{E})\,\alpha_0.
\]

Jordan-frame static spherical metric:

\[
g_{00}^{\rm J}=-1+\frac{2GM}{r}\left(1+\alpha_0\,\alpha_{\rm eff}\right),
\qquad
g_{rr}^{\rm J}=1+\frac{2GM}{r}\left(1-\alpha_0\,\alpha_{\rm eff}\right).
\]

\[
\gamma_{\rm PPN}-1=-\frac{2\,\alpha_0\,\alpha_{\rm eff}}{1+\alpha_0\,\alpha_{\rm eff}}
=-\frac{4\beta_A^2\,S_\Sigma}{1+2\beta_A^2\,S_\Sigma},
\qquad
\beta_{\rm PPN}-1=\frac{\alpha_0\,\alpha_{\rm eff}}{2\left(1+\alpha_0\,\alpha_{\rm eff}\right)^2}\,\frac{d\alpha_{\rm eff}}{d\phi_0},
\]

with $\alpha_{\rm eff}=S_\Sigma\alpha_0$ and the DEF normalization $\alpha_0=\sqrt{2}\,\beta_A$. The $\gamma$ deviation is linear in the screened source charge because the photon probe is unscreened.

\(\mathcal{S}_\Sigma\to 0\) \(\Rightarrow\) \(\gamma_{\rm PPN}=\beta_{\rm PPN}=1\) independently of cosmological \(\alpha_0\). Screening does *not* set \(d\ln A/d\phi=0\); it suppresses the exterior charge.

---

## 7. Mount Wilson Equivalence Theorem (cosmological redshift)

**Isochrony Axiom (what is being dropped).** After GR gravitational and kinematic clock effects from a *single* metric, no independent field may rescale matter proper time across epochs.

**Ansatz.** Homogeneous, spatially flat, static gravitational geometry \(g_{ij}\) (\(a=1\)), \(\phi=\phi(t)\), universal coupling to \(\tilde g\), \(A^2-B\dot\phi^2>0\).

Effective lapse \(N^2=A^2-B\dot\phi^2\). Spatial homogeneity \(\Rightarrow\) three translational Killing vectors. The covariant spatial wavevector

\[
q_i=A^2\delta_{ij}k^j
\]

is conserved along the photon. Null condition:

\[
N^2(k^0)^2=q^2/A^2.
\]

Measured frequency:

\[
\omega=-\tilde g_{\mu\nu}\tilde u^\mu k^\nu=Nk^0=\frac{q}{A}.
\]

\(q\) conserved \(\Rightarrow\)

\[
1+z=\frac{A_0}{A_{\rm em}}.
\]

\(B\) enters \(N\) and cancels from the homogeneous *endpoint* ratio. Holonomy is not this ratio; it needs inhomogeneous non-exact transport.

On this branch the redshift observable is degenerate with FLRW \(1+z=a_0/a_{\rm em}\). That is the theorem. It is not a proof that \(g_{ij}\) is static.

### Frame obstruction (why this is not “just a conformal rewrite”)

Let \(g'_{\mu\nu}=\Omega^2(\phi)g_{\mu\nu}\). Then

\[
\tilde g_{\mu\nu}=(A^2/\Omega^2)\,g'_{\mu\nu}+B\nabla_\mu\phi\nabla_\nu\phi.
\]

- \(\Omega=A\) removes the conformal factor from matter but makes \(g'\) time-dependent.
- Constant \(\Omega\) keeps \(g\) static but leaves matter \(A(\phi)\)-dependent.

No single conformal frame both keeps gravitational geometry static and strips temporal calibration from matter. Inter-sector observables (sirens: tensors on \(g\), clocks on \(\tilde g\)) are the degeneracy breakers.

---

## 8. Temporal horizon regularity

Identify the observational clock map with the redshift definition:

\[
A_{\rm clock}(z)=(1+z)^{-1}.
\]

Near the past boundary write, in the temporal-horizon conformal coordinate \(\eta\),

\[
A_{\rm clock}(\eta)=C\,\eta^{-p}.
\]

**Proposition (TH).** For \(0<p\le 1/2\):

- polynomial curvature invariants vanish as \(A_{\rm clock}\to 0\);
- timelike proper time diverges;
- null geodesics have divergent affine parameter;
- the Strong Energy Condition is violated on this branch (needed for regularity).

\(A_{\rm clock}(z)=(1+z)^{-1}\) is the redshift definition. \(0<p\le 1/2\) is an independent regularity condition on the boundary. The apparent \(a\to 0\) is rewritten as \(A_{\rm clock}\to 0\), a temporal horizon \(\mathscr{T}^-\), not a curvature singularity *in this branch*.

This does not by itself derive recombination chemistry or the CMB blackbody. Those are delegated to the BBN/opacity papers and are not part of this theorem.

---

## 9. Pure-conformal EFT closure (HC)

Bellini–Sawicki functions, production branch \(B=0\), \(\alpha_T=0\):

\[
\alpha_A\equiv\frac{d\ln A}{d\ln\tilde a}=-\frac{d\ln A}{d\ln(1+z)},
\]

\[
\alpha_M=\frac{d\ln M_{\rm eff}^2}{d\ln a}=-\frac{d\ln A^2}{d\ln a}=-2\alpha_A
\quad\text{(adopted closure)},
\]

\[
\alpha_B=2\alpha_A,
\qquad
\alpha_K=-5\alpha_A^2,
\qquad
\alpha_T=0.
\]

No-ghost discriminant:

\[
D=\alpha_K+\frac32\alpha_B^2=\alpha_A^2\ge 0.
\]

Matter-frame Hubble under the HC reference-density convention \(\tilde\rho=\rho_{\Lambda\mathrm{CDM}}\), \(H_E=A^2 H_{\Lambda\mathrm{CDM}}\):

\[
\tilde H(z)=\frac{A(z)}{1-\alpha_A(z)}\,H_{\Lambda\mathrm{CDM}}(z),
\qquad
M_{\rm exact}(z)=\frac{A}{1-\alpha_A}.
\]

Acoustic integrals are preserved at the part-per-million level *because* they are conformal images of the same sound horizon. That is consistency of the map, not a selection of static \(g_{ij}\).

Growth in the unscreened Yukawa window (Paper 0):

\[
G_{\rm eff}(k,a)=G\left[1+\frac{2\alpha_0^2}{1+m_\phi^2 a^2/k^2}\right].
\]

---

## 10. Saturation radius (UCD) — derived *if* one accepts the soliton identification

If the Temporal Topology saturates at a characteristic density \(\rho_T\) for a spherical mass \(M\),

\[
R_T(M)=\left(\frac{3M}{4\pi\rho_T}\right)^{1/3}.
\]

If the GNSS covariance length is identified with that geometric scale for Earth,

\[
L_c\equiv R_T(M_\oplus)
\quad\Longrightarrow\quad
\rho_T=\frac{3M_\oplus}{4\pi L_c^3}.
\]

With the programme’s \(L_c\approx 4200\,\mathrm{km}\) this returns \(\rho_T\approx 20\,\mathrm{g\,cm^{-3}}\).

This is a **Level-2 identification**, not a theorem of §§1–7. It stands or falls with the MGEX per-constellation test already written in Paper 6.

---

## 11. Quantum tangent limit (QF) — endpoint only

Matter-frame Klein–Gordon from the causal-metric kinetic term:

\[
\bigl(\tilde\Box+m^2c^2/\hbar^2\bigr)\Psi=0.
\]

Eikonal limit \(\Rightarrow\) Hamilton–Jacobi equation of \(\tilde g\), not of \(g\).

Dirac operator recovered when \(\Sigma_\mu\to 0\) and \(B\to 0\):

\[
\bigl(i\gamma^\mu\partial_\mu-mc/\hbar\bigr)\psi=0,
\]

as the local Clifford/tetrad representation on an isochronous background. Spin / antiparticle structure is *proposed* as temporal-orientation holonomy on that bundle. Born rule, interactions, and a lab number are **not** derived.

---

## 12. What is *not* a final derivation

Do not treat these as consequences of the action alone:

- any specific \(\kappa_X\) (flyby \(2.56\times 10^{-3}\), LLR \(\eta_{\rm resid}\), pulsar \(\Gamma\), Cepheid bias, SN \(\epsilon_T^{\rm los}\));
- \(\lambda_T\) numerical values and the MGEX shift;
- J0437 \(\psi\) as holonomy, until Maxwell-on-\(\tilde g\) gives the transfer map \(\omega_{\rm eff}\);
- phantom-mass lensing \(\Gamma_t\) as a theorem rather than a shear projection;
- BBN abundances as an attractor;
- black holes as temporal wells without a solved \(f(\phi)\mathcal{G}\) spacetime;
- \(V(\phi)\) itself — the action writes it; no unique potential is derived.

---

## The closed chain, in one line

\[
\begin{aligned}
&\tilde g=A^2 g+B\,d\phi\otimes d\phi
\\
&\Rightarrow
\tilde\nabla\cdot\tilde T=0,\;
\Box\phi-V_{,\phi}=-\mathcal Q,\;
c_{\rm loc}\text{ invariant}
\\
&\Rightarrow
\oint d\ln A=0,\;
H_{\rm resid}=\oint(\tilde\sigma-\sigma_{\rm GR})\text{ needs }d(\Delta\sigma)\neq 0
\\
&\Rightarrow
1+z=A_0/A_{\rm em}\text{ on static homogeneous }g_{ij}
\\
&\Rightarrow
\text{PPN recovered by }\alpha_{\rm eff}=\mathcal S_\Sigma\alpha_0,\;
\text{not by }\alpha(\phi)=0.
\end{aligned}
\]

That is the mathematical skeleton. Everything else in the corpus is a map from this skeleton onto a dataset. Strengthening the programme means deriving those maps from a solved \(\phi\) with one frozen \(\mathcal{S}_\Sigma\), not adding more endpoint formulae.

The disformal connection is the mixed time–space piece of \(\tilde g_{\mu\nu}\) after the conformal background has been factored out. Everything below is the Paper 0 Appendix A3 argument written as an explicit expansion.

---

## 1. What “connection” means here

The matter metric is

\[
\tilde g_{\mu\nu}=A^2(\phi)\,g_{\mu\nu}+B(\phi)\,\nabla_\mu\phi\nabla_\nu\phi.
\]

Thread it with respect to a physical clock congruence (the network that actually holds the clocks):

\[
d\tilde s^2
=\tilde g_{00}\bigl(dt+\tilde\sigma_i\,dx^i\bigr)^2
+\tilde h_{ij}\,dx^i dx^j.
\]

The coordinate synchronization one-form is the mixed threading component

\[
\tilde\sigma_i=\frac{\tilde g_{0i}}{\tilde g_{00}}
\qquad
\bigl(\text{opposite sign convention}\Rightarrow\text{overall minus; loops unchanged}\bigr).
\]

This is not an affine connection on spacetime. It is the one-form whose line integral is the synchronization transport assigned to a spatial leg.

---

## 2. Gauge: why only the residual is physical

Re-label simultaneity on the same worldlines,

\[
t'=t+\chi(x^i),\qquad dt=dt'-\partial_i\chi\,dx^i.
\]

Then

\[
dt+\tilde\sigma_i dx^i
=dt'+(\tilde\sigma_i-\partial_i\chi)\,dx^i,
\]

so

\[
\tilde\sigma\longrightarrow\tilde\sigma-d\chi.
\]

The GR reference threading \(\sigma_{\rm GR}\) of the *same* congruence transforms identically. Therefore

\[
\Delta\sigma\equiv\tilde\sigma-\sigma_{\rm GR}
\]

is invariant, and

\[
H_{\rm resid}(C)=\oint_C\Delta\sigma
\]

is the loop class of \(\Delta\sigma\) modulo exact forms. Kinematic pieces that live in both connections (Sagnac, Thomas, coordinate rotation) cancel in \(\Delta\sigma\).

The loop is a spatial loop at a fixed epoch. Time dependence of \(\phi\) changes \(\Delta\sigma(t)\); it does not reopen the gauge argument.

---

## 3. Split the metric

Write

\[
\tilde g_{\mu\nu}
=\underbrace{A^2 g_{\mu\nu}}_{\text{conformal}}
+\underbrace{B\nabla_\mu\phi\nabla_\nu\phi}_{\delta_B\tilde g_{\mu\nu}}.
\]

Let \(u^\mu\) be the clock four-velocity, \(u^\mu u_\mu=-1\), and

\[
P_\mu{}^\nu=\delta_\mu{}^\nu+u_\mu u^\nu
\]

the projector into the local rest space. Decompose the scalar gradient:

\[
\nabla_\mu\phi
=D_\mu\phi-(u\cdot\nabla\phi)\,u_\mu,
\qquad
D_\mu\phi\equiv P_\mu{}^\nu\nabla_\nu\phi.
\]

The mixed projection of the disformal piece onto the congruence is then exact:

\[
\begin{aligned}
P_\mu{}^\alpha u^\beta\,\delta_B\tilde g_{\alpha\beta}
&=B\,(P_\mu{}^\alpha\nabla_\alpha\phi)\,(u^\beta\nabla_\beta\phi)\\
&=B\,(u\cdot\nabla\phi)\,D_\mu\phi.
\end{aligned}
\]

That scalar–vector object is the source of the disformal connection. If \(u\cdot\nabla\phi=0\) and \(D_\mu\phi=0\) fail to hold simultaneously in a way that varies around a loop, there is nothing to integrate.

---

## 4. Leading-order \(\delta\tilde\sigma\) in the clock rest frame

Choose local rest-frame coordinates of the congruence so that, at the event of interest,

\[
u^\mu=(1,0,0,0),\qquad
\tilde g_{0i}^{(A)}=0,\qquad
\tilde g_{00}^{(A)}=-A^2.
\]

Then the only mixed component at linear order in \(B\) is disformal:

\[
\delta_B\tilde g_{0i}=B\,(\nabla_0\phi)\,(\nabla_i\phi)=B\,(u\cdot\nabla\phi)\,D_i\phi.
\]

The disformal correction to \(\tilde g_{00}\) is

\[
\delta_B\tilde g_{00}=B(u\cdot\nabla\phi)^2,
\]

but it multiplies a vanishing conformal numerator \(\tilde g_{0i}^{(A)}=0\). Expanding

\[
\tilde\sigma_i
=\frac{\tilde g_{0i}^{(A)}+\delta_B\tilde g_{0i}}{\tilde g_{00}^{(A)}+\delta_B\tilde g_{00}}
\]

therefore gives, to first order in \(B\),

\[
\delta\tilde\sigma_i
=\frac{\delta_B\tilde g_{0i}}{\tilde g_{00}^{(A)}}
+O\bigl(B^2,\,B(u\cdot\nabla\phi)^2\tilde g_{0i}^{(A)}\bigr).
\]

Hence

\[
\delta\tilde\sigma_i
\simeq\frac{B(u\cdot\nabla\phi)D_i\phi}{-A^2}
=-\frac{B}{A^2}(u\cdot\nabla\phi)\,D_i\phi.
\]

Covariantly:

\[
\boxed{
\delta\tilde\sigma_\mu
\simeq
-\frac{B}{A^2}(u\cdot\nabla\phi)\,P_\mu{}^\nu\nabla_\nu\phi
=
-\frac{B}{A^2}(u\cdot\nabla\phi)\,D_\mu\phi.
}
\]

This is the disformal synchronization connection. It is linear in \(B\), linear in the clock-measured time derivative \(u\cdot\nabla\phi\), and linear in the clock-measured spatial gradient \(D\phi\).

Higher-order terms: \(\delta_B\tilde g_{00}\) in the denominator, lapse–shift mixing, and \(O\bigl((B X/A^2)^2\bigr)\).

---

## 5. ADM form: where the extra \(1/N\) comes from

In \(3+1\) variables of the Einstein-frame metric,

\[
ds^2=-N^2 dt^2+h_{ij}(dx^i+N^i dt)(dx^j+N^j dt),
\]

so

\[
g_{00}=-N^2+N_k N^k,\qquad g_{0i}=N_i,\qquad g_{ij}=h_{ij}.
\]

The matter-frame components are

\[
\begin{aligned}
\tilde g_{00}
&=A^2(-N^2+N_k N^k)+B(\partial_0\phi)^2,\\
\tilde g_{0i}
&=A^2 N_i+B(\partial_0\phi)(\partial_i\phi),\\
\tilde g_{ij}
&=A^2 h_{ij}+B(\partial_i\phi)(\partial_j\phi).
\end{aligned}
\]

The unit normal is \(n^\mu=N^{-1}(1,-N^i)\), so

\[
\partial_0\phi
=N\,(n\cdot\partial\phi)+N^i\partial_i\phi.
\]

Take small shift \(N_i\to 0\) and hypersurface-orthogonal clocks \(u^\mu=n^\mu\). Then

\[
\tilde g_{00}^{(A)}\simeq -A^2 N^2,
\qquad
\delta_B\tilde g_{0i}\simeq B\,N\,(n\cdot\partial\phi)\,\partial_i\phi,
\]

and

\[
\delta\tilde\sigma_i
=\frac{\delta_B\tilde g_{0i}}{\tilde g_{00}^{(A)}}
\simeq
\frac{B N(n\cdot\partial\phi)\,\partial_i\phi}{-A^2 N^2}
=
-\frac{B}{A^2 N}(\partial_i\phi)\,(n\cdot\partial\phi).
\]

So

\[
\boxed{
\delta\tilde\sigma_i
\approx
-\frac{B}{A^2 N}(\partial_i\phi)\,(n\cdot\partial\phi)
\qquad(u=n,\;N_i\simeq 0).
}
\]

The extra \(1/N\) is coordinate: \(\partial_0\phi=N(n\cdot\partial\phi)\) while the threading denominator is \(\tilde g_{00}\propto N^2\). In the local rest frame of §4 one has \(N=1\) and the \(1/N\) disappears. Paper 0’s “proper-time-normalized representative differs by a lapse factor” is this distinction. Closed-loop residuals use the same representative for \(\tilde\sigma\) and \(\sigma_{\rm GR}\), so the convention cancels.

---

## 6. When the connection is exact, and when it is not

Write the leading one-form as

\[
\delta\tilde\sigma
=
f\,D\phi,
\qquad
f\equiv-\frac{B}{A^2}(u\cdot\nabla\phi).
\]

Its exterior derivative on a spatial slice is

\[
d(\delta\tilde\sigma)
=df\wedge D\phi+f\,d(D\phi).
\]

- If \(f\) is a function of \(\phi\) alone and \(D\phi=d\phi\) is exact on the slice, then \(\delta\tilde\sigma\propto d\Phi(\phi)\) and \(\oint\delta\tilde\sigma=0\) on a contractible loop. That is the conformal situation: \(d\ln A\) is exact.
- \(d(\delta\tilde\sigma)\neq 0\) when \(f\) varies independently of \(\phi\) around the loop: \(\partial_t\phi\neq 0\), inhomogeneous \(B(\phi)\), spatially varying screening of \(B\), lapse/shift structure, or a congruence \(u^\mu\) whose factor \(u\cdot\nabla\phi\) is not a gradient.

Leading holonomy on a surface \(\Sigma\) bounded by \(C\):

\[
H_{\rm resid}^{(B)}(C)
=\iint_\Sigma d(\delta\tilde\sigma)
+O(B^2,\text{GR-subtraction error}).
\]

### Static field, moving clocks

If \(\phi\) is static and the slicing is adapted to \(n\) with \(N_i=0\), then \(n\cdot\partial\phi=0\) and the ADM formula vanishes. That is Theorem 2 in connection language: a static gradient plus hypersurface-orthogonal clocks produces no linear mixed term.

A GNSS or satellite congruence is *not* \(u=n\). A clock with spatial velocity \(v^i\) relative to the rest frame of \(\nabla\phi\) has

\[
u\cdot\nabla\phi
=u^i\partial_i\phi\neq 0
\]

even if \(\partial_t\phi=0\). Then

\[
\delta\tilde\sigma_i
\simeq
-\frac{B}{A^2}(u^j\partial_j\phi)\,P_i{}^k\partial_k\phi
\]

is present at linear order in \(B\) and linear order in velocity-through-the-gradient. That is the kinematic window in which a closed Earth–MEO triangle can see disformal holonomy while a static same-path one-way test cannot.

---

## 7. Relation to the photon cone (same \(B\), different observable)

The inverse (Appendix A4) is

\[
\tilde g^{\mu\nu}
=A^{-2}\left[
g^{\mu\nu}
-\frac{(B/A^2)q^\mu q^\nu}{1+(B/A^2)(q\cdot q)}
\right],
\quad q_\mu\equiv\partial_\mu\phi.
\]

Photons satisfy \(\tilde g^{\mu\nu}k_\mu k_\nu=0\). In a local inertial frame \(g_{\mu\nu}=\eta_{\mu\nu}\),

\[
\eta^{\mu\nu}k_\mu k_\nu
=\frac{(B/A^2)(q\cdot k)^2}{1+(B/A^2)(q\cdot q)}.
\]

For \(k^\mu=(\omega/c)(1,\hat n)\) and static \(q_\mu=(0,\nabla\phi)\),

\[
v_{\rm ph}
=c\left[
1+\frac{(B/A^2)(\partial_{\hat n}\phi)^2}{1+(B/A^2)|\nabla\phi|^2}
\right]^{-1/2}
\approx
c\left[
1-\frac12\frac{(B/A^2)(\partial_{\hat n}\phi)^2}{1+(B/A^2)|\nabla\phi|^2}
\right].
\]

Gravitons stay on \(g\), so \(c_g=c\), and

\[
\frac{|c_\gamma-c_g|}{c}
\lesssim
\frac{B}{2A^2}|\nabla\phi|^2
\quad(B X/A^2\ll 1).
\]

GW170817 therefore bounds the *same* combination \(B|\nabla\phi|^2/A^2\) that multiplies \(\delta\tilde\sigma\), but along an open common path and as an even cone split. Holonomy is the odd, closed-loop integral of \(\delta\tilde\sigma\). The two observables are not interchangeable.

---

## 8. Residual used in experiment

The object that is compared with data is not \(\delta\tilde\sigma\) raw. It is

\[
\Delta\sigma=\tilde\sigma-\sigma_{\rm GR},
\qquad
\tilde\sigma=\sigma^{(A)}+\delta\tilde\sigma+O(B^2).
\]

\(\sigma^{(A)}\) contains \(d\ln A\) plus the GR-like threading of \(A^2 g_{\mu\nu}\). After a complete GR subtraction on the same congruence, the conformal exact piece drops from the loop and

\[
\oint_C\Delta\sigma
=\oint_C\delta\tilde\sigma
+O(B^2,\text{model error in }\sigma_{\rm GR}).
\]

That is the disformal connection term in the form the theory actually predicts.

