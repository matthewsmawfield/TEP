> **PARTIALLY SUPERSEDED 2026-09-13 (v5.9).** The statement that
> $S_\Sigma$, $B$, the transfer map and several amplitudes are "not yet
> derived from one frozen action" is resolved by the cuscuton gradient floor
> $\mu^2/(2K)$ (CORPUS_PLAN.md §7.8 Target 1). The authoritative account is
> `CORPUS_PLAN.md` §7 (v5.9).

Yes. I have now pushed this beyond the roadmap and produced a **full analytic closure plus an explicit microscopic completion candidate** for TEP v0.11.

[Download the full derivation and canonical completion](sandbox:/mnt/data/TEP_v0.11_FULL_DERIVATION_AND_CANONICAL_COMPLETION.md)

The important distinction is that I did **not** simply fill the gaps with whatever functions make the existing results work. The corpus itself admits that \(S_\Sigma\), \(B\), the transfer map and several amplitudes are not yet derived from one frozen action.  The previous closure derivation also established that the corpus does not uniquely determine \(V\), \(K\), or the disformal function.  So I separated:

* what follows rigorously from existing TEP;
* what is an unavoidable no-go or consistency condition;
* what must be introduced as **new canonical microphysics**.

## What is now actually derived

The document now contains the complete low-energy variational structure:

$$
S=
\int d^4x\sqrt{-g}
\left[
\frac{M_{\rm Pl}^{2}}{2}R
+K(\Phi,X)-V(\Phi)
\right]
+
S_m[\psi,\tilde g_{\mu\nu}],
$$

with

$$
X=-\frac12g^{\mu\nu}\nabla_\mu\Phi\nabla_\nu\Phi
$$

and

$$
\boxed{
\tilde g_{\mu\nu}
=
A^2(\Phi)g_{\mu\nu}
+
E(\Phi,X)\nabla_\mu\Phi\nabla_\nu\Phi
},
\qquad
E=\frac{D}{M_D^4}.
$$

I derived the inverse metric, determinant and Lorentzianity condition,

$$
\Delta=1+\frac{Eq^2}{A^2}>0,
$$

the Einstein equation, and—importantly—the **full scalar variation including \(E_{,X}\)**. That latter term was still absent from the earlier closure draft.

The general matter source is now

$$
\boxed{
Q=
AA_{,\Phi}g_{\mu\nu}\mathcal T^{\mu\nu}
+
\frac12E_{,\Phi}\mathcal S
-
\nabla_\mu
\left[
E\mathcal T^{\mu\nu}q_\nu
-\frac12E_{,X}\mathcal S q^\mu
\right]
}
$$

with

$$
\mathcal S=\mathcal T^{\mu\nu}q_\mu q_\nu ,
$$

so

$$
\boxed{
\nabla_\mu(K_{,X}\nabla^\mu\Phi)
+K_{,\Phi}-V_{,\Phi}
=-Q.
}
$$

That closes the actual Euler–Lagrange problem for the proposed conformal-disformal theory.

---

# A much better screening architecture

I do **not** think the old idea of letting a fundamental \(\beta_A\) arbitrarily “run with environment” should survive.

Instead I constructed a bounded least-coupling law:

$$
\boxed{
\alpha(\Phi)
\equiv
M_{\rm Pl}\frac{d\ln A}{d\Phi}
=
\tanh
\left(
\frac{\Phi-\Phi_*}{\mu_A}
\right)
}
$$

which integrates exactly to

$$
\boxed{
A(\Phi)
=
A_*
\left[
\cosh
\left(
\frac{\Phi-\Phi_*}{\mu_A}
\right)
\right]^{\mu_A/M_{\rm Pl}}.
}
$$

This does something extremely useful.

At

$$
\Phi=\Phi_*,
$$

we have

$$
\alpha=0.
$$

Dense environments can therefore dynamically drive the field to a **real least-coupling point**, rather than saying

$$
\beta_A=-1
$$

fundamentally and then multiplying it by an empirical \(S_\Sigma\).

But on the negative unscreened branch,

$$
\Phi-\Phi_*\ll-\mu_A,
$$

we recover

$$
\boxed{\alpha\rightarrow-1}.
$$

So the corpus value \(-1\) survives, but as the **unscreened asymptotic coupling**, rather than an impossible universal local coupling.

That is, in my view, a major improvement.

---

# A concrete potential

I chose the deliberately minimal runaway form

$$
\boxed{
V(\Phi)=
V_0
\exp
\left[
\lambda
\frac{\Phi-\Phi_*}{M_{\rm Pl}}
\right]
+V_{\rm off}.
}
$$

For non-relativistic matter,

$$
V_{\rm eff}\simeq V+\rho A.
$$

The minimum condition gives

$$
\lambda V+\rho A\alpha=0,
$$

hence

$$
\boxed{
\alpha_{\rm min}
=
-\frac{\lambda V}{\rho A}.
}
$$

This is important because screening now follows directly.

For

$$
\rho A\gg\lambda V,
$$

$$
|\alpha_{\rm min}|\ll1.
$$

So high-density screening is no longer inserted by hand.

It is derived from the same \(A\) and \(V\).

---

# The effective mass is also derived

The exact local curvature is

$$
\boxed{
m_{\rm eff}^2
=
\frac{\lambda^2V}{M_{\rm Pl}^2}
+
\rho A
\left[
\frac{\alpha^2}{M_{\rm Pl}^2}
+
\frac{\operatorname{sech}^2y}
{M_{\rm Pl}\mu_A}
\right]
},
$$

where

$$
y=(\Phi-\Phi_*)/\mu_A.
$$

Near the screened point,

$$
m_{\rm eff}^2
\simeq
\frac{\lambda^2V}{M_{\rm Pl}^2}
+
\frac{\rho A}{M_{\rm Pl}\mu_A}.
$$

Far onto the unscreened branch,

$$
m_{\rm eff}^2
\simeq
\frac{(\lambda+1)\rho A}{M_{\rm Pl}^2}.
$$

That gives TEP something it badly needed: **the same scalar can have radically different local stiffness in screened and unscreened environments without changing its fundamental coupling by hand**.

---

# I also added a stable kinetic screening sector

The proposed C1 action uses

$$
\boxed{
K(X)=X+\frac{X^3}{\Lambda_K^8}.
}
$$

This particular choice has a nice mathematical property:

$$
K_{,X}
=
1+\frac{3X^2}{\Lambda_K^8}>0
$$

and

$$
K_{,X}+2XK_{,XX}
=
1+\frac{15X^2}{\Lambda_K^8}>0.
$$

So the pure k-essence principal part does not acquire the obvious ghost/gradient sign failure that some simpler nonlinear choices would.

For a static spherical source,

$$
r^2\Phi'
\left(
1+
\frac{3\Phi'^4}{4\Lambda_K^8}
\right)=C.
$$

That gives a second, genuinely dynamical way of suppressing source gradients.

This is useful because the previous derivation proved that the \(\rho_T\sim20\) g/cm³ / \(M^{1/3}\) scale cannot simultaneously be the sole Cassini screening scale. 

So I now explicitly separate:

**Temporal-topology/correlation scale**

from

**Solar-System scalar-charge suppression.**

That resolves an important conceptual conflation.

---

# The disformal sector is now much better defined

Rather than blindly importing the black-hole quartic-Gaussian \(B(\phi)\), I constructed

$$
\boxed{
D(\Phi,X)=D_0\,S_\Phi(\Phi)S_X(X)
}
$$

with

$$
\boxed{
S_X(X)
=
\frac12
\left[
1-\tanh(X/X_D)
\right].
}
$$

Because

$$
X>0
$$

for a predominantly homogeneous timelike cosmological gradient while

$$
X<0
$$

for a quasi-static spatial gradient, this provides a principled way for the disformal sector to be strongly suppressed on cosmological multimessenger propagation while remaining available in local spatial-gradient configurations.

That addresses the actual GW170817 issue more correctly than simply saying “\(B\neq0\) but is somehow screened.”

The bound is imposed on the solved combination

$$
\boxed{
\left|
\frac{D}{M_D^4A^2}
(\partial_t\Phi)^2
\right|
\lesssim10^{-15},
}
$$

rather than on \(D_0\) by itself. Multimessenger observations constrain photon/graviton propagation extraordinarily strongly, which is why disformal scalar models require precisely this sort of care. ([arXiv][1])

---

# Holonomy is now properly derived

The synchronization connection is

$$
\boxed{
\delta\sigma_\mu
=
-\frac{E}{A^2}
(u\cdot q)
P_\mu{}^\nu q_\nu.
}
$$

Therefore

$$
\boxed{
F_{\mu\nu}^{(D)}
=
-2\nabla_{[\mu}
\left[
\frac{E}{A^2}
(u\cdot q)
P_{\nu]}{}^\alpha q_\alpha
\right]
}
$$

and

$$
\boxed{
H_{\rm resid}[C]
=
\oint_C\delta\sigma
=
\int_\Sigma F^{(D)}+\cdots .
}
$$

This confirms and sharpens one of the strongest theoretical ideas in the corpus:

$$
\oint d\ln A=0.
$$

Pure conformal TEP cannot make a closed-loop signal.

But an equally important result is:

$$
B\neq0
$$

alone is **still not enough**.

You need genuinely non-exact spacetime structure.

That means observer motion, evolving gradients, anisotropic boundaries, a nontrivial shift/lapse structure or \(D(\Phi,X)\).

That makes the optical holonomy experiment considerably more meaningful.

---

# The transfer-map problem is formally solved

Instead of dozens of unexplained \(\kappa_X\)'s, the fundamental mapping becomes

$$
\boxed{
\delta O_X
=
\int d^4x\,
K_X^{\mu\nu}
\delta\tilde g_{\mu\nu}.
}
$$

And the actual prediction chain is

$$
\boxed{
\text{action}
\rightarrow
\Phi,g
\rightarrow
K_X
\rightarrow
P_X
\rightarrow
O_X^{\rm predicted}.
}
$$

This gives us three clean layers:

1. **field solution**;
2. **physical measurement transfer**;
3. **data/estimator transfer**.

For GNSS specifically,

$$
\boxed{
C_{\rm obs}
=
P_P C_{\rm field}P_P^T+C_{\rm noise}.
}
$$

That is the mathematically correct foundation for the product-operator idea in your finalisation plan. The plan was right that the product dependence must be *computed*, rather than narrated retrospectively. 

---

# There is one major result I would not hide

The full derivation confirms the earlier cosmological no-go.

For homogeneous TEP,

$$
\tilde a=Aa_g.
$$

Therefore

$$
\boxed{
\tilde H=
\frac{\dot A/A+\dot a_g/a_g}
{\sqrt{A^2N^2-E\dot\Phi^2}}.
}
$$

If you impose

$$
a_g=\mathrm{constant}
$$

but let \(A\) evolve,

$$
\tilde H\neq0.
$$

So under the existing universal matter metric, matter-frame rulers live in an effectively FLRW geometry.

The Mount Wilson result remains important:

$$
1+z=\frac{A_0}{A_e}
$$

for static \(a_g\).

But mathematically this is a **frame degeneracy**, rather than proof of static matter-frame distances. The previous closure draft established the same result. 

And the thermodynamic result follows:

$$
\boxed{
T_\gamma(z)=T_0(1+z)
}
$$

if radiation is separately conserved.

So the current minimal theory cannot simultaneously have:

* universal matter coupling;
* \(1+z=A_0/A_e\);
* and “the high-redshift universe was not physically hotter.”

The earlier derivation explicitly reached this conclusion too. 

This is now the **single biggest remaining foundational choice**.

It cannot be repaired honestly by another screening function.

You either:

**A. keep the minimal theory**, in which case TEP is Einstein-frame static/background-degenerate but retains ordinary matter-frame hot-redshift thermodynamics;

or

**B. retain the stronger eternal/cold ontology**, in which case we need to construct a new non-equilibrium radiation/matter interaction sector or change the universal matter coupling.

That would be a genuine TEP v0.12-level theoretical addition.

---

## What remains now

The important point is that most of what remains is no longer “missing derivation.”

It is **numerical forward calculation**.

The new document reduces the next stage to ten concrete calculations:

1. solve the Sun/Earth nonlinear boundary-value problem for C1;
2. see whether its **exterior** field really produces the 4,200 km GNSS covariance;
3. calculate Cassini/MICROSCOPE from that exact same solution;
4. put its covariance through CODE/GFZ/MGEX/PPP \(P_P\);
5. put the acceleration kernel directly into INPOP/DE for LLR;
6. calculate Cepheid stellar-response and the direct core–disk \(q_i\);
7. calculate the J0437 EM/scintillation kernel;
8. calculate signed optical-loop \(H_{\rm resid}\);
9. solve C1 cosmological background + perturbations + recombination/high-\(\ell\);
10. only afterwards build the global BH strong-curvature completion.

That is considerably better than the original plan's idea that Paper 0 v0.11 was merely assembling some existing functions.

The corpus audit itself correctly said the work remaining is substantive derivation and calculation, not errata. 

### Where I think TEP now stands

We now have a **coherent candidate mathematical theory to try to kill**.

That is the right milestone.

C1 is not retroactively “what TEP always was.” It is a concrete proposed completion generated from the requirements of the existing framework. Its coupling, screening, kinetic and disformal functions are now explicit enough that they cannot be changed quietly every time a dataset disagrees.

The next phase should therefore be the **numerical C1 closure programme**, beginning with the Earth/Sun boundary-value equations. That calculation will immediately tell us whether this proposed completion can simultaneously survive Cassini and produce anything resembling the GNSS field. If it fails there, there is no reason to propagate C1 into another twenty papers.

[1]: https://arxiv.org/abs/1710.05893?utm_source=chatgpt.com "Implications of the Neutron Star Merger GW170817 for Cosmological Scalar-Tensor Theories"

# TEP v0.11 — Full Closure Derivation and Canonical Microscopic Completion

**Date:** 5 September 2026  
**Status:** Proposed canonical completion for falsification, not a claim that the microscopic functions were already derived by the earlier corpus.

---

## 0. Executive conclusion

The existing TEP corpus is sufficient to derive the universal two-metric variational structure, exact scalar source equation, matter-frame conservation law, weak-field transfer operator, PPN/source-charge mapping, disformal photon-cone correction, and closed-loop holonomy invariant.

It is **not** sufficient to derive unique microscopic functions \(A(\Phi)\), \(K(X)\), \(V(\Phi)\), and \(D(\Phi,X)\). Any attempt to fill those remaining functions is a new theory choice.

This document therefore does two things:

1. completes everything that follows uniquely from the present TEP axioms; and  
2. proposes one explicit **Canonical Completion C1** designed to satisfy the constraints identified by the corpus audit.

C1 is deliberately rigid. Once its parameters are fixed on a small theory/constraint set, its downstream predictions must be computed without refitting. If it fails, C1 fails and the failed prediction remains in the evidence ledger.

The largest remaining no-go is cosmological: under universal matter coupling to the current conformal-disformal metric, a homogeneous evolving conformal factor produces a matter-frame FLRW scale factor and the ordinary adiabatic \(T_\gamma\propto1+z\) law. A genuinely static *matter-frame* universe with a non-hot high-redshift radiation bath cannot be obtained from the minimal universal action alone. It requires a new non-equilibrium matter/radiation sector or a different matter coupling. This is not repaired by relabelling frames.

---

# Part I — Formal TEP theory

## 1. Field content and conventions

Use signature \((-+++)\) and a dimension-one scalar \(\Phi\).

\[
X\equiv-\frac12 g^{\mu\nu}\nabla_\mu\Phi\nabla_\nu\Phi.
\]

The universal matter metric is

\[
\boxed{
\tilde g_{\mu\nu}
=
A^2(\Phi)g_{\mu\nu}
+
E(\Phi,X)\nabla_\mu\Phi\nabla_\nu\Phi
}
\]

with

\[
E(\Phi,X)\equiv \frac{D(\Phi,X)}{M_D^4}.
\]

The explicit scale \(M_D\) removes the dimensional ambiguity in earlier \(B(\phi)\) notation.

Define

\[
q_\mu=\nabla_\mu\Phi,\qquad
q^2=g^{\mu\nu}q_\mu q_\nu=-2X.
\]

For \(E=E(\Phi)\), the inverse metric is exactly

\[
\boxed{
\tilde g^{\mu\nu}
=
A^{-2}
\left[
g^{\mu\nu}
-
\frac{E q^\mu q^\nu}
{A^2+Eq^2}
\right]
}
\]

and

\[
\boxed{
\sqrt{-\tilde g}
=
A^4\sqrt{-g}\sqrt{\Delta},
\qquad
\Delta=1+\frac{Eq^2}{A^2}.
}
\]

The necessary matter-metric signature/invertibility conditions are

\[
A>0,\qquad \Delta>0.
\]

For \(E(\Phi,X)\), the algebraic inverse remains of disformal form but the variational equations acquire \(E_{,X}\) contributions.

---

## 2. Canonical low-energy action class

The minimal low-energy theory is

\[
\boxed{
S_{\rm LE}
=
\int d^4x\sqrt{-g}
\left[
\frac{M_{\rm Pl}^2}{2}R
+
K(\Phi,X)
-
V(\Phi)
\right]
+
S_m[\psi,\tilde g_{\mu\nu}].
}
\]

The canonical limit is \(K=X\).

Strong-curvature operators, such as scalar-Gauss-Bonnet terms, are **not** required by the low-energy GNSS/LLR/cosmology sector and should not be inserted into the universal freeze as though they were already derived. They belong to a separate strong-curvature completion,

\[
S_{\rm SC}
=
\int d^4x\sqrt{-g}\,
\alpha_{\rm GB}f(\Phi)\mathcal G\,W(\mathcal K/\mathcal K_*)
+\cdots ,
\]

with \(W\to0\) in the weak-curvature EFT regime.

---

## 3. Einstein equation

Define the matter-frame stress tensor by

\[
\tilde T^{\mu\nu}
\equiv
\frac{2}{\sqrt{-\tilde g}}
\frac{\delta S_m}{\delta\tilde g_{\mu\nu}}.
\]

For \(E=E(\Phi)\), variation with respect to \(g_{\mu\nu}\) gives

\[
\delta \tilde g_{\alpha\beta}
=
A^2\delta g_{\alpha\beta},
\]

so the Einstein-frame matter tensor is

\[
\mathcal T^{\mu\nu}
=
A^2\frac{\sqrt{-\tilde g}}{\sqrt{-g}}
\tilde T^{\mu\nu}
=
A^6\sqrt{\Delta}\,\tilde T^{\mu\nu}.
\]

For general \(K(\Phi,X)\),

\[
T_{\mu\nu}^{(\Phi)}
=
K_{,X}\nabla_\mu\Phi\nabla_\nu\Phi
+
g_{\mu\nu}(K-V).
\]

Thus

\[
\boxed{
M_{\rm Pl}^2G_{\mu\nu}
=
T_{\mu\nu}^{(\Phi)}
+
\mathcal T_{\mu\nu}
}
\]

for \(E=E(\Phi)\). If \(E\) depends on \(X\), extra metric-variation terms proportional to \(E_{,X}\) must be included; C1 below therefore treats the \(X\)-dependent disformal piece as part of the matter-sector variation and requires a full symbolic variation before numerical deployment.

---

## 4. Scalar equation for \(E=E(\Phi)\)

Matter variation gives

\[
\delta S_m
=
\frac12
\int d^4x\,\sqrt{-\tilde g}\,
\tilde T^{\mu\nu}\,
\delta\tilde g_{\mu\nu}.
\]

Using

\[
\delta\tilde g_{\mu\nu}
=
2AA_{,\Phi}\delta\Phi\,g_{\mu\nu}
+
E_{,\Phi}\delta\Phi\,q_\mu q_\nu
+
E(\nabla_\mu\delta\Phi\,q_\nu+q_\mu\nabla_\nu\delta\Phi)
\]

and defining

\[
\mathcal T^{\mu\nu}
\equiv
\frac{\sqrt{-\tilde g}}{\sqrt{-g}}\tilde T^{\mu\nu},
\]

one obtains the matter source

\[
Q
=
AA_{,\Phi}g_{\mu\nu}\mathcal T^{\mu\nu}
+
\frac12E_{,\Phi}\mathcal T^{\mu\nu}q_\mu q_\nu
-
\nabla_\mu(E\mathcal T^{\mu\nu}q_\nu).
\]

For canonical \(K=X\),

\[
\boxed{
\Box\Phi - V_{,\Phi}=-Q.
}
\]

For general \(K\),

\[
\boxed{
\nabla_\mu(K_{,X}\nabla^\mu\Phi)
+
K_{,\Phi}
-
V_{,\Phi}
=
-Q.
}
\]

Matter is covariantly conserved in the matter frame:

\[
\boxed{
\tilde\nabla_\mu\tilde T^{\mu\nu}=0.
}
\]

---

## 5. Weak-field perturbations and the GNSS covariance theorem

Let

\[
\Phi=\bar\Phi+\varphi .
\]

For a locally linear perturbation operator,

\[
Z_{\rm eff}\Box\varphi-m_{\rm loc}^2\varphi=J.
\]

After canonical normalization,

\[
(\Box-m_{\rm eff}^2)\varphi=J,
\qquad
m_{\rm eff}^2=\frac{m_{\rm loc}^2}{Z_{\rm eff}}.
\]

For static perturbations,

\[
(\nabla^2-m_{\rm eff}^2)\varphi=J,
\]

with Yukawa Green function

\[
G(r)=-\frac{e^{-m_{\rm eff}r}}{4\pi r}.
\]

If the source fluctuations are short-correlated,

\[
\langle J(\mathbf k)J(\mathbf k')\rangle
=
(2\pi)^3P_J\delta^3(\mathbf k+\mathbf k'),
\]

then

\[
P_\varphi(k)=\frac{P_J}{(k^2+m_{\rm eff}^2)^2},
\]

and

\[
\boxed{
C_\varphi(r)
=
\frac{P_J}{8\pi m_{\rm eff}}e^{-m_{\rm eff}r}.
}
\]

For a locally exponential conformal slope
\(\alpha=M_{\rm Pl}\,d\ln A/d\Phi\),

\[
\delta\ln A=\frac{\alpha}{M_{\rm Pl}}\varphi
\]

and therefore

\[
\boxed{
C_A(r)
=
\frac{\alpha^2P_J}{8\pi M_{\rm Pl}^2m_{\rm eff}}
e^{-r/\lambda_\Phi},
\qquad
\lambda_\Phi=m_{\rm eff}^{-1}.
}
\]

If, and only if, the observed GNSS exponential is identified with this local propagator,

\[
\lambda_T=4200\,{\rm km}
\]

corresponds to

\[
\boxed{
m_Tc^2=\frac{\hbar c}{\lambda_T}
\simeq4.70\times10^{-14}\ {\rm eV}.
}
\]

This interpretation is **conditional**. GNSS clocks lie in the exterior terrestrial environment, so the actual Earth boundary-value solution must demonstrate that the relevant exterior perturbation operator has this mass. A bulk-Earth density inserted directly into \(m_{\rm eff}\) is not a derivation.

---

## 6. \(R_T\), \(\rho_T\), and the independence issue

The UCD relation is

\[
R_T(M)=
\left(\frac{3M}{4\pi\rho_T}\right)^{1/3}.
\]

If one identifies \(R_T(M_\oplus)=\lambda_T\),

\[
\boxed{
\rho_T=
\frac{3M_\oplus}{4\pi\lambda_T^3}.
}
\]

For \(\lambda_T=4200\) km,

\[
\rho_T\simeq19.24\ {\rm g\,cm^{-3}},
\]

rounded in the corpus to \(\sim20\ {\rm g\,cm^{-3}}\).

Thus \(\rho_T\) is presently a transform of \(M_\oplus\) and \(\lambda_T\), not an independent empirical constant. It becomes an independent cross-domain constant only after a second instrument or source determines it without importing the GNSS length.

---


## 4A. Exact matter variation for \(E(\Phi,X)\)

For the C1 disformal sector \(E=E(\Phi,X)\), the \(X\)-dependence must be varied explicitly.

Since

\[
X=-\frac12g^{\alpha\beta}q_\alpha q_\beta,
\]

the scalar variation gives

\[
\delta X=-q^\mu\nabla_\mu\delta\Phi.
\]

Define

\[
\mathcal S\equiv\mathcal T^{\mu\nu}q_\mu q_\nu .
\]

The derivative part of the matter variation is then

\[
\left[
E\mathcal T^{\mu\nu}q_\nu
-\frac12E_{,X}\mathcal S\,q^\mu
\right]\nabla_\mu\delta\Phi .
\]

Therefore the exact matter source generalises to

\[
\boxed{
Q
=
AA_{,\Phi}g_{\mu\nu}\mathcal T^{\mu\nu}
+
\frac12E_{,\Phi}\mathcal S
-
\nabla_\mu
\left[
E\mathcal T^{\mu\nu}q_\nu
-\frac12E_{,X}\mathcal S q^\mu
\right].
}
\]

The scalar equation is consequently

\[
\boxed{
\nabla_\mu(K_{,X}\nabla^\mu\Phi)
+
K_{,\Phi}
-
V_{,\Phi}
=
-Q.
}
\]

The \(X\)-dependence also modifies the Einstein-frame matter stress. Under variation of the covariant Einstein metric,

\[
\delta X
=
\frac12q^\mu q^\nu\delta g_{\mu\nu},
\]

so

\[
\delta\tilde g_{\alpha\beta}
=
A^2\delta g_{\alpha\beta}
+
\frac12E_{,X}q_\alpha q_\beta q^\mu q^\nu
\delta g_{\mu\nu}.
\]

Hence

\[
\boxed{
\mathcal T_g^{\mu\nu}
=
\frac{\sqrt{-\tilde g}}{\sqrt{-g}}
\left[
A^2\tilde T^{\mu\nu}
+
\frac12E_{,X}
\left(
\tilde T^{\alpha\beta}q_\alpha q_\beta
\right)
q^\mu q^\nu
\right].
}
\]

The Einstein equation for C1 is therefore

\[
\boxed{
M_{\rm Pl}^2G_{\mu\nu}
=
T_{\mu\nu}^{(\Phi)}
+
\mathcal T^{(g)}_{\mu\nu}.
}
\]

This closes the variational step introduced by the \(X\)-selective disformal switch. Any numerical C1 solver must use these equations rather than the simpler \(E(\Phi)\) formulas.

---

# Part II — Correct screening concept

## 7. Source charge, not a fitted algebraic multiplier

Define the microscopic logarithmic coupling

\[
\alpha(\Phi)
=
M_{\rm Pl}\frac{d\ln A}{d\Phi}.
\]

For a compact source the exterior solution is

\[
\Phi(r)
=
\Phi_\infty
-
\frac{Q_\Phi}{4\pi r}e^{-m_\infty r}
+\cdots .
\]

The physically relevant screened charge is \(Q_\Phi\). Define

\[
\boxed{
\alpha_{\rm eff}
\equiv
\frac{Q_\Phi}{M_{\rm source}}
}
\]

up to the chosen scalar normalization.

The corpus screening operator must therefore become

\[
\boxed{
S_\Sigma(\mathcal E)
\equiv
\frac{\alpha_{\rm eff}(\mathcal E)}
{\alpha(\Phi_\infty)}.
}
\]

It is an **output** of the nonlinear field solution, not an input function fitted separately in each domain.

For a light scalar on the measurement scale,

\[
\gamma_{\rm PPN}-1
\simeq
-2\alpha_{\rm source}\alpha_{\rm test}.
\]

A Yukawa range adds the appropriate \(e^{-r/\lambda}\) factor.

---

## 8. Why the UCD radius cannot be the sole PPN screening radius

For the cubic derivative model

\[
\mathcal L
=
X
-
\frac{1}{\Lambda_3^3}X\Box\Phi
+
\frac{\beta\Phi}{M_{\rm Pl}}T,
\]

the static spherical equation integrates to

\[
r^2\Phi'
+
\frac{2}{\Lambda_3^3}r(\Phi')^2
=
C,
\qquad
C=\frac{\beta M}{4\pi M_{\rm Pl}}.
\]

The physical branch is

\[
\Phi'
=
\frac{\Lambda_3^3r}{4}
\left[
\sqrt{1+(r_V/r)^3}-1
\right],
\]

with

\[
r_V^3=\frac{8C}{\Lambda_3^3}.
\]

The suppression relative to the linear field is

\[
S_V(r)
=
\frac{2}
{1+\sqrt{1+(r_V/r)^3}}.
\]

Thus \(r_V\propto M^{1/3}\). But setting \(r_V=R_T\) with
\(\rho_T=20\ {\rm g\,cm^{-3}}\) gives

\[
r_{V,\odot}\simeq2.87\times10^5\,{\rm km},
\]

inside the solar photosphere. It cannot suppress an order-unity bare conformal coupling sufficiently at a grazing Cassini ray.

Therefore the corpus-wide theory must distinguish:

1. a mesoscopic/topological/covariance scale; and
2. source-charge/PPN suppression.

They may come from the same action but they cannot be assumed to be the same radius.

---

# Part III — Canonical Completion C1

## 9. Design principles

C1 is a new microscopic branch chosen to satisfy these requirements:

- no arbitrary environment-dependent \(\beta_A\);
- bounded microscopic coupling;
- a true least-coupling point for dense objects;
- asymptotic recovery of the corpus sign and magnitude \(\alpha\to-1\);
- a stable nonlinear kinetic sector for additional source screening;
- a disformal term that is suppressed on homogeneous timelike cosmological propagation but can activate on quasi-static spatial gradients;
- one frozen action to be used in every downstream domain.

---

## 10. C1 conformal coupling: bounded least-coupling law

Define

\[
y\equiv\frac{\Phi-\Phi_*}{\mu_A}.
\]

Choose

\[
\boxed{
\alpha(\Phi)
=
M_{\rm Pl}\frac{d\ln A}{d\Phi}
=
\tanh y.
}
\]

Integrating,

\[
\boxed{
A(\Phi)
=
A_*
\left[
\cosh\left(
\frac{\Phi-\Phi_*}{\mu_A}
\right)
\right]^{\mu_A/M_{\rm Pl}}.
}
\]

Properties:

- at \(\Phi=\Phi_*\), \(\alpha=0\): dense objects can be driven to a genuine least-coupling point;
- for \(\Phi-\Phi_*\ll-\mu_A\), \(\alpha\to-1\): the canonical TEP sign and order-unity unscreened slope are recovered;
- for \(\Phi-\Phi_*\gg\mu_A\), \(\alpha\to+1\);
- the coupling is bounded: \(|\alpha|\le1\);
- there is no arbitrary \(\beta_A(\mathcal E)\).

The old registry value \(\beta_A=-1\) is therefore reclassified as the **negative unscreened asymptote**, not a constant slope at every field value.

This change is necessary if TEP wants genuine density screening from one matter coupling rather than fitting a separate multiplicative suppression factor.

---

## 11. C1 potential

Choose the minimal monotonic potential

\[
\boxed{
V(\Phi)
=
V_0
\exp\left[
\lambda\frac{\Phi-\Phi_*}{M_{\rm Pl}}
\right]
+
V_{\rm off}.
}
\]

\(V_{\rm off}\) is an irrelevant constant for the scalar equation but contributes to the gravitational background and must be fixed by the chosen cosmological branch.

For non-relativistic matter and negligible disformal corrections,

\[
V_{\rm eff}(\Phi)
\simeq
V(\Phi)+\rho A(\Phi).
\]

The stationary condition is

\[
V_{,\Phi}
+
\rho A_{,\Phi}
=0,
\]

or

\[
\boxed{
\lambda V+\rho A\,\alpha=0.
}
\]

Therefore

\[
\boxed{
\alpha_{\rm min}
=
-\frac{\lambda V}{\rho A}.
}
\]

This is the key screening result.

In a dense environment with

\[
\rho A\gg\lambda V,
\]

one obtains

\[
|\alpha_{\rm min}|\ll1,
\]

and the field sits near the least-coupling point.

In a sufficiently dilute environment, the field moves onto the negative branch and \(\alpha\to-1\).

No phenomenological screening factor has been inserted.

---

## 12. C1 effective mass

Because

\[
\frac{A_{,\Phi}}{A}
=
\frac{\alpha}{M_{\rm Pl}},
\]

and

\[
\frac{d\alpha}{d\Phi}
=
\frac{1}{\mu_A}\operatorname{sech}^2 y,
\]

we have

\[
\frac{A_{,\Phi\Phi}}{A}
=
\frac{\alpha^2}{M_{\rm Pl}^2}
+
\frac{\operatorname{sech}^2y}
{M_{\rm Pl}\mu_A}.
\]

Also

\[
V_{,\Phi\Phi}
=
\frac{\lambda^2}{M_{\rm Pl}^2}V.
\]

Thus at an effective-potential minimum

\[
\boxed{
m_{\rm eff}^2
=
\frac{\lambda^2V}{M_{\rm Pl}^2}
+
\rho A
\left[
\frac{\alpha^2}{M_{\rm Pl}^2}
+
\frac{\operatorname{sech}^2y}
{M_{\rm Pl}\mu_A}
\right].
}
\]

Dense least-coupling limit \(y\simeq0\):

\[
\boxed{
m_{\rm eff}^2
\simeq
\frac{\lambda^2V}{M_{\rm Pl}^2}
+
\frac{\rho A}{M_{\rm Pl}\mu_A}.
}
\]

Dilute asymptotic limit \(y\ll-1\), \(\alpha\simeq-1\), \(\operatorname{sech}^2y\simeq0\):

\[
m_{\rm eff}^2
\simeq
\frac{\lambda^2V+\rho A}{M_{\rm Pl}^2}.
\]

Using \(\lambda V\simeq\rho A\) at the minimum,

\[
\boxed{
m_{\rm eff}^2
\simeq
\frac{(\lambda+1)\rho A}{M_{\rm Pl}^2}
}
\]

in the dilute tracking branch.

This has the desired qualitative multi-scale structure: very high local curvature is possible near the least-coupling point through the \(1/(M_{\rm Pl}\mu_A)\) term, while the dilute branch reverts to a Planck-suppressed Hubble-like mass scale.

The actual Earth exterior solution must still be solved before identifying \(m_{\rm eff}^{-1}\) with the GNSS length.

---

## 13. C1 nonlinear kinetic sector

Choose

\[
\boxed{
K(X)
=
X+\frac{X^3}{\Lambda_K^8}.
}
\]

Then

\[
K_{,X}=1+\frac{3X^2}{\Lambda_K^8}>0,
\]

and

\[
K_{,X}+2XK_{,XX}
=
1+\frac{15X^2}{\Lambda_K^8}>0.
\]

Thus the pure k-essence principal part has positive kinetic and radial sound coefficients for all real \(X\).

For a quasi-static spherical profile outside a source,

\[
\frac1{r^2}
\frac{d}{dr}
\left[
r^2K_{,X}\Phi'
\right]
\simeq
V_{,\Phi}
+
\frac{\alpha A}{M_{\rm Pl}}\rho.
\]

Outside the material source and over a region in which the potential is negligible,

\[
r^2K_{,X}\Phi'=C.
\]

With \(X=-\Phi'^2/2\),

\[
K_{,X}
=
1+\frac{3\Phi'^4}{4\Lambda_K^8}.
\]

Hence

\[
\boxed{
r^2\Phi'
\left(
1+\frac{3\Phi'^4}{4\Lambda_K^8}
\right)
=C.
}
\]

This supplies an additional source-gradient suppression independently of the least-coupling mechanism. Its scale is not identified with \(\rho_T\).

---

## 14. C1 disformal sector

The GW170817 constraint is on the realized combination
\(E(\nabla\Phi)^2/A^2\), not on a bare dimensionless number.

To suppress the disformal sector on a homogeneous timelike cosmological background while allowing activation in quasi-static spatial-gradient regions, define

\[
\boxed{
D(\Phi,X)
=
D_0\,S_\Phi(\Phi)\,S_X(X),
}
\]

where

\[
S_X(X)
=
\frac12
\left[
1-\tanh\left(\frac{X}{X_D}\right)
\right].
\]

Because

- homogeneous rolling backgrounds have \(X>0\),
- quasi-static spatial gradients have \(X<0\),

one obtains

\[
S_X\to0
\quad (X\gg X_D),
\]

and

\[
S_X\to1
\quad (X\ll-X_D).
\]

A bounded field-space window may be chosen, for example

\[
S_\Phi(\Phi)
=
\operatorname{sech}^2
\left(
\frac{\Phi-\Phi_D}{\mu_D}
\right),
\]

if the data require a localized shear zone. This is optional and must not be introduced unless demanded by a frozen fit/constraint stage.

The matter metric is then

\[
\boxed{
\tilde g_{\mu\nu}
=
A^2g_{\mu\nu}
+
\frac{D_0S_\Phi S_X}{M_D^4}
\nabla_\mu\Phi\nabla_\nu\Phi.
}
\]

The signature condition is

\[
1+
\frac{D_0S_\Phi S_X}{M_D^4A^2}q^2>0.
\]

The multimessenger bound is imposed directly on the solved background:

\[
\boxed{
\left|
\frac{D(\Phi,X)}{M_D^4A^2}
(\partial_t\Phi)^2
\right|
\lesssim10^{-15}
}
\]

along the relevant late-time astrophysical propagation path.

This construction does not prove that a measurable holonomy exists. It makes it possible for the disformal sector to be negligible on the GW path and non-negligible in quasi-static gradient regions without simply setting \(D=0\) everywhere.

---

# Part IV — Disformal observables

## 15. Photon cone

In a local Einstein-frame inertial patch,

\[
\tilde g_{\mu\nu}
=
A^2\eta_{\mu\nu}
+
E q_\mu q_\nu.
\]

For propagation in spatial direction \(\mathbf n\), define

\[
q_t=\partial_t\Phi,\qquad
q_n=\nabla\Phi\cdot\mathbf n.
\]

The matter-null condition is

\[
A^2(-dt^2+d\ell^2)
+
E(q_tdt+q_nd\ell)^2
=0.
\]

At first order in \(Eq^2/A^2\),

\[
\boxed{
\frac{v_\gamma(\mathbf n)}{c}
=
1-
\frac{E}{2A^2}(q_t+q_n)^2
+\mathcal O(E^2).
}
\]

Reversal gives

\[
\frac{v_\gamma(-\mathbf n)}{c}
=
1-
\frac{E}{2A^2}(q_t-q_n)^2
+\cdots ,
\]

so

\[
\boxed{
\frac{v_\gamma(\mathbf n)-v_\gamma(-\mathbf n)}{c}
=
-\frac{2E}{A^2}q_tq_n
+\mathcal O(E^2).
}
\]

---

## 16. Closed-loop synchronization curvature

Let \(u^\mu\) be the clock-network congruence and

\[
P_\mu{}^\nu=\delta_\mu{}^\nu+u_\mu u^\nu.
\]

The leading disformal synchronization one-form is

\[
\boxed{
\delta\sigma_\mu
=
-\frac{E}{A^2}
(u\cdot q)
P_\mu{}^\nu q_\nu.
}
\]

Its curvature is

\[
\boxed{
F^{(D)}_{\mu\nu}
=
2\nabla_{[\mu}\delta\sigma_{\nu]}
=
-2\nabla_{[\mu}
\left[
\frac{E}{A^2}
(u\cdot q)
P_{\nu]}{}^\alpha q_\alpha
\right].
}
\]

The GR-subtracted holonomy is

\[
\boxed{
H_{\rm resid}[C]
=
\oint_C\delta\sigma
=
\int_\Sigma F^{(D)}
+\cdots .
}
\]

Consequences:

1. pure conformal transport gives \(d\ln A\) and vanishes on a smooth closed loop;
2. \(D\neq0\) is not sufficient;
3. if \(u\cdot q=0\), the leading local holonomy vanishes;
4. if the full one-form is exact, the loop still vanishes;
5. a measurable local holonomy requires independent spacetime structure: observer motion, evolving gradients, anisotropic boundaries, nontrivial shift/lapse, \(D(\Phi,X)\), or topology.

Therefore the optical-triangle amplitude must be calculated from the solved Earth/spacetime field. A universal \(10^{-19}\) s value cannot be assigned before that solution exists.

---

# Part V — Universal transfer map

## 17. Field → instrument → estimator

For an observable \(O_X\),

\[
K_X^{\mu\nu}(x)
\equiv
\frac{\delta O_X}
{\delta\tilde g_{\mu\nu}(x)}.
\]

Then

\[
\boxed{
\delta O_X
=
\int d^4x\,
K_X^{\mu\nu}(x)\,
\delta\tilde g_{\mu\nu}(x).
}
\]

The complete prediction chain is

\[
\boxed{
\text{frozen action}
\rightarrow
\Phi(x),g_{\mu\nu}(x)
\rightarrow
K_X
\rightarrow
P_X
\rightarrow
\widehat O_X^{\rm pred}.
}
\]

Here \(P_X\) is the actual data-reduction/estimator operator.

A published coefficient \(\kappa_X\) is a **prediction** only after this entire chain is closed.

### Direct clocks

For an Einstein-frame worldline \(u^\mu\),

\[
\frac{d\tilde\tau}{d\tau_g}
=
A
\sqrt{
1-\frac{E(u\cdot q)^2}{A^2c^2}
}.
\]

Hence

\[
\boxed{
\delta\ln\nu_{\rm clock}
=
\delta\ln A
-
\frac12
\delta
\left[
\frac{E(u\cdot q)^2}{A^2c^2}
\right].
}
\]

### Slow-body dynamics

At leading conformal order,

\[
\boxed{
\mathbf a_\Phi
=
-c^2\nabla\ln A
=
-\frac{c^2\alpha(\Phi)}{M_{\rm Pl}}
\nabla\Phi,
}
\]

with source charge supplied by the nonlinear solution.

LLR, flyby, and wide-binary claims must be generated by this acceleration inside their dynamical integrators.

### GNSS

If the field covariance is \(C_{\rm field}\),

\[
\boxed{
C_{\rm obs}
=
P_P C_{\rm field}P_P^T+C_{\rm noise}.
}
\]

This is the product-operator closure for CODE, GFZ, MGEX, broadcast/SPP and PPP.

### Cepheids

The scalar theory does not by itself yield a magnitude coefficient. It must be coupled to a stellar-pulsation response operator. The direct target is the predicted systemic/core versus Cepheid-field timing/spectroscopic differential \(q_i\), not another fitted ladder coefficient.

### J0437

The action does not by itself predict \(\bar\psi\simeq1\) rad. One must compute the EM/scintillation transfer operator and then apply the closure-phase estimator.

---

# Part VI — Cosmological closure and no-go

## 18. Background theorem

Take

\[
ds_g^2
=
-N(t)^2dt^2
+
a_g(t)^2h_{ij}dx^idx^j
\]

and homogeneous \(\Phi(t)\).

Then

\[
ds_{\tilde g}^2
=
-[A^2N^2-E\dot\Phi^2]dt^2
+
A^2a_g^2h_{ij}dx^idx^j.
\]

Matter proper time is

\[
d\tilde\tau
=
\sqrt{A^2N^2-E\dot\Phi^2}\,dt
\]

and the matter-frame scale factor is

\[
\boxed{
\tilde a=Aa_g.
}
\]

Therefore

\[
\boxed{
\tilde H
=
\frac1{\tilde a}
\frac{d\tilde a}{d\tilde\tau}
=
\frac{
\dot A/A+\dot a_g/a_g
}{
\sqrt{A^2N^2-E\dot\Phi^2}
}.
}
\]

If \(a_g\) is static,

\[
\tilde H
=
\frac{\dot A/A}
{\sqrt{A^2N^2-E\dot\Phi^2}}.
\]

Thus an evolving \(A\) produces an evolving physical matter-frame spatial metric.

In the conformal limit,

\[
1+z
=
\frac{\tilde a_0}{\tilde a_e}
=
\frac{A_0a_{g0}}{A_ea_{ge}}.
\]

For static \(a_g\),

\[
1+z=\frac{A_0}{A_e}.
\]

This is a frame/background degeneracy with FLRW.

It is **not** a proof that matter-measured physical spatial distances remain static.

---

## 19. Thermodynamics

Matter-frame conservation gives

\[
\frac{d\tilde\rho}{d\tilde\tau}
+
3\tilde H(\tilde\rho+\tilde p)
=0
\]

for a separately conserved component.

For radiation,

\[
\tilde\rho_\gamma\propto\tilde a^{-4},
\qquad
T_\gamma\propto\tilde a^{-1}.
\]

Therefore

\[
\boxed{
T_\gamma(z)=T_0(1+z)
}
\]

under the minimal universal matter theory.

A non-hot high-redshift branch requires

\[
\dot\rho_\gamma+4\tilde H\rho_\gamma=Q_\gamma
\]

with an explicit interaction model. To keep \(\rho_\gamma\) constant would require

\[
Q_\gamma=4\tilde H\rho_\gamma.
\]

No such term is generated by the present minimal universal scalar/conformal coupling for classical radiation.

Hence:

\[
\boxed{
\text{cold eternal matter-frame cosmology is not closed by C1 or the existing minimal TEP action.}
}
\]

It requires a new, separately derived radiative/chemical non-equilibrium sector or a foundational change to the matter coupling.

---

# Part VII — Strong curvature

## 20. Black-hole completion

The low-energy functions must be frozen first.

Then:

1. add the strong-curvature operator;
2. solve the coupled \(g_{\mu\nu},\Phi\) equations with one boundary-value problem;
3. verify \(\Delta>0\);
4. verify tensor/scalar hyperbolicity and absence of ghosts;
5. verify finite curvature invariants and geodesic completeness;
6. derive the causal horizon from the solved matter metric;
7. compute EHT and QNM observables from the same global solution;
8. compare against Kerr using the same likelihood.

The present prescribed BH profiles cannot be declared a global solution until this is done.

---

# Part VIII — Parameter freeze

## 21. Proposed registry classes

### Fundamental / theory-choice parameters

- \(M_{\rm Pl}\)
- \(A_*,\Phi_*,\mu_A\)
- \(V_0,\lambda,V_{\rm off}\)
- \(\Lambda_K\)
- \(D_0,M_D,X_D\)
- optional \(\Phi_D,\mu_D\) only if the field-space disformal window is retained

### Derived field quantities

- \(\Phi_{\rm env}\)
- \(\alpha_{\rm env}=\tanh[(\Phi_{\rm env}-\Phi_*)/\mu_A]\)
- source charge \(Q_\Phi\)
- \(S_\Sigma=Q_\Phi/[\alpha_\infty M]\)
- local \(m_{\rm eff}\)
- field covariance \(C_\Phi(x,x',t,t')\)
- holonomy curvature \(F^{(D)}_{\mu\nu}\)

### Calibrated but not independent

- \(\lambda_T\) until predicted on a held-out product
- \(\rho_T\) while defined from \(M_\oplus,\lambda_T\)

### Channel outputs

Every \(\kappa_X\) must be tagged as:

- fitted
- inherited
- proxy
- predicted_through_transfer

---

# Part IX — Mandatory execution order

## 22. Theory-only fit/freeze

C1 must first be constrained using only theory consistency and genuinely external constraints:

1. matter-metric Lorentzianity/invertibility;
2. hyperbolicity of the scalar principal part;
3. local PPN/Cassini/MICROSCOPE constraints;
4. GW170817 path constraint;
5. absence of pathological cosmological background evolution.

Do **not** use GNSS, wide binaries, H0, J0437, LLR, BBN or JWST to tune the remaining C1 parameters during this stage.

---

## 23. Publish the forward table

Before opening new held-out data, compute and publish:

- Earth scalar profile and source charge;
- predicted GNSS field covariance;
- predicted CODE/GFZ/PPP/MGEX projections through \(P_P\);
- wide-binary force/clock response;
- flyby amplitude/sign from the trajectory integrator;
- LLR \(\eta\) from the source-level model;
- H0 host tracer differential \(q_i\);
- J0437 predicted phase/chromaticity after EM transfer;
- optical triangle signed \(H_{\rm resid}\);
- cosmological \(a_g,\tilde a,\Phi\) background;
- CMB/recombination prediction for the chosen thermodynamic branch.

A miss after this freeze is a failure of C1. A later C2 model may be developed, but the data used to repair C1 cannot then count as a held-out confirmation of C2.

---

# Part X — What has now been completed

## 24. Closed analytically

The following gaps are now closed at formal level:

- exact universal matter metric dimensions;
- inverse and determinant;
- Einstein-frame matter stress mapping for \(E(\Phi)\);
- exact scalar matter-source equation;
- general \(K(\Phi,X)\) scalar equation;
- matter-frame conservation;
- weak-field propagator;
- exponential covariance theorem;
- conditional \(\lambda_T\leftrightarrow m_T\) mapping;
- \(\lambda_T\leftrightarrow\rho_T\) dependence;
- source-charge definition of screening;
- proof that UCD radius cannot be the sole Cassini/Vainshtein radius;
- explicit bounded least-coupling completion;
- analytic high-/low-density screening limits;
- stable nonlinear kinetic completion;
- explicit \(X\)-selective disformal completion;
- local photon-cone anisotropy;
- synchronization-curvature two-form;
- field→measurement→estimator transfer map;
- matter-frame cosmological background theorem;
- ordinary radiation-temperature theorem;
- exact statement of the cold-eternal-cosmology gap;
- correct strong-curvature completion sequence.

---

## 25. Still requires numerical solution rather than more algebra

These are no longer undefined conceptual gaps; they are computations:

1. solve the Earth/Sun boundary-value problem for C1;
2. determine whether the Earth exterior field actually produces the measured 4,200 km covariance scale;
3. compute Cassini and MICROSCOPE observables from the same solution;
4. propagate the solved field through GNSS product operators;
5. source-level LLR/INPOP or DE refit;
6. stellar pulsation/atmosphere transfer for Cepheids;
7. J0437 EM/scintillation transfer;
8. optical triangle holonomy;
9. full cosmological perturbations/high-\(\ell\) likelihood;
10. strong-curvature global solution.

These tasks can fail. That is a feature of the freeze.

---

# 26. Canonical scientific claim after this derivation

The strongest defensible statement is:

> TEP can be formulated as a universal conformal-disformal scalar theory in which environmental screening is the solved source charge of the field rather than a fitted multiplier. The formal variational and transfer structure is closed. A concrete canonical microscopic completion can be chosen with a bounded least-coupling conformal law, stable nonlinear kinetic screening and a path-selective disformal sector. That completion is now sufficiently rigid to generate held-out predictions across GNSS, Solar-System dynamics, astrophysics and cosmology. The existing corpus does not yet establish that this particular completion is correct; the next stage is numerical forward solution and pre-registered falsification. Under the present universal matter metric, the Mount Wilson mapping is a frame/background degeneracy and ordinary matter-frame radiation still obeys \(T_\gamma\propto1+z\); a genuinely non-hot eternal matter-frame cosmology requires additional derived microphysics.

---

## 27. Source basis

This derivation is built from the current TEP corpus architecture, especially:

- TEP Foundation / Jakarta v0.10
- TEP-UCD / New Delhi
- TEP-EFA / Yogyakarta
- TEP-HC / Cambridge
- TEP-C0 / Athens
- TEP-BH / Bahrain
- TEP-H0 / Kingston upon Hull
- TEP-VOID / Valencia
- TEP-J0437 / Sintra
- the 5 September 2026 corpus analysis and finalisation plan
- the prior `TEP_v0.11_CLOSURE_DERIVATION.md`

The least-coupling \(A(\Phi)\), nonlinear \(K(X)\), and \(X\)-selective \(D(\Phi,X)\) in **Canonical Completion C1** are new explicit model choices introduced here. They must not be represented as previously derived corpus results.

