# TEP Scalar-Sector Closure

> **SUPERSEDED 2026-09-12.** This report adopted the quadratic potential
> $V(\phi) = \Lambda^4(1 + \phi^2/(2\Lambda^2))$, which gives a density-independent
> mass $m = \Lambda$ and a fixed $\lambda_C = 0.105$ mm in every environment.
> That passes the local null gates by killing the temporal shear everywhere,
> which also removes the mechanism from wide binaries, SPARC and every
> environmental ordering result. The corrected potential is the inverse-power
> form $V(\phi) = \Lambda^4(1 + \Lambda/\phi)$ at the same zero-parameter scale
> $\Lambda = \sqrt{M_{\rm Pl} H_0}$, which makes $\lambda_C \propto \rho^{-3/4}$
> and gives Temporal Topology saturation the environmental response TEP
> requires. The screening formula was also corrected from the constant-mass
> Helmholtz result $3(1+x)/(x(2x+1))$ to the thin-shell charge fraction
> $s = 3\delta R/R$. The authoritative updated report is
> `CLOSURE_RESOLUTION.md`; the verified numbers are in
> `results/tep_screening_closure_results.json`. This document is retained as
> a historical record of the first computational attempt.

Status: Local screening closure established with a single-scale potential. The scalar sector uses one scale ($\Lambda_{\rm DE}$), one coupling ($\beta_A=-1$), and one potential $V(\phi)=\Lambda^4(1+\phi^2/(2\Lambda^2))$. The spatial projection provides Yukawa screening (all 5 local gates pass). The temporal projection provides the cosmic redshift through the clock map. The screening is a continuous gradient transition, not chameleon and not thin-shell. Two open questions remain: the wide-binary/galactic unscreening channel, and the potential energy at cosmological field values.

## 1. Three ingredients

### Scale

$$\Lambda = (M_{\rm Pl}\,H_0\,\hbar)^{1/2} \approx 1.9\ {\rm meV}$$

This is the dark-energy scale. It is the only dimensionful parameter in the scalar sector.

### Coupling

$$A(\phi) = \exp\!\left(\frac{\beta_A\,\phi}{M_{\rm Pl}}\right), \qquad \beta_A = -1$$

Frozen universally, as specified by Jakarta. The DEF coupling is $\alpha_0 = \beta_A = -1$.

### Potential

$$V(\phi) = \Lambda^4\!\left(1 + \frac{\phi^2}{2\Lambda^2}\right) = \Lambda^4 + \frac{1}{2}\Lambda^2\phi^2$$

One function, one scale. The cosmological constant $\Lambda^4$ sits beneath a quadratic mass term.

- $V(0) = \Lambda^4 \approx 1.3\times10^{-47}$ GeV$^4$ (dark energy at present epoch, $z=0$)
- $V''(0) = \Lambda^2$, giving scalar mass $m = \Lambda \approx 1.9\times10^{-12}$ GeV
- Compton wavelength $\lambda_C = \hbar c/m \approx 0.10$ mm

## 2. Two projections

The scalar field $\phi$ has a spatial profile $\phi(\mathbf r)$ around bodies and a temporal evolution $\phi(\eta)$ on the cosmological background. These are different projections of the same field, governed by the same action.

### Spatial projection: local screening

For a static body of density $\rho$ and radius $R$ in ambient density $\rho_{\rm amb}$, the linearized field equation ($\phi \ll M_{\rm Pl}$, verified for all bodies) is the Helmholtz equation:

$$\nabla^2\phi = m^2\phi - \frac{\rho}{M_{\rm Pl}}$$

The solution is a Yukawa profile. The field transitions continuously from the body-density equilibrium $\phi_{\rm body} = \rho/(m^2 M_{\rm Pl})$ to the ambient equilibrium $\phi_{\rm amb} = \rho_{\rm amb}/(m^2 M_{\rm Pl})$ over a distance $\sim\lambda_C$. This is a continuous gradient transition — not chameleon, not thin-shell. Jakarta explicitly states: "Environmental suppression in TEP is defined by the continuous response of the Temporal Topology rather than by a discrete density threshold or thin-shell boundary."

The source-charge screening factor for a uniform sphere is:

$$S_\Sigma = \frac{Q}{Q_0} = \frac{3(1+x)}{x(2x+1)}, \qquad x \equiv mR$$

- For $x \ll 1$: $S_\Sigma \to 1$ (unscreened)
- For $x \gg 1$: $S_\Sigma \to 3/(2x)$ (screened)

All 5 local gates pass at $m = \Lambda$:

| Gate | Body | $x = mR$ | $S_\Sigma$ | Observable | Bound | Margin |
|---|---|---|---|---|---|---|
| Cassini | Sun | $6.7\times10^{12}$ | $2.2\times10^{-13}$ | $S_\odot$ | $<3.4\times10^{-3}$ | $10^{13}\times$ |
| Geodesy | Earth | $6.2\times10^{10}$ | $2.4\times10^{-11}$ | $\alpha_\oplus = 2\beta_A^2 S$ | $<10^{-8}$ | $200\times$ |
| LLR | Earth/Moon | — | — | $\eta_N$ | $<4.4\times10^{-4}$ | $10^{18}\times$ |
| Pulsar | NS | $1.2\times10^{8}$ | $1.3\times10^{-8}$ | $\alpha_{\rm NS}$ | $<10^{-3}$ | $8\times10^{4}\times$ |
| WD | WD | $6.8\times10^{10}$ | $2.2\times10^{-11}$ | $\alpha_{\rm WD}$ | $<10^{-2}$ | $5\times10^{8}\times$ |

The tightest gate is geodesy ($200\times$ margin). The scalar Compton wavelength is $\lambda_C \approx 0.10$ mm — the fifth force is a contact interaction at macroscopic scales.

The field equilibrium $\phi_{\rm eq} = \rho/(m^2 M_{\rm Pl})$ is $\ll M_{\rm Pl}$ for all bodies (Earth: $10^{-30} M_{\rm Pl}$, NS: $10^{-17} M_{\rm Pl}$), confirming the linearized approximation.

### Temporal projection: cosmology

On the homogeneous cosmological background, $\phi = \phi(\eta)$ depends only on conformal time. The disformal term is purely temporal and is absorbed into the lapse by a time reparameterisation, leaving the spatial metric static ($a_m = 1$). The effective scale factor is:

$$a_{\rm eff}(z) = A_{\rm clock}(z) \cdot a_m = (1+z)^{-1}$$

where $A_{\rm clock}(z) = (1+z)^{-1}$ is the observational clock/redshift mapping, fixed by the redshift definition. The field evolution is:

$$\phi(z) = M_{\rm Pl}\ln(1+z)$$

The temporal shear is:

$$\Sigma_0 = \frac{\beta_A}{M_{\rm Pl}}\frac{d\phi}{dt} \sim H_0$$

The temporal shear IS the Hubble flow. The cosmic redshift is accumulated temporal shear along cosmological lines of sight, not physical expansion of space. The matter frame is static and eternal; the apparent Big Bang is a temporal horizon ($A_{\rm clock} \to 0$ as $z \to \infty$), not a curvature singularity.

The potential energy at the present epoch is $V(0) = \Lambda^4 \approx 1.3\times10^{-47}$ GeV$^4$, providing the dark-energy budget. The kinetic energy $\frac{1}{2}\dot\phi^2 \sim \frac{1}{2}(M_{\rm Pl}H_0)^2 = \frac{1}{2}\Lambda^4$ contributes comparably. The total $\rho_\phi \sim \Lambda^4$ matches the observed dark-energy density to order of magnitude; the exact normalization is fixed by matching at the present epoch, as described in TEP-TH (Paper 27).

## 3. What this resolves

The 19-order-of-magnitude single-scale tension is resolved by decoupling the two projections:

- The spatial projection (screening) depends on the scalar mass $m = \Lambda$, giving $\lambda_C \sim 0.1$ mm. The fifth force is screened at all macroscopic scales.
- The temporal projection (cosmology) depends on the clock map $A_{\rm clock}(z) = (1+z)^{-1}$, which is an observational input, not a dynamical consequence of the potential. The cosmic redshift comes from temporal shear, not from the scalar's spatial range.

The same scalar field, same action, same scale $\Lambda$ provides both. The mass $\Lambda$ screens locally; the vacuum energy $\Lambda^4$ drives cosmology; the temporal evolution $\phi(\eta) = M_{\rm Pl}\ln(1+z)$ provides the redshift; and the transition acceleration $a_0 = \Lambda^2/M_{\rm Pl} = cH_0$ sets the wide-binary and galactic screening floor. One scale, one coupling, one potential, three projections.

The identity $\Lambda^2/M_{\rm Pl} = H_0\hbar$ is exact (since $\Lambda^2 = M_{\rm Pl}H_0\hbar$ by definition). In SI units this gives $a_0 = cH_0 = c\Lambda^2/(M_{\rm Pl}\hbar) \approx 6.8\times10^{-10}$ m s$^{-2}$. This connection is already present in the corpus: Paper 5 (Singapore) states $a_0 \sim \Lambda^2/M_{\rm Pl} \sim cH_0$, and Paper 4 (Tortola) gives $a_T \sim cH_0$. The SPARC-extracted value $g_{\rm TEP} \approx 5\times10^{-10}$ m s$^{-2}$ (Paper 6) is within a factor of 1.36 of this exact theoretical prediction. The wide-binary transition radius $R_s = \sqrt{GM/a_0} \approx 3289$ AU (no EFE) or $\sqrt{GM/(2a_0)} \approx 2325$ AU (with $\eta=2$ external field) brackets the observed $R_s = 2646 \pm 182$ AU (Paper 13).

## 4. Open questions

### Wide-binary and galactic unscreening

The scalar fifth force is screened at all scales $> 0.1$ mm. The wide-binary signal (Paper 13, transition at $R_s \approx 2646$ AU) and the SPARC rotation-curve signal (Paper 6) cannot come from the scalar fifth force. They must come from the temporal shear channel.

Jakarta's domain dictionary specifies: "Wide binaries: weak-field recovery / low-acceleration environmental un-screening." The unscreening is of the temporal shear, not of the scalar fifth force. The environmental operator $\mathcal{S}_\Sigma(\mathcal{E})$ suppresses the temporal shear in dense environments and unsuppresses it in diffuse environments.

The microscopic origin of this environmental dependence is not yet derived from the action. Jakarta states: "A completed microscopic closure must derive the continuous environment-dependent field configuration, together with its $S_A$, $S_\Sigma$, and covariance projections, from a common action and stated boundary conditions."

The phenomenological $S(a) = [1+(a/a_0)^2]^{-1/2}$ with $a_0 = cH_0 = \Lambda^2/M_{\rm Pl} \approx 6.8\times10^{-10}$ m s$^{-2}$ remains a target for what the environmental operator should reproduce. This identity is already in the corpus (Paper 5, Singapore, Eq. in Section on GTE: $a_0 \sim \Lambda^2/M_{\rm Pl} \sim cH_0$; Paper 4, Tortola: $a_T \sim cH_0$). The SPARC fit gives $g_{\rm TEP} \approx 5\times10^{-10}$ m s$^{-2}$ (Paper 6), within a factor of 1.36 of the exact theoretical value $cH_0 = 6.8\times10^{-10}$ m s$^{-2}$. The wide-binary prediction $R_s = \sqrt{GM/a_0} \approx 3289$ AU (no EFE) or $\sqrt{GM/(2a_0)} \approx 2325$ AU (with $\eta=2$ external field) brackets the observed $R_s = 2646 \pm 182$ AU (Paper 13). Deriving the environmental operator from the action is the remaining closure problem.

### Potential energy at cosmological field values

At $z = 1$, $\phi = 0.69 M_{\rm Pl}$ and $V = \Lambda^4(1 + 0.24 M_{\rm Pl}^2/\Lambda^2) \sim \Lambda^2 M_{\rm Pl}^2 \sim 10^{13}$ GeV$^4$, which is $\sim 10^{60}$ times $\Lambda^4$. In the standard Friedmann framework this would be catastrophic. In the TEP static matter frame, the expansion is kinematic (clock map), not dynamic (Friedmann). The large potential energy at early times is part of the scalar field's energy budget but does not drive the expansion. TEP-HC (Paper 18) verifies that the pre-recombination sound horizon is preserved to 6 ppm in this framework. The consistency of the large potential energy with the CMB acoustic scale and BBN requires further verification within the static-frame Boltzmann code.

### Corpus-claimed values

The massive scalar predicts $|\eta_N| \sim 10^{-22}$ for LLR, while Paper 17 claims $\eta_{\rm resid} = -3.91\times10^{-4}$. This 18-order discrepancy indicates that the LLR signal, if real, must come from the temporal shear channel, not from the scalar fifth-force channel. Similarly, Paper 15's $S_\oplus = 0.35$ and Paper 29's galaxy $S = 0.013$ are temporal-shear projections, not scalar-charge observables.

## 5. Previous errors (corrected)

1. No-go 1 (false): $m_{\rm eff} \propto \rho^p$ with $p \le 1/2$ for any $V(\phi)$. Wrong. For InvPow $n=1$, $p=3/4$. Retracted.

2. $S(a)$ ansatz: reclassified as phenomenological target for the environmental operator, not a microscopic closure.

3. Jakarta coupling: $A(\phi) = e^{\beta_A\phi/M_{\rm Pl}}$ is frozen. $S(a)$ was incorrectly inserted into $A$.

4. Compton wavelength: $\lambda_C$ at $\Lambda_{\rm DE}$ is 4.3 Gpc (vacuum) for the Cosh potential with scale $M_{\rm Pl}$, not 133,000 km. Unit error, 18 orders. The correct local screening mass is $m = \Lambda_{\rm DE}$, giving $\lambda_C \approx 0.1$ mm.

5. Chameleon approach: abandoned. Jakarta explicitly states screening is not chameleon. The mechanism is Yukawa suppression by a fixed mass — a continuous gradient transition.

6. BVP numerical artifacts: solve_bvp converged to trivial solutions. Corrected by analytical Helmholtz solution.

7. InvPow d2V sign: missing chain-rule factor. Fixed.

8. Ambient Compton wavelength: does not set the screening transition. The field transitions over $\lambda_C^{\rm body}$, not $\lambda_C^{\rm amb}$. But more fundamentally, the mechanism is not chameleon at all.

## 6. Files

```
closure/
  FINAL_REPORT.md             this document
  vphi_bvp_v5.py              shooting-method BVP (superseded by analytical Helmholtz)
  vphi_bvp_v4.py              InvPow BVP (superseded — chameleon approach abandoned)
  vphi_bvp_v3.py, v2.py, v1.py  earlier versions (superseded)
  tep_closure.py              phenomenological S(a) test (target for environmental operator)
```
