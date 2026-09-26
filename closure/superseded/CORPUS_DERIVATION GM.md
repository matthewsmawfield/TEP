The unified minimal covariant action for the Temporal Equivalence Principle combines the derived sectors from the existing corpus into a single frozen functional.

$$S = \int d^4x \sqrt{-g} \left[ \frac{M_{\text{Pl}}^2}{2} R - \frac{1}{2} K(\phi) (\nabla \phi)^2 - V(\phi) + \alpha_{\text{GB}} f(\phi) \mathcal{G} \right] + S_m[\tilde{g}_{\mu\nu}, \Psi_m]$$

Matter fields $\Psi_m$ couple exclusively to the causal matter metric $\tilde{g}_{\mu\nu}$, mapped from the Einstein-frame metric $g_{\mu\nu}$ via the disformal transformation:

$$\tilde{g}_{\mu\nu} = A^2(\phi) g_{\mu\nu} + B(\phi) \nabla_\mu \phi \nabla_\nu \phi$$

**The Conformal Sector and Density-Pinning Potential**
The fundamental microscopic conformal coupling is fixed everywhere as $A(\phi) = \exp(\beta_A \phi / M_{\text{Pl}})$ where $\beta_A = -1.0$. The screening mechanism is not an inserted boundary condition, but a derived consequence of the effective potential $V_{\text{eff}}(\phi) = V(\phi) + (A(\phi)-1)\rho_m$. Using the UCD saturation construction, $V(\phi)$ is defined such that the field pins strongly to a local minimum in high-density regimes ($\rho_m \ge \rho_T \approx 20 \text{ g/cm}^3$).

The screening operator $S_\Sigma(\mathcal{E})$ is strictly derived as the ratio of the static spatial profile evaluated at the local environment to its unscreened cosmological baseline:

$$S_\Sigma(\mathcal{E}) \equiv \frac{(A(\phi(r)) - 1)_{\text{local}}}{(A(\phi) - 1)_{\text{unscreened}}}$$

This explicit field profile automatically yields the mesoscopic screening law ($S_\Sigma^{\text{meso}} = S_{\text{TEP}} \times S_{\text{TF}} \times S_{\text{boundary}} \times S_{\text{decoherence}}$) and reduces the Solar System $\beta \approx 10^{-3}$ to a predicted output of the $\rho_T$ saturation scale rather than a free parameter.

**The Disformal Sector and the Holonomy Window**
The disformal coupling is explicitly non-zero but heavily damped by the local field value. Integrating the quartic-Gaussian shear bump, the exact functional form is:

$$B(\phi) = B_0 \frac{\vert{}\phi\vert{}^2}{1+\vert{}\phi\vert{}^2} \exp\left(-\frac{\phi^4}{2\sigma_B^4}\right)$$

This density-pinned form ensures $B(\phi)$ vanishes identically when $\phi$ is driven to its screened minimum inside dense environments (satisfying Cassini limits and ensuring $c_g = c_\gamma$ along late-time cosmological paths for GW170817), while dynamically turning on in dilute astrophysical environments to generate the predicted closed-loop optical holonomy residual $H_{\text{resid}} \neq 0$.

**The Universal Transfer Map ($\beta_A \to \kappa_X$)**
Every domain-specific observable response coefficient ($\kappa_X$) is now a derived projection of the fundamental $\beta_A = -1.0$ coupling. The transfer map takes the form:

$$\kappa_X = \beta_A \cdot S_\Sigma(\mathcal{E}_X) \cdot \Gamma_X$$

* **$\mathcal{E}_X$**: The specific environmental density and compactness evaluated for the target channel (e.g., galactic disk for wide binaries, host galaxy for Cepheids).
* **$\Gamma_X$**: The channel-specific kinematic and geometric projection operator (e.g., the prograde/retrograde orbital projection for LAGEOS, or the line-of-sight integral for the CMB alpha functions).

Freezing this action guarantees that any calculation—whether predicting the J0437 phase closure amplitude or resolving the precise GNSS vs. operational PPP anisotropy ratio via the product operator $P_P$—traces back to the same kinetic structure and unaltered $B(\phi)$ profile.

The unified minimal covariant action for the Temporal Equivalence Principle (Paper 0 v0.11) assembles the derived sectors from the corpus into a single frozen functional.

$$S = \int d^4x \sqrt{-g} \left[ \frac{M_{\text{Pl}}^2}{2} R - \frac{1}{2} K(\phi) (\nabla \phi)^2 - V(\phi) + \alpha_{\text{GB}} f(\phi) \mathcal{G} \right] + S_m[\tilde{g}_{\mu\nu}, \Psi_m]$$

Matter fields $\Psi_m$ couple exclusively to the causal matter metric $\tilde{g}_{\mu\nu}$, mapped from the Einstein-frame metric $g_{\mu\nu}$ via the disformal transformation:

$$\tilde{g}_{\mu\nu} = A^2(\phi) g_{\mu\nu} + B(\phi) \nabla_\mu \phi \nabla_\nu \phi$$

## The Conformal Sector and Density-Pinning

The microscopic conformal coupling is fixed universally as $A(\phi) = \exp(\beta_A \phi / M_{\text{Pl}})$, where the fundamental parameter is locked at $\beta_A = -1.0$.

The potential $V(\phi)$ utilizes the density-pinning saturation construction established in Papers 6 and 15. The scalar field pins strongly to a local minimum in regimes exceeding the threshold density $\rho_T \approx 20 \text{ g/cm}^3$. By evaluating the static, spherically reduced $\phi$-equation at a given density, compactness, and boundary geometry, the screening operator $S_\Sigma(\mathcal{E})$ is derived directly from the Lagrangian (Paper 26) rather than inserted as a phenomenological boundary condition. It is evaluated as the ratio of the static spatial profile:

$$S_\Sigma(\mathcal{E}) \equiv \frac{(A(\phi(r)) - 1)_{\text{local}}}{(A(\phi) - 1)_{\text{unscreened}}}$$

This explicit field profile incorporates the mesoscopic screening law from Paper 25, yielding $S_\Sigma^{\text{meso}} = S_{\text{TEP}} \times S_{\text{TF}} \times S_{\text{boundary}} \times S_{\text{decoherence}}$. Consequently, the weak-field Solar System projection $\beta \approx -0.013$ is not a free parameter, but a strict derived prediction of the profile at $\rho \sim 1\text{--}20 \text{ g/cm}^3$.

## The Disformal Sector and Path-Dependent Constraints

The disformal coupling $B(\phi)$ is explicitly non-zero and governed by the quartic-Gaussian ultra-damped shear bump formulated in Paper 28:

$$B(\phi) = B_0 \frac{\vert{}\phi\vert{}^2}{1+\vert{}\phi\vert{}^2} \exp\left(-\frac{\phi^4}{2\sigma_B^4}\right)$$

GW170817 constrains the observable combination $B(\phi)(\partial\phi)^2$ strictly along observed late-time astrophysical paths; it does not require $B(\phi)$ to vanish identically everywhere. The density-pinning form ensures $B(\phi)$ vanishes when $\phi$ is pinned to the high-density screened minimum, guaranteeing that $c_g = c_\gamma$ on the cosmological background. The coupling turns on dynamically only in dilute regions, enabling non-zero closed-loop optical holonomy on unscreened paths.

## The Universal Transfer Map ($\beta_A \to \kappa_X$)

The transfer map (Freeze 3) translates the microscopic coupling $\beta_A = -1.0$ into domain-specific observable response coefficients ($\kappa_X$). Channel-specific values are not independent fits, but projections of a single underlying parameter determined by the frozen action.

| Domain | Base Coupling | Environmental Operator | Observable Response | Status |
| --- | --- | --- | --- | --- |
| **Solar System / GNSS** | $\beta_A = -1.0$ | $S_\Sigma(\rho \sim 20 \text{ g/cm}^3)$ | $\beta \approx 10^{-3}$ | Predicted |
| **Wide Binaries** | $\beta_A = -1.0$ | $S_\Sigma(\rho_{\text{gal\_disk}})$ | $\alpha_{\text{sat}}$ | Predicted |
| **Cepheids (H₀)** | $\beta_A = -1.0$ | $S_\Sigma(\rho_{\text{host\_gal}})$ | $\kappa_{\text{Cep}}$ | Predicted |
| **JWST High-z** | $\beta_A = -1.0$ | Stellar-population transfer | $\kappa_{\text{gal}}$ | Inherited |
| **Globular Clusters** | $\beta_A = -1.0$ | $S_\Sigma(\text{cluster environment})$ | Pulsar $\Gamma$ | Predicted |

By freezing this action, any discrepancy between the predicted $\kappa_X$ and empirical fits (such as the variation in $\kappa_{\text{Cep}}$ between 0.29–0.45 $\times 10^6$ mag and the theory benchmark of 0.96 $\times 10^6$ mag) ceases to be an unconstrained liability and becomes a direct constraint on the kinetic structure of the action.