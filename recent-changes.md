**Release Notes: Temporal Equivalence Principle v0.9 (Jakarta)**
*Date: June 5, 2026*

Version 0.9 (Jakarta) adds the Screening Dictionary and Domain Mapping, clarifies PPN recovery, updates cosmology to match tested implementations, and resolves cross-paper screening inconsistencies identified in the full theoretical consistency audit.

**1. Screening Dictionary and Domain Mapping (§7)**
* **Universal definition:** $\Sigma_\mu^{\text{obs}} = \mathcal S_\Sigma(\mathcal E) \, \Sigma_\mu$ with $\mathcal E = \{\rho, \Phi/c^2, \text{source structure}, \text{ambient environment}, \text{boundary conditions}, z, \text{measurement channel}\}$
* **Universal warning:** $\rho_T \approx 20$ g/cm$^3$ is not a binary on/off switch; it is the cross-scale saturation scale where the non-linear Temporal Topology response saturates
* **Domain table:** 8-row canonical map linking each TEP observational domain (GNSS, Cepheids, pulsars, wide binaries, LLR, JWST, cosmology, LHC) to its screening variable, observable response, and paper number
* **Explicit cross-scale statement:** $\rho_{\text{core}}$, $\rho_c$, and $\rho_{\text{half}}$ are different effective projections of the same non-linear response; the first-principles transfer relation remains an open derivation

**2. PPN Mapping Clarification (§7)**
* Removed incorrect claim that screening makes $d(\ln A)/d\phi \to 0$
* For $A(\phi) = \exp(\beta_A\phi/M_{\rm Pl})$, the bare derivative $\beta_A/M_{\rm Pl}$ is constant
* Screening suppresses the effective exterior charge: $\alpha_{\rm eff} = \mathcal S_\Sigma(\mathcal E) \, \alpha_0$, with $\mathcal S_\Sigma \to 0$ in dense environments
* PPN recovery to GR ($\gamma = \beta_{\text{PPN}} = 1$) follows from $\alpha_{\rm eff} \to 0$, not from a vanishing bare derivative

**3. Cosmology Update (§8)**
* Replaced "early dark energy shifts sound horizon by 1–2%" framing with CMB-preserving realization
* Early homogeneous temporal-shear mode frozen/suppressed before recombination by $S(z) = \exp[-(z/z_T)^{n_T}]$
* Sound horizon $r_s$ and acoustic morphology preserved to sub-percent precision
* Dominant observables arise in late universe through environmental and line-of-sight temporal transport
* Matches tested implementations: native hi_class MCMC (Paper 18) and Cobaya joint-likelihood (Paper 26)

**4. Cross-Paper Consistency**
* Fixed LHC (Paper 20) "high-density = maximally coupled" inconsistency — now framed as channel-specific response $\kappa_{\text{LHC}}$, not bulk-density screening
* SPIN (Paper 24) already correctly labels $\rho_c$ as many-body; no change needed

---

**Release Notes: Temporal Equivalence Principle v0.8 (Jakarta)**
*Date: April 29, 2026*

Version 0.8 (Jakarta) revises the core TEP manuscript with tightened geometric definitions and canonical nomenclature, removing residual chameleon/binary-threshold language.

**1. Canonical Definition Block (§7)**
* **Temporal Topology:** Now defined as the spatial and covariance structure of the conformal-factor field $\ln A(\phi(x))$, represented by the field configuration and its covariance $C_A(x,x') = \langle \delta\ln A(x)\,\delta\ln A(x') \rangle$. For compactness one may write $\Theta \equiv \ln A$ and $C_\Theta \equiv C_A$, but the canonical series notation remains $\ln A$, $\Sigma_\mu = \nabla_\mu \ln A$, and $C_A$. Describes how matter-clock proper-time accumulation varies across extended regions.
* **Temporal Shear:** Defined as the locally active gradient sector $\Sigma_\mu \equiv \nabla_\mu \ln A = \alpha(\phi)\nabla_\mu\phi$. In equilibrium regions where $\ln A$ is approximately constant, $\Sigma_\mu \approx 0$. Observable shear is further reduced by source/environment screening.
* **Screening Clarification:** The saturation scale $\rho_T$ is **not** a binary local-density threshold of the form $\rho > \rho_T \Rightarrow$ GR, $\rho < \rho_T \Rightarrow$ active. Recovery of GR in local tests is controlled by suppression of the observable shear/source-charge sector: $\Sigma_\mu^{\text{obs}} = \mathcal S_\Sigma(\mathcal E)\Sigma_\mu$ with $\mathcal S_\Sigma \to 0$ in screened regimes.

**2. Axiom A4 (Screening and Universality)**
* Removed "chameleon-like potential" language
* Now states: "Environmental suppression of the locally observable Temporal Shear/source-charge sector reconciles local tests with cosmological dynamics. Chameleon, Vainshtein, Galileon, DBI, and symmetron mechanisms are treated as candidate microscopic completions, not as the defining ontology of TEP."

**3. Synchronization Holonomy (§6)**
* Clarified distinction between raw loop integral $H = \oint_C \tilde{\sigma}$ and physical TEP observable $H_{\rm resid} = \oint_C (\tilde{\sigma} - \sigma_{\rm GR})$
* Emphasized $H_{\rm resid}$ is the GR-subtracted residual, not raw holonomy

**4. Response Coefficients**
* Removed $\alpha_{\text{eff}} = \alpha(\phi) \cdot f(\rho_T)$ formula
* Replaced with channel-level observable: $\Delta O_X = \kappa_X \cdot \mathcal S_X(\mathcal E) \cdot \mathcal F_X[\Delta\ln A, \Sigma_\mu, C_A; \Phi, \rho, z]$ (compactly denoted $\Delta\Theta$, $C_\Theta$ where convenient)
* Clarified $\kappa_X$ is an observable response coefficient, not the microscopic coupling $\beta$ and not a PPN coupling

**5. Glossary Updates**
* $\ln A(\phi)$: "conformal factor controlling matter-frame proper time; compact notation $\Theta \equiv \ln A$"
* $C_A$: "covariance of conformal-factor fluctuations, $C_A(x,x') = \langle \delta\ln A(x) \, \delta\ln A(x') \rangle$; also denoted $C_\Theta$ when using $\Theta \equiv \ln A$"

**Read the updated manuscript:** [https://mlsmawfield.com/tep/theory/](https://mlsmawfield.com/tep/theory/)