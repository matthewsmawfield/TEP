# TEP Corpus Review: Complete Derivations

> **PARTIALLY SUPERSEDED 2026-09-13 (v5.9).** The BBN amplification mechanism
> (D11, "no candidate mechanism") is resolved by the cuscuton gradient floor
> `μ²/(2K)` (CORPUS_PLAN.md §7.8 Target 1). The "Gate 4B remains open" statement
> is withdrawn. The authoritative account is `CORPUS_PLAN.md` §7 (v5.9).

Every quantitative claim in the corrections document, derived from stated inputs. Nothing is copied from the papers except the inputs, which are labelled with paper and line number. Where a paper's own figure is quoted it appears next to the recomputed value so the two can be compared.

The accompanying script `tep_derivations.py` reproduces all of it in one run.

**Constants used throughout**

| Symbol | Value | Note |
|---|---|---|
| G | 6.674 × 10⁻⁸ cm³ g⁻¹ s⁻² | |
| c | 2.99792458 × 10¹⁰ cm/s | |
| M☉ | 1.989 × 10³³ g | |
| pc | 3.0857 × 10¹⁸ cm | |
| M⊕ | 6.0 × 10²⁷ g | Paper 6, line 221 |
| ρ_T | 20 g/cm³ | Paper 6, line 29 |

---

## D1 — ρ_T from the GNSS correlation length

**Setup.** Paper 6 Eq. (rt_definition):

$$R_T = \left(\frac{3M}{4\pi\rho_T}\right)^{1/3}$$

Inverting for the density:

$$\rho_T = \frac{3M}{4\pi R_T^3}$$

Paper 6 §2 identifies $R_T(M_\oplus)$ with the GNSS correlation length $L_c$, so

$$\rho_T(L_c) = \frac{3M_\oplus}{4\pi L_c^3}$$

**Sensitivity.** Taking logs and differentiating:

$$\frac{d\ln\rho_T}{d\ln L_c} = -3$$

A factor $f$ error in $L_c$ is a factor $f^3$ error in $\rho_T$. This is the reason the whole calibration is fragile.

**Evaluation.**

| L_c (km) | Source | ρ_T (g/cm³) | ×20 |
|---|---|---|---|
| 4,200 | P6 adopted (line 259) | 19.33 | 0.97 |
| 4,152.8 | exact R_T(M⊕) at ρ_T = 20 | 20.00 | 1.00 |
| 4,700 | P6 stated extreme (line 233) | 13.80 | 0.69 |
| 3,700 | P6 stated extreme (line 233) | 28.28 | 1.41 |
| 4,549 | P1 upper (line 15) | 15.22 | 0.76 |
| 3,330 | P1 lower (line 15) | 38.79 | 1.94 |
| 4,201 | P2 central (line 15) | 19.32 | 0.97 |
| 6,168 | P2 +1σ | 6.10 | 0.31 |
| 2,234 | P2 −1σ | 128.47 | 6.42 |
| 4,767 | P3 ionofree 2024 (line 55) | 13.22 | 0.66 |
| **1,396** | **P14 MGEX (line 34)** | **526.51** | **26.3** |
| 1,100 | P3 MSC upper (line 435) | 1,076.18 | 53.8 |
| 1,069 | P3 ionofree best (line 51) | 1,172.55 | 58.6 |
| 700 | P3 MSC lower (line 435) | 4,176.08 | 208.8 |

**Checks against the paper.** R_T(M⊕) at ρ_T = 20 gives 4,152.8 km; P6 line 279 states 4,146 km. P6's stated extremes reproduce: L = 3,700 → 28.3 (paper says ~30), L = 4,700 → 13.8 (paper says ~14). Paper 6's arithmetic is correct.

**Uncertainty propagation.** P6 line 227 adopts ±12% on L_c. Since 1.12³ = 1.40, that gives ρ_T = 20 ± 7 g/cm³ (≈ ±36%), which is internally consistent.

But the underlying measurement is P2's 4,201 ± 1,967 km, i.e. ±47%. Propagating that instead:

$$L_c \in [2{,}234,\ 6{,}168]\ \text{km} \;\Longrightarrow\; \rho_T \in [6.1,\ 128.5]\ \text{g/cm}^3$$

a dynamic range of **21×**, not ±36%.

P6 justifies preferring ±12% on the grounds that the centres "share largely overlapping underlying data." That argument establishes ±12% as a *floor* on the systematic. It does not license replacing the ±47% statistical error with it.

**The ρ_T ≈ 20 g/cm³ ≈ "plausible terrestrial-material density" argument does not survive any L_c below about 2,500 km.**

---

## D2 — Paper 0's actual forecast versus the band cited to it

**Jakarta v0.10, line 688:** λ_T ~ 2,000–3,000 km for viable screening parameters.

**Cited to Jakarta as "1,000–10,000 km" by:** P1 (lines 42, 312, 343, 364, 1405, 2773), P2 (15, 330), P5 (309, 409, 949, 1834), P9 (120). P1 line 343 cites Jakarta's DOI directly. The band appears nowhere in Jakarta v0.10.

**Test against the actual forecast:**

| Value (km) | Source | In 2,000–3,000? |
|---|---|---|
| 3,330 | P1 lower | OUT |
| 4,549 | P1 upper | OUT |
| 4,201 | P2 central | OUT |
| 700 | P3 MSC lower | OUT |
| 1,100 | P3 MSC upper | OUT |
| 1,069 | P3 ionofree | OUT |
| 3,485 | P3 phase-alignment upper | OUT |
| 1,396 | P14 MGEX | OUT |
| 4,200 | P6 adopted | OUT |

Nine of nine central values fall outside.

**The second forecast, Jakarta line 690:** "<5% variation in fitted parameters" across multi-centre cross-validation.

Observed, P1 line 15: 3,330–4,549 km.

$$\text{midpoint} = \frac{3330+4549}{2} = 3{,}940\ \text{km}$$
$$\text{half-range} = \frac{4549-3330}{2 \times 3940} = \pm 15.5\%, \qquad \text{full range} = 30.9\%$$
$$\frac{15.5}{5} = 3.1\times \text{ the stated tolerance}$$

This criterion is never cited or tested in P1, P2, P3, P5 or P14.

---

## D3 — "r = −0.888, 5.1σ"

**Where the 5.1σ comes from.** P2 ran 5,000,000 phase-randomised surrogates and found 0 exceedances. With 0 out of N exceedances the tightest attainable bound is

$$p < \frac{1}{N} = \frac{1}{5\times10^6} = 2\times10^{-7}$$

Converting the bound to a sigma:

$$\Phi^{-1}(1 - 2\times10^{-7}) = 5.07\sigma \ \text{(one-sided)}, \qquad 5.20\sigma \ \text{(two-sided)}$$

P5 line 29 quotes 5.1σ. **This is the Monte Carlo resolution floor.** It states how many surrogates were run, not how strong the signal is. With 50M surrogates and the same zero exceedances it would "become" 5.49σ, with no new data.

**P2's own parametric test (line 805).** With r = −0.888 and N_eff ≈ 11:

$$SE = \sqrt{\frac{1-r^2}{N_{\rm eff}-2}} = \sqrt{\frac{1-0.7885}{9}} = 0.1533$$
$$t = \frac{r}{SE} = -5.793 \ \text{on } 9 \ \text{dof} \;\Longrightarrow\; p = 2.62\times10^{-4} = \mathbf{3.65\sigma}$$

**How many points would 5.07σ actually require?** Solving $|r|/\sqrt{(1-r^2)/(n-2)}$ against the t-distribution:

$$n = 19$$

P2 has N_eff ≈ 11.

P2 v0.19 handles this correctly — it reports both p-values, states they test different nulls, and calls the surrogate result a bound. P5, the synthesis paper, drops the caveats and presents the bound as a measurement.

---

## D4 — The CMB-alignment variance ratio

**Inputs (P2 line 1380):** R²(best-fit direction) = 55.7%, R²(exact CMB dipole) = 35.7%, R²(Solar Apex) = 0.0054%. P5 rounds the last to 0.01%.

$$\frac{55.7}{0.0054} = 10{,}315 \qquad \text{P2 quotes} \sim 10{,}300$$
$$\frac{35.7}{0.0054} = 6{,}611 \qquad \text{P2 quotes} \sim 6{,}600$$
$$\frac{55.7}{0.01} = 5{,}570 \qquad \text{P5 quotes } 5{,}570$$

The 1.85× spread between P2 and P5 is entirely the rounding of the denominator. A ratio whose denominator is statistically indistinguishable from zero is unbounded and is not an evidence measure.

P2's substantive argument survives deletion of the ratio: the Solar Apex vector (+30° Dec) predicts N–S anisotropy while the observation is E–W (EW/NS = 2.16). That is a geometric incompatibility and needs no ratio.

---

## D5 — κ_gal: the benchmark against every measurement of it

**Provenance (P12 line 311, the paper's own words).** κ_gal = 0.960 × 10⁶ mag is the value that produces "a ~0.3 mag characteristic modulus shift" through the Leavitt slope b ≈ −3.30. That is the shift required to close the Hubble tension. The benchmark is reverse-derived from the target.

**Failure threshold (P12 line 770).** At κ ≤ 0.56 × 10⁶ the corrected Red Monster SFE rises to ~0.26, above the ΛCDM limit of 0.20.

**Tension of each estimate against the benchmark**, computed as (0.960 − v)/σ_v:

| Estimate (×10⁶ mag) | Value | σ below 0.960 | Clears 0.56? |
|---|---|---|---|
| P11 Step 44, σ_v = 182.1 | 0.326 ± 0.206 | 3.08 | **No** |
| P11 Step 44, σ_v = 150 | 0.290 ± 0.227 | 2.95 | **No** |
| P11 Step 42 restricted | 0.369 ± 0.310 | 1.91 | **No** |
| P11 redshift-only, σ_v = 150 | 0.452 ± 0.220 | 2.31 | **No** |
| P11 redshift-only, σ_v = 182 | 0.417 ± 0.249 | 2.18 | **No** |
| P11 mixed model (both free) | 0.116 ± 0.377 | 2.24 | **No** |
| P12 JWST-side recovery | 0.600 ± 0.380 | 0.95 | Yes, marginally |
| P12 cited as P11 S44 (**stale**) | 0.400 ± 0.270 | 2.07 | **No** |

Every independent estimate fails the paper's own threshold. The single value that clears it is recovered from the JWST data being explained, so it is not independent evidence.

**Stale citation.** P12 lines 311, 752 and 1156 attribute 0.400 ± 0.270 to P11 Step 44. P11 v0.9 lines 17 and 304 give Step 44 as 0.326 ± 0.206 (σ_v = 182.1) or 0.290 ± 0.227 (σ_v = 150). Neither matches.

---

## D6 — Paper 7: the R_T = R_S crossover

**Scaling.** The two radii scale differently in mass:

$$R_T(M) = \left(\frac{3M}{4\pi\rho_T}\right)^{1/3} \propto M^{1/3}, \qquad R_S(M) = \frac{2GM}{c^2} \propto M$$
$$\Longrightarrow \frac{R_T}{R_S} \propto M^{-2/3}$$

So above some mass the TEP saturation radius falls inside the event horizon.

**Solving for the crossover.** Set R_T = R_S and cube both sides:

$$\frac{3M}{4\pi\rho_T} = \frac{8G^3M^3}{c^6}$$
$$\frac{3c^6}{4\pi\rho_T} = 8G^3M^2$$
$$\boxed{M_{\rm crit} = c^3\sqrt{\frac{3}{32\pi G^3 \rho_T}}}$$

Evaluating at ρ_T = 20 g/cm³:

$$M_{\rm crit} = 6.036\times10^{40}\ \text{g} = \mathbf{3.03\times10^7\ M_\odot}$$

**Check against Paper 7's own numbers.** For RBH-1 at M = 2 × 10⁷ M☉:

| Quantity | Computed | P7 abstract |
|---|---|---|
| R_T | 7.80 × 10⁷ km | 7.8 × 10⁷ km ✓ |
| R_S | 5.91 × 10⁷ km | — |
| R_T/R_S | 1.32 | ~1.3 ✓ |

The observed wake is 62 kpc = 1.91 × 10¹⁸ km, so

$$\frac{\text{wake}}{R_T} = 2.45\times10^{10}$$

The calibration lands 24 billion times short of the structure it is invoked to explain. P7's abstract says so correctly ("not the full 62 kpc wake length"); P7 line 591 contradicts it.

**The consequence, stated nowhere in the corpus.** Above M ≈ 3.03 × 10⁷ M☉ the TEP saturation radius lies inside the event horizon and the framework predicts nothing observable. RBH-1 sits a factor of 1.52 below the crossover, so R_T ≈ 1.3 R_S is a coincidence of that particular mass rather than a prediction. Sgr A* (4 × 10⁶ M☉) is below; M87* (6.5 × 10⁹ M☉) is 214× above.

This is a genuine falsifiable boundary of the framework and should be stated in Papers 7 and 28.

---

## D7 — Paper 26's master operator against Paper 6's ρ_T

**The master form (P26 line 96):**

$$\mathcal{S}_\Sigma(\mathcal{E}) = \left[1 + \left(\frac{\Sigma_\mu\Sigma^\mu}{g_t^2}\right)^n + \left(\frac{\rho}{\rho_{\rm half}}\right)^2\right]^{-1}, \qquad \rho_{\rm half} \approx 0.5\ M_\odot/\text{pc}^3$$

**Unit conversion:**

$$\rho_{\rm half} = \frac{0.5 \times 1.989\times10^{33}\ \text{g}}{(3.0857\times10^{18}\ \text{cm})^3} = 3.385\times10^{-23}\ \text{g/cm}^3$$

**Against Paper 6:**

$$\frac{\rho_T}{\rho_{\rm half}} = \frac{20}{3.385\times10^{-23}} = 5.91\times10^{23} \qquad (\mathbf{23.8\ \text{orders of magnitude}})$$

**Evaluating the density term at ρ = ρ_T:**

$$\left(\frac{\rho_T}{\rho_{\rm half}}\right)^2 = 3.49\times10^{47} \;\Longrightarrow\; \mathcal{S}_\Sigma(\rho_T) = 2.86\times10^{-48}$$

**ρ_T plays no role whatsoever in the corpus's own master screening operator**, despite anchoring Papers 6, 7, 12, 13, 15 and 23.

Physically the two are different kinds of object: ρ_half sits at a galactic ambient density (solar neighbourhood ≈ 0.1 M☉/pc³), ρ_T at a condensed-matter density. No paper in the corpus reconciles them.

**Note on a superseded finding.** My earlier assessment that no master functional form existed anywhere is wrong and should not be repeated. Paper 26 §2 constructs one, and line 91 states the intent correctly: "a complete theory must supply a single covariant realization of this operator, not a patchwork of scale-specific proxies." The problem is not absence — it is that the form is not cited by Jakarta, is contradicted within Paper 29 line 159 ("No master functional form of ℰ has been written down," in the same paragraph that quotes it), and has not been shown to reproduce any of the domain projections built on ρ_T.

---

## D8 — Paper 17: the blind hold-out cannot yield 14.3σ

**Inputs.**

| | η_resid | SE | σ | N |
|---|---|---|---|---|
| Headline (line 15) | −3.91 × 10⁻⁴ | 5.63 × 10⁻⁵ | 6.94 | 25,445 |
| Hold-out (line 908) | −3.66 × 10⁻⁴ | 2.56 × 10⁻⁵ | 14.30 | 20 × 20% splits |

The hold-out SE is **2.20× smaller** than the full-sample SE while each split uses 20% of the data. Under any correct combination this is impossible.

**Reconstruction.** Scaling a single 20% split as 1/√N from the full-sample SE:

$$SE_{\rm single} = \frac{5.63\times10^{-5}}{\sqrt{0.20}} = 1.259\times10^{-4} \;\Longrightarrow\; 2.91\sigma$$

Inverse-variance combining k such splits **as though independent** gives $SE_k = SE_{\rm single}/\sqrt{k}$:

| k | SE | σ | |
|---|---|---|---|
| 5 | 5.63 × 10⁻⁵ | **6.50** | maximum number of *disjoint* 20% partitions |
| 10 | 3.98 × 10⁻⁵ | 9.19 | |
| 20 | 2.82 × 10⁻⁵ | **13.00** | matches the reported 14.3σ |

**The reuse factor.** Each calendar year appears in 20 × 0.20 = 4 of the twenty hold-out sets. The number of genuinely independent replicates available at a 20% hold-out fraction is 1/0.20 = 5. Treating 20 overlapping splits as independent inflates the significance by

$$\sqrt{\frac{20}{5}} = 2.00\times$$

**The sanity check it fails.** Correctly combining over 5 disjoint partitions gives 6.50σ, which recovers the full-sample 6.94σ. Any hold-out procedure that returns *more* significance than the full sample it was drawn from has double-counted data.

**The honest number is in the same paragraph.** Line 908 reports 15 of 20 replicates negative, 5 positive.

$$\text{binomial } p(\geq 15 \mid n=20, p_0=0.5) = 0.0207 \;\Longrightarrow\; \mathbf{2.04\sigma}$$

A genuine 6.94σ effect should give 20/20. Separately, the era × station × lunation grid (line 2226) gives 13 of 26 sign-consistent, binomial p = 1.00 — exactly chance.

---

## D9 — Paper 17: the GP absorption tests

| Variant | η_resid | Result |
|---|---|---|
| Full 2D GP (elongation × time) | +5.30 × 10⁻⁵, p = 0.36 | sign flip, signal gone |
| Cross-validated GP | +2.45 × 10⁻⁵, 0.42σ | sign flip, signal gone |
| Time-only GP | −3.26 × 10⁻⁴, 5.50σ | survives |
| Adversarial PCA (20 of 82) | −2.93 × 10⁻⁴, 2.95σ | attenuated |

**Why the passing test is not a test.** A Gaussian process with the elongation dimension *removed* cannot absorb an elongation-locked signal. It is structurally incapable of failing. Two of the three GP variants kill the signal, and the one that preserves it has no discriminating power.

**PCA attenuation:**

$$\frac{2.93\times10^{-4}}{3.91\times10^{-4}} = 75\% \text{ of headline amplitude}$$
$$\frac{6.94}{2.95} = 2.35\times \text{ loss in significance}$$

**Variance model.** χ²_red = 0.0038 (Birge R_B = 0.062) on the headline, against χ²_red ≈ 0.48 at line 706 on published σ_m. The two variance definitions differ by

$$\frac{0.48}{0.0038} = 126\times$$

The Birge policy of applying no downward rescaling is conservative and correct in direction. But a precision-weighted estimator whose weights are wrong by ~16× in σ cannot deliver a calibrated σ_η, and the bootstrap and LOSO intervals should be primary.

**Defensible range.** Taking the adversarial PCA as the binding constraint, the supportable significance is approximately **3σ**, not 6.94σ.

---

## D10 — Paper 16: formal error against the paper's own bootstrap

**Inputs.** Weighted circular mean ψ̄ = 0.984 ± 0.046 (formal circular SE). Reported epoch-level bootstrap 95% CI = [0.737, 1.235].

$$SE_{\rm boot} = \frac{1.235 - 0.737}{2 \times 1.96} = 0.127 = 2.76 \times SE_{\rm formal}$$

| Basis | σ |
|---|---|
| Formal circular SE (0.046) | 21.4 |
| Bootstrap SE (0.127) | 7.75 |
| V-test, p = 2.04 × 10⁻⁵ | 4.26 |

The V-test agrees with the bootstrap, not the formal SE. The formal error understates the uncertainty by a factor of 2.76.

**Rayleigh consistency check.** For large n, $p \approx \exp(-n\bar{R}^2)$, so $n_{\rm eff} = -\ln p / \bar{R}^2$. With R̄ = 0.308:

| Statistic | p | Z = −ln p | n_eff |
|---|---|---|---|
| Unweighted | 1.34 × 10⁻⁴⁴ | 101.0 | 1,065 |
| Weighted | 1.39 × 10⁻¹³ | 29.6 | 312 |

**This confirms the pseudo-replication problem is fixed.** The unweighted statistic uses all 1,093 epochs (n_eff ≈ 1,065), not the 19,167 triplets. v0.3's adoption of epoch-level independence is real and my earlier criticism no longer applies. The weighted statistic loses effective sample size to the inverse-variance weights, which is expected.

---

## D11 — Paper 29: the screening fit and the amplification gap

**The three-anchor calibration.**

$$\mathcal{S}_\Sigma(|\Phi|/c^2) = \frac{1}{1 + (|\Phi|/c^2 / \Phi_{\rm half})^n}$$

Free parameters: Φ_half and n — **two**. Anchors: absorber, Cassini, galactic — **three**. But P29 line 159 states that the Cassini value (S_Σ ≈ 10⁻⁴) is "chosen to satisfy the Cassini bound," so this is two parameters against two real constraints.

$$\frac{|\Phi_\odot|/c^2}{|\Phi_{\rm gal}|/c^2} = \frac{2.12\times10^{-6}}{5.38\times10^{-7}} = 3.941 \qquad \text{(paper: 3.94)}$$
$$3.941^{3.56} = 131.9 \qquad \text{(paper: 132× contrast)}$$

The 132× is the arithmetic output of the fit, reproducing the number the fit was constructed to produce. n = 3.56 is fitted, not derived. Gate 10 is declared "conditionally closed" on this basis in the same paragraph that states cross-domain continuity is unproven.

**The required amplification.**

$$\mathcal{A}_{\rm env} = \frac{|\Delta v_T|}{|\Delta v_{\rm conf}|} = \frac{81.6\ \text{km/s}}{0.009\ \text{km/s}} = 9.07\times10^3 \qquad \text{(paper: } 9.3\times10^3\text{)}$$

Falsification floor stated by the paper: 𝒜_env < 10³.

**The mechanism assessment (P29 line 284, the paper's own calculation).**

| Quantity | Value |
|---|---|
| Edge-transition channel requires | B̃_eff ~ 10³⁴ |
| Physically motivated m_φ (10⁻³⁰ to 10⁻¹⁵ eV) supplies | B̃_eff ≲ 10⁵ |
| **Shortfall** | **10²⁹** |
| Even at ΔR ~ 1 AU, still requires | B̃_eff ~ 10²⁵ |

The only identified mechanism is excluded by 29 orders of magnitude. "Gate 4B remains open" understates this: the amplification currently has **no candidate mechanism**, only the uncomputed global holonomy of 𝒞_T,∥.

---

## D12 — Paper 15: what an n = 3 monotonic ordering is worth

P15's abstract calls the monotonic ordering of anomaly magnitudes with the trajectory-asymmetry factor "the primary systematic discriminator," on a sample of three.

$$\text{orderings of 3 items} = 3! = 6$$
$$P(\text{one specific monotonic order}) = \tfrac{1}{6} = 0.167 \;\Longrightarrow\; 0.97\sigma$$
$$P(\text{monotonic in either direction}) = \tfrac{2}{6} = 0.333 \;\Longrightarrow\; 0.43\sigma$$

**A perfect result at n = 3 cannot exceed about 1σ.** It cannot serve as a primary discriminator for any model.

**Model budget (line 271).** Seven fixed geometry-envelope coefficients plus one fitted amplitude, against four retained anomalies.

**Failed prediction (line 297).** Cassini: predicted −0.023 mm/s, observed +0.11 mm/s. The sign is wrong on one of the four retained anomalies and is attributed to a pre-existing literature mismatch rather than counted as a miss.

---

## D13 — Paper 26: Δχ² on scale, and the failed cross-checks

**Δχ² to sigma, for equal-parameter non-nested models:**

$$\Delta\chi^2 = -3.4 \;\Longrightarrow\; \sqrt{3.4} \approx 1.84\sigma \quad (z_{\rm los}=5)$$
$$\Delta\chi^2 = -7.5 \;\Longrightarrow\; \sqrt{7.5} \approx 2.74\sigma \quad (z_{\rm los}=100)$$

Bayes factors of 4.6 and 61.8 are attached to preferences of roughly 1.8σ and 2.7σ in the underlying likelihood.

**BAO distance-duality (line 400).**

$$\bar\eta = 0.866 \pm 0.020 \;\Longrightarrow\; \frac{1 - 0.866}{0.020} = \mathbf{6.7\sigma \text{ from unity}}$$

Both ΛCDM and TEP predict η ≡ 1 exactly. A 6.7σ failure of a quantity that both models fix means the distance pipeline is wrong — and it is the same pipeline that produces BF = 4.6.

**Redshift-cut robustness (line 428).** Removing the z < 0.023 anchors drops BF to ≈ 0.8, which slightly favours ΛCDM. The headline evidence lives in the low-redshift tail where peculiar velocities and calibration systematics dominate.

**Joint SNe+CMB (line 305).**

$$\epsilon_T^{\rm CMB} = -0.0015 \pm 0.0037 \;\Longrightarrow\; 0.41\sigma, \text{ consistent with zero}$$

Gelman–Rubin R−1 = 0.0276 against a convention of < 0.01, so the joint chain is not converged.

---

## D14 — Corpus-level duplicate quantities

**κ_Cep / κ_gal (× 10⁶ mag)** — ten values for one quantity:

| Value | Source |
|---|---|
| 0.116 | P11 mixed model (0.31σ, consistent with zero) |
| 0.290 | P11 Step 44, σ_v = 150 |
| 0.326 | P11 Step 44, σ_v = 182.1 |
| 0.369 | P11 Step 42 restricted |
| 0.400 | P12 cites as P11 Step 44 — **stale** |
| 0.417 | P11 redshift-only, σ_v = 182 |
| 0.452 | P11 redshift-only, σ_v = 150 |
| 0.600 | P12 JWST-side recovery |
| 0.960 | P12 canonical benchmark |
| 1.05 | P17 citation |

**κ_MSP:**

| Value | Source |
|---|---|
| (2.9 ± 4.5) × 10⁴ | P10's own value — **0.64σ, consistent with zero** |
| 2.9 × 10⁴ (no error) | P13 line 190, cited as a point value |
| ~10⁶–10⁷ | P17 |
| ~10⁶ unscreened | P13 / P10 geometric factor |

A quantity at 0.64σ cannot anchor a cross-scale hierarchy.

**R_s (AU), all within Paper 13:** 2,646 ± 182 (headline, below all five of its own |Z| sub-bins); 3,447 / 4,912 (joint profile likelihood, which the body says the inference rests on); 4,662 / 7,131 (two-bin fixed α, which the abstract quotes); 8,670 → 4,326 (five-bin); 7,709 ± 3,222 (RV subset); 3,551 / 10,307 (distance split).

**Chameleon index n, all within Paper 13:** 1.02 ± 0.14 (five-bin direct fit, abstract); 1.5 ± 0.9 (two-bin inversion); 1.60 ± 0.89 (log-linear).

**H₀ (km s⁻¹ Mpc⁻¹):** 66.65 (P11 abstract), 66.70 (P26 joint MCMC), 68.36 (P11 Hubble-flow intercept), 68.61 (P11 Step 44), 71.77 (P11 full 3,490-row SH0ES propagation).

**M31 P–L gradient:** 2.6σ (P11 plain), 3.24σ (P11 PHAT-matched), 3.65σ (P31).

---

## Summary table: claimed versus supportable

| Claim | As published | Recomputed | Basis |
|---|---|---|---|
| Orbital-velocity coupling (P5) | 5.1σ | **3.65σ** | P2's own N_eff = 11 test; 5.1σ is the MC floor |
| LLR Nordtvedt (P17) | 6.94σ | **~3σ** | adversarial PCA, 20 components |
| LLR blind hold-out (P17) | 14.3σ | **6.5σ** | 5 disjoint partitions, not 20 overlapping |
| LLR sign consistency (P17) | implied ~20/20 | **2.04σ** | 15/20 binomial |
| Scintillation closure (P16) | 21.4σ formal | **4.3–7.8σ** | paper's own bootstrap and V-test |
| Flyby ordering (P15) | "primary discriminator" | **≤1σ** | n = 3 has a 1/6 floor |
| SNe temporal shear (P26) | BF 4.6 "substantial" | **1.84σ** | Δχ² = 3.4 at equal k |
| ρ_T (P6) | 20 ± 7 g/cm³ | **6–128 g/cm³** | P2's ±47% propagated through L_c⁻³ |
| CMB variance ratio (P2/P5) | 5,570–10,300× | **not an evidence measure** | denominator consistent with zero |
| κ_gal (P12) | 0.960 benchmark | **0.29–0.60 measured** | all P11 estimates below P12's own 0.56 failure point |

**Results that survive unchanged:** Paper 31's R_H = 1.009 ± 0.006 rejecting KBC/MOND at 8.4σ; Paper 18's sound-horizon ratio at <6 ppm; Paper 19's blind sign-coherence result; Paper 30's equivalence theorem; Paper 23's δm/m ≲ 1.5 × 10⁻¹⁹ falsifiability floor.
