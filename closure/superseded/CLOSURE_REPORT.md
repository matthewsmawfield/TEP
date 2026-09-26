# TEP Scalar Field Closure Programme: Full Report

> **SUPERSEDED 2026-09-12.** This report's candidate potential was the
> quadratic $V(\phi) = \Lambda^4(1 + \phi^2/(2\Lambda^2))$, which gives a
> density-independent mass and therefore no environmental screening response.
> The corrected potential is the inverse-power form
> $V(\phi) = \Lambda^4(1 + \Lambda/\phi)$ at the same zero-parameter scale
> $\Lambda = \sqrt{M_{\rm Pl} H_0}$. The authoritative updated report is
> `CLOSURE_RESOLUTION.md`; verified numbers are in
> `results/tep_screening_closure_results.json`. This document is retained as
> a historical record of the first computational attempt.

## 1. What was done and why

### 1.1 Motivation

The Jakarta manuscript (Paper 0) defines the TEP scalar field architecture with frozen conformal coupling $A(\phi) = \exp(\beta_A \phi/M_{\rm Pl})$, $\beta_A = -1$, and an unresolved scalar self-interaction $V(\phi)$. The disformal coupling $B(\phi)$ is presented as a prescribed strong-field realization from Paper 28, not as a derived law. The five Jakarta edits completed in the previous pass established the correct theoretical framing:

- The common-gradient identity $\mathcal{D}_{\mu\nu} = \frac{B M_{\rm Pl}^2}{A^2 \beta_A^2} \Sigma_\mu \Sigma_\nu$
- The holonomy kernel in Temporal Shear language
- The demotion of the quartic-Gaussian to a prescribed strong-field realization
- The correction of the late-time cone statement to bound the observable deformation, not $B$ alone
- The downgrade of experimental forecasts to benchmark sensitivity windows

What remained was the forward-closure calculation: does a single $V(\phi)$ with $K=1$ generate the screened profiles required by the corpus? This report describes the first computational attempt.

### 1.2 Methodology

The programme was designed as a code-first closure test with predeclared pass/fail gates, following the principle that generators get frozen after they generate, not before.

The scalar equation solved is:

$$\phi'' + \frac{2}{r}\phi' = V_{,\phi}(\phi) + \rho_*(r) A_{,\phi}(\phi)$$

with $A_{,\phi} = (\beta_A/M_{\rm Pl}) A$, $\beta_A = -1$, solved as a two-point boundary value problem using scipy's `solve_bvp` collocation solver on a logarithmic radial mesh.

The calculation was performed in dimensionless variables:
- $u = \phi/M_{\rm Pl}$ (dimensionless field)
- $x = r/R$ (dimensionless radius)
- $s = \rho_*/\rho_{\rm ref}$ (dimensionless density)
- $\lambda = \rho_{\rm ref} R^2/M_{\rm Pl}^2$ (single dimensionless group)

This avoids the unit-confusion problem that affected earlier cross-paper comparisons.

The candidate potential is:

$$V(\phi) = \Lambda^4 \left(\cosh\frac{\phi}{M_{\rm Pl}} - 1\right)$$

This is the simplest viable candidate because:
- Minimum at $\phi = 0$ (consistent vacuum, $V = 0$, $V_{,\phi} = 0$)
- $V_{,\phi} > 0$ for $\phi > 0$ (satisfies the equilibrium condition $V_{,\phi} = \rho_* A/M_{\rm Pl} > 0$)
- $V_{,\phi\phi} > 0$ everywhere (stable)
- No singularity at any finite $\phi$
- Connects to $\phi \to +\infty$: $V \to (\Lambda^4/2) e^{\phi/M_{\rm Pl}} \to \infty$ (temporal horizon, $A \to 0$)
- One parameter ($\Lambda$, or equivalently $\eta = \Lambda^4 R^2/M_{\rm Pl}^2$)

The pure exponential $V = \Lambda^4 e^{\phi/M_{\rm Pl}}$ was rejected because it has no finite vacuum ($V \to 0$ only as $\phi \to -\infty$), making $u_{\rm ambient} = 0$ inconsistent.

No thin-shell estimates, Khoury-Weltman analytic solutions, or chameleon mechanism assumptions were used. The Temporal Topology is the ontology; $V(\phi)$ generates $\phi(x)$; $S_\Sigma$ is read off the solved exterior profile, not estimated from a geometric formula.

### 1.3 Validation stages

Four validation tests were run before any production calculation:

**Stage 0a — Linear Yukawa:** $V = m^2\phi^2/2$ with a uniform sphere source. The solver converges, the exterior decays as $\exp(-\mu x)/x$ with fitted slope $-2.97$ vs expected $-3.00$ (1% error), and the field is regular at the center. PASS.

**Stage 0b — Conformal exactness (Theorem 2):** The numerical integral of the conformal connection $\beta_A u' dx$ matches the analytic $\beta_A \Delta u$ to $6 \times 10^{-7}$. This confirms $d(\ln A)$ is an exact differential and $A(\phi)$ cannot source synchronization holonomy. PASS.

**Stage 0c — Weak source:** At $\lambda = 10^{-9}$ (Earth-like), the solver converges to tol=$10^{-12}$ and the field response is linear in $\lambda$ and bounded by the infinite-medium particular solution. PASS.

**Stage 0d — Nonlinear (cosh):** With an informed initial guess, the solver converges for the cosh potential and produces a physically sensible profile (field rises toward equilibrium inside, relaxes to ambient outside). PASS.

### 1.4 Production calculations

**Stage 1 — Earth and Sun profiles:** An $\eta$ scan from $10^{-4}$ to $10^{14}$ was run, solving both Earth and Sun with the same potential parameter ($\eta_{\rm sun} = \eta_{\rm earth} \times (R_\odot/R_\oplus)^2 = \eta_{\rm earth} \times 11{,}934$). Key findings:

- $\lambda_{\rm earth} = 4.17 \times 10^{-9}$, $\lambda_{\rm sun} = 1.27 \times 10^{-5}$
- The field excursion is tiny everywhere: $u \sim 10^{-9}$ (Earth) to $10^{-5}$ (Sun)
- Cassini ($S_\Sigma^{(\odot)} < 3.4 \times 10^{-3}$) is satisfied for all 19 $\eta$ values tested
- The Sun is always heavily screened because $\eta_{\rm sun}$ is 4 orders of magnitude larger than $\eta_{\rm earth}$
- The interesting regime is $\log \eta_{\rm earth} = -1$ to 0, where $S_{\rm earth} \sim 10^{-3}$ to $10^{-7}$

**Stage 2 — $B_0$ window scan:** For each $\eta$, the disformal deformation shape $\epsilon_B/B_0$ was computed from the solved profiles. The quartic Gaussian envelope is confirmed numerically inert: $\sigma_B$ is irrelevant in the weak-field regime (the envelope equals 1 to $>30$ decimal places when $|\phi|/M_{\rm Pl} < 10^{-4}$). The scan is one-dimensional in $B_0$.

### 1.5 Critical correction

The initial "optimal $B_0$" was computed as the geometric mean of $B_0^{\rm min}$ (triangle detection) and $B_0^{\rm max}$ (Cassini limit). This was wrong. At the geometric mean, $\epsilon_B \sim 10^{29}$, which completely violates EFT validity ($\epsilon_B \ll 1$).

The correct operating point is near $B_0^{\rm min}$, where:
- $H_{\rm resid} \sim 10^{-18}$ (at detection threshold)
- $\epsilon_B \sim 10^{-18} \ll 1$ (EFT valid)
- $\epsilon_{\rm Cassini} \sim 10^{-108}$ to $10^{-170}$ (trivially safe)
- $M_B = B_0^{-1/4} \sim 30$–$200$ eV (light but not absurd EFT scale)

At $H_{\rm resid} = 10^{-16}$ (100× above threshold): $\epsilon_B = 10^{-16}$, $M_B \sim 10$–$60$ eV. Still EFT-valid.

---

## 2. Results

### 2.1 What passes

| Gate | Status | Details |
|---|---|---|
| A1. Cassini ($S_\Sigma < 3.4 \times 10^{-3}$) | PASS | All 19 $\eta$ values |
| A2. Earth profile smooth | PASS | Field flattens inside, recovers outside |
| A5. Stability ($m_{\rm eff}^2 > 0$) | PASS | Cosh is convex everywhere |
| A6. No environment-specific retuning | PASS | One $\eta$ for Earth and Sun |
| C1. $\epsilon_B$ below Cassini/GW170817 | PASS | $10^{-108}$ to $10^{-170}$ at Cassini |
| C2. Triangle holonomy detectable | PASS | $H \sim 10^{-18}$ at $B_0^{\rm min}$ |
| C3. $B_0$ window exists | PASS | $B_0^{\rm min}$ to $B_0^{\rm max}$ spans $\sim 90$ orders |
| EFT validity ($\epsilon_B \ll 1$) | PASS | $\epsilon_B \sim 10^{-18}$ at $B_0^{\rm min}$ |
| Theorem 2 (conformal exactness) | PASS | Numerical verification to $6 \times 10^{-7}$ |
| B-envelope inertness | PASS | $\sigma_B$ irrelevant for $|\phi|/M_{\rm Pl} < 10^{-4}$ |

### 2.2 What has not been tested

| Gate | Status | Details |
|---|---|---|
| A3. $\lambda_T$ correlation length | NOT TESTED | Requires stochastic source model |
| A4. Ambient unscreened branch | NOT TESTED | Requires cosmological solution |
| B1. Temporal Horizon continuation | NOT TESTED | Cosmological evolution with cosh V not solved |
| B2. No singularity before TH | NOT TESTED | Requires cosmological integration |
| B3. Acceptable stress-energy | NOT TESTED | Requires cosmological integration |
| B4. Noncanonical $K(\phi)$ compatibility | NOT TESTED | $K=1$ assumed throughout |
| Data calibration | NOT TESTED | $\eta$ not fitted to GNSS/LLR/wide-binary data |
| Held-out prediction | NOT TESTED | No environment withheld from calibration |

### 2.3 What does not work across the corpus

The calculation is a **local weak-field closure test**. It does not automatically extend to all TEP domains:

**Paper 28 (Bahrain/BH):** Uses $\phi/M_{\rm Pl} \sim 2.0$, a completely different field regime. The cosh $V$ at this field value gives $V = \Lambda^4(\cosh 2 - 1) = 2.76 \Lambda^4$, which is $O(\Lambda^4)$ — the right order for the strong-field construction. But Paper 28 uses a different scalar equation (with scalar-Gauss-Bonnet coupling) and a different field profile. The cosh $V$ has not been tested in the BH context. The $B(\phi)$ envelope is active at $|\phi|/M_{\rm Pl} \sim 2$ (the quartic Gaussian suppresses it by $\sim 20\%$), while in the weak field it is completely inert. These are disjoint regimes.

**Paper 29 (Dubai/BBN):** Uses constant $B = B_0/M_{\rm Pl}^2$ or inverse-field $B = B_0 M_{\rm Pl}^2/\phi^2$ in the diffuse absorber regime. The weak-field form $B(\phi) \sim B_0 (\phi/M_{\rm Pl})^2$ from this calculation differs from the constant form by a factor of $(M_{\rm Pl}/\phi)^2 \sim 10^{18}$. The open reconciliation target identified in Paper 29 is not resolved by this calculation. The two forms operate in disjoint field regimes, which is why Paper 29 tags this as an open reconciliation target rather than a structural inconsistency.

**Paper 31 (Temporal Horizon):** Requires $\phi \to +\infty$ ($A \to 0$). The cosh $V$ does connect to this asymptotic ($V \to \infty$), but the cosmological evolution has not been solved. The density-mass scaling $m_{\rm eff} \propto \rho^{1/4}$ (for the cosh potential in the small-field limit) is too weak to span the 31 orders of magnitude in density from $\rho_T$ to cosmological mean density while keeping the field both locally screened and cosmologically light. This is the known chameleon tension, but framed in TEP language: the same $V$ must work locally and cosmologically, and the cosh $V$ has not been shown to do both.

**Papers 2-7 (GNSS/clock):** The $S_A$ values from the cosh $V$ have not been compared to actual GNSS clock-rate data. The $\eta$ parameter is not calibrated.

**Papers 14-15 (wide binaries, UCD):** The screening profiles for these environments have not been computed. A held-out prediction against wide binaries or UCDs would be the most convincing validation.

**Paper 10 (cosmology):** The cosmological evolution with the cosh $V$ has not been solved. The $c_T = 1$ constraint, BBN bounds, CMB constraints, and structure formation all require a cosmological integration with the same $V$.

---

## 3. Assessment

### 3.1 The cosh potential is a candidate, not a result

The calculation shows that the minimal $K=1$ closure with $V = \Lambda^4(\cosh(\phi/M_{\rm Pl}) - 1)$ passes the local gates. This is a necessary condition for the closure programme, not a sufficient one. The potential is a legitimate candidate canonical $V$ for the local weak-field EFT, but it is not yet frozen.

What would justify freezing it:
1. Calibrate $\eta$ against actual GNSS/LLR data
2. Predict a held-out environment (wide binaries or UCDs)
3. Solve the cosmological continuation (Gate B)
4. Show the same $V$ is compatible with Paper 28's strong-field branch

What would justify not freezing it:
1. If no $\eta$ simultaneously fits GNSS and LLR
2. If the held-out prediction fails
3. If the cosmological continuation is singular or unstable
4. If the required $B_0$ is incompatible with Paper 29's absorber analysis

### 3.2 The $B_0$ window is real but narrow in practice

The window technically spans $\sim 90$ orders of magnitude, but the EFT-valid region is much smaller. At $B_0^{\rm min}$ (triangle detection threshold), $\epsilon_B \sim 10^{-18}$ and the EFT is valid. At $B_0^{\rm max}$ (Cassini limit), $\epsilon_B \sim 10^{-15}$ and the EFT is still valid. But the "optimal" geometric-mean $B_0$ gives $\epsilon_B \sim 10^{29}$, which is completely unphysical.

The correct interpretation is: the disformal signal is at the very edge of detectability ($10^{-18}$) with a coupling scale $M_B \sim 30$–$200$ eV. This is not a strong prediction — it is a sensitivity statement. The Jakarta edits correctly labelled these as "benchmark sensitivity windows," and the calculation confirms that label.

### 3.3 The fundamental challenge

The disformal deformation in the weak field scales as:

$$\epsilon_B \sim B_0 \left(\frac{\phi}{M_{\rm Pl}}\right)^2 \left(\frac{\phi'}{M_{\rm Pl}/R}\right)^2 \sim B_0 \times 10^{-36}$$

because the field excursion is $\phi/M_{\rm Pl} \sim 10^{-9}$ (Earth) to $10^{-6}$ (Sun). The double suppression by the tiny field means the disformal signal is inherently small in the weak-field regime. This is consistent with the common-gradient identity: the disformal deformation is quadratic in the shear, and the shear is already suppressed.

This is not a failure of the theory — it is the reason GW170817 is satisfied. But it means the disformal sector is only detectable at the $10^{-18}$ level with next-generation metrology, which is exactly what the Jakarta experiments target.

### 3.4 The cross-paper reconciliation remains open

The calculation does not resolve the Paper 28/Paper 29 $B(\phi)$ mismatch. The weak-field form $B \sim B_0 (\phi/M_{\rm Pl})^2$ and the Paper 29 constant form $B = B_0/M_{\rm Pl}^2$ differ by $\sim 10^{18}$. A single composite $B(\phi)$ that reduces appropriately in both regimes has not been derived. This remains an open reconciliation target, as stated in Paper 29.

---

## 4. Recommendations

### 4.1 Immediate (do not edit manuscripts yet)

1. **Do not freeze $V(\phi)$ in Jakarta.** The cosh potential passes local gates but has not been calibrated against data or tested cosmologically. The current Jakarta language ("the unique corpus-wide microscopic form and normalization of $B(\phi)$ therefore remain part of the action-closure problem") is correct and should remain.

2. **Do not add the cosh $V$ to the parameter registry.** It is a candidate, not a frozen parameter. Adding it prematurely would repeat the density-pinning error from the previous pass.

3. **Keep the benchmark sensitivity window language.** The calculation confirms that the experimental forecasts are sensitivity targets, not derived predictions. The Jakarta edits are correct.

4. **Add one quantitative sentence to the Jakarta §2.2 disformal coupling paragraph.** The B-envelope inertness is now a numerical result, not just an interpretation: the quartic Gaussian envelope equals 1 to $>30$ decimal places for $|\phi|/M_{\rm Pl} < 10^{-4}$, confirming that $\sigma_B$ is irrelevant in the weak-field regime. This strengthens the existing demotion language.

### 4.2 Near-term computational programme

1. **Calibrate $\eta$ against GNSS data.** Use the solved Earth profile to compute $S_A^{(\oplus)}$ and compare with GNSS clock-rate residuals. This is the first real data comparison.

2. **Compute a held-out prediction for wide binaries.** Solve the scalar profile for a wide-binary system ($\rho \sim 10^{-24}$ g/cm$^3$, $R \sim 10^4$ AU) using the same $\eta$ and predict the screening signature. Compare with Paper 14's data.

3. **Solve the cosmological continuation.** Integrate the cosh $V$ in a Friedmann background and check whether it permits the temporal horizon asymptotic ($\phi \to +\infty$, $A \to 0$) without singularity. This is Gate B.

4. **Test the cosh $V$ in the Paper 28 BH context.** At $\phi/M_{\rm Pl} \sim 2$, the cosh gives $V = 2.76 \Lambda^4$. Check whether this is compatible with Paper 28's strong-field construction. This does not require redoing Paper 28's calculations — it requires checking whether the cosh $V$ is in the admissible class for the BH branch.

### 4.3 Medium-term

1. **If $\eta$ calibration succeeds:** Upgrade the Jakarta language from "candidate $V$" to "calibrated local weak-field EFT" with the specific $\eta$ value and the data it was fitted to.

2. **If the cosmological continuation succeeds:** Add a sentence stating that the same $V$ permits the temporal horizon asymptotic.

3. **If the held-out prediction succeeds:** This is the strongest result. State that the same $V$ calibrated on Earth/Sun predicts the wide-binary (or UCD) screening signature without retuning.

4. **If any gate fails:** Do not manufacture a potential. Report the negative result. A failed closure is itself a finding about the framework.

### 4.4 What not to do

- Do not introduce $P(X, \phi)$ or noncanonical $K$ without first exhausting the $K=1$ cosh candidate
- Do not freeze $V$ before data calibration
- Do not claim the $B_0$ window is a prediction — it is a sensitivity statement
- Do not resolve the Paper 28/Paper 29 $B(\phi)$ mismatch by proclamation — it requires a composite derivation
- Do not edit Papers 28, 29, or 31 based on this local calculation

---

## 5. Files produced

```
closure/
  radial_bvp.py              — original dimensional solver (superseded)
  radial_bvp_dimless.py      — dimensionless solver, observable extraction, loop test
  stage0_validation.py       — four validation tests (all PASS)
  stage1_earth_sun.py        — eta scan for Earth and Sun (19 values, all Cassini PASS)
  stage2_b0_scan.py          — B0 window scan (4 eta values, all windows exist)
```

All code uses the correct TEP framing: BVP solve, direct observable extraction, no thin-shell estimates, no chameleon mechanism. The Temporal Topology is the ontology; $V(\phi)$ generates $\phi(x)$; $S_\Sigma$ is read off the solved profile.

---

## 6. Conclusion

The minimal $K=1$ closure with $V = \Lambda^4(\cosh(\phi/M_{\rm Pl}) - 1)$ passes the local weak-field gates. The cosh potential is a viable candidate for the local EFT, with one parameter ($\eta$) that generates screened profiles for both Earth and Sun without environment-specific retuning. The $B_0$ window exists and is EFT-valid at the detection threshold. The conformal exactness theorem (Theorem 2) is verified numerically. The B-envelope inertness is confirmed quantitatively.

The calculation does not extend to the BH regime (Paper 28), the absorber regime (Paper 29), or the cosmological regime (Paper 31). The cross-paper $B(\phi)$ reconciliation remains open. The $\eta$ parameter is not calibrated against data. No held-out prediction has been made.

The correct next step is data calibration and held-out validation, not manuscript editing. The Jakarta language is already correct. The calculation supports it.
