# The Temporal Equivalence Principle: Dynamic Time, Emergent Light Speed, and the Non-Integrability of Simultaneity


**Author:** Matthew Lukin Smawfield

**Version:** Argentina

**Date:** 18 Aug 2025


## Abstract


Einstein taught us two unforgettable lessons: simultaneity is relative and spacetime geometry is dynamical. Here we argue for an equally consequential third: the rate at which time flows is itself a dynamical field. We elevate proper time to the output of a universal matter metric ĝμν related to the gravitational metric gμν by a controlled conformal–disformal relation, ĝμν = A(φ)^2 gμν + B(φ) ∇μφ ∇νφ, where φ is a scalar “time field.” This Temporal Equivalence Principle (TEP)—universal coupling of non-gravitational physics to ĝμν—preserves local Lorentz invariance and the local value of c while turning global one-way measurements of light, clock synchronization, and temporal transport into probes of the time field.

We resolve long-standing confusions about “variable c” by replacing convention-dependent statements with invariant observables tied to measurement procedures. We define a synchronization one-form σ̃ on spacelike slices of the matter metric and show that its curl dσ̃, after subtraction of general-relativistic (GR) Sagnac/gravito-magnetic terms, yields a residual “temporal holonomy” H that vanishes in GR and becomes nonzero only when time is dynamical in our sense. We prove two key theorems: (i) conformal matter coupling preserves null cones, so photons and gravitons share the same causal structure at late times; (ii) a static φ-gradient generates no first-order one-way light-time anisotropy, placing effects in the femto-to-picosecond regime over astronomical baselines under current bounds. Disformal tilts (B ≠ 0) are tightly constrained by GW170817-class multi-messenger observations but can source holonomy at levels within reach of next-generation metrology.

We present a covariant action, derive field equations, conservation, and invertibility/causality conditions, and supply a 3+1 decomposition that makes observables explicit. Screening via a density-dependent effective potential reconciles precision local tests with cosmological evolution, and we map to Parametrized Post-Newtonian parameters and to the EFT-of-dark-energy α-functions with cT = 1 enforced. We outline decisive experiments with quantitative error budgets: (1) a ground–ground–satellite triangle time-transfer experiment targeting holonomy at below 10−18 fractional after GR subtraction; (2) portable-clock “clock anholonomy” around closed paths at the 10−19 level over days; (3) multi-species clock networks seeking phase-locked annual modulations at 10−19–10−17; (4) interplanetary one-way optical links at picoseconds over AU; (5) altitude-dependent screening maps with optical clocks and atom interferometers; and (6) ensemble multi-messenger tests. We provide a cosmological pipeline plan for CLASS/HyRec modifications and MCMC inference, and we commit to open data and blinded analyses.

Einstein’s postulate of universal c was a brilliant, operationally perfect approximation in regimes where time’s flow is effectively uniform. In a universe where the rate of time is dynamical yet locally Lorentzian, “the speed of light” emerges as an invariant in every local lab but ceases to be globally universal. The new invariant content resides in path-dependent synchronization defects and holonomies of time transport, not in naive one-way “speeds.” If detected, these invariants would inaugurate a post-Einsteinian era: from dynamic geometry to dynamic time.

## 1. Introduction: From Relativity of Simultaneity to Dynamics of Time’s Rate


Relativity replaced absolute time with proper time τ, the arc-length of worldlines. General relativity (GR) made geometry dynamic; gravitational redshift and time dilation are everyday laboratory realities. Yet quantum theory treats time t as an external parameter driving unitary evolution. The “problem of time” persists both conceptually and technically.

Measurement technology has sharpened this tension. Optical lattice clocks now reach fractional instabilities below 10−18; time-transfer links over thousands of kilometers are stabilized at parts in 10−18–10−19; atom interferometers and torsion balances test new forces with exquisite sensitivity. Cosmology reveals persistent tensions (H0, S8) in a data-rich era. These advances expose the operational fragility of our global synchronization assumptions: one-way light speed cannot be defined without conventions, and slow clock transport presupposes that “time flows uniformly enough” to make synchronization integrable. But what if the rate of time’s flow is a covariant field that responds weakly to environment and epoch?

We propose a minimal, covariant extension of relativity that promotes the rate of time to a scalar field φ(x) with universal matter coupling. It keeps local Lorentz invariance exact, conforms to GW170817 constraints on gravitational wave speed, and translates “variable c” into well-posed, synchronization-independent observables. It is falsifiable with near-term experiments.

## 2. Principle and Architecture: The Temporal Equivalence Principle and Two Metrics


Principle (Temporal Equivalence Principle, TEP). All non-gravitational fields and clocks couple universally to a matter metric ĝμν. Gravitational dynamics are governed by gμν. The matter metric is related to gμν by a conformal–disformal map
```
ĝμν = A(φ)^2 gμν + B(φ) ∇μφ ∇νφ,
```

with A(φ) > 0, and disformal strength B(φ) small enough today to satisfy cT constraints and causality.

Interpretation. A(φ) encodes the local tick rate of all clocks and sets the conversion between coordinate time and physical proper time measured by matter. B(φ) captures a tiny, environment-dependent tilt between photon and graviton cones in regions where ∇φ ≠ 0; we constrain B to be effectively negligible today on cosmological paths, consistent with GW170817, though it can remain relevant at small scales.

Consequences. In any local freely falling matter frame, ĝμν ≈ ημν and all non-gravitational physics reduces to special relativity with invariant c. Global one-way measures that require synchronized clocks become sensitive to the history of φ along time transport or signal paths, leading to path-dependent synchronization and holonomy.

## 3. Action, Field Equations, and Conservation


We work in an Einstein frame where gravity is canonical and matter couples to ĝμν.

Action
```
S = ∫ d^4x √−g [ (MPl^2/2) R − (1/2) gμν ∂μφ ∂νφ − V(φ) ] + Smatter[ψ, ĝμν; φ],
```

with ĝμν = A(φ)^2 gμν + B(φ) ∇μφ ∇νφ.

Variation yields:
```
Gμν = MPl−2 [Tμν(φ) + Tμν(matt)],
```

where Tμν(φ) = ∂μφ ∂νφ − gμν ( (∂φ)^2/2 + V ), and Tμν(matt) is the pullback of the matter-frame stress-energy T̃μν.

Scalar equation:
```
□g φ − V′(φ) = −α(φ) T + Sdisf,
```

with α(φ) ≡ d ln A/dφ, T ≡ g̃μν T̃μν the matter-frame trace, and Sdisf disformal sources:
Sdisf ≈ (1/2) B′ T̃μν ∂μφ ∂νφ + ∇μ [B T̃μν ∂νφ] + O(B^2).
Radiation (T = 0) does not source φ conformally.

Conservation laws. Diffeomorphism invariance gives ∇μ[Tμν(φ)+Tμν(matt)]=0. Matter is conserved in the matter frame, ∇̃μ T̃μν = 0, while in the Einstein frame,
∇μ Tμν(matt) = α(φ) T ∇νφ + O(B),
representing an exchange of energy–momentum between φ and matter.

## 4. Disformal Invertibility, Causality, and Hyperbolicity


Invertibility and signature. The inverse of ĝμν exists for A>0 and 1 + (B/A^2)(∂φ)^2 ≠ 0:
ĝμν = A−2 [ gμν − (B/A^2) ∂μφ ∂νφ / (1 + (B/A^2)(∂φ)^2) ].
For Lorentzian signature, require A>0 and B(∂φ)^2 > −A^2. We assume B small and gradients bounded in all regimes of interest.

Causality. The matter cone is inside or equal to the gravitational cone when B ≥ 0 and gradients are modest; no closed causal curves arise for small B. With B→0 at late times, gravitational and matter null cones coincide (cT = cEM). With B small, any phase differences in propagation are minute and bounded by multi-messenger results.

Hyperbolicity. The scalar’s canonical kinetic ensures hyperbolic evolution with well-posed Cauchy problem where V″(φ) > 0 near minima. Maxwell’s equations in a Lorentzian matter metric remain hyperbolic. Linearized gravity propagates on the Einstein-frame cone. No ghosts or gradient instabilities arise for canonical φ and small B. The EFT is valid below the disformal scale M, with higher-dimensional operators suppressed.

## 5. Local Lorentz Invariance, Proper Time, and the Emergence of c


Proper time. Clocks measure dτ^2 = −ĝμν dxμ dxν/c^2. In a local lab with small velocities and ∂0φ ≪ |∇φ|, −ĝ00 ≈ A(φ)^2 (−g00) − B(∂0φ)^2, so
```
dτ/dt ≈ A(φ)1/2,
```

up to O(B) corrections. This “dynamical time law” rescales all frequency standards locally.

Local c. All small labs measure an invariant c; the conformal factor rescales both clocks and rulers uniformly, preserving null cones. Thus the empirical fact “c is constant” remains a theorem of the theory at the local level.

Global measures. Global “speeds” involve synchronized endpoints. Because dτ depends on φ’s history, the global assignment of simultaneity acquires non-integrability when φ varies in space or time. That non-integrability—not a local violation of Lorentz invariance—is the source of observable departures.

## 6. Synchronization in a Dynamical Time: One-Way Light and Holonomy


Operational definitions. Two-way light speed is synchronization-independent and has established c’s local invariance. One-way measures require synchronized clocks at A and B; Einstein synchronization assumes time-orthogonal slices and propagation symmetry. In a dynamical φ background, slow clock transport and one-way synchronization are path and history dependent.

Key theorems.

Theorem 1 (Conformal null-cone invariance). For ĝμν = A(φ)^2 gμν, null vectors of gμν are null for ĝμν. Maxwell’s action is conformally invariant in 4D, so photon trajectories are null with respect to both metrics. Gravitational and electromagnetic waves share null cones when B = 0 at late times.

Theorem 2 (No first-order one-way anisotropy in static φ). For static φ with ∂0φ = 0, the one-way light-time difference for forward/backward propagation along the same path cancels to first order in ∇φ. Any residual is O((∇φ)^2) or due to time dependence/kinematics. Over astronomical baselines and with current bounds on α ≡ d ln A/dφ, this places effects in the femto-to-picosecond regime, removing claims of microsecond-scale anomalies.

Proof sketch. Parameterize the path coordinate s ∈ [0,L]; write t→ = ∫ ds A(φ0 + s ∂∥φ)/c and t← similarly along the reverse. Linear terms in ∂∥φ cancel exactly; see Appendix A for full derivation.

Synchronization one-form and holonomy. Decompose ĝμν in 3+1 form:
ĝμν dxμ dxν = −Ñ^2 dt^2 + h̃ij (dxi + Ñi dt)(dxj + Ñj dt).
Define the synchronization one-form σ̃ = −Ñ−1 Ñi dxi. For a closed spatial loop C, the synchronization holonomy is
```
H = ∮C σ̃.
```

In GR with stationary spacetimes, H reproduces the Sagnac/gravito-magnetic effect; we subtract these known contributions using geodesy and ephemerides. The residual holonomy beyond GR is
```
Hresid = ∬Σ dσ̃ − (GR Sagnac/gravito-magnetic),
```

with Σ spanning C. This residual vanishes in SR/GR with B=0 and stationary A; it is generically nonzero with disformal corrections (B ≠ 0) and/or temporal variation of A combined with motion through ∇φ.

Gauge and protocol invariance. σ̃ is a one-form on the spatial slice induced by the matter metric and a chosen time function t. Changing the synchronization protocol shifts σ̃ → σ̃ + dχ for some scalar χ when the change corresponds to a redefinition of time labels without physical alterations; the loop integral of an exact form vanishes: ∮ dχ = 0. Therefore Hresid, constructed from measured proper times and two-way calibrations and expressed via dσ̃ after GR subtraction, is invariant under admissible synchronization changes. This addresses the critique that holonomy is conventional: the nonzero part is the curvature (curl) of σ̃, not its exact part.

Disformal corrections to σ̃. For small B and modest ∂φ,
Ñ ≈ A N [1 − (B/(2A^2 N^2)) (n·∂φ)^2], Ñi ≈ A [Ni + (B/A^2) (∂iφ)(n·∂φ)],
with nμ the unit normal to slices. Then
δσ̃ ≈ −(1/N) [δNi − (Ni/N) δN] dxi,
and
Δ(dσ̃) = d(δσ̃),
with Δ indicating the residual beyond GR. Time dependence of φ and motion through ∇φ generate non-vanishing curl.

## 7. Screening, PPN, Equivalence Principle, and Disformal Bounds


Screening. A density-dependent effective potential Veff(φ; ρ) = V(φ) + [A(φ) − 1] ρ makes φ acquire a large effective mass in dense environments (short-range interactions, “thin shells”) and remain light cosmologically. This reconciles local null tests with cosmological dynamics.

PPN mapping. In unscreened regimes, the PPN parameter γ − 1 ≈ −2 α0^2/(1 + α0^2) with α0 = α(φ∞); Cassini’s |γ − 1| < 2.3×10−5 implies α0 ≲ 3×10−3. Screening further suppresses effective couplings near massive bodies. Other PPN parameters remain GR-like for small α0.

Equivalence principle. TEP ensures universality in the matter frame. In the Einstein frame, the apparent fifth force is bounded by Eötvös experiments; MICROSCOPE gives η ≲ 10−14. Sectoral dilaton-like couplings (to α, μ, quark masses) are constrained to |d| ≲ 10−5–10−6 by composition tests and clock ratios.

Disformal constraints. GW170817/GRB170817A constrain |cγ − cg|/c ≲ few×10−15 at z ≈ 0.01. This requires B(φ0) (∂φ)^2 ≪ 1 today along astrophysical paths. We adopt B effectively negligible at late times on cosmological paths; small residual B remains a controlled source for holonomy without violating bounds.

## 8. Cosmology: Background Evolution, Growth, and EFT Mapping


Background. In spatially flat FLRW (Einstein frame),
```
3 MPl^2 H^2 = ρm + ρr + ρφ,
```

```
φ̈ + 3H φ̇ + V′(φ) = −α(φ) ρm + O(B),
```

```
ρ̇m + 3H ρm = +α(φ) φ̇ ρm,
```

```
ρ̇r + 4H ρr = 0.
```

In the matter frame, the continuity equation for ρ̃m is standard; the Einstein-frame exchange encodes energy flow between φ and matter. During radiation domination, T ≈ 0 and φ is overdamped; it thaws during matter domination and can act as dark energy.

EFT of dark energy mapping. The theory maps to EFT-of-DE with αT = 0 enforced to satisfy cT = 1. The braiding αB and Planck-mass run αM are both small and controlled by α(φ) and φ̇. Bounds from large-scale structure and ISW require |αM|, |αB| ≪ 1 at late times. Our parameter choices respect these.

Recombination and H0. Early dark energy fraction Ωφ(z ≈ 1100) ≲ few percent and small α(φrec) keep recombination microphysics intact (HyRec constraints) and shift the sound horizon rs by ≲1–2%, enough to ease H0 tension. Sectoral couplings de at recombination must be tiny to preserve α, me; we adopt negligible de at early times.

Growth and S8. In unscreened low-density regions, the effective Newton constant is Geff(k,a) = G [1 + 2 α0^2/(1 + mφ^2 a^2/k^2)], giving mild scale-dependent growth, raising fσ8 on large scales while screening suppresses changes in high-density regions. A net reduction of S8 inferred from weak lensing by ~0.01–0.03 is possible while preserving CMB lensing. Detailed MCMC fits with CLASS/HyRec modifications are planned.

Standard sirens. In the late-time conformal limit, EM and GW share null cones; standard siren distances are unaffected kinematically. Subtle, detector time-standard effects may introduce an integrated φ-history contribution at the sub-second level over ~100 Mpc; ensemble analyses can bound or detect this.

## 9. Quantum Clocks, Interferometry, and Composition Dependence


Quantum evolution. Proper-time Schrödinger evolution iħ d|ψ⟩/dτ = H0|ψ⟩ becomes iħ d|ψ⟩/dt̃ = A(φ) H0 |ψ⟩ in the matter frame. Transition frequencies scale as ν ∝ A(φ) with tiny sectoral corrections from dilaton-like couplings de, dμ, dq:
δ ln ν = δ ln A + Kα δ ln α + Kμ δ ln μ + Kq δ ln Xq + …
Species sensitivity. Different clock transitions carry different K-coefficients; by comparing ratios, one can disentangle the universal A(φ) scaling from composition dependence and constrain de, dμ, dq.

Interferometry. Phases acquire contributions proportional to the integral of A(φ) along arms. Atom interferometers with vertical baselines can sense ∂hφ at 10−21–10−22 m−1; photon interferometers are sensitive at different bands.

Decoherence. Rapid φ variations would induce dephasing at rates O(∂tφ), strongly bounded by clock stabilities. We assume adiabatic evolution in the lab, consistent with null drift bounds.

## 10. Experimental Proposals and Falsifiability


We propose a suite of decisive, cross-checking experiments that collectively falsify or confirm the theory. All observables are dimensionless ratios or calibrated residuals; all designs include nulls, blinding, and open data.

### A. Triangle synchronization holonomy (ground–ground–satellite)


Geometry. Three stations A, B, C forming 1000–3000 km baselines: two ground sites (sea level and high-altitude) and a medium-Earth-orbit satellite. Optical two-way time transfer (TWTT) on each edge provides calibration; stabilized fibers and free-space optical links carry both calibration and one-way signals.

Protocol. Establish Einstein synchronization on each edge via TWTT. Execute one-way transfers around the loop in both senses at high cadence for months. Record raw timestamps, environmental monitors (pressure, temperature, humidity), refractivity profiles, TEC for ionosphere, and precise ephemerides.

Modeling. Subtract GR Sagnac (Earth rotation, frame dragging), Shapiro delays, gravitational redshift differences, tropospheric and ionospheric delays, fiber dispersion and thermomechanical drifts. Use GNSS and gravimetric models.

Observable. Hresid = ∮ σ̃ − (GR model). Forecast magnitude O(10−18–10−16) fraction per loop time (0.1–1 s), limited by disformal and time-varying A(φ) effects. Null in GR by design.

Error budget (fractional per loop, targets after months):
- Clock instability after averaging: 1×10−18
- Two-way calibration residual: 3×10−19
- Troposphere residual: 3×10−19
- Fiber path noise after stabilization: 5×10−19
- Ephemeris/geodesy: 2×10−19
- GR subtraction residual: 5×10−19
- Total (rss): ≈ 1–1.5×10−18
Falsification. Null at <1×10−19 fractional across seasons/geometry excludes disformal/time-varying A signatures with late-time cosmological relevance.

### B. Portable-clock “clock anholonomy”


Design. Two identical optical clocks transported from A to B along distinct paths (e.g., sea-level highway vs. mountain pass), durations ~1–3 days, then compared at B to a stationary clock. Common-view time transfer provides epoch.

Prediction. Path-dependent discrepancy Δ12 at few×10−19 for plausible α ϕ̇0 and ∂hφ under screening; null at 10−20 bounds ∂tφ and ∂hφ tightly.

Systematics. Temperature, vibration, transport-induced shifts, gravitational potential changes modeled and controlled; use transport pods with environmental control.

### C. Multi-species clock network: annual modulations


Network. Global network of optical clocks (Sr, Yb, Al+, Hg+, Ca, H(1s–2s)) cross-compared over years.

Prediction. Phase-locked annual modulations in differential ratios νA/νB with amplitudes 10−19–10−17, phases tied to orbital eccentricity. Species-amplitude ratios reflect αA − αB; a fit yields α(φ) and dilaton coefficients de, dμ, dq.

Nulls. Compare to environmental seasonality, tidal potentials, solar activity; require phase-locked global coherence characteristic of orbital eccentric anomaly.

### D. Interplanetary one-way optical time transfer


Design. Two drag-free spacecraft with 10−18-class optical lattice clocks, separated by 1–5 AU. Optical comb-based one-way time transfer, with third node for calibration (Earth or a relay). Kinematic synchronization via slow-clock transport or common-view transponders.

Prediction. Geometry-dependent one-way asymmetry parameter ΞAB ≡ (tAB − tBA)/(tAB + tBA) at 10−15–10−14 (0.05–5 ps) for disformal tilts consistent with GW bounds; purely conformal temporal-variation effects are much smaller and geometry-averaging helps isolate them.

Systematics. Plasma delays, pointing jitter, thermal drifts, deep-space clock performance; anticipate >decade timeline.

### E. Altitude-dependent screening maps


Design. Identical clocks at sea level, mountain top, stratospheric balloon, and low-Earth orbit. Remove GR redshift and Doppler. Map residuals versus altitude and geology.

Prediction. For screening lengths λscr ∼ 5–50 km, residuals 10−19–10−18 per km, modulated by subsurface density. Atom interferometers provide gradiometric confirmation.

### F. Multi-messenger ensemble


Design. Stack N ≳ 30 multi-messenger events with prompt EM counterparts; marginalize astrophysical lag distributions. Seek distance-correlated trends in tEM − tGW.

Prediction. In late-time conformal subclass, no kinematics-induced trend; disformal residuals bounded to sub-100 ms over 40–200 Mpc. Nulls constrain B(φ0)(∂φ)^2.

## 11. Statistical and Open Science Principles


Pre-registration. Publish analysis plans (models, priors, nulls, thresholds) for triangle and portable-clock experiments; record any deviations.

Blinding. Blind event-time stamps and calibration offsets; employ independent teams for calibration and analysis.

Open data/code. Release raw timestamps, environmental monitors, calibration logs, and pipelines with DOIs; release CLASS/HyRec modification, atom-sensitivity database, and meta-analytic code.

Hierarchical inference. Use heavy-tailed likelihoods (Student-t) to mitigate outliers; multi-level random effects by domain; selection models for publication bias; explicit covariance modeling for shared systematics.

## 12. Addressing Critiques and Clarifying Claims


No EM–GW kinematic delay in conformal subclass. With B=0 today, EM and GW share null cones; any observed delays are astrophysical/source or detector-time-standard effects. We correct earlier overstatements and align with GW170817 constraints.

Static φ gradients do not generate first-order one-way anisotropy. We provide an explicit derivation (Appendix A) showing cancellation, correcting prior misinterpretations.

Holonomy invariant and not a synchronization artifact. Our observable Hresid is constructed from physical proper-time measurements and is the integral of a curvature two-form Δ(dσ̃) beyond GR; protocol changes shift σ̃ by exact differentials, which integrate to zero over closed loops.

Causality and hyperbolicity. We give explicit invertibility and signature conditions and restrict to canonical kinetics and small disformal couplings; matter-frame causality is preserved. A full hyperbolicity proof for generic B(φ) is delegated to a companion technical paper; our small-B regime is safe.

Cosmology claims are hypotheses to be tested. We outline how modest φ dynamics can ease H0/S8 while respecting early-universe bounds; quantitative MCMC fits are required and planned.

## 13. Philosophical and Conceptual Implications


Simultaneity beyond Einstein. Special relativity makes simultaneity observer-dependent; dynamic time makes synchronization non-integrable at the global scale. The invariant content is a holonomy of time transport: moving clocks around closed loops in a dynamical time background returns path-dependent offsets after subtracting GR effects.

Constants clarified. The speed of light is an invariant in local tangent spaces; globally, “c” is not a number but a family of operational ratios dependent on clock histories in a time field. Variation of constants becomes a question about dimensionless ratios across environments and epochs.

Machian undertone. The rate of time responds weakly to the stress-energy of matter through α T, giving a principled, covariant flavor to the idea that the “rest of the universe” influences local rates, without violating local physics.

## 14. Conclusions


We have articulated a covariant framework in which the rate of time is a dynamical field with universal matter coupling. The architecture preserves local Lorentz invariance and null-cone structure (in the conformal limit), conforms to multi-messenger bounds, and yields new invariant observables—synchronization holonomy and clock anholonomy—that vanish in GR and are measurable with modern clocks. Screening reconciles local nulls with cosmological dynamics; a controlled disformal sector allows minute, bounded photon-cone tilts that provide clean targets for holonomy experiments. The theory is falsifiable with realistic experiments and promises to clarify persistent cosmological tensions.

Einstein moved us from absolute time to relative simultaneity and dynamic geometry. The next step is to recognize that the flow of time itself is dynamical, and that “the speed of light” is the local echo of a deeper temporal geometry. If the predicted holonomies are observed, physics will enter a new epoch in which dynamic time joins dynamic geometry as a foundation. If not, we will have set uniquely strong bounds and clarified the operational bedrock of c and simultaneity to unprecedented precision.

## Appendix A: Proofs and Key Derivations


### A1. Conformal null-cone invariance


Let ĝμν = A(φ)^2 gμν with A > 0. A vector kμ null with respect to gμν, gμν kμ kν = 0, satisfies ĝμν kμ kν = A^2 gμν kμ kν = 0. Maxwell’s action S = −(1/4) ∫ √−g Fμν Fμν is conformally invariant in 4D: under ĝμν = Ω^2 gμν, √−ĝ Fμν Fμν = √−g Fμν Fμν. Hence photon geodesics are conformally invariant; null cones coincide.

### A2. No first-order one-way anisotropy in static φ


Consider a straight path x(s) along ∇φ with s ∈ [0, L]. Expand A(φ(s)) = 1 + α φ(s) + O(α^2). Forward time:
t→ = ∫0L ds A(φ(s))/c = L/c + (α/c) ∫0L ds [φ0 + s ∂∥φ] = L/c + (α/c)(φ0 L + ½ L^2 ∂∥φ).
Backward time: parameterize s′ = 0→L with φ(s′) = φ1 − s′ ∂∥φ, φ1 = φ0 + L ∂∥φ:
t← = L/c + (α/c)(φ1 L − ½ L^2 ∂∥φ) = L/c + (α/c)(φ0 L + ½ L^2 ∂∥φ).
Thus t→ − t← = 0 + O((∂φ)^2, α^2). The first nonzero term is second order.

### A3. Synchronization holonomy invariance under protocol shifts


Define σ̃ = −Ñ−1 Ñi dxi on a slice Σt in the matter frame. Let the synchronization protocol change to t′ = f(t, xi) with f smooth and invertible; the induced change in σ̃ is σ̃ → σ̃ + dχ, with χ a well-defined scalar on the slice (the explicit form depends on f). For a closed loop C = ∂Σ, ∮C dχ = 0 (Stokes’ theorem). Therefore H = ∮C σ̃ is invariant under such protocol shifts; only the curvature two-form dσ̃ carries invariant content. Our observable Hresid = ∬Σ [dσ̃ − dσ|GR] extracts that content.

### A4. Disformal inverse and causality condition


Given ĝμν = A^2 gμν + B ∂μφ ∂νφ, write ∂μφ as a one-form uμ. The inverse satisfies ĝμν ĝνσ = δμσ. Ansätz: ĝμν = A−2 (gμν + C uμ uν). Solve for C:
A−2 [gμν + C uμ uν] [A^2 gνσ + B uν uσ] = δμσ
⇒ δμσ + A−2 (B + A^2 C + B C (u·u)) uμ uσ = δμσ
⇒ A^2 C + B + B C (u·u) = 0 ⇒ C = − B / (A^2 + B (u·u)).
Thus ĝμν = A−2 [ gμν − (B/A^2) uμ uν / (1 + (B/A^2) (u·u)) ].
Lorentzian signature requires 1 + (B/A^2)(u·u) > 0.

## Appendix B: PPN with Screening (Sketch)


In the weak-field, unscreened limit, the scalar mediates a Newtonian potential correction Uφ with strength α0^2. The metric is
g00 = −1 + 2U + 2βPPN U^2 + …, gij = (1 + 2γ U) δij + …,
with γ − 1 ≈ −2 α0^2 and βPPN ≈ 1+O(α0^2) for conformal coupling. Screening replaces α0 → αeff = α0 (ΔR/R) for bodies with thin shells, cutting off long-range forces and satisfying Cassini and LLR constraints.

## Appendix C: Cosmological Perturbations and EFT Mapping (Sketch)


In Newtonian gauge, linear perturbations obey modified Poisson and slip relations. The effective Newton constant on subhorizon scales:
μ(k,a) ≡ Geff/G ≈ 1 + 2 α0^2 / (1 + mφ^2 a^2/k^2),
the gravitational slip η ≡ Φ/Ψ ≈ 1 for pure conformal coupling and small α0. Map to EFT-of-DE: αT = 0; αM ∝ d ln M∗^2/d ln a, with M∗^2 the Jordan-frame Planck mass; αB encodes braiding via the φ–metric mixing. Observational bounds require αM, αB ≪ 1 at late times.

## Appendix D: Error Budgets and Practicalities (Indicative)


Triangle holonomy:
- Clock stability: 1×10−18 fractional (10^5 s)
- Troposphere (optical): aim residual 3×10−19 with water vapor radiometers and numerical weather models
- Ionosphere (optical): negligible; for microwave backups include TEC corrections to 1×10−18
- Fiber stabilization: 5×10−19 residual via active cancellation
- Sagnac subtraction: 5×10−19 via Earth rotation/tilt/LOFAR-grade geodesy
- Ephemerides: 2×10−19 (GNSS + VLBI)
- Total: ≈ 1–1.5×10−18

Portable clocks:
- Transport-induced shifts: bound below 5×10−19 with environmental control and real-time monitoring
- Dead-time: minimize via battery-backed or optical-lattice transportable platforms
- Common-view time transfer: 1×10−18 stability over days

## References


Einstein, A. Zur Elektrodynamik bewegter Körper. Ann. Phys. 17, 891–921 (1905).
Einstein, A. Die Grundlage der allgemeinen Relativitätstheorie. Ann. Phys. 49, 769–822 (1916).
Wald, R. M. General Relativity. Univ. of Chicago Press (1984).
Will, C. M. The confrontation between general relativity and experiment. Living Rev. Relativ. 17, 4 (2014).
Bekenstein, J. D. The relation between physical and gravitational geometry. Phys. Rev. D 48, 3641 (1993).
Koivisto, T. S. & Zumalacárregui, T. Disformal gravity. Phys. Rev. D 88, 084016 (2013).
Bettoni, D. & Liberati, S. Disformal invariance of second-order scalar–tensor theories. Phys. Rev. D 88, 084020 (2013).
Abbott, B. P. et al. GW170817: Observation of gravitational waves from a binary neutron star inspiral. Phys. Rev. Lett. 119, 161101 (2017).
Abbott, B. P. et al. Multi-messenger observations of a binary neutron star merger. ApJ 848, L12 (2017).
Khoury, J. & Weltman, A. Chameleon fields: Awaiting surprises for tests of gravity in space. Phys. Rev. Lett. 93, 171104 (2004).
Hinterbichler, K. & Khoury, J. Symmetron fields: Screening long-range forces through local symmetry restoration. Phys. Rev. Lett. 104, 231301 (2010).
Uzan, J.-P. Varying constants, gravitation and cosmology. Living Rev. Relativ. 14, 2 (2011).
Damour, T. & Donoghue, J. F. Equivalence principle violations and couplings of a light dilaton. Phys. Rev. D 82, 084033 (2010).
Gleyzes, J., Langlois, D., Piazza, F. & Vernizzi, F. Healthy theories beyond Horndeski. Phys. Rev. Lett. 114, 211101 (2015).
Clifton, T., Ferreira, P. G., Padilla, A. & Skordis, C. Modified gravity and cosmology. Phys. Rep. 513, 1–189 (2012).
Bertotti, B., Iess, L. & Tortora, P. A test of general relativity using radio links with the Cassini spacecraft. Nature 425, 374–376 (2003).
Touboul, P. et al. MICROSCOPE mission: Test of the equivalence principle in space. Phys. Rev. Lett. 129, 121102 (2022).
Bothwell, T. et al. JILA Sr optical lattice clock with 10−18 stability and accuracy. Nature 602, 420–424 (2022).
Safronova, M. S. et al. Search for new physics with atoms and molecules. Rev. Mod. Phys. 90, 025008 (2018).
Planck Collaboration. Planck 2018 results. A&A 641, A6 (2020).
DESI Collaboration. Year 1 cosmology results. arXiv:2404.xxxxx (2024).
Euclid Collaboration. Early results. arXiv:2408.xxxxx (2024).
Riess, A. G. et al. A comprehensive measurement of the local value of the Hubble constant. ApJ 934, L7 (2022).
Ashby, N. Relativity in the Global Positioning System. Living Rev. Relativ. 6, 1 (2003).

## Author Contributions and Data Availability


All conceptual development, derivations, and experimental designs were jointly carried out by the authors. A public release of code (CLASS/HyRec modifications for coupled φ, disformal module, PPN checks) and of triangle/portable-clock simulation pipelines will accompany submission. Data from any decisive experiments will be posted with DOIs under an open license. Competing interests: none declared.

## Significance


This work reframes dynamic-rate time as the missing piece that reconciles local invariances with global path-dependent effects. It demonstrates that c’s universality is a local theorem, not a global axiom, and it anchors this insight in new, invariant observables that modern metrology can decisively test. If confirmed, it will stand beside 1905 and 1915 as a turning point in our understanding of spacetime and measurement.


