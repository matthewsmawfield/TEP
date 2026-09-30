# Master-Sector Constraint Ledger (AUD-0 / AUD-1)

Plan references: `notes/plan.md` — finding AUD-0 (master sector), attack
protocol T0.1–T0.4; finding AUD-1, task T1.0.  Doctrine D7: functional
forms are outputs of a constraint search, not inputs.

Numerical admissibility results cited below are produced by
`scripts/steps/step_54_master_sector_admissibility.py` →
`results/step_54_master_sector_admissibility.json`.

---

## 1. Fixed inputs (not candidates)

| Object | Value | Source |
|---|---|---|
| Metric map | g̃ = A²(φ)g + B(φ)∇φ∇φ | Rule 1; Paper 0 §1 |
| Conformal coupling | A(φ) = e^{β_A φ/M_Pl}, β_A = −1 | Rule 3 |
| Sign convention | φ > 0 in wells; φ = 0 ambient today | Rule 4 |
| Matter-frame transport | proper time of g̃ | Rule 2 |
| Linear response | δφ/M_Pl = 2β_A Φ_N; ∇ln Ñ = (1+2β_A²𝒮_Σ)∇Φ_N | step_26 |

## 2. Ledger A — conditions on P(X, φ)

| ID | Source | Condition |
|---|---|---|
| P-1 | Paper 0 §4 | Hyperbolicity / well-posedness on both X signs |
| P-2 | Paper 0 §4 | No ghosts: perturbation kinetic coefficient > 0 |
| P-3 | Paper 0 §4 | Stability: Legendre factor P_X + 2X P_XX > 0 |
| P-4 | Paper 0 §7 | Large-|X| screening: 𝒮_Σ = 1/P_X → 0 in wells |
| P-5 | Recovery F1 | Shear channel must act as gravity at low g |
| P-6 | Recovery F4 | Solar-System benchmarks recover the incumbent 𝒮_Σ = (1+ξ)⁻¹ predictions (Cassini δa ≲ 5×10⁻¹³ m s⁻² class) |
| P-7 | plan AUD-1 | Small-ξ flat tail: point-source a_φ ∝ 1/r (deep-MOND/RAR asymptote g_obs ≈ √(g_N a_eff)) |
| P-8 | plan AUD-1 | a_eff must equal the measured a₀ ≈ 1.2×10⁻¹⁰ m s⁻² = 0.366 g_t (WB plateau α_sat) with no free scale |
| P-9 | Rule 12 | Same P(X) on the cosmological branch X > 0 — ambient drift 𝒮_Σ,cosmo(z) |
| P-10 | Paper 0 App E | Single kinetic scale: Λ⁴ = M_Pl²H₀², Λ ≈ 1.91 meV |
| P-11 | plan AUD-8 | Transition-region exponent must reproduce the observed pairwise profile (s^{4/3}-class law, step_30) |
| P-12 | Rule 20 | Interior roll of deep wells must remain admissible (R3 conditions propagate through P_X) |

### Admissibility results (step_54, section A)

Candidate space: `baseline 1+2ξ`; `two_branch k√ξ+2ξ`;
`twobranch_u 1+k√ξ+2ξ`; `exp_interp 1−e^{−k√ξ}+2ξ`; `aqual_only k√ξ`;
power-law family `kξ^p+2ξ`, p ∈ {1/3, 1/2, 2/3, 1}.

**Tail-exponent uniqueness (partial uniqueness theorem).**
For `P_X ~ kξ^p`, flux conservation `r²P_Xφ' = C` gives
`φ' ∝ r^{−2/(2p+1)}`. A flat rotation asymptote (`a_φ ∝ 1/r`)
requires **p = 1/2 uniquely**.  Numeric scan confirms:
p = 1/3 → −1.200, p = 1/2 → −1.000, p = 2/3 → −0.857, p = 1 → −0.667.

**Deep-MOND coefficient (derived, supersedes plan's k² ≈ 11 estimate).**
```
a_φ² = a_eff g_N,   a_eff = 2√2 |β_A|³ Λ²/(k M_Pl) = 4√2 |β_A| g_t / k
k = 4√2 g_t/a₀ ≈ 16.0
```

**Per-candidate verdicts:**

| Candidate | Flat tail | P_X>0 | Legendre>0 | c_s² | Verdict |
|---|---|---|---|---|---|
| baseline 1+2ξ | ✗ (slope −2) | ✓ | ✓ | [1/3, 1] | REJECTED — Keplerian tail, fails RAR (P-7); this is AUD-1's bug |
| two_branch | ✓ (−1.000) | ✓ | ✓ | [1/3, 1/2] | ADMISSIBLE |
| twobranch_u | ✗ (−2.000) | ✓ | ✓ | [1/3, 1] | REJECTED — unit term restores linear small-ξ limit (P-7) |
| exp_interp | ✓ (−1.000) | ✓ | ✓ | [1/3, 0.78] | ADMISSIBLE — transition-region rival to two_branch; discriminant at ξ ~ 1 |
| aqual_only | ✓ (−1.000) | ✓ | ✓ | 1/2 | REJECTED — no large-ξ screening (P-4) |

**Branch crossover.** `ξ* = (k/2)² ≈ 64`.  Environments:

| Environment | ξ | branch | benchmark shift |
|---|---|---|---|
| Cassini limb / Earth / Saturn / LLR / 1 AU | 4×10¹⁰–1×10²³ | large-X | ≈1.000 (P-6 safe by construction) |
| GC core | 4.2×10⁴ | large-X | 1.04 |
| Wide binary 2646 AU | 9.5 | **small-X** | 3.4 |
| MW solar circle (ambient) | 0.35 | **small-X** | 6.0 |
| Void/web | 8.6×10⁻⁶ | small-X | 𝒮_Σ ≈ 21 → anti-screening (Rule 10 consistent: faster in voids) |

The pair-transition region sits at the branch crossover by construction
of k; the s^{4/3} interior profile (mutual-dominated, ξ_pair > ξ*) and
the saturated plateau both survive, but the vertex factors
𝒮_Σ(X_env)² are now `≈0.107²` rather than `0.59²` — the external-field
suppression emerges naturally (this is the AQUAL-type EFE).  T-W3/T1.3
must recompute the plateau normalization and fitted transition radius.

**Cosmological branch (ξ_cosmo = H²(z)/(2H₀²)).**

| z | 𝒮_Σ baseline | G_eff/G baseline | 𝒮_Σ two-branch | G_eff/G two-branch |
|---|---|---|---|---|
| 0 | 0.500 | 2.00 | 0.081 | **1.16** |
| 1 | 0.244 | 1.49 | 0.043 | 1.09 |
| 3 | 0.048 | 1.10 | 0.014 | 1.03 |
| 10 | 2.5×10⁻³ | 1.005 | 1.6×10⁻³ | 1.003 |
| 100 | 3.2×10⁻⁶ | 1.000 | 3.2×10⁻⁶ | 1.000 |

The two-branch completion weakens the late-time growth enhancement from
G_eff ≈ 2G to ≈ 1.16G — direction consistent with structure-growth
measurements sitting below ΛCDM; this is the AUD-3 discriminant.
`exp_interp` retains the baseline ≈2G boost — a clean discriminator
between the two admissible candidates (transition-region shape vs
cosmological boost).

## 3. Ledger B — conditions on V(φ)

| ID | Source | Condition |
|---|---|---|
| V-1 | Recovery F1 | Admissible minimum on the β_A = −1 branch: V_eff = V + ρe^{−u} |
| V-2 | Rule 11 | Density-dependent m_eff (environmental amplitude screening) |
| V-3 | R3 / closure_status | Locally flat floor at deep u — constant V_,u + Q_m rolling-floor self-consistency (7×10⁻¹⁴) |
| V-4 | R4 / step_12 | Closed-static benchmark V → 2M_Pl²/a² floor identification |
| V-5 | Rule 9 | Cosmological range Δu ~ ln(1+z); ambient drift A(η) from the same V |
| V-6 | Rule 20 | Floor V_0 supplies scale-invariant well-core density; N_min > 0 |
| V-7 | step_20 | Void sector must not run secularly (|Δu| < 10⁻²) |
| V-8 | Paper 0 §4 | V'' > 0 at matter-domain minima |
| V-9 | Rule 23 | Potential energy density accommodates Planck-scale well interiors |

### Candidate results (step_54, section B)

| Candidate | Interior min (ρ range) | m_eff ∝ ρ^p | Verdict |
|---|---|---|---|
| quadratic m²u²/2 | 4/4 | p ≈ 0.015 ≈ 0 | REJECTED — fixed λ_C kills environmental screening (V-2) |
| inverse-power Λ⁴(1+1/u) | 0/4 (all runaway) | — | REJECTED — V-1 fails on β_A = −1 branch, confirming the registry note numerically |
| quartic λu⁴/4 | 4/4 | p ≈ 1/3 matter-domain | VIABLE — plus emergent feature: for ρ ≳ λ-scale densities the minima saturate **logarithmically**, u_min ~ ln(ρ/λ) — self-limiting well depth |
| exponential Λ_V⁴e^u | 4/4 | p ≈ 0.4 (transition-mixed) | VIABLE runner-up — emergent feature: ambient-pinned (u → 0) for ρ ≲ Λ_V⁴, interior min for denser — a natural φ_ambient = 0 pin |
| symmetron −μ²u²/2 + λu⁴/4 | 4/4 | p ≈ 0.15 | REJECTED — vacuum (ρ→0) sits at u = v ≠ 0: breaks the φ_ambient = 0 convention and puts dense media at the symmetric point — wrong environmental direction (V-5/Rule 4) |
| **master family λu⁴/4·e^{−(u/u_s)⁴} + V_0e^{−(u_s/u)⁴}** | 4/4 | p ≈ 0.5 | **SELECTED** — only candidate meeting V-1–V-9 jointly |
| constant V_0 | 0/4 | — | REJECTED alone (no matter-domain min) — but is the u→∞ limit of the master family |

**Master-family structure (computed, step_54):** the quartic term governs
u < u_c and the floor term V_0e^{−(u_s/u)⁴} governs u > u_c, where the
crossover solves λu_c⁴/4·e^{−(u_c/u_s)⁴} = V_0e^{−(u_s/u_c)⁴}:
u_c ≈ 2.82 for V_0 ~ M_Pl⁴ (u_c ≈ 3.26 for V_0 ~ 10⁻³⁰).  Matter-domain
minima (u_min ~ 5×10⁻⁴ at Earth density) sit deep in the quartic regime;
the well-interior field (u ~ u_s ≈ 9–10, step_16/R3) lives on the floor —
which is precisely the Rule-20 role.

**Global-floor question (AUD-0 open structural point) — now sharpened.**
In the master family, V(u→∞) = V_0 is the *same* limit the ambient
power-law branch approaches (A_clock → 0 ⟺ u → ∞): the floor is global,
not merely effective.  On the flat floor V_,u → 0 and the roll
asymptotes — A → 0 is approached over unbounded coordinate time and
never completed, while well cores pin at finite u through **source
starvation** (matter coupling ρA_,φ ~ −ρe^{−u} → 0 at u ≫ 1 starves the
drive; verified structurally — the quartic's own u_min saturation is the
same mechanism).  Whether this yields Rule 20's N_min > 0 verbatim
(finite minimum) or an asymptotic approach is the T-B1 numerical target;
note the divergence between `tep-rules.md` (N_min > 0 finite) and
`.devin/rules/tep.md` (asymptotic approach) must be resolved against the
BVP result — the master family favors the asymptotic reading, with the
core *value* of A finite at finite radius.

## 4. Ledger C — conditions on B(φ)

| ID | Source | Condition |
|---|---|---|
| B-1 | Paper 0 §4 | Signature: B(∂φ)² > −A² pointwise on realized profiles |
| B-2 | Paper 0 §4 / step_52 | Hyperbolicity: Z_t, Z_s > 0; B ≥ 0 unconditional, B < 0 pointwise — all realized negative-branch reconstructions fail (Q → −1) |
| B-3 | GW170817 / step_52 | B → 0 at φ ~ 0; path integral of B(∂φ)² below the \|c_g − c_γ\| ~ 10⁻¹⁵ bound |
| B-4 | Rule 13 | Holonomy admissibility: nonzero H_resid requires spatially varying (B/A²)(n·∂φ) — conformal-only sector gives identically zero |
| B-5 | Rule 16 | g-sector wave stays a pure tensor ripple; B must not inject scalar polarization at leading order |

### Candidate gate table (step_54, section C)

| Candidate | B(0)=0 | Verdict |
|---|---|---|
| B ∝ φ² | ✓ | REJECTED — Solar-System cone tightness (Paper 29) |
| B ∝ φ⁵ | ✓ | REJECTED — signature loss on z ≳ 3.4 paths (Paper 29) |
| bump B_0φ²e^{−φ²/φ_c²} | ✓ | ADMISSIBLE — B_0 window [2.1×10²¹, 6.1×10¹⁶⁴]; fixed-action normalization open (R6) |
| Paper-28 envelope B_0u²/(1+u²)·e^{−u⁴/2σ_B⁴} | ✓ | ADMISSIBLE — prescribed strong-field realization; weak-field ≈ B_0u² |
| B < 0 branch | — | CONDITIONALLY EXCLUDED — own null cone fails on every realized reconstruction (R5, step_52) |

The surviving class is the **bump/envelope family** — B(0)=0 enforced by
GW170817, mid-field support required by Rule 13, large-field damping by
signature/hyperbolicity margins.  Constrained, not free; amplitude
carrier not yet derived from the action.

## 5. Ledger D — A(η) ambient drift (AUD-4 junction)

| ID | Source | Condition |
|---|---|---|
| A-1 | Rule 9 | Endpoint clock ratio A_clock = (1+z)⁻¹; power-law benchmark A_clock ∝ η^{−p}, p ≈ 1 (paper 18: ln(1+z) ≈ 0.96 ln η) |
| A-2 | Paper 0 §8 | Ambient branch sees the u→∞ limit — the V_0 floor domain |
| A-3 | AUD-4 | Junction between drift and the screened local sector must be one junction, not a per-channel match |

### T4.1 result (step_64) — the junction window and the unique reconciler

`step_64_ambient_junction.py` → `results/step_64_ambient_junction.json`.

Continuity of the field and its flux (C¹ in A_clock) fixes the
power-law tail index to the local slope of the conformal image at the
junction:

  p_eff(z) = [E(z)/(1+z)] · χ(z),   χ(z) = ∫₀^z dz'/E(z')
           → 0.68 at z=1,  → 2 in the matter era,  → 67 at z=1100.

Proposition 1 requires the tail index 0 < p ≤ 1/2.  Therefore:

- **z_j ≤ 0.687** — the latest admissible C¹ junction sits at p = 1/2
  exactly (boundary of the regularity window).  Reference choice in
  the JSON: z_j = 0.6 → p = 0.447, C_tail = 0.4646.
- **Both legacy conventions are irregular**: Paper 18's full
  conformal identification has p_eff → 67 at recombination; the
  p ≈ 1 power-law benchmark sits outside the window.  The piecewise
  A(η) is the unique construction satisfying Rule 9 *and* temporal-
  horizon regularity — this resolves the Papers-18/27
  incompatibility structurally, not by preference.
- The junction carries a genuine curvature kink (ΔA'' ≈ 1.03 at the
  reference junction — a localized pulse at finite η, not a tail
  term).  Proposition-1's s ≥ 2p condition is satisfied on the tail
  (the clock field's own kinetic curvature supplies s = 2).
- Tail redshift map: recombination lands at η ≈ 1.14×10⁶ H₀⁻¹ on the
  p = 0.447 tail vs χ_ΛCDM ≈ 3.18 — the eternal-past reading is
  quantitatively realized; the ambient departs the ΛCDM image at
  z ≈ 0.6–0.7, and ALL earlier physics (CMB, BBN, structure) is
  carried by the well network, per plan condition (c).
- Consistent with the master family: V₀ is the u→∞ floor the
  η^{−p} roll approaches — the 𝒯⁻ asymptote is the floor, approached
  never completed (the T-B1 asymptotic Rule-20 reading propagates to
  the ambient branch).

### 5b. T4.2 — acoustic-side bound on z_j — RESOLVED (no acoustic bound)

`step_65_acoustic_junction.py` → `results/step_65_acoustic_junction.json`.

The acoustic sector does NOT bound z_j from below.  Three established
facts decide it:

1. **Deep-ambient well starvation (measured).**  Re-running the step_58
   well solve with u(r_max) = u_amb > 0 gives an exponential starvation
   law for the local well depth:

      u_loc ≈ u_loc(0) · e^{−0.45 u_amb}

   Wells nested on deep ambient are almost entirely ambient-depth: a
   source observed at z = 1100 (u_emit = ln 1101 = 7.00 — fixed by the
   Rule-9 endpoint ratio, independent of the ambient map) decomposes as
   u_amb ≈ 6.96 + u_loc ≈ 0.05, i.e. ~99.3% ambient.  This is the
   quantitative realization of Rule 23 (high-z emitters sit deep inside
   temporal wells on the ambient landscape).

2. **Acoustic observables are relational.**  r_s in matter units is
   invariant under a uniform clock shift of the emitter's environment;
   the ambient map enters only through well-interior structure
   (starvation) and transport — not the sound horizon itself.  A naive
   FLRW-image reading, where η_lb(z) maps directly to observed
   distances, would fail catastrophically (the tail gives
   η_lb(z = 1100) = 1.14×10⁶ H₀⁻¹ vs the ΛCDM-image 3.18) — but that
   reading is precisely the pre-derivation ansatz the Projection
   Dictionary retires: distances come from the constraint-slice
   geometry, not the ambient bookkeeping map.

3. **Binding constraint on z_j.**  Regularity only: z_j ≤ 0.687
   (step_64).  The physical z_j is where the well network hands the
   ambient its drift (the junction-kink reading already in step_64).
   The remaining acoustic dependency — well-frame r_s and peak
   morphology under deep-ambient nesting — is the AUD-6/T6.3
   constraint-slice computation, not a bound on the ambient map.

## 6. T0.4 status — uniqueness statement

- **P_X:** unique power-law tail (p = 1/2) + incumbent large-X term +
  single calibrated coefficient k ≈ 16.0.  Within the surveyed power-law
  space this is a **unique minimal completion**; `exp_interp` survives as
  a one-parameter rival discriminated by (i) ξ ~ 1 transition shape and
  (ii) G_eff(z=0) ≈ 2.0 vs 1.16.  Wider K(X) uniqueness (cuscuton-class
  members carry no propagating scalar dof) remains open per Paper 0 §2.2.
- **V:** master family is the unique surveyed member meeting V-1–V-9;
  exponential is the viable runner-up to test against the ambient drift.
- **B:** bump/envelope class survives; amplitude normalization is the
  single open number.
- **A(η):** power-law benchmark on the floor domain; junction to the
  screened sector is AUD-4.

## 7. T1.1 result — SPARC RAR under the two-branch sector

`TEP-UCD/scripts/steps/step_10_two_branch_sparc_rar.py` →
`TEP-UCD/results/outputs/step_10_two_branch_sparc_rar.json`
(migrated from Paper 0 step_55; SPARC Table 2 mass models, TEP-UCD: 3,389 points / 175 galaxies,
Υ_d = 0.5, Υ_b = 0.7).

Derived prediction (zero free parameters — k = 16.03 fixed by a₀):

  g_obs = g_N + 2g_t·ū,   (k/√2)ū² + ū³ = g_N/g_t
  deep tail: g_obs² = a_eff g_N;   V_flat⁴ = GM a_eff   (BTFR)

| Statistic | two-branch | baseline 1+2ξ | empirical MOND ν ref |
|---|---|---|---|
| RMS log₁₀ resid | 0.219 | 0.281 | 0.199 |
| median log₁₀ resid | −0.086 | −0.032 | −0.014 |
| frac within 0.1 dex | 0.447 | 0.235 | — |

Free-k scan: `k_fit ≈ 30.9` → `a_eff ≈ 6.2×10⁻¹¹ ≈ 0.52 a₀`
(the fitted normalization is ~a₀/2).

**Binned residuals are uniform (ratio ≈ 0.80 ± 0.05 over −12.0 < log₁₀g_bar
< −8.4)** — the derived *shape* tracks the RAR across four decades; the
residual is an ~20% normalization offset, not a shape failure.  The deep
tail (g_bar < 0.3a₀, n = 2061) against the pure `√(a₀g_N)` form has
median offset +0.031 dex — the additive `g_N + a_φ` term is what places
the transition bins ~20% high.

**Identified candidate resolutions for the normalization gap** (in order
of expected size):

1. *Ambient-field suppression (EFE).*  The algebraic law is the
   isolated-source channel.  Real galaxies sit in a cosmic-web ambient
   ξ_env; the nested operator adds vertex factors
   𝒮_Σ(ξ_env)² ~ (0.1)² — a uniform suppression in exactly the right
   direction (data want less boost).  Requires the embedded-propagator
   treatment (extend step_32 to the two-branch flux law) — this is the
   T1.3 target.
2. *Disk geometry.*  The spherical flux `r²P_Xφ' = β_AM(r)/4πM_Pl`
   underestimates the planar field gradient for a thin disk (the flux
   does not spread over 4πr²) — an order-unity downward correction.
3. *Υ systematics.*  ±20% mass-to-light shifts translate directly into
   g_bar.

The 4651 vs 2646 AU wide-binary discrepancy is recast by the branch
split: the WB transition sits on the small-X branch (ξ ≈ 9.5 < ξ* = 64),
so the transition profile and plateau normalization must be recomputed
under the new vertex factors — the old number is superseded pending
T-W3.

### 7b. T1.5 result — BTFR under the two-branch sector

`TEP-UCD/scripts/steps/step_13_two_branch_btfr.py` →
`TEP-UCD/results/outputs/step_13_two_branch_btfr.json` (migrated from
Paper 0 step_66; SPARC Tables 1+2, 135 galaxies with measured Vflat).

Two zero-parameter tests:

- **Outer-edge (honest) test** — predict v at each galaxy's last 3
  measured radii from the mass-model g_N through the exact flux law:
  median v_obs/v_pred = 0.926, rms log₁₀ residual 0.089 dex, median
  g_N/g_t = 0.022 — zero-parameter velocities within ~7% at the outer
  edge, well inside the mass-model Υ systematics.
- **Deep-limit BTFR** — naive v_flat⁴ = G M_b a_eff
  (a_eff = 4√2 g_t/k = 1.08×10⁻¹⁰ m/s² = 0.90 a₀_MOND) under-predicts
  v by ~9% median because SPARC outer radii sit at ū ~ 0.1–0.3, where
  the finite-ū correction to g_obs²/(g_N a_eff) runs to 1.9–3.1 — the
  deep-limit formula is not the right comparison at those depths.
  Empirical M_b–v fit slope 3.19 (shallower than 4 — expected since
  low-mass curves have not reached their asymptote); median
  per-galaxy implied a_eff = 1.52×10⁻¹⁰.

Interpretation: the sector's *derived* acceleration scale
a_eff/g_t = 4√2/k = 0.353 ≈ a₀/g_t lands the BTFR normalization in
the right neighbourhood with no fitting; the ~7–9% residual is
consistent with the RAR normalization offset already characterised
(same ū regime) and with Υ mass-model systematics.  Note the RAR
itself is already a fully blind prediction (k was derived from
g_t/a₀, never fit to SPARC), so the whole 3,389-point curve is the
T1.6 "blind curves" content; a held-out-galaxy split adds no new
information beyond the existing scatter.

## 8. T1.3 result — embedded propagator under the two-branch sector

`step_56_embedded_two_branch.py` → `results/step_56_embedded_two_branch.json`
(port of step_32's nonlinear elliptic solve, cylindrical (s,z) grid,
finite-volume flux `J_i = a_i f(|a|)`, Newton + sparse Jacobian).

**Port validation:** baseline `f(q) = 1+q²` control reproduces step_32's
embedded far-field suppression exactly (y_far = 0.558 vs 0.56).

**Two-branch sweep** (orientation-mean far-field response y_far, and the
pure-propagator plateau `α_sat = √(1+2y_far)−1`):

| ambient u0 (g_t units) | P_X(ambient) | y_far | α_sat |
|---|---|---|---|
| 0.05 | — | 1.22 | 0.85 |
| 0.12 | — | 0.566 | 0.46 |
| **0.15** | 1.72 | **0.456** | **0.38** |
| 0.30 | 3.5 | 0.227 | 0.21 |
| 0.570 (solar circle) | 6.8 | 0.116 | 0.11 |
| 0.721 (X_GAL-calibrated) | 8.7 | 0.091 | 0.086 |
| 1.10 | — | 0.057 | 0.054 |

At u0 = 0.15 the half-excess radius is 0.785 r_* — the transition
compresses (baseline embedded gave ~0.6–0.8 r_* scale ratios; the fitted
R_s then lands in the ~3.3–3.7 kAU band against the observed 2646 ± 609,
comparable to the incumbent's post-propagator ratio ~1.26–1.4).

**Estimator-level verdict (step_019, same estimator as TEP-WB
step_017 — same data, mass draw, normalisation, canonical fit):**

JOINT (α, R_s) CLOSURE FAILS AT EVERY AMBIENT, under every reading:

| u0 | pure (α, R_s) | consistent-vertex (α, R_s) |
|---|---|---|
| 0.12 | 0.518, 7232 | 0.302, 7669 |
| 0.15 | 0.423, 5844 | 0.161, 6258 |
| 0.30 | 0.222, 3138 | 0.020, 3283 |
| 0.57 | 0.118, 1950 | 0.003, 1994 |
| 0.72 | 0.093, 1732 | 0.001, 1762 |
| 1.10 | 0.060, 1572 | 0.000, 1589 |

The two observables pull the ambient in **opposite directions**:
α = 0.366 wants u0 ≈ 0.19 (where R_s ≈ 5,100 AU — ~4σ high);
R_s = 2,646 wants u0 ≈ 0.37 (where α ≈ 0.17 — ~16σ low).  The
small-X branch gives a slow algebraic far-field tail that the
canonical exponential fit reads as a stretched R_s through the mass
convolution.

**Honest incumbent comparison.**  Running the BASELINE sector through
the same estimator at u_sun = 0.570:

- pure-propagator: α = 0.567 (17σ overshoot), R_s = 3,602 (1.6σ)
- consistent-vertex: α = 0.354 ✓, R_s = 3,748 (1.8σ) — the corpus's
  published quasi-pass.

*Re-reading the vertex (correction to the earlier "double-counted"
characterisation).*  The factorisation 𝓡 = 𝒮_env² · y(s) is NOT a
double count: the propagator y(s) measures the mutual field
transmission, while 𝒮_env per star is the object's effective scalar
CHARGE suppression — the same physics the starved two-centre solve
computes explicitly as ρ e^{−φ} (step_59).  They are distinct
channels (source vs propagator) and the corpus convention is
defensible.  What the convention omits is that the objects' charge
suppression is set by the TOTAL local field — ambient level, ambient
gradient, and companion well — which step_59 now computes rather than
approximating by a uniform 1/P_X(u0).

**Consequence — REVERSED by the ambient-semantics audit (step_020,
TEP-WB).**  The joint-closure failure above was manufactured by
three estimator errors, not by the sector:

1. **u0 was the wrong variable** (suspect 1 — reversion).  The
   embedded solve's boundary gradient is the scalar field's own
   gradient ū (the Temporal Shear), not the total solar-circle
   acceleration.  The corpus's u_sun = 0.57 set u0 = a_sun/g_t =
   1.94e-10 m s⁻² — the FULL force including the g-sector Newtonian
   part.  The self-consistent ambient solves
   ū(P_X(ū²/2) + 2β²) = g_obs,sun/g_t = 0.631, giving u0 = 0.163
   (two-branch), 0.207 (baseline), 0.213 (exp_interp).  The implied
   baryonic share of v_c² at R☉ is a falsifiable side output:
   two-branch 48% (v_b ≈ 162 km/s — consistent with the MW baryonic
   inventory), baseline 34% (136 km/s — low).
2. **The consistent-vertex q² = 𝒮_Σ(u0)² double-counts under the
   resolved machinery** (suspect 2 — pre-derivation ansatz).  y_emb
   already carries the ambient suppression once; a per-star charge
   vertex is set by each star's own well depth, not by the ambient
   stiffness.  Pure propagator is the honest reading of the resolved
   solve; all three readings are still reported.
3. **Deprojection was missing** — the observed ṽ(s) bins on projected
   separation while the propagator lives on d ≥ s (⟨d⟩ ≈ 1.57 s);
   convolving y(d/r*) over the isotropic kernel reduces the fitted
   R_s by ~30%.

**Corrected joint test (step_020, derived ambient, pure propagator,
deprojected):**

| sector | u0_derived | α_pred | σ_α | R_s (AU) | σ_Rs | joint χ² |
|---|---|---|---|---|---|---|
| **two-branch** | **0.163** | **0.371** | **0.44σ** | **3676** | **1.69σ** | **3.06** |
| baseline | 0.207 | 0.719 | 29σ | 3202 | 0.91σ | 867 |
| exp_interp | 0.213 | 0.731 | 30σ | 3263 | 1.01σ | 929 |

**Verdict: PASS for two-branch — the WB channel now SELECTS the
derived sector.**  The amplitude lands at 0.4σ with zero free
parameters, close to the plan's predicted identity
α_sat ≈ a_eff/g_t = 4√2/k = 0.353; the scale sits 1.7σ inside the
systematic band.  Baseline and exp_interp are excluded at ~30σ on
amplitude — the channel that was the sector's liability is now its
sharpest discriminator in the right direction.  The earlier baseline
"quasi-pass" (α = 0.354 at u0 = 0.57, corpus vertex) was the two
errors partially cancelling: the inflated ambient over-suppressed the
propagator, and the squared vertex suppressed it again.

**Refined estimator (step_020b/020c) — the scale residual closes.**
The two-centre solve at the derived ambient (020b) shows the parallel
pair is pessimistic mid-transition (each star sits inside the other's
compressed axial wake: y_pair(d=1) = 0.15), while the perpendicular
channel sits beside the deep wake and recovers by s ~ 2r*
(y_perp ~ 0.44).  The isotropic mutual-force estimate
⟨y⟩ = (2/3)y_perp + (1/3)y_par, marginalised over the real per-bin
component masses with per-star r*ᵢ = √(Gmᵢ/g_t) — each star's own
well sets its response scale, shifting the transition √2 earlier
than the total-mass convention — gives:

| run | α_pred | σ_α | R_s (AU) | σ_Rs |
|---|---|---|---|---|
| component masses (m_corr), orient-avg pair | 0.347 | 1.57σ | 2940 | 0.48σ |

versus observed α = 0.366 ± 0.012, R_s = 2646 ± 609 AU (systematic
band).  The outer plateau itself is 1.35–1.36 predicted vs 1.37
observed.  Baseline-window sensitivity on the DATA side sharpens the
agreement further: the innermost bin (59 AU, N = 128, residual −0.10)
carries documented unresolved-companion dilution that deflates the
normalisation; excluding it from the baseline window refits the
observed profile to α = 0.346–0.350, R_s = 3071–3168 AU — within
~0.3σ of the prediction on BOTH parameters.

The canonical-fit χ² is no longer the honest test (the exponential
family cannot hold both the predicted transition shape and plateau;
direct profile comparison at the outer bins gives 1.356 vs 1.366).
The joint two-parameter test passes at ~1.5σ / ~0.5σ, and ~0.1σ /
~0.5σ under the contamination-aware normalisation — zero free
parameters throughout: u0 derived from the solar-circle flux law, the
profile from the embedded solve, the masses from the catalogue.

**Environmental ordering under corrected semantics (step_020c/d/e).**
The step_005 strata — midplane (z_med = 47 pc): α = 0.244 ± 0.006,
R_s = 2821 ± 214 AU; halo (z_med = 248 pc): α = 0.401 ± 0.018,
R_s = 4681 ± 428 AU — each close with a SINGLE ambient under the
per-stratum forward model (own mass marginals, own two-centre pair
profile solved at the stratum ambient, isotropic orientation
average): midplane u0 = 0.20 predicts (0.236, 2612) — 1.3σ / 1.0σ;
halo u0 = 0.11 predicts (0.414, 4514) — 0.7σ / 0.4σ.  The pair
suppression deepening at higher ambient steepens the α(u0) map, so
the required environmental contrast is ~1.8 (ambient law
u_eff ∝ ρ^0.67), down from the monopole-map estimate of 2.0.  The
physical channels bracket it: global Galactic shear has the wrong
sign (0.95), local kpc-scale generated field ~1.1, k-branch planar
projection ρ^{1/2} gives 1.56, local-column ρ^1 gives 2.44.

Distribution-averaging sharpens the required law: evaluating
ρ(u_eff)^n over the actual |Z| distributions (halo stratum p84 = 404 pc,
max ~1 kpc) rather than at the stratum medians gives contrast 1.63
(n = 0.5), 1.90 (n = 0.67), 2.10 (n = 0.78) — i.e. the required
exponent drops to n ≈ 0.6, midway between the k-branch planar
projection (n = 0.5, contrast 1.63 — within ~10%) and the local-column
projection (n = 1.0, contrast 2.55).

**Field-magnitude channels excluded definitively — including the
exact nonlinear superposition.**  In the k-branch the flux
F = (k/√2)|∇u|∇u obeys a LINEAR Gauss law, so F(x) is exactly the
Newtonian-like field of the source set and |∇u| = √(√2|F|/k) is
the exact superposition — no approximation.  Computed over the
disk well network plus smooth components (step_020h field_census):

| component | midplane | halo (z=347) | verdict |
|---|---|---|---|
| vertical F_z-derived | 0.030 | 0.086 | rises — wrong sign |
| discreteness (nearest wells) | 0.0041 | 0.0030 | ~50× short |
| Galactic in-plane F_R | 0.222 | 0.222 | right amplitude, z-flat |
| total \|∇u\| | 0.222 | 0.223 | flat — no ordering |

The measured u_eff (0.205 → 0.116) declines while the true field
is flat-to-rising; no real-space field quantity at the pair's
location reproduces either the slope or a declining trend at the
pair's halo scale (~0.02 pc ≪ the ~2 pc well spacing).  Likewise
point-source neighbour wells fail the normalisation by ~30×
(u(1 pc, 1 M☉) = 0.006), the disk's integrated field VALUE depth
is ≪ the starvation scale, and the inter-well saddle amplitude
u_saddle ≈ u_c − A ln(D/2r*) is ≲ 0 (clamped) with the many-well
accumulation ~0.14 density-independent by construction
(n·D³ ∝ const).  The ambient that produces the ordering must
track the local density/density-column — a projection to be
derived from the constraint slice, not the local |∇φ| at the
pair.

**Density-coupled channel identified.**  Solving V_eff(u) =
V_master(u) + ρe^{−u} at physical disk densities
(ρ_mid = 6.3×10⁻¹¹⁵, ρ_halo = 2.6×10⁻¹¹⁵ M_Pl⁴): the minimum is
ambient-pinned (u_min ≈ 4×10⁻³⁷) and the effective mass is
source-curvature dominated, m_eff² = V''(u_min) + ρe^{−u_min} ≈ ρ —
i.e. m_eff ∝ ρ^{1/2}, the same exponent family the ordering
requires.  The Compton range at disk density,
λ_C ≈ 3.3 Mpc, is far too long to truncate the ~0.01-pc pair halo
directly — but the ρ^{1/2} density coupling is the correct mechanism
CLASS (environment enters through the field's in-medium response,
not the local field strength).

**[step_67 correction — the quoted densities were wrong.]**  The
in-medium propagator computation (step_67, tep_model machinery)
finds ρ_mid = 6.3×10⁻¹¹⁵ M_Pl⁴ corresponds to ~5×10⁻²⁷ g/cm³ —
about 1100× BELOW the real disk midplane (0.085 M_sun/pc³ =
5.8×10⁻²⁴ g/cm³).  At physical densities the matter curvature
contributes only ~10⁻¹⁶ of m_eff²: the propagator is
POTENTIAL-CURVATURE dominated, V'' ~ u_eq² ∝ ρ^{2/3}, giving
m_eff ∝ ρ^{1/3} exactly (numerical fit 0.333) and
λ_C = 0.08–0.13 pc ≈ 6–10× the pair-halo r* — the pair sits
INSIDE one Compton length.  The measured u_eff ∝ ρ^{0.4–0.6}
therefore sits between the k→0 propagator exponent (1/3) and the
leading pair-scale correction (m_eff/k_pair)² ∝ ρ^{2/3}: the
density law is consistent with the propagator evaluated AT the
halo scale (k = 1/r*), not in the infrared — the operative
projection is now the in-well response function
D⁻¹(k_pair; ρ), which supplies both the absolute normalisation
and the in-well X_env for the fσ8 channel.  An exact
k-branch Monte-Carlo landscape census (linear-Gauss flux
superposition over Poisson well networks, |du| =
sqrt(sqrt(2)|F|/K)) confirms every real-space channel is
amplitude-eliminated: the discreteness shear reaches only
|du| ~ 0.004-0.007 (~40x short of u_eff ~ 0.16-0.18) with
exponent n^{0.26} ~ nearest-well n^{1/3}.  

**[step_67 embedded-response pass — D^-1(k_pair) now computed.]**
The interior well machinery required already exists in the
corpus (TEP-BH steps 32-34/47 solve the deep interior well
from the master action; tep_model.solve_sphere implements the
nested phi_total = phi_env + delta_phi embedded construction),
so the remaining piece was the response function itself.  The
pair-scale source (1.4 M_sun at r* = 4.9 kAU) solved embedded
on the constraint-slice background gives

    charge_ratio(rho_bg, u_amb): rho_bg sensitivity < 0.01%
      across the full disk range (1.4-5.8e-24 g/cm^3) at fixed
      ambient; ambient dependence steep,
      charge ~ u_amb^{-2.00} over u_amb = 0.10-0.40
      (9.3e-3 at u_amb = 0 -> 1.6e-8 at 0.40).
    [Numbers quoted at the operative Cassini branch
     lambda_Cassini = 1e5 * lambda_ref; the reference branch
     (lambda_ref = 7.526e-71, excluded by the corrected
     linear-in-S_Sigma Cassini evaluation) gives slope -1.58
     and <0.4% sensitivity — retained in the output as
     embedded_response_lambda_ref for comparison.]

The environmental ordering is therefore carried ENTIRELY by
the ambient boundary field value at the vertex — local
density does not screen the pair directly, confirming the
R = S_Sigma(X_env)^2 y(s) vertex structure at response-
function level.  u_eff(rho) is the landscape ambient-map
value, and D^-1(k_pair; rho) = charge_ratio(u_amb(rho)).
Unit bookkeeping fixed: the well-solve u is charge-normalised
(u = delta-varphi/psi_uns), so the ambient enters as
phi_env = u_amb * psi_uns ~ 1e-12 varphi, not ~0.1 varphi
(which would collapse lambda_C to ~1 m).  A further
convention note: the WB machinery's ambient u0 is the
gradient/tilt variable of the step_56 construction
(phi_tot = PHI0 + u0 z + psi, xi = u0^2/2), not a field-value
offset — both live in the same solve family, applied at
different boundaries.

**Ambient-map candidate census (step_67, completed).**  With
D^-1 shown to be ambient-controlled, the remaining question
is which landscape functional sets u_amb(rho).  Evaluated
numerically across the stratum densities:

| channel | exponent | amplitude verdict |
|---|---|---|
| propagator m_eff (k->0) | 1/3 | rate ~1e11 g_t — wrong scale |
| discreteness shear MC | ~1/3 | ~40x short (020h consistent) |
| packing / quartic minimum / tail-at-Compton / halo-truncation | 1/3 | wrong exponent or amplitude |
| pair-scale correction (m_eff/k)^2 | 2/3 | bracketing bound, only ~2.4% amplitude at k_pair |
| **planar/local k-branch u_amb ~ rho^{1/2}** | **1/2** | contrast 2.44 (stratum-mean rho) / 2.01 (median) vs required 1.78-2.0 — MATCHES |

**Ambient map identified.**  The local k-branch density
response — the disk landscape's shear amplitude to its own
source, u_amb ~ (rho L)^1/2 at the disk's coherence scale —
is the only channel that sits inside the measured band AND
reaches the required amplitude: J^{1/2} at stratum-averaged
densities gives contrast 2.44 (median mapping 2.01) against
the required 1.78-2.0, and the endpoint law
u_eff ~ rho^{0.46} matches rho^{1/2} within the estimator
and e-law systematics.  The physical reading: in the
nonlinear k-branch the landscape's shear responds to local
density as sqrt(source x coherence-length) — a square-root
law, not the infrared propagator scaling; the midplane's
deeper local field suppresses the pair vertex more, in the
direction and amplitude measured.

**Coherence length derived (ambient_closure, step_67).**  The
k-branch planar balance P_X u = J with P_X ~ K u/sqrt(2) and
local source J = 4 pi G rho L/g_t gives
u_amb = sqrt(4 sqrt(2) pi G rho L/(K g_t)).  Calibrating L
once at the midplane inversion (u = 0.18 at rho_mean =
0.0852 M_sun/pc^3) yields L = 837 pc — a physical disk
coherence scale between the thin (300 pc) and thick (900 pc)
scale heights — and the single-parameter law then predicts
the remaining four strata to within 10.1%: u_pred =
0.180/0.161/0.141/0.117/0.074 vs measured
0.18/0.17/0.15/0.13/0.08.  The ambient map u_amb(rho) is
closed at coefficient level.  020d's smaller 1.56 contrast
used a two-bin density split, not a different mechanism.

**Closed (f sigma8 X_env, step_74).**  The step_62 growth channel
used xi_cosmo = H^2/(2 H0^2) from the temporal kinetic variable
alone; step_74 (`scripts/steps/step_74_fs8_inwell_convolution.py` ->
`results/step_74_fs8_inwell_convolution.json`) performs the remaining
convolution: the collapsed fraction of the mode support
(Press-Schechter, sigma_M = sigma8 (M/M8)^{-0.28}) reverts to
G_eff = 1 through the x_env -> large limit of the in-well regime
map, diluting the ambient boost.  Mass-weighted (M > 1e12 Msun) and
field-weighted (M > 1e10 Msun, lognormal-delta resolved) dilutions
agree: two-branch chi2 = 27.0 / 26.8 against the LCDM control floor
26.0, fs8(0.61) = 0.485 / 0.483 vs data 0.436 — the residual is the
shared S8 tension, not a TEP-specific excess (baseline remains at
33.5-35.8).  The scalar sector is closed at coefficient level.

**Empirical ambient map measured (step_020g, corrected).**  Free
(α, R_s) refits of the five |Z|-quintile strata give a clean
monotone ordering — α: 0.231 → 0.446 across |Z| = 22 → 347 pc.
Audit found and fixed a real bias: the first inversion marginalised
over the GLOBAL mass kernel, but stratum median mass runs
1.11 → 1.38 M_sun with |Z| (selection), so each stratum's forward
map must use its own mass distribution.  With per-stratum mass
kernels and density-weighted ρ, the α channel yields

    u_eff(|Z|) = 0.205 → 0.116,   u_eff ∝ ρ^{0.40 ± 0.05},
    contrast mid/halo = 1.78.

The R_s channel agrees at the endpoints (0.196 / 0.080) but
diverges mid-strata (z = 68 pc: u_α = 0.177 vs u_Rs = 0.241) — a
single constant ambient per stratum does not jointly reproduce
both parameters there; the stratum-internal ρ spread and the
z = 196 pc profile-shape bump remain visible in the two-channel
residuals.  Range sensitivity: refits restricted to s < 7000 AU
move the tension but do not remove it (z = 121 closes at
0.147/0.142 while z = 22 opens to 0.189/0.149) — the canonical
saturating-exponential family is being distorted by real
bin-level structure (z = 68 turnover beyond ~9 kAU, z = 196 late
bump), so part of the mid-strata disagreement is a fit-family
artifact rather than an ambient failure.  [020k later shows
that structure is itself selection-manufactured: the
R_chance < 0.01 cut truncates the high-v_tilde tail at
s ≳ 9 kAU — see the step_020k resolution below.]  The
corrected exponent
is consistent with both the ρ^{1/2} density-coupled family
(m_eff² ≈ ρ; k-branch planar u ∝ J^{1/2} over the local column)
and the ρ^{1/3} well-packing family within ~1.5σ each.

**Depth channel tested and rejected (step_020h).**  The
integrated field depth u_depth(z) = ∫|u_z| to a fixed outer edge
reproduces the measured exponent (0.52 at z_max = 600 pc) only as
a boundary artifact — the exponent sweeps 1.30 → 0.11 over
z_max = 400 → 2000 pc because depth is dominated by remaining path
length, and the apparent ρ^{0.5} match is the numerical
coincidence (z_max − z) ∝ ρ^{~0.5} over the sampled range.
Local-window depth (the pair's own screening column, W ≤ 300 pc)
has the wrong sign entirely (exponent −0.25 → −0.02: |u_z| rises
with z), and the depth's field-value amplitude is ~5×10⁻⁹ in
well-solve units — far too shallow to act as a starvation
ambient.  The whole-disk stochastic rms (~0.12, nearly flat in z)
is likewise excluded.

**Ambient-amplitude starvation tested — wrong sign (step_56
PHI0 scan).**  Enabling the dormant STARVE source channel
ρe^{−(PHI0+u0z+ψ)} and scanning the ambient amplitude PHI0 at
fixed u0 = 0.163 (proper STARVE_LAM continuation) gives
y_far = 0.28 → 0.82 over PHI0 = 0 → 0.2: deeper ambient
amplitude STRENGTHENS the pair response (starving the deep
interior relieves the self-screening shell more than the
source).  The inter-well saddle amplitude u_ISM, which is
HIGHER in denser packing (u_ISM ≈ u_c − A ln(D/2r*), A ≈ 0.20,
D ∝ n^{−1/3}), therefore drives the ordering the wrong way —
the channel is excluded.

**Internal-mixture forward model (2026-09-28).**  Rather than one
ambient per stratum, each pair was assigned a per-pair effective
ambient u_eff(z_i) ∝ ρ(z_i)^{0.4} normalised to the stratum's
α-channel inversion, and the full forward model (per-pair masses,
per-star r*, isotropic deprojection, ambient-resolved response
map) was evaluated over the stratum's actual |Z|/mass/separation
distribution.  Result: the α channel reproduces the observed
ordering essentially exactly (predicted 0.240/0.273/0.292/0.355/
0.467 vs observed 0.227/0.257/0.267/0.334/0.443 — within ~5%
across all strata), while the predicted R_s under-runs the
observed high-z values by ~25% (3835 vs 5173 AU at z = 347 pc;
3409 vs 4412 at z = 196 pc).  The two-parameter tension is
therefore concentrated in the R_s channel: the observed high-z
profiles keep rising through the 10–25 kAU bins where the
saturating-exponential family has already plateaued — a profile-
morphology mismatch, not purely an ambient mismatch.  Whether a
non-exponential response family (the true embedded-solve profile
shape, or a second environmental parameter) absorbs it is the
remaining question; the α ordering itself needs nothing beyond
the density-scaled ambient.  [020k resolution: the outer-bin
morphology is dominated by the R_chance selection function —
fits and inversions should be restricted to s ≲ 9 kAU or carry
the selection in the likelihood.]

**Flagged residual systematic — PM-noise floor in ṽ
(step_020i).**  The observable ṽ = |Δμ|·d/v_circ carries no
noise subtraction, and the per-pair PM error translates into a
ṽ noise floor that grows with distance d (hence with |Z|
stratum) and with s (v_circ ∝ s^{−1/2}): median wide-bin
(8–25 kAU) noise floor rises 0.44 → 0.63 in v_circ units from
z = 22 → 347 pc.  Since |Δμ| is a magnitude, Gaussian PM noise
rectifies upward (Rice bias) and the pipeline's bin medians are
inflated — the step_012 Newtonian null is noiseless, so the
published comparison baseline underestimates the expected
profile at wide bins.

step_020i adds the measured per-pair errors to the step_012
Newtonian sampler and reports (a) the noise-rectified null and
(b) an unbiased low-noise subsample test (cut on the absolute
noise floor σ_ṽ < 0.10/0.15/0.20, which does not condition on
the observed ṽ).  Results:

* The noise floor genuinely inflates the raw profile at the
  widest bins — a noise-aware thermal null reaches outer median
  ~1.27 vs the noiseless ~1.0 baseline — so raw α and R_s at
  s ≳ 8 kAU carry a real systematic inflation that must be
  reported alongside the signal.
* However, the low-noise subsample shows the boost and the
  environmental ordering survive where the floor is provably
  negligible: at σ_ṽ < 0.15 the z = 22 pc stratum sits at
  ~1.04–1.11 in the outer bins while z = 347 pc sits at
  ~1.16–1.33.  The ordering is stronger at high |Z| — the
  direction TEP predicts — and the ~1.10–1.30 residual boost is
  consistent with the published WB literature.
* Net: the noise floor inflates the *raw* α/R_s at the widest
  bins (a real systematic, now quantified and reproducible in
  020i), but it is NOT the origin of the signal — the detection
  and ordering persist in the clean subsample.
* Joint-likelihood resolution (implemented): a boosted-Newtonian
  model (α, R_s) forward-modeled through the SAME per-pair noise
  kernel is compared against the noise-aware pure-Newtonian null
  per stratum, per eccentricity law, on common random numbers.
  The boosted model wins decisively under every e-law and every
  stratum except a single (z4, superthermal) tie:
    Δχ² ranges 34–2249 (uniform), 179–1922 (circular),
    189–1054 (thermal), 0–1217 (superthermal; the z5 stratum
    rejects the null hardest at Δχ² = 1217 with α_fit = 0.40).
  Fitted residual boosts α_fit ≈ 0.05–0.40 with R_s ≈ 500–6600 AU
  — the corrected amplitude is smaller than the raw profile fit
  (α ≈ 0.44 at high z), as expected once the Rice floor is in
  the null, but the detection is unambiguous: the noise-aware
  Newtonian null is rejected at Δχ² ≳ 34–2249 in 19 of 20
  stratum×e-law cells.
* Remaining caveat: the σ_ṽ cut selects nearer/tighter/
  higher-mass pairs within each stratum, and the joint grid is
  coarse (17×20); a refined boost-grid or free eccentricity
  marginalisation is the follow-up.  The detection, the
  environmental ordering, and a residual boost consistent with
  the published literature all survive the noise-aware null.
* Noise-corrected ambient inversion (step_020j): the two-branch
  y(d/r*, u0) forward model pushed through the per-pair noise
  kernel gives u_eff = 0.60 (z~22 pc), 0.27 (68), ≥1.5 (121, 195;
  grid edge), 0.18 (346).  The endpoint ordering is the right
  direction (deep screening at the midplane, least screening in
  the halo) but the mid-strata are non-monotone and chi²_min is
  large (103–370 on ~12 dof): the single-u0-per-stratum model
  is inadequate for mid-strata profile shapes even after noise
  correction — the same tension the raw-profile inversion showed,
  now confirmed not to be the noise floor.
* Internal-mixing test (same 020j machinery): assigning each
  pair its own ambient u0(z) = 0.25 (rho(z)/rho_med)^p with the
  disk column rho ~ e^{-|z|/260 pc} and scanning p = 0–1 leaves
  the z~121 pc stratum at chi² ~ 1160 for every p — internal
  density mixing does NOT produce the mid-strata shape either.
* Morphology diagnosis — RESOLVED (step_020k): the observed
  per-stratum profiles are bump-shaped because the purity cut
  manufactures the downturn.  Rebuilding the pool from the raw
  catalog with all step_001 cuts except the R_chance_align
  threshold shows the per-stratum profiles are IDENTICAL at
  every threshold for s ≲ 9 kAU (the cut removes <1.3% of
  pairs there) and diverge sharply above it: under the
  baseline R_chance < 0.01 the profiles turn over (z1 peaks
  1.23 @ 9k AU -> 1.19 @ 25k; z5 peaks 1.27 @ 18k -> 1.22),
  while under R_chance < 0.05–1.0 they keep rising to
  1.4–1.8.  At s = 8–30 kAU the baseline cut removes 20% of
  the otherwise-qualified pool, preferentially the
  high-v_tilde tail (the residual R_chance ≳ 1% systems are
  unbound interlopers with median v_tilde ~ 2.4x the kept
  population), so the 1% threshold truncates the real outer
  tail and depresses the outer-bin medians.  Neither extreme
  is the true bound-binary profile: the 1% cut truncates, the
  relaxed cut admits unbound contaminants.  Consequences:
  (i) the mid-strata inversion tension (u0 grid-edge at
  z3/z4, chi² ~ 1000) is explained — the canonical family was
  being forced to fit a selection-manufactured downturn, not
  a physical turnover; (ii) the reliable region for fitting
  and environmental tests is s ≲ 9 kAU, where the profile is
  selection-insensitive, or equivalently the selection cut
  must be forward-modelled in the likelihood for wider bins;
  (iii) in the cut-insensitive region the environmental
  ordering persists — median v_tilde over s = 2–10 kAU rises
  0.718 -> 0.728 -> 0.741 -> 0.764 from z2 to z5 (+6.4%),
  and the peak-amplitude ordering z5 > z3 > z1 remains the
  TEP direction on the unaffected profile; (iv) canonical
  (alpha, R_s) refits restricted to s < 9 kAU show the low-z
  strata saturating at alpha = 0.25–0.28 while the high-z
  profiles are still rising at 9 kAU and drive the fit to
  the alpha = 1, R_s ~ 50 kAU bounds — i.e. the dilute-
  environment response may be less saturated than the
  headline fit implied, consistent with a stronger and
  later-rising high-z signal.  The ordering direction is
  robust; a clean corrected exponent requires the selection
  function in the likelihood (or an honest outer-bin window).
* Selection-safe coupling quantification (step_020l,
  `results/outputs/020l_selection_safe_environment.json`):
  restricting to the cut-insensitive window s < 9 kAU, the
  per-pair slope of v_tilde on |Z| with separation x distance
  fixed effects and a mass covariate is +0.0046/100 pc at
  +6.65 sigma against a within-cell permutation null
  (R_chance < 0.01, N=305,778), rising to +11.7 sigma
  (R<0.005) and +17.2 sigma (R<0.002); against the disk-column
  proxy log rho the slope is -0.0277 (same z).  Robustness
  anatomy: (i) the OLS slope collapses to +0.3 sigma when the
  ~1,964 boundary pairs in R_chance in [0.01,0.05) are
  re-admitted -- but those pairs are independently identified
  as unbound (median v_tilde 2.07 vs 0.685, median |Z| 86 pc,
  i.e. plane-concentrated chance alignments); (ii) the
  MEDIAN-based z-quintile ordering is threshold-independent
  -- quintiles 0.674/0.668/0.680/0.696/0.727 with NO purity
  cut at all, identical to baseline -- so the coupling lives
  in the bulk bound population; (iii) a v_tilde<1 boundness
  proxy with no R_chance cut returns a null slope, but that
  selection truncates the signal channel (boosted pairs sit
  at v_tilde ~ 0.8-1.5), so it bounds fragility rather than
  refutes the coupling.  The per-cell ordering matrix shows
  the z5-z1 contrast monotone in z in EVERY separation bin
  and growing with s (+0.018 to +0.084), reaching
  +0.043 +/- 0.0063 (6.9 sigma, 100% bootstrap-positive) at
  s > 3 kAU.  Consequence for the manuscript: the
  environmental ORDERING (dilute environment -> stronger,
  later-rising response) is robust and now demonstrated
  selection-independently at ~7 sigma; the fitted R_s VALUES
  and the p = 0.50 exponent derived from full-range profiles
  remain contaminated by the outer-bin truncation and require
  either the restricted-window refits or forward-modelled
  selection.

* Window-restricted noise-aware inversion — RESOLVED the
  mid-strata tension (step_020l).  Two residual artifacts
  beyond the outer-bin truncation were isolated and removed:
  (a) the outer bins s > 9 kAU carry the R_chance
  selection function (020k), so chi2 is restricted to the
  selection-insensitive window; (b) the pipeline's
  per-stratum self-normalisation anchors each stratum on a
  different separation range — the high-|Z| anchor bins are
  themselves noise-inflated — which manufactured the
  apparent mid-strata "Newtonian preference".  On the
  ABSOLUTE median-v_tilde comparison (mass-calibration bias
  cancels because the forward model uses the same
  mass_total column), through the same per-pair noise
  kernel, the inversion is clean and monotone in the TEP
  direction under every realistic eccentricity law:

    thermal:      u_eff = 0.18, 0.17, 0.15, 0.13, <0.08
                  (dchi2 vs Newtonian = 1589–2787; slope
                  vs ln rho_disk = 0.50, with z5 a lower-edge
                  bound — consistent with rho^{1/2})
    superthermal: u_eff = 0.11, 0.12, 0.09, <0.08, <0.08
                  (dchi2 = 2200–3340; slope ~0.55)
    uniform:      u_eff = 0.45, 0.45, 0.35, 0.27, 0.15
                  (dchi2 = 571–1794; slope ~0.93)
    circular:     Newtonian-favoured at z1–z4 (the same
                  eccentricity-prior caveat the published WB
                  anomaly carries; circular orbits
                  under-predict the observed velocity
                  distribution outright)

  The noise-aware Newtonian null is rejected at
  dchi2 ~ 1600–2800 per stratum under the reference thermal
  law, and the corrected environmental exponent comes out
  ~0.5 — now identified with the derived k-branch ambient law
  u_amb ~ rho^{1/2} (step_67 ambient_closure; the earlier
  m_eff^2 ~ rho source-curvature claim was corrected: the
  propagator regime is potential-dominated, m_eff ~ rho^{1/3}),
  from a defensible estimator (selection-clean window + noise
  kernel + absolute-profile comparison).
  The z5 stratum reaches the grid floor with chi2 flat for
  u0 <~ 0.08 — a genuine bound: the dilute environment is
  effectively unscreened.  Corpus consistency: the midplane
  stratum (|Z| ~ 22 pc, the solar-neighbourhood column)
  inverts to u_eff = 0.18, within ~10% of the ambient
  u0 = 0.163 derived independently at the solar circle from
  g_obs = g_t·u_bar·(P_X+2) — the WB environment chain is
  anchored to the same scalar configuration, not a free
  parameter.  Residual structure at best fit: the observed
  profile sits +2–10 sigma above the model at 0.3–2.4 kAU
  and slightly below at the window edge — a steeper mid-
  transition than the canonical y-family shape; the u_eff
  ordering and endpoints are stable across eccentricity
  laws, so this is a response-shape refinement (candidate:
  orientation-weighting or family shape), not an ordering
  issue.

* Covariate-hardened environmental coupling (step_020l
  selection-safe version).  The pool-level |Z| slope is
  distance-confounded — high-|Z| pairs are more distant and
  carry a larger PM-noise floor, and adding the per-pair
  noise floor or fine distance cells absorbs the naive
  bulk-median ordering (the raw 0.669 -> 0.726 quintile
  medians are mostly this effect).  The decisive test fixes
  separation, distance decile, noise-floor decile and mass
  decile in fine cells and asks whether residual v_tilde
  still varies with |Z|.  The coupling is SEPARATION-
  LOCALISED on the metallicity-corrected sample: z-slopes
  per 100 pc of +0.002 (s=0.5–1.5 kAU) rising monotonically
  to +0.025 at s=5–9 kAU (t = +13.0, +12.2 with the feh
  covariate).  A uniform distance/noise artifact cannot
  switch on at a separation scale — and the observed
  turn-on tracks the model's own prediction: the z5−z1
  model contrast grows 0.001 -> 0.073 across the window
  while the observed contrast grows 0.022 -> 0.110, the
  s-dependent part matching within ~40%.  Direct bootstrap
  contrast at s>3 kAU: median z5−z1 = +0.0433 ± 0.0063
  (6.9 sigma, positive in 100% of resamples).  Slice-
  consistency audit inside the 5–9 kAU window: the |Z|
  slope is positive in 10/10 distance deciles, 10/10
  noise-floor deciles and 10/10 mass deciles (30/30
  slices, median slopes +0.019 to +0.025) — the coupling
  is not carried by any confounded subset.  Numerical
  hardening: OLS replaced by QR/lstsq + pinv with HC1
  heteroscedasticity-consistent SEs (t = +12.7, +11.9
  with feh; the earlier warnings were stale LAPACK FP
  flags, not data).  Remaining degenerate alternatives,
  stated plainly: thick-disk orbital-population
  differences (eccentricity priors) and the Galactic-tide
  survival argument both predict directionally similar
  effects and are not fully excluded by the covariate
  battery; the feh covariate reduces the slope only
  +12.7 -> +11.9 sigma.

* Free eccentricity-mixture fit (step_020m).  Letting the
  eccentricity distribution float as a mixture over
  thermal/superthermal/uniform components on the safe window,
  every stratum selects the uniform-dominated component and the
  absolute profiles fit at chi2/dof ~ 1-8 — the thermal-law
  mid-window residual structure was eccentricity-prior shape,
  not response-family failure.  The u_eff ordering stays
  monotone (0.40 -> 0.15); under the data-preferred uniform
  mixture the density exponent is ~0.7, so the honest exponent
  band across e-law systematics is 0.5-0.7, centred on
  rho^{1/2} but not precision-pinned.  The noise-aware
  Newtonian null is rejected at chi2 ~ 800-2000 per stratum
  under the same mixture.  The fitted e-law
  weights are |Z|-independent (uniform ~1.0 in every stratum):
  a thick-disk eccentricity-population explanation of the
  ordering would have shown up as stratum-dependent weights —
  partially closing that degenerate alternative.

* Jacobi-radius bound on the tide alternative (step_020o, real
  per-pair numbers).  With Mamajek masses from the pipeline pool
  (median M_bin = 1.34 M_sun) and the solar-circle tidal tensor
  (M_enc = 1e11 M_sun, R_GC = 8.2 kpc): r_J percentiles 10/50/90 =
  221/279/323 kAU; response-window pairs (s = 5–9 kAU) sit at
  median s/r_J = 0.023 with 0.0% beyond 0.1; the widest sample
  bin (9–30 kAU) contributes 3.5% beyond 0.1.  Beyond ~30 kAU the
  tide does become the dominant control (median s/r_J = 0.13,
  92% beyond 0.1) — bounding the reach of pair-intrinsic
  inference rather than challenging the response window.
  Differential tidal truncation cannot manufacture an effect 40×
  inside the tidal radius; the remaining tide channel
  (formation-era survival demographics) does not act on the
  measured v_tilde ordering of surviving bound pairs.  Of the two
  alternatives, the eccentricity-population explanation is now
  the only one that remains materially open, and it is
  constrained by the |Z|-independent mixture weights.

* Purity-weighted profiles (step_020n).  The hard R_chance cut
  is replaced by the selection-consistent estimator: every pair
  carries bound-probability weight w = 1 - R_chance_align, so
  interlopers are down-weighted smoothly with no
  separation-dependent truncation anywhere.  On the full pool
  (N = 358,359; <w> = 0.995, N_eff = 357,863) the weighted
  profiles rise MONOTONICALLY to 1.6-1.8 at 30 kAU in every
  stratum — the downturn is removed entirely (weighted outer
  decline = 0.0 vs hard-cut 0.04-0.08), confirming the 020k
  diagnosis with a continuous estimator rather than threshold
  variants.  Absolute (un-normalised) weighted medians show the
  environmental ordering growing with separation: the z5-z1
  bootstrap contrast is +0.018 +/- 0.006 (3.3 sigma) at
  s = 2-5 kAU, +0.042 +/- 0.009 (4.7 sigma) at 5-9 kAU,
  +0.066 +/- 0.014 (4.8 sigma) at 9-15 kAU, and
  +0.094 +/- 0.018 (5.3 sigma) at 15-30 kAU — the
  separation-growing coupling now measured without any hard
  purity cut.  A smaller inner-window contrast exists as well
  (+0.040 +/- 0.005, 8.4 sigma at s = 0.5-2 kAU): a
  baseline component consistent with the stellar-population
  residual already identified in the 020l tight-bin slopes,
  over which the separation-growing response sits.  Canonical saturating-exponential fits on the
  weighted profiles push R_s to the 60 kAU bound in the
  dilute strata (z3-z5; the restricted <9 kAU fits do the same):
  the response is unsaturated at the widest measured
  separations, so hard-cut transition radii are lower-bounded
  effective values and the high-|Z| response is stronger and
  later-rising than any truncated-profile fit indicated.

* Bound-component deconvolution (step_020n, pessimistic
  interloper sweep).  CDF-level mixture subtraction on the
  20-50 kAU bin, with the interloper shape from clearly unbound
  pairs (R_chance > 0.2) rescaled by sqrt(s).  Because
  R_chance is a prior that plausibly underestimates the true
  interloper fraction at wide s — the v_tilde > 1.5 tail is
  0.34-0.38 there vs 0.06-0.09 at small s — the bound
  component is evaluated under f_unbound = 0.11 (the R_chance
  estimate), 0.20 and 0.30.  The bound-component median stays
  0.86-1.11 across all five strata vs the s < 2 kAU baseline
  0.65-0.69 — elevated ~25-70% even under a 3x-prior
  interloper fraction.  The wide-separation elevation is a
  bulk property of the retained distribution, not a few
  residual interlopers pulling a tail.

**Manuscript-number audit (p = 0.50 ± 0.03 claim).**  The
published exponent traces to step_009's "direct" fit, which is a
one-parameter ANCHORED estimator (amplitude fixed at bin 1), not a
free fit: p ≈ ln(R_s,1/R_s,5)/ln(ρ_1/ρ_5) = ln 2.00/ln 4.07 = 0.49
by construction.  Under a free two-parameter fit the same five
fixed-α radii give p = 0.36–0.39; the quoted ±0.03 is the raw
Δχ² = 1 error on a fit with χ² = 10.7 on 3 dof (rescaling gives
±0.06).  The radii themselves carry the fixed-α = 0.4 convention
compression.  Removing the predicted mass covariate (heavier halo
pairs would *raise* the halo R_s, opposite the observed dip)
steepens the exponent to ~0.47.  Net: the honest data-side
exponent is ~0.4 ± 0.1 under free fits, consistent with both
ρ^{1/3} and ρ^{1/2}; the manuscript's precision claim should be
revisited at the next revision (either free-fit 0.39 ± 0.13 from
the same points, or the corrected ambient-inversion 0.40 ± 0.05).

**Surviving statement.**  The corrected measured law
u_eff ∝ ρ^{0.40 ± 0.05} stands as the data-side target, and
the step_020l window-restricted noise-aware inversion now
returns an independent estimate centred on ρ^{1/2}
(0.50 non-edge slope under thermal, 0.55 superthermal),
lifting the target range to ρ^{0.4–0.6}.  Every literal real-space ambient
candidate is now tested and
rejected (gradient channels wrong sign; rms flat; depth
boundary-artifact; starvation amplitude wrong sign; packing
amplitude ~40× short), so u_eff is an effective parameter of
the in-medium response — not a literal boundary field.
Step_67 corrects the earlier propagator estimate: at physical
disk densities m_eff ∝ ρ^{1/3} (potential-curvature dominated,
NOT m_eff² ≈ ρ — that evaluation ran ~1100× below the real
midplane density), and λ_C ~ 0.1 pc encloses the pair halo.
The measured exponent then sits between the infrared
propagator (1/3) and the pair-scale correction (2/3): the
operative projection is the in-well response function
D⁻¹(k = 1/r*; ρ) on the constraint slice — the AUD-3
remainder; no phenomenological exponent is inserted.

**Mass-stratified transition (step_005 vs corrected estimator).**
The same forward model split into mass terciles predicts
R_s = (2456, 3322, 4173) AU at M_med = (0.74, 1.20, 1.70) M☉ —
a predicted mass exponent n = 0.64 — versus the observed
(2778, 2966, 4847) AU and n = 0.635 (bootstrap [0.40, 4.88]).  The
pipeline-shaped scaling (per-star r* ∝ √M modulated by per-bin mass
marginalisation and deprojection) reproduces the measured exponent
to ~0.5% — a third independent signature with no free parameters.

## 9. T1.2 result — PPN under the two-branch sector (Rule 17)

`step_57_ppn_two_branch.py` → `results/step_57_ppn_two_branch.json`.

Isolated-source response `y(x)`, x = r/r_*, solves the flux law

  baseline:    y (1 + y²/x⁴) = 1
  two-branch:  k y²/(√2 x²) + y³/x⁴ = 1

Both have the screened interior asymptote y ~ x^{4/3} (verified:
numerical exponents 1.3333 and 1.3211).  The two branches inside the
Sun's own field cross over at x_c = (√2/k³)^{1/2} ≈ 0.0186 r_*
≈ 78 AU — far outside every bound-constraining probe:

| probe | x / r_* | y (two-branch) | γ − 1 |
|---|---|---|---|
| Cassini conj. (1.6 R_☉) | 1.8e-6 | 2.16e-8 | −8.63e-8 |
| Mercury orbit | 9.3e-5 | 4.16e-6 | −1.66e-5 |
| Saturn orbit | 2.3e-3 | 2.8e-4 | −1.12e-3 (diagnostic) |

Cassini bound |γ−1| < 2.3e-5 → required S_eff < 5.8e-6; the Sun's
exterior response at conjunction is 2.16e-8 → **margin ~270×, PASS**.
Two-branch shift at conjunction: −0.06%.  The divergence between
sectors begins at x ~ O(1) — the wide-binary regime — where the
isolated profile is ambient-clamped and step_56's embedded solve is
the correct estimator.  Sign-definiteness preserved: γ < 1 always
(small-𝒮 expansion of a DEF γ formula with α_eff = 𝒮_Σ α_0).

### Follow-up: interacting-fields (two-centre) solve — step_59

`step_59_two_center.py` → `results/step_59_two_center.json`.
Two Gaussian sources (Q = 4π each) at separation d on-axis, ambient
gradient parallel (only axisymmetric configuration); mutual force via
momentum-flux surface integral, B-only solution subtracted to isolate
the interaction; linear-sector reference F_lin ∝ 1/d² verified to ~4%.

**The mechanism is real but directionally wrong.**  Each star's well
raises the companion's local ambient, so the mutual interior is
suppressed MORE than the monopole solve (y_pair(1 r*) = 0.15 vs
monopole ~0.27 at u0 = 0.15, two-branch) and the rise to the plateau
is a slow algebraic approach.  Through the estimator:

| sector | u0 | α | R_s (AU) |
|---|---|---|---|
| two-branch | 0.15 | 0.503 | 15,331 (20.8σ) |
| two-branch | 0.57 | 0.152 | 6,085 (5.6σ) |
| baseline | 0.57 | 0.681 | 8,945 (10.3σ) |

All three candidate escapes are now tested and closed:
(i) heterogeneous ambient — convex combination interpolates along the
    failing trade-off curve, never reaches (0.366, 2646);
(ii) exp_interp rival — baseline-identical at Galactic ambient,
    surrenders the SPARC gain;
(iii) interacting two-centre fields — deepens interior suppression,
    inflates R_s.

**Extension — object-level physics (starved source + K(phi)).**
Following the objects-as-wells correction, the two-centre solve was
extended with (a) starved sources ρ → ρ e^{−φ_tot} (continuation-Newton;
partial convergence only — the e^{−λφ} source shell is sub-cell at the
fiducial grid, residual floor ~1e2) and (b) a φ-dependent kinetic
stiffness J_i = K(φ)·f(|a|)·a_i with K = e^{cψ} (Picard-frozen
coefficients, converged dmax ~ 7e-7).  Forces extracted by the direct
matter integral F = |∫ρ_B w ∂_zψ dV| (cross-validated against the
stress integral to ~1–3% on the unstarved sector).

Findings:

- Starvation (φ0 = 1): uniform amplitude suppression, no shape rescue —
  at u0 = 0.57 the plateau falls further (y(5r*) ≈ 0.06 vs 0.15
  unstarved); at u0 = 0.15 the profile shape is unchanged.
- Sign convention check (audited 2025-09-28): the pair solve uses
  ψ < 0 inside wells (div J = ρ Poisson sign), so solver `c > 0` in
  K = e^{cψ} corresponds to the PHYSICAL coupling K(φ) = e^{−cφ}
  with φ > 0 in wells — i.e. K < 1, a SOFTENED bridge.
- K(ψ) pair channel, CORRECTED bookkeeping (2025-09-28): the direct
  matter integral requires the B-only subtraction when K ≠ 1 —
  B's own modulated well refracts ambient flux and exerts a net
  force on B (dielectric self-term), which dominated the earlier
  unsubtracted numbers (spurious y > 1 "focussing" and the apparent
  large-d turnover).  With the self-term removed the mutual
  response is a clean monotone function of d and the enhancement
  is modest:

  physical e^{−cφ} (soft bridge, solver c > 0):
    u0 = 0.15, c = 1.0:  y = 0.29 → 0.71  (vs plain 0.15 → 0.44)
    u0 = 0.30, c = 1.0:  y = 0.21 → 0.33  (vs plain 0.13 → 0.26)
    u0 = 0.57, c = 1.0:  y = 0.14 → 0.11  (vs plain 0.10 → 0.14)

  physical e^{+cφ} (stiff wells, the SPARC-required sign, solver
  c < 0):
    u0 = 0.15, c_phys = 0.5:  y = 0.12 → 0.38  (suppressed ~15%)
    u0 = 0.57, c_phys = 0.5:  y = 0.09 → 0.14  (≈ unchanged)

  i.e. at the solar-circle ambient u0 = 0.57 the modulation is
  ineffective in EITHER direction; where it does act (low u0) it
  changes amplitude, not shape — the profile still rises on the
  wrong (broad) scale.
- **Sign tension between the two channels (the operative
  exclusion).**  SPARC requires the OPPOSITE sign: with physical
  K(φ) = e^{+cφ} (K > 1 in wells) the RAR improves monotonically,
  rms 0.219 (c = 0) → 0.199 (c = 3, matching the empirical ν-function
  reference 0.199), median offset −0.086 → −0.020; with K < 1 in
  wells (c < 0) SPARC degrades (0.239 at c = −1, 1.42 at c = −3).
  The WB pair instead improves only under K < 1 in the merged well.
  Since galactic interiors (ψ_med ≈ 0.35) and the pair bridge
  (|ψ| ~ 0.3–5) sample the same field-depth range, no monotone
  K(φ) — and a fortiori no exponential — can satisfy both.
- solver bugs found and fixed during this audit: (i) `solve_linear`
  in step_59 patched f → 1 but left KPHI_C/KPHI_PERT/STARVE active,
  silently modulating any reference computed after a coupled run —
  fixed by save/restore of all coupling flags (step_59's main()
  computes references before enabling couplings and is unaffected);
  (ii) `pair_force_direct` skipped the B-only subtraction on the
  argument that the self-field is symmetric — true only for K = 1;
  under K(ψ) ≠ 1 the ambient refracted through B's own modulated
  well produces a net dielectric self-force, now subtracted;
  (iii) W_of_q's flux-potential table was capped at q ≤ 1 while
  the ambient field can exceed that — table now extends to the
  realized q range.
- Corrected status of the earlier claim: the "~9-decade field range"
  structural exclusion was an artifact of the mass-unit bug in what is
  now TEP-UCD step_11 (former Paper 0 step_60)
  (M = gR²/G missing a factor of R, inflating interior ψ to ~10⁹;
  corrected ψ_med ≈ 0.35).  The φ-level coupling is now excluded on
  the stronger, correct ground of the sign tension above.

**Closure of the remaining live options.**

(a) ṽ mapping — CLOSED analytically.  The observable
`ṽ = √(1 + 2β²·q²·⟨y⟩)` is sector-general: the factor 2β² derives
from the fixed conformal coupling A(φ) = e^{−φ} (Rule 3), not from
P(X); the only sector-dependent inputs are the vertex q² and the
propagator y(s), and step_019 already sweeps all three vertex
readings (pure-propagator, legacy 0.433, consistent P_X(u0)−²).  The
inner-bin normalization cancels ⟨y⟩_inner which is small in every
sector (~x^{4/3}).  No mapping freedom remains.

(b) Perpendicular orientation — CLOSED by proxy.  step_56 now stores
the directional decomposition: the ambient wake is AXIS-focused, so
the transverse monopole channel (linearized proxy for the ⊥ pair
mutual force, stars on the s-axis at z = 0) is MORE suppressed than
the axial one.  Two-branch at u0 = 0.57: y_far_par ≈ 0.21,
y_far_perp ≈ 0.10; at u0 = 0.15: 0.78 vs 0.41.  Orientation
averaging lowers the plateau relative to the ∥-only estimate — the
∥ solve was the optimistic bound.  CONFIRMED by the full 3D solve:
step_75 (Cartesian box, parameters matched to step_59) finds
y_pair(θ) = 0.270, 0.258, 0.242, 0.219, 0.158 at θ = 0, 30, 45, 60,
90 deg for the two-branch sector (baseline 0.460 → 0.376): the
perpendicular pair is ~1.7× MORE suppressed than parallel, as the
proxy predicted.  θ = 0 lands at 0.270 vs the axisymmetric 0.152 —
the same-side result, offset by the coarse 3D source representation
(ratio stable N=24 vs N=64; absolute drift ~20%).  Orientation is a
directional, bounded ~2× correction downward for orientation-averaged
populations — the ∥ plateau remains the optimistic bound.

**Structural exclusions — revised (final form).**  φ-level couplings
of monotone exponential form are excluded by the SIGN TENSION: the
RAR requires K(φ) > 1 in wells (stiffening suppresses the predicted
response toward the data), while the pair channel requires K(φ) < 1
in the merged well (softening to reach the observed amplitude at the
solar-circle ambient).  Both environments sample φ ~ 0.3, so the
field value alone cannot distinguish them.  The earlier 9-decade
exclusion argument is superseded (it rested on a unit bug).  What
survives as the residual live option is a non-monotone or
multi-argument K(X, φ) — e.g. keyed on the AMBIENT gradient rather
than well depth — but no concrete form has been identified, and the
burden is on any proposal to show it does not simply reintroduce the
two-channel conflict.

**Verdict — SUPERSEDED (step_020).**  The structural exclusion of
φ-level couplings above stands as physics (no monotone K(φ) can serve
both channels), but the premise that the WB channel needed rescuing
was wrong: the corrected estimator (scalar-shear ambient u0 = 0.163,
resolved propagator, deprojection) closes the joint (α, R_s) test at
0.44σ / 1.69σ for the two-branch sector while excluding baseline and
exp_interp at ~30σ.  The pair-channel investigations above remain
valid as mechanism surveys — none is needed, and the sign tension
now records why the K(φ) direction is unavailable as a sector
modification rather than as a rescue requirement.

## 10. T-B1 result — well-family solve settles Rule 20

`step_58_well_family.py` → `results/step_58_well_family.json`.

Static spherical well family: fixed total charge Q = 4π, uniform
sphere radius R swept over 8 decades toward the compact limit.
Field equation ∇·[f(|∇u|)∇u] = ρ_shape e^{−u} − C_V V'(u), solved by
lagged-flux Picard iteration (each step exact: f·q = F(r)/r² inverted
algebraically per the step_55 law).

**Result: asymptotic approach — the .devin Rule 20 reading is
dynamically correct.**

- u_c(R) grows unbounded but only ∝ ln ln(1/R): the lnln fit beats the
  ln fit in both sectors (baseline residual 0.024 vs 0.142; two-branch
  0.038 vs 0.051).  Analytic self-consistency: the exterior tail sees
  only the starved charge e^{−u_c}, giving u_c·e^{u_c/2} ~ C ln(1/R).
- N_min = e^{−u_c} → 0 as a power of 1/ln(1/R): an approach, never a
  completed halt.  The file version's "finite N_min > 0" is NOT
  supported by the master-sector dynamics.
- The V₀ floor does NOT cap the field: master-potential restoring
  changes u_c by <0.3% at the deepest compaction (both branches of
  V' die at large u — floor tail ∝ u^{−5}, quartic-envelope
  super-exponentially — while the starved source e^{−u} decays but
  the 1/r² flux-focusing keeps integrating).  What V₀ supplies is
  the *regular interior geometry* (V → const = de-Sitter-like
  scale-invariant core in g) — that part of the file's Rule 20
  survives, as a statement about geometry, not about a clock floor.
- Physical bound: growth is so slow that at Planckian compaction
  (ln(1/R) ~ 100) u_c ~ O(5–10), N ~ e^{−5..10} — nonzero for every
  realizable source; the halt is only formal at R = 0.
- Two-branch vs baseline: same asymptote family (lnln), two-branch
  systematically shallower (u_c ≈ 1.41 vs 1.57 at R = 10⁻³) — the
  stiff small-X branch starves faster.

**Rule-file resolution:** sync `TEP/tep-rules.md` Rule 20 to the
.devin asymptotic wording for the clock rate; retain the V₀
scale-invariant-core language re-scoped to the g-sector interior.

## 11. T1.4 result — stability / hyperbolicity certificate (PASS)

`step_61_stability_certificate.py` → `results/step_61_stability_certificate.json`.

Scanned ξ ∈ [10⁻¹², 10⁶] on the static branch and
ξ_cosmo = H²/(2H₀²) on the timelike branch for all three candidate
sectors.  Result:

- **P_X > 0, Z_∥ = P_X + 2ξP_XX > 0, Z_⊥ = P_X > 0** for all three
  sectors on both branches — strictly elliptic, no ghosts
  (Z_t > 0 at every sampled z up to z = 100).
- Sound-speed index c_s² = P_X/Z_∥ ∈ [1/3, 1/2] for two-branch —
  **subluminal everywhere**; the ξ → 0 vacuum limit is the marginal
  ellipticity point (P_X → 0) but the fluctuation equation stays
  hyperbolic with c_s² → 1/2.  Baseline's vacuum limit is c_s² → 1.
- Two-branch is the stiffest sector on the cosmological branch
  (P_X ≈ 12.3 at z = 0 vs baseline 2.0) — consistent with the
  relaxed G_eff found in step_62.

## 12. AUD-3 first pass — growth under G_eff(z) (step_62)

`step_62_growth_fsigma8.py` → `results/step_62_growth_fsigma8.json`.

Quasi-static growth on the ΛCDM conformal image (AUD-4 branch a) with
G_eff(z) = 1 + 2β²/P_X(H²/2H₀²).  Reference data: BOSS DR12/eBOSS
fσ₈ (5 pts) + DES Y3/KiDS S₈ (2 pts), ΛCDM-normalized σ₈ = 0.834.

| sector | G_eff(0) | fσ₈(0.61) | γ | χ² (7 pts) |
|---|---|---|---|---|
| baseline | 2.00 | 0.532 | 0.282 | 44.2 |
| **two-branch** | **1.16** | **0.489** | **0.505** | **27.7** |
| exp_interp | 2.00 | 0.532 | 0.282 | 44.2 |
| data | — | 0.436 ± 0.034 | ~0.55 | — |

The same small-X branch that improves SPARC halves the growth-channel
discrepancy: baseline's G_eff ≈ 2 badly overproduces late structure;
two-branch remains ~0.05 above the fσ₈ points — the residual excess is
the plan's flagged debug signal (X_env on linear scales is set by the
landscape constraint slice, AUD-6/T6.3, not by H(z)).  Recorded as
directional improvement, not closure — the full growth equation on
static g is still owed.

**Ambient-projection inversion (R3.2 check).**  Fitting
u_env(z) = u0(1+z)^p through the same growth ODE and reference data:
at u0 = 1 the data strongly prefer p > 0 (χ² = 14.4 at p = 1.0 vs
19.3 at p = −0.5) — i.e. ambient GROWING with z, more screening in
the past.  This is the direction the drift-linked ambient
(u_env ~ E(z)) supplies and the OPPOSITE of a naive
structure-deepening ambient (u_env falling with z).  The E-image
projection is therefore not obviously mis-projected — the residual
χ² = 27.7 sits at the known S8-tension floor (control 26.0), and the
constraint-slice derivation remains the closure task, now with a
measured target: u_env(z) rising into the past, roughly E-like.

**Post-junction clarification (AUD-4 completion).**  Under the
piecewise A(η) (step_64/65), the fσ₈ data range z < 1 lies entirely
on the ΛCDM-image branch, where the ambient drift rate IS the image
H(z) — so the ξ = E(z)²/2 projection is the *correct* ambient-drift
rate within the measured range, not merely a provisional proxy.
The residual fσ₈ excess must therefore come from the well-interior
projection (perturbations feel the local field inside their host
wells, not the ambient) — i.e. the X_env that enters the growth
equation is the in-well kinetic variable, which the constraint slice
(T6.3) must supply; the ambient-map part of the projection is
settled.

## 13. AUD-2 / R2.1 — zero-parameter onset-radius law (TEP-UCD step_12)

`TEP-UCD/scripts/steps/step_12_onset_radius_law.py` →
`TEP-UCD/results/outputs/step_12_onset_radius_law.json`
(migrated from Paper 0 step_63).

54 SPARC galaxies bracket the crossing g_bar = a_eff
(a_eff = 4√2 g_t/k = 1.20×10⁻¹⁰ m/s²); M_bar from the outermost
baryonic rotation component.

- Fixed slope 1/2: fitted intercept −4.445 vs zero-parameter
  prediction −4.467 → **normalization reproduced to 0.022 dex**,
  rms 0.276 dex.  (Audit fix 2025-09-28: the crossing interpolation
  fed np.interp a descending xp, pinning r_onset to the bin's outer
  edge; corrected values supersede the earlier 0.08-dex estimate.)
- Zero-parameter curve R_onset = √(GM_bar/a_eff): rms 0.277 dex,
  bias +0.022 dex — onset radii sit ~5% above the pure law.
- Free fit prefers slope 0.272 (rms 0.256): the known F5
  size–mass degeneracy — the slope test is degenerate with the disk
  size–mass relation; the normalization test is not.

## 13b. Cross-corpus significance audit — GNSS + J0437 (2026 audit)

**TEP-GNSS (step_4_4b_autocorrelation_audit.py ->
results/outputs/step_4_4b_autocorrelation_audit.json).**  The
manuscript's gravitational-temporal section quoted raw p-values
(10^-29..10^-59) for all but one smoothing window while §2.3's own
Bretherton N_eff protocol was applied only at 227 d.  The audit
replicates the step_4_4 computation on the pipeline-exported daily
series (site/data/step_4_4) and computes corrected values for every
window: N_eff = 46.7 throughout; r = -0.362 (31 d, p_corr = 1.3e-2)
strengthening monotonically to r = -0.503 (227 d, p_corr = 3.3e-4);
unsmoothed r = -0.164, p_corr = 2.5e-5.  Because the reported window
is the maximum |r| of a 7-window scan, the honest headline is the
selection-adjusted p ~ 2.3e-3 (~3 sigma) — the seasonal-timescale
anti-correlation survives consistent correction at ~3-4 sigma rather
than the raw 10^-59.  Individual planetary values previously quoted
(Jupiter -0.256, Sun -0.230, Mars -0.111) did not trace to current
outputs; corrected traced values on coherence_std: Jupiter r = -0.175,
p_corr = 7.2e-6 (strongest); Saturn -0.101, 1.0e-2; Sun -0.065, 0.10;
Mars -0.072, 0.068.  Manuscript §3.4.1 updated to corrected values;
step_4_4 code now records the full window_scan and per-planet
corrected p-values natively for future runs.

**TEP-J0437 (step_060_mcl_epoch_sem_audit.py ->
results/step_060_mcl_epoch_sem_audit.json).**  The delay-domain
M_cl excess t = 12.7 sigma used H_sem = 1/sqrt(sum w_i) = 0.102 ns —
a measurement-noise IVW SEM, not the declared epoch independence
unit (triplets within an epoch share arclets/legs).  Exact
replication of the pipeline estimator (H_mean = 8.1002 ns) under
epoch-level bootstrap gives SEM = 0.260 ns (iid) for the raw
magnitude; under contiguous 15-240-epoch blocks the raw-mean SEM
rises 0.39 -> 0.55 ns without reaching a plateau (epoch-mean ACF
~0.15-0.2 out to lag 200), so the blocked value is a lower
envelope, and resampling collapses for blocks ~n/3.  The first-pass
step_058 audit divided the noise-subtracted excess by that raw-mean
SEM while holding the correlated folded-normal floor fixed
(corr(m_i, f_i) = 0.91) — the correct statistic is the paired
difference m_i - f_i resampled jointly, giving excess SEM = 0.106 ns
(iid) and 0.12-0.17 ns blocked: the excess is ~12.2 sigma (iid),
~7.6-10.8 sigma (blocked, 15-240 epochs).  Per calendar half the
paired excess is positive in both (1.59 ns, t = 8.9 Jan-Jun; 1.14
ns, t = 8.7 Jul-Dec).  A std-based floor variant bounds the
floor-convention sensitivity (~4-5 sigma).  J1603: paired excess
0.225 ns, SEM = 0.051 ns -> ~4.5 sigma iid, ~3.7-4.1 sigma under
5-30-epoch blocks (raw-mean epoch SEM 0.225 ns).  Signed mean
-0.184 ns: -1.80 -> -1.17 sigma.  The primary Phase Closure
statistic is unaffected — it is already epoch-level (Z = 59.2,
epoch bootstrap CI [+0.737,+1.235] rad; leg-permutation null
p = 0.0025 saturating at 400 reps; the trimmed-excess channel
(step_052) gives 19.1 sigma under iid epoch resampling and
3.8-9.0 sigma under genuine contiguous 15-240-epoch blocks —
step_052 was found to implement iid resampling despite the
'block' name and now runs real MJD-ordered moving blocks; the
trimmed channel carries more campaign correlation than the raw
paired excess (SEM 0.44 iid -> 2.24 ns at 240 epochs, no
plateau)).  Manuscript updated across 0_abstract,
2_theory, 3_data, 4_results, 5_discussion, 6_conclusion,
8_reproducibility; step_060 audit JSON records the replicated
estimator and all SEMs, and supersedes the step_058 audit
(step_058 is the registered scintillation forward-null step).

**TEP-SLR.**  The resampling-artifact resolution was already
propagated: step_2_8 sampling-matched coloured-noise nulls reproduce
the 14.12x TEP-band ratio (white noise through the identical
resample-interpolate chain gives 13.97; 0/46 stations exceed
coloured-noise nulls; real-epoch Lomb-Scargle ratio 1.04 vs null
1.05).  Manuscript correctly presents the spectral channel as a
sampling-kernel property; the surviving claims (pass-bin spatial
test p = 0.0025 nominal / 0.046 label-swap; sub-diurnal coherence
binomial 0.088) are already appropriately caveated.

**GNSS auxiliary-analysis circularity (2026 audit).**
`scripts/null_model_covariance.py` silently regenerated synthetic
distance bins from the exponential best fit whenever the pair-level
`step_2_0_correlation_data_*.csv` was absent — the existing
`null_model_covariance_summary.json` therefore showed r2 = 1.0 and
"exponential wins" by construction.  The script now returns a
`synthetic` flag, prints a CIRCULAR warning, and the regenerated
JSON marks every centre `data_source =
synthetic_bins_from_exponential_best_fit, circular_if_synthetic =
true`; its AIC ordering is non-evidentiary.  The manuscript's
model-comparison claims are unaffected — they trace to step_2_0's
fit on the real 28 bins.  Clustered inference already exists in
the pipeline: pair-level bootstrap CIs exclude zero amplitude in
all three centres ([0.073,0.123] CODE, [0.164,0.222] IGS,
[0.242,0.283] ESA; clustered half-width ~4x the formal error) and
LOSO/LODO/station-block CV leave lambda stable.

**J0437 leg-permutation floor (step_059).**  The step_055 xepoch
null at 400 reps reported p <= 2.5e-3, its resolution floor.
Extended to 5,000 replicates on the cached leg phases: 0
replicates reach the real Rbar = 0.3076 (null mean 0.066, sd
0.032, max 0.190).  Empirical bound p < 2.0e-4; parametric
z = 7.6 sigma (Gaussian-tail p ~ 1.8e-14, reported as indicative
only).  Manuscript §4 updated.

## 14. Next actions (plan order)

1. ~~T1.3 / T-W3~~ — done (steps 56, 59, 020/020c; **PASS** —
   α = 0.347 vs 0.366 (1.6σ), R_s = 2940 vs 2646 AU (0.5σ), ~0.1σ/0.4σ
   under contamination-aware normalization; baseline and exp_interp
   miss amplitude by ~30σ — WB uniquely selects the two-branch sector)
2. ~~T1.2~~ — done (step_57; PASS §9)
3. ~~T1.4~~ — done (step_61; PASS §11)
4. ~~T-B1~~ — done (step_58; asymptotic, §10)
5. ~~T1.5~~ — done (step_66; §7b — outer-edge v within ~7%,
   0.089 dex; RAR is already the blind-curve content of T1.6)
6. ~~AUD-4~~ — done (steps 64/65; junction z_j ≤ 0.687 regularity-
   bounded; T4.2 resolved: no acoustic bound under endpoint-ratio
   transport + starvation law — Papers 18/27 reconciled structurally)
7. **AUD-3 remainder** — ambient map CLOSED at coefficient level
   (step_67 ambient_closure): u_amb = √(4√2πGρL/Kg_t), L = 837 pc
   calibrated once at the midplane inversion, all five strata
   predicted within 10.1%.  Remaining piece: the fσ8 in-well X_env
   volume weighting.  Regime map computed (step_62): extended
   regions carry spatial ξ ~ 1e-3 (below ambient ξ_t = 0.5);
   compact objects ξ ~ 1e20–1e23; the operative in-well variable
   is the positive-operator vertex map x = (1−y)/y — y→0 in wells
   gives x_env→large (S_Σ→0, G_eff→1), collapsed regions grow
   Newtonianly.  Sign question (X = X_t − X_s) resolved in favour
   of the vertex response.  The collapsed-fraction convolution is
   bounded below measurement precision: in-well weighting can move
   fσ8(0.61) by at most the ambient boost (0.489 − 0.477 = 0.012 ≪
   the 0.034 datum error), and the two-branch fit already sits at
   the ΛCDM χ² floor — no observable consequence until the data
   tighten ~3×.  AUD-3 remainder closed in practice; the formal
   convolution remains as a refinement.
8. ~~AUD-5~~ — done (step_68): GW170817 lies at z = 0.0098 ≪ z_j —
   the path integral is unchanged pointwise.  Under the corrected
   ambient the tail drift u_dot = p(1+z)/η slows ~10× at z ~ 3,
   raising the realized lapse cap from 0.03 (LCDM-image) to 0.744
   at z = 1.23 — margins narrow but remain excluding: B0_GW = 77.7
   exceeds the cap 104×, B0_LENS = 422 exceeds it 567×.  The
   ambient-gate resolution stays mandatory.
9. **AUD-6** — well-BBN / constraint slice: gate machinery exists
   (TEP-BBN gate10b, saturating G ~ 1e-3 in wells); the regime map
   above supplies the in-well variable structure; the BBN-epoch
   ambient inputs are computed (step_68 bbn_probe): u_bar = 20.7,
   u_dot = 1.9e-11 H0 (asymptotically static under the tail),
   envelope shape = 0 — the BBN ambient is gate-closed, and
   proto-wells at δ ~ 1e-5 are ambient-dominated.  Remaining:
   the collapsed-fraction convolution in #7.
10. **WB residuals** — 3D off-axis pair solve DONE (step_75):
    y_pair falls 0.270 → 0.158 parallel → perpendicular (two-branch),
    confirming the ∥-solve was the optimistic bound; estimator-level
    orientation averaging remains the residual (~2× downward knob,
    same direction).  Ambient z-profile mapping done (negative result
    for |∇φ| channels — recorded §8).

================================================================================
MANUSCRIPT HARMONIZATION PASS (post-AUD closure)
================================================================================
The derivations closed above have now been propagated into the prose
layer; several stale statements predating the closure were corrected
rather than defended.

1. **lambda branch harmonization.**  Paper 0 Sec 3 previously quoted
   only lambda_ref = 7.526e-71 while TEP-BH quoted 7.5e-66 with no
   stated relation.  Both sites now state the canonical relationship:
   lambda_Cassini = 1e5 * lambda_ref is the operative weak-field value
   (corrected linear-in-S_Sigma Cassini evaluation fails at lambda_ref
   by ~250x, gamma_PPN-1 ~ -5.8e-3 vs bound 2.3e-5; step_03); both are
   fixed entries of the shared corpus constants
   (LAMBDA_QUARTIC_REF / LAMBDA_QUARTIC_CASSINI), not per-paper choices.

2. **Environmental ordering direction — real inversion fixed.**
   Paper 0 Sec 3 (and the older bookkeeping it summarized) stated the
   midplane required a SMALLER effective ambient (u_eff 0.52 midplane
   vs 0.71 halo) with R_s and alpha "rising toward midplane".  Paper
   13's data statement is the opposite: midplane is MORE suppressed
   (alpha = 0.244 vs 0.401 at |Z|<0.1 vs >0.15 kpc), requiring the
   ambient ~36% STRONGER at midplane (u_eff 0.71 vs 0.52).  Paper 0 now
   agrees.  Corollary: the linear estimator u = |a|/g_t fails on SIGN
   as well as magnitude (~4% vs required ~36%) — it grows with |Z|
   while the required ambient strengthens toward the midplane.  The
   ordering is produced by the k-branch density response
   u_amb = sqrt(4 sqrt(2) pi G rho L / (K g_t)), L ~ 0.84 kpc,
   five strata within ~10%, p ~ 0.46-0.5 — an operator output, not a
   proxy.  Historical inverted statement superseded.

3. **Kinetic-action prose.**  TEP-BH 04_coupled_solution.html's stale
   description "K(X) = 1 + (X/Lambda^4)^n" replaced by the action the
   step_47 well solve actually uses: P = X - V + X|X|/Lambda^4,
   P_X = 1 + 2|X|/Lambda^4 on the spacelike branch.  TEP-UCD
   appendix_c_eft.html no longer calls a candidate kinetic completion
   "the corpus's primary screening vehicle"; completions are now
   labeled candidate/benchmark/canonical correctly.

4. **Gate-A propagation.**  TEP-WB 2_theory.html now carries the
   k-branch ambient-map identification (same bracketing honesty as
   6_discussion: two-bin check brackets rho^{1/3}-rho^{1/2},
   five-stratum inversion fixes p ~ 0.46-0.5) and replaces the "open
   Gate-A question" ending with the embedded-response result
   (D^{-1} ~ u_amb^{-2.0}, rho_bg sensitivity <0.01% at the
   operative lambda_Cassini branch).  Field-value vs
   shear-variable distinction made explicit (u_min ~ rho^{1/3} is the
   field-value equilibrium; the pair ambient is the shear variable).
   The TEP-H0/TEP-VOID "Gate-A" references are the Cepheid/ladder
   transfer sector (S_Cep * Gamma_Cep) — a different, still-open item,
   correctly labeled.

5. **Sites rebuilt:** TEP (21 sections), TEP-WB (10), TEP-BH,
   TEP-UCD — all clean.

Remaining open items: B(phi) microscopic carrier, Gamma_WB projector
(partially constrained), and the Cepheid-transfer Gate-A sector.
The f sigma8 in-well convolution is closed by step_74 (two-branch
chi2 = 26.8-27.0 at the LCDM floor); the 3D off-axis pair solve is
closed by step_75 (orientation dependence confirmed at the ~2x
level, downward for averaged populations, theta=0 same-side vs the
axisymmetric reference).


---

## Step 69 — resolved scalar-dipole radiation (added 2025-09-30)

`scripts/steps/step_69_dipole_transmission.py` -> `results/step_69_dipole_transmission.json`.

The step_31 residual debt ("exact perturbation propagator on the nonlinear shell") is closed from the action, not by selecting a reading:

- **Charge conserved.** The static scalar equation is a Gauss law for the current P_X grad phi: the shell suppresses grad phi inside r* but cannot change the enclosed charge. Each compact body carries alpha_A = alpha_0(1-2 s_A); the beta_0 = 0 "identical bare charges" premise in step_31 was wrong (it removed scalarization only, not sensitivities) and is corrected there and in Paper 0 SS8 and Paper 28 SS8.6.
- **Suppression in propagation only.** The l=1 radial equation with the fluctuation metric Z_par = 1+6u, Z_perp = 1+2u (matching the step_32 embedded propagator) is solved with an outgoing-wave condition at x_w = 5 sqrt(2)/Om and an absorbing inner boundary; the transmitted flux is compared to the canonical solve. T_amp ~ a/r*: J1738+0333 1.4e-5, J0348+0432 8.4e-6, J0737-3039 8.3e-6, J0437-4715 4.4e-5. Effective asymmetries Delta_alpha_eff ~ 5e-7 - 2e-5 against bounds 4e-3 - 1e-1: margins ~5e2-6e3 under the full EOS sensitivity bracket.
- **J0337+1715 SEP channel** (static, so propagation suppression does not apply): eta_pred <= 2.5e-7 (weakest single-vertex reading) vs the 2.6e-6 bound.
- **Pre-committed exclusion criterion** recorded in the step docstring: Delta_alpha_eff above the bound under the EOS bracket excludes the specified P(X) completion.
- **Causality ledger:** in-shell radial scalar sound speed c_r^2 -> ~3 c^2 for u >> 1 (standard k-essence; the matter metric carries no superluminal signal). Logged alongside the step_52 disformal margins.

Superseded: the step_31 three-reading bracket (static-vertex conflates the screened static response with the charge; canonical sqrtZ is a point evaluation; Vainshtein (a/r*)^{9/2} is a cubic-Galileon impedance the specified P(X) operator does not generate).

## Step 70 — standard-siren invariant check (added 2025-09-30)

`scripts/steps/step_70_siren_invariant.py` -> `results/step_70_siren_invariant.json`.

The reviewer question on the R10 siren chain — "if the source is locally standard, does the standard (1+z) mass map return?" — is settled symbolically. With locally-measured constants the source observer infers M_loc = M_tilde (1+alpha_A alpha_B)/(1+alpha_e^2) -> M_tilde: the binary IS standard in its own units, so no hidden source-side assumption exists. The end-to-end dimensionless invariant I = G_loc M_chirp f/c^3 obeys I_det/I_src = A_o^2(1+alpha_o^2)/(1+alpha_A alpha_B) -> 1 today. The (1+z)^{-1} map survives because the nonstandard link is transport — coordinate frequency conserved on the static background — not source physics.

Exposed cost, stated in the manuscript (SS9): the local Cavendish constant scales as G_loc ∝ A^2(phi)(1+alpha_amb^2), so Gdot/G = 2 alpha_amb phidot/M_Pl. Solar-screened: ~2e-17/yr, orders below the 4e-13/yr LLR bound. Unscreened-ambient rate ~1.4e-10/yr marks the scale confronting the LVK mass-spectrum population channel — the item that decides whether the epoch-varying G_loc is falsified.

## Step 71 — intrinsic-mass backreaction on the siren discriminant (added 2025-09-30)

`scripts/steps/step_71_intrinsic_mass_drift.py` -> `results/step_71_intrinsic_mass_drift.json`.

The step_21 chirp identity was also corrected to the locally-measured-G form: G_loc,o M_det = A_o G_dyn M_e, i.e. M_det = (A_e/A_o)(1+a_A a_B)/(1+a_o^2) M_tilde — identical at A_o=1, and reducing to the source-local mass at A_o=A_e (the old A_o A_e form violated that locality check). No headline numbers change.

Backreaction: interior k-mouflage tracking makes G_loc(z) ∝ (1+z)^{-2 eps} inside source galaxies; characteristic mass scales ∝ G^{-3/2} then drift as (1+z)^{3 eps}. M_det ∝ (1+z)^{3 eps - 1} vs GR's (1+z): eps = 2/3 nullifies the discriminant, eps = 1 (naive additive tracking) flips it to (1+z)^2 — 'over-massive' high-z BBH as a signature — and eps = 0 restores the naive map. Cross-sector bound, resolved by step_73 (`step_73_sn_luminosity_bound.py` -> `results/step_73_sn_luminosity_bound.json`): the Chandrasekhar channel gives dm = -7.5 eps log10(1+z), unabsorbable by ladder normalization; the SN residual floor (~0.10-0.15 mag) bounds eps <= ~0.05-0.10. Full tracking predicts -2.3 mag at z=1 — excluded by an order of magnitude — so interiors are effectively pinned against ambient drift over the SN baseline (consistent with shallower early wells and path-sector drift). The discriminant survives at near-naive strength: M_det ~ (1+z)^{-0.70..-0.85} vs GR's (1+z), separation (1+z)^{1.7-1.85}. Deriving the residual eps from the interior solution remains the shared item with the compact-star solver.

## Step 72 — unequal two-centre screening diagnostic

`scripts/steps/step_72_mutual_screening_topology.py` -> `results/step_72_mutual_screening_topology.json`.

The Section 2.2 habitat/probe assignment is a leading-order test-probe projection, not a field-level decoupling. Step 72 reuses the step_56 nonlinear finite-volume solver for two unequal Gaussian sources (charge ratio 1:81, separation 3 r*, parallel ambient 0.15 g_t); the specified two-branch example has a sampled bridge minimum |grad phi|/g_t = 0.203 and P_X = 2.34 there versus 1.72 in the ambient. The kinetic current P_X dphi can have nonzero curl while dphi and d ln A remain exact. This curl is not the GR-subtracted disformal photon synchronization connection; it cannot by itself restore conformal propagation delays or lensing. The two-branch P_X tends to the regulator at zero gradient, not 1 as in the baseline sector. No Earth-Moon-Sun boundary-value solve, finite-body scalar charges, flyby tracking calculation, LLR range fit or convergence sweep is claimed. The full finite-mass calculation is needed before changing the Cassini/LLR/flyby inference.

## Paper 17 Step 084 — physical-scale Earth–Moon bridge

The physical-scale bridge diagnostic, tests, and JSON live in `TEP-LLR/scripts/steps/step_084_earth_moon_solar_bridge.py` and `TEP-LLR/results/outputs/step_084_earth_moon_solar_bridge.json`. Paper 0 retains the general kinetic-sector two-centre diagnostic (step_72); Paper 17 owns the matched-profile solar-orientation scan and its comparison to the residual-channel fit $\eta_{\rm resid}$. The collinear local saddle is a concrete target for the coupled boundary-value problem, not an integrated force or range prediction. No Cassini/LLR bound is relaxed from a local minimum alone.


## Step 74 — f sigma8 in-well X_env convolution (scalar sector closed at coefficient level)

`scripts/steps/step_74_fs8_inwell_convolution.py` -> `results/step_74_fs8_inwell_convolution.json`.

Closes the residual flagged in step_62: the mode-weighted collapsed fraction (Press-Schechter, sigma_M = sigma8 (M/M8)^{-0.28}) reverts to G_eff = 1 through the x_env -> large limit of the in-well regime map. Mass-weighted (M > 1e12) and field-weighted (M > 1e10, lognormal-delta resolved) dilutions agree: two-branch chi2 = 27.0/26.8 vs the LCDM control floor 26.0, fs8(0.61) = 0.485/0.483 vs data 0.436 — the residual is the shared S8 tension, not a TEP excess; baseline remains at 33.5-35.8.

## Paper 17 Step 086 — coupled Earth–Moon–Sun nonlinear solve at syzygy

`TEP-LLR/scripts/steps/step_086_earth_moon_coupled_pde.py` -> `TEP-LLR/results/outputs/step_086_earth_moon_coupled_pde.json`.

Extends the step_84 matched-profile ansatz to the real nonlinear three-source PDE at the one axisymmetric configuration (syzygy: Sun–Earth–Moon collinear, pair axis parallel to the solar gradient). Two Gaussian sources (Q_E = 4 pi a_E, Q_M = 4 pi a_M in g_t/d_EM units) solved on the step_56 cylinder grid with u0-continuation. Result: the bridge minimum is |grad phi| ~ 335 g_t at z/d ~ 0.93 in both sectors — the solar ambient gradient dominates the pair bridge, P_X ~ 1.2e5 there, screened volume fraction ~ 0. The matched-ansatz near-exact-cancellation point (bridge min ~ 1e-13) does NOT survive the nonlinear solve: the sources' flux-law saturation prevents the fields from annihilating to the ansatz floor. No P_X < 1 region forms along the pair axis at syzygy; off-syzygy orientations remain bounded by the step_84 scan, and a full 3D solve plus LLR range-model integration remains the residual.

## Step 75 — 3D off-axis two-centre solve (orientation coverage closed)

`scripts/steps/step_75_offaxis_pair_3d.py` -> `results/step_75_offaxis_pair_3d.json`.

Cartesian 3D generalisation of the step_59 axisymmetric solve: div J = rho with the same flux law on a 64^3 box (half-length 4 d), pair separation vector at theta = 0, 30, 45, 60, 90 deg to the ambient gradient u0 z_hat; parameters matched to step_59 (sigma = 0.15, u0 = 0.15, Q = 4 pi per star, d = 1). Newton solver with the exact 19-point anisotropic Jacobian (finite-difference verified to <0.1% interior); mutual force along n_hat by the step_59 stress convention T_ij = f a_i a_j - delta_ij W through a sphere with B-only subtraction; linear-sector reference at each orientation.

Results: theta = 0 lands at y_pair = 0.270 vs the step_59 axisymmetric 0.152 (same-side; coarse 3D source representation inflates the absolute forces ~20% but the ratio is stable between N = 24 and N = 64). The response is orientation-dependent in the direction the step_56 transverse-monopole proxy predicted: y_pair falls monotonically 0.270 -> 0.158 parallel -> perpendicular for two-branch (baseline 0.460 -> 0.376), bridge minimum |a| rising 0.21 -> 0.32 g_t with bridge P_X 2.4 -> 3.8. A perpendicular pair is ~1.7x more suppressed than parallel — the parallel solve is the optimistic bound for the orientation-averaged population, so the WB plateau estimate is conservative in amplitude. Residual: this is an orientation diagnostic at fixed (u0, d, sigma), not the estimator-level integration over the observed orientation distribution, selection function and noise — that mapping belongs to the WB pipeline (TEP-WB step_017 family).

## Step 76 — ambient-excursion transmission into interior clocks

`scripts/steps/step_76_ambient_transmission.py` -> `results/step_76_ambient_transmission.json`.

Measures the clock-channel transmission T = du_int/du_amb directly on the solved unified profile (reuses the step_53 solver): the outer Dirichlet value u_amb is perturbed and the interior-median response is read off for Moon, Earth, Sun, Cepheid envelope and a retained-branch toy control. Motivation: the canonical law S_A = min[1,(rho_bar/rho_T)^{1/3}] is the equilibrium-level ratio u_eq(rho)/u_eq(rho_T) — correct only if interiors sit at the matter-coupled minimum; step_53 placed real bodies on the leaking (flux-following) branch, where TEP's kinetic screening does not require the interior to reach u_eq at all.

Results: leaking-branch bodies transmit ambient field-value excursions ~additively — T = 0.997 (Moon), 0.996 (Earth) at the 2.5e-11 LLR excursion scale, 0.99+ at 1e-9 for all four real bodies; departures are of order the body's own plateau depth (roughly T ~ delta/(delta + u_int): Sun 0.81, Cepheid 0.64 at 2.5e-11). Retained control pins (T ~ 0 for delta << u_eq). Verdict: excursion transmission ~ unity, not (rho/rho_T)^{1/3}.

Propagated: TEP-LLR step_083 regenerated with delta_phi_eff = T_Moon delta_phi_sun = 2.53e-11 (projection ceiling 9.7 mm vs measured 5.1 mm; required ephemeris transfer ~0.5, still order-unity). step_70 corrected: Gdot/G ~ 2(A-dot/A) — the A^2 factor dominates and the screened shear charge alpha_0 S_Sigma ~ 1e-7 enters only its own term — so the LLR bound |Gdot/G| < 4e-13/yr pins the ambient drift at the Earth-Moon location to |A-dot/A| <~ 2e-13/yr = 0.0028 of H0 (output `local_drift_fraction_of_H0`). Paper 17 gains a secular-channel subsection (the same clock-sector projection carries drift; common mode enters clock residuals only at the Earth-Moon differential, the binding observable is orbital) and the p=0.64 residual null is qualified as post-fit, not a drift bound. Consequence recorded for the cosmology programme: interior pinning is NOT available as a drift escape — the local ambient itself must be ~static, and the redshift trend must live in the emission-epoch contrast / spatial landscape structure / non-exact transport. Remaining item: the depth bookkeeping S_A (equilibrium-level ratio) is distinct from the leaking-branch plateau depth (u_int ~ tail-matched, ~ 4 sqrt(Q_bare g_t)/c scaling); which of the two enters each channel's static clock depth (GNSS, Cepheid, kappa coefficients) is the next propagation check.

Scope note (non-integrability, added after review): T~1 is a single-body Dirichlet statement — "interior rides the local ambient value at its boundary." In multi-centre configurations the local ambient itself is supplied by the non-exact solution: on the static branch J = P_,X du can have dJ = dP_,X wedge du != 0 (two-branch current, Appendix A; the pair/bridge solves of steps 59/75 and TEP-LLR 84-86). T therefore applies per body to whatever ambient the solved multi-centre field supplies — which need not equal the superposed-analytic ambient — rather than fixing that ambient's value. Separately, a non-exact scalar current does NOT by itself generate synchronization holonomy: d ln A remains exact; the disformal/metric connection must carry the loop residual. The clock-amplitude transmission result neither creates nor removes holonomy — different connection, different observable.
