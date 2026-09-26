# Pre-registration — §7.8 Target 3 closure computation (shared environmental clock distortion)

Written 2026-09-13, before any fit of this computation was run.
Script: `closure/tep_shared_distortion_closure.py`. Output: `closure/results/tep_shared_distortion_closure.json`.

## Mechanism (CORPUS_PLAN §7.4)

Under universal coupling, the host environment distorts every clock-based distance indicator in a
host together. Distances to a host with potential coordinate `X_i = (S_i σ_i² − U_ref)/c²` are
underestimated by `κ(X_i − X_anchor)` magnitudes. For Hubble-flow SNe this gives a Hubble residual
`r_i = M − κ X_i` (deeper hosts appear brighter at fixed redshift), and the ladder satisfies
`H₀,ladder = H₀,true × 10^{κ(X̄_HF − X̄_anchor)/5}`.

## Step 1 — κ_SN from Pantheon+ Hubble-flow SNe

- **Sample.** `USED_IN_SH0ES_HF = 1` (0.023 < z < 0.15), `HOST_LOGMASS > 5`.
- **Standardisation without a host-mass step.** `mB − μ_ΛCDM(zHD, zHEL; Ω_m = 0.334) = M − α x1 + β c − κ_SN X_i`.
  α, β, M and κ_SN are fitted jointly.
- **Errors.** `mBERR`, `x1ERR` and `cERR` with their covariances (converted from x0), peculiar velocity
  σ_v = 250 km/s, and σ_int tuned so that χ²/dof = 1.
- **Coordinate.** `σ_i = 10^{−0.16 + 0.211 logM}`, the Paper 11 calibrator relation; `S_i = 1` (field hosts).
- **Recovery test.** Inject κ = 4 × 10⁵ into mock residuals built from the fitted model and recover it
  before the result is read.
- **Variants.** `biasCor_m_b` subtracted; published `m_b_corr`; 0.01 < z < 0.15 non-calibrators;
  `logM > 9` only; σ–mass slope ×0.5 and ×1.5; `S = 0.867` (median calibrator S); joint host-mass step
  (reported to show absorption, not primary).

## Step 2 — agreement with Paper 11's velocity-space κ

`κ_vel = Γ_X / ((ln10/5) H_app)` from Step 42, with the canonical σ_v = 250 km/s as primary:
`(3.69 ± 3.10) × 10⁵`. The 182.1 and 150 km/s values are reported alongside.

**Agreement** means `|κ_SN − κ_vel| < 2σ_combined`. **Detection** means `κ_SN > 0` at ≥ 2σ one-sided.

## Step 3 — zero-free-parameter H₀ prediction

`H₀,true = 73.04 × 10^{−κ(X̄_HF − X̄_anchor)/5}`, where `X̄_HF` is the mean over SH0ES Hubble-flow SNe
and `X̄_anchor = U_ref/c²`.

- Primary: κ = κ_vel (σ_v = 250) with both anchor references, `U_ref = 87.165²` (unscreened) and
  `30.507²` (screened).
- Also reported: κ = κ_SN, κ = 0.96 × 10⁶, and the `X̄_HF` required to reach 67.4.

**Consistency check.** The predicted uncorrected calibrator-host mean `cz/d`, computed from the Step 04
calibrator `Sσ²` values, is compared with the observed 67.72.

## Addendum — TEP versus dust discriminator (written after Steps 1–3, before this test was run)

Step 1 detects `κ_SN > 0`. The same host-potential dependence is conventionally attributed to dust,
following the Brout & Scolnic (2021) model, which the `biasCor_m_b` correction implements. The two
explanations make different predictions.

**Predictions.**
- **TEP.** The distortion is a clock effect and therefore achromatic. `κ` is the same for blue
  (`c < 0`) and red (`c ≥ 0`) SNe, and is independent of stretch.
- **Dust.** The dependence is concentrated in red SNe, with `κ_blue ≈ 0` and `κ_red > 0`.

**Test.** Refit the Step 1 primary model separately for `c < 0` and `c ≥ 0`, and separately for
`x1 < 0` and `x1 ≥ 0`. Also fit an interaction term `κ_c · (c − c̄) · X` on the full sample, using the
same sample, coordinate and errors as Step 1.

**Decision.**
- **Supports TEP (achromatic):** `κ_blue` and `κ_red` agree within 2σ, `κ_blue > 0` at ≥ 1.5σ, and `κ_c` is consistent with 0.
- **Supports dust (chromatic):** `κ_blue` is consistent with 0 and `κ_red` differs from it by more than 2σ, or `κ_c > 0` at ≥ 2σ.
- **Otherwise:** undetermined.
