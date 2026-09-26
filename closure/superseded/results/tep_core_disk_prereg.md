# Pre-registration — TEP core–disk observing-chain test on Pantheon+

Written 2026-09-13, before any fit of this test was run.
Script: `closure/tep_core_disk_pantheon_test.py`. Output: `closure/results/tep_core_disk_pantheon_test.json`.

## Mechanism under test (Papers 30, 31; Paper 11 v0.9)

TEP gives `1 + z = A_obs / A_em`, with `A < 1` in deeper potentials. A host's spectroscopic redshift is weighted
toward its core. A supernova sits in the disk, where the potential is shallower, so its clock factor differs.
The core redshift therefore exceeds the redshift appropriate to the SN by `Δ_i > 0`.

The model is `Δ_i = k · D_i`. Here `D_i = (σ_i² / c²) · ln(1 + R_i² / R_s,i²)` is the isothermal
potential difference between the core and the SN's projected radius `R_i`. The constant `k > 0` is the
clock amplification factor, the analogue of κ.

Two effects follow, and both push the same way:

1. **Distance.** The SN is placed at `z_true + Δ`, so its model distance modulus is too large. Its Hubble
   residual shifts by `−(dμ/dz) · Δ`, which grows as `1/z` at low redshift.
2. **Time dilation.** Rest-frame timescales are underestimated by `Δ`, which lowers `x1`. Taking `α ≈ 0.15`
   and `x1 ≈ 10(s − 1)`, the Hubble residual shifts by about `−1.5 Δ` mag.

**Prediction.** Hubble residual `r_i = −k · T_i`, with `T_i = D_i · (dμ/dz|_{z_i} + 1.5)` and **k > 0**.
The TEP-unique signature is the potential × offset × `1/z` structure. Host-mass steps and local-environment
effects do not scale as `1/z`.

## Primary analysis (fixed)

- **Sample.** Pantheon+SH0ES, `IS_CALIBRATOR = 0`, `HOST_ANGSEP > 0`, `HOST_LOGMASS > 5`, `0.01 < zHD < 0.8`.
- **Residual.** `m_b_corr − μ_ΛCDM(zHD, zHEL; Ω_m = 0.334)`. The intercept absorbs `H₀` and `M`.
- **Potential depth.** `σ_i` from `log10 σ = −0.16 + 0.211 · logM`, fitted on the Paper 11 calibrator hosts
  (r = 0.68).
- **Offset.** `R_i = HOST_ANGSEP × d_A`, with `R_s = 3 kpc × 10^{0.25(logM − 10.5)}`.
- **Model.** GLS with the STAT+SYS covariance: `r = a + γ·[logM > 10] + b·log10(R/R_s) + η·T`, and `k = −η`.
- **Decision.**
  - **Support:** `k > 0` at one-sided p < 0.01 in the primary analysis, the same sign in ≥ 80% of robustness variants, and a permutation p < 0.01.
  - **Excluded:** `k` at or below zero at the primary level.
  - **Amplitude check.** The required `k_req` (below) is compared with the fitted `k`. Resolution needs `k_req` inside the fitted 95% interval.

## Robustness (fixed)

- Redshift windows: 0.01–0.15, 0.023–0.15, and all z.
- `R_s` fixed at 1, 3 or 5 kpc.
- σ–mass slope ×0.5 and ×1.5.
- Diagonal errors instead of the full covariance.
- No host-mass or offset controls.
- Duplicate SNe removed.
- Leave-one-survey-out.
- Permutation of `D` across SNe (2,000 draws).

## Amplitude required to resolve the tension

On the SH0ES Hubble-flow rows present in the sample, `ln(73.04/67.4) = (ln10/5) · k_req · ⟨T⟩_HF`, minus the
calibrator time-dilation term where calibrator offsets exist.
