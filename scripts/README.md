# TEP Scripts

This directory contains utility scripts for the Temporal Equivalence Principle (TEP) publication pipeline.

## Scripts

- `generate_site_pdf.py` — Converts the built `index.html` into a publication-quality PDF with metadata embedding.
- `bump_papers.py` — Helper for versioning and cross-paper consistency checks.
- `utils/` — Shared utility modules.

## Pipeline Steps

Step-specific scripts live in `steps/`.

- `step_01_radial_ode.py` — Radial scalar field BVP solver (Cassini, wide-binary, nested hierarchy). Outputs `results/step_01_radial_ode.json`.
- `step_02_earth_topology_closure.py` — Cross-scale BVP test: injects the reference quartic coupling into an Earth PREM-simplified BVP. Outputs `results/step_02_earth_topology_closure.json`.
- `step_03_rho_T_derivation.py` — Saturation scale derivation: tests whether the reference quartic predicts the universal ratio λ_c/R_T and derives ρ_T from self-consistency. Outputs `results/step_03_rho_t_derivation.json`.
- `step_04_gamma_derivation.py` — Geometric projector derivation: derives Γ_GNSS from the linearized scalar fluctuation Green's function in the ground-receiver geometry; emits the covariance 1/e crossing under source-spectrum, resolution and coupling scans. Outputs `results/step_04_gamma_derivation.json`.
- `step_05_gamma_all_channels.py` — All-channel geometric projector derivation: derives the dimensionless projectors Γ_X for GNSS (|β_A|, quadratic covariance), Cepheids (5/ln10), MSP (1), galactic (Cepheid + redshift transfer), and LLR/flyby (unit range projector with the corpus S_Σ(g) suppression, κ_LLR ≈ 1.6e-14); the wide-binary saturation amplitude is flagged partial. Outputs `results/step_05_gamma_all_channels.json`.
- `step_06_strong_field_hyperbolicity.py` — Hyperbolicity analysis: conformal subsystem strongly hyperbolic; GW170817 bounds the propagation-sector cone split only; necessary conditions Z_t > 0, Z_s > 0 with counterexamples showing metric regularity is insufficient; scalar-Gauss-Bonnet has second-order Horndeski EOM requiring separate characteristic analysis. Outputs `results/step_06_strong_field_hyperbolicity.json`.
- `step_07_quantum_formalization.py` — Quantum formalization: distinguishes established kinematic results (from the proper-time postulate) from candidate dynamic results (spin-1/2, antiparticles, Born rule, renormalizability, Standard Model). Outputs `results/step_07_quantum_formalization.json`.
- `step_08_empirical_tensions.py` — Empirical tension analysis: MGEX λ = 1,862 ± 112 km replication scale (axis 21.4° from CMB dipole), JWST mixed Bayesian evidence across comparison spaces, and LLR/flyby post-fit residual circularity. Outputs `results/step_08_empirical_tensions.json`.
- `step_09_transport_geometry.py` — Symbolic ADM, clock-transport and cosmological consistency identities. Outputs `results/step_09_transport_geometry.json`.
- `step_10_cassini_quartic.py` — Quartic-branch wide-binary test. Outputs `results/step_10_wb_quartic.json`.
- `step_11_derivative_screening.py` — Analytical derivation of derivative screening. Outputs `results/step_11_derivative_screening.json`.
- `step_12_cosmology_solution.py` — Static closed Einstein-frame background check (manuscript §8 benchmark). Outputs `results/step_12_cosmology_solution.json`.
- `step_13_disformal_closure.py` — Single-action admissibility test for the disformal completion; reconstructs B(u). Outputs `results/step_13_disformal_closure.json`.
- `step_14_horizon_balance.py` — Residual-geometry test: horizon saturation of the flow/drift balance. Outputs `results/step_14_horizon_balance.json`.
- `step_15_holonomy_amplitude.py` — Derived synchronization-holonomy amplitude for loop experiments (manuscript §10). Outputs `results/step_15_holonomy_amplitude.json`.
- `step_16_interior_roll.py` — Interior rolling-floor consistency solve on the Painlevé–Gullstrand background (temporal-well roll, A → 0). Outputs `results/step_16_interior_roll.json`.
- `step_17_inhomogeneous_closure.py` — Algebraic closure of the Lichnerowicz Hamiltonian constraint for the universal potential on a non-expanding slice. Outputs `results/step_17_inhomogeneous_closure.json`.
- `step_18_landscape_accommodation.py` — Numerical Lichnerowicz BVP: Planck-scale curvature quarantined inside the dense emitting regions with asymptotically flat exterior. Outputs `results/step_18_landscape_accommodation.json`.
- `step_20_landscape_existence.py` — Inhomogeneous non-compact slice existence demonstration: bordered-Newton sparse-direct solve of the Lichnerowicz constraint on a periodic cell (universal cover R^3) over a well/void temporal landscape; derives the required sign-changing total-energy structure (negative-energy void sector supplied by V(u<0)<0) and the matter-frame flatness ledger |Omega_k^eff| ~ <|grad u|^2>/H_drift^2. Outputs `results/step_20_landscape_existence.json`.
- `step_19_operator_evaluation.py` — Two-body kinetic operator evaluation under the resolved nested operator R = S_Σ(X_env)²·y(s) (flux-conserving profile, embedding-ambient vertex factors; step_30/issue 0-27): Cassini, Saturn, LLR, GP-B suppression factors and the nested-hierarchy field decomposition (manuscript §3). Outputs `results/step_19_operator_evaluation.json`.
- `step_26_weak_field_lapse_ledger.py` — Canonical weak-field normalization ledger: derives `delta_phi/M_Pl = 2 beta_A Phi_N`, the unscreened matter lapse `N_tilde/N_tilde_infinity = exp[(1 + 2 beta_A^2) Phi_N]`, `G_eff/G = 3` at `beta_A = -1`, and screened recovery as `S_Sigma -> 0`. Outputs `results/step_26_weak_field_lapse_ledger.json`.
- `step_28_momentum_shear_sector.py` — Momentum/shear (A_ij != 0) sector of the landscape-slice construction: solves the momentum constraint via the conformal Killing potential for the drift-sourced source `j_i = c_pi d_i u`, then re-solves the full Hamiltonian constraint including the `A^2 psi^{-7}/8` shear term; scans drift amplitude, reports positivity, momentum-constraint residuals, and the linearized-operator spectrum (local uniqueness). Outputs `results/step_28_momentum_shear.json`.
- `step_30_nested_two_body_operator.py` — Bi-directional nested evaluation of the corrected pairwise operator R = q1*q2*y(s) with the exact flux-conserving profile y(1+y^2(r*/r)^4)=1 (transition-region exponent 4/3, mutual-dominated interior asymptote s^4) and per-body ambient response factors q_i = S_Sigma(X_amb,i); evaluates Cassini, Saturn, Mercury, Earth-Sun, Earth-Moon, LLR differential, and the wide-binary vertex-structure readings. Outputs `results/step_30_nested_two_body_operator.json`.

## Archive

`steps/archive/` retains the v0.13→v0.14 closure campaign (former steps 11, 15–27, 29, 33, 33b, the symbolic static-cosmology check, the toy loop example, and the step_34 potential audit) — the systematic exclusion of transport-, volume- and partition-based redshift mechanisms that motivated the endpoint-drift solution — together with their result JSONs under `archive/results/`. These are kept as the documented exclusion record; they are no longer part of the active pipeline. `step_33_interior_backreaction.py` remains cited by Paper 28 (TEP-BH) as the conventional-evolution comparison computation.
