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
- `step_05_gamma_all_channels.py` — All-channel geometric projector derivation: derives Γ_X for GNSS, Cepheids, MSP, wide binaries, galactic, and LLR/flyby channels from the linearized scalar fluctuation equation. Outputs `results/step_05_gamma_all_channels.json`.
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
- `step_19_operator_evaluation.py` — Two-body kinetic operator evaluation: Cassini, Saturn, LLR, GP-B suppression factors and the nested-hierarchy field decomposition (manuscript §3). Outputs `results/step_19_operator_evaluation.json`.

## Archive

`steps/archive/` retains the v0.13→v0.14 closure campaign (former steps 11, 15–27, 29, 33, 33b, the symbolic static-cosmology check, the toy loop example, and the step_34 potential audit) — the systematic exclusion of transport-, volume- and partition-based redshift mechanisms that motivated the endpoint-drift solution — together with their result JSONs under `archive/results/`. These are kept as the documented exclusion record; they are no longer part of the active pipeline. `step_33_interior_backreaction.py` remains cited by Paper 28 (TEP-BH) as the conventional-evolution comparison computation.
