# Closure directory

Working space for the TEP action-level closure programme, relocated to the TEP
project so that the definitive closure work lives with the foundational paper.

## Canonical artifacts

- `CORPUS_DERIVATION_GR.md` — the corpus-level derivation skeleton
  (two-metric structure, field equations, Theorems 1–2, holonomy gauge
  invariance, redshift endpoint identity, PPN screening map, EFT closure).
- `closure_status.yaml` — the closure status registry.
- **Definitive screening closure**: `../scripts/steps/step_27_master_action_screening_closure.py`
  → `../results/step_27_master_action_screening_closure.json`, cited as
  Appendix E R11 and presented in §2.2(ii). The kinetic completion
  `P(X,phi) = X - V(phi) + X|X|/Lambda^4` with `Lambda^4 = M_Pl^2 H_0^2`
  derives `S_Sigma = [1 + (g/g_t)^2]^-1`, `g_t = cH_0/(2 beta_A^2)`, and the
  pairwise projection `[1 + (R_s/s)^4]^-1`, `R_s = sqrt(GM/g_t)`, at zero free
  parameters; the quartic `V = lambda phi^4/4` carries the clock-amplitude
  sector `S_A = min[1, (rho_bar/rho_T)^{1/3}]`.

## Supporting analyses

Kinetic-operator checks and field solvers that remain valid:
`recovery_operator_*.py` (non-analytic recovery structure, k >= 4),
`stage6_kinetic_closure.py` (corpus-wide kinetic-class test),
`stage5_corrected.py`, `radial_bvp*.py`, `benchmark_bvp.py`,
`corpus_consistency_theorems.py`, `independent_foundations_checks.py`,
`tep_local_force_ppn3.py`, `vphi_bvp_v5.py` — outputs in `results/`.

## `superseded/`

Withdrawn closure attempts kept for provenance only: inverse-power and
quadratic potential branches (F1-excluded on beta_A = -1), the stage0–stage5
potential-class pipeline, pre-GR derivation docs (CL, GM, GP), and exploratory
one-off probes. Nothing in `superseded/` is cited by any current manuscript.
