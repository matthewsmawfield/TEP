# Closure directory

Working space for the TEP action-level closure programme, relocated to the TEP
project so that the definitive closure work lives with the foundational paper.

## Canonical artifacts

- `CORPUS_DERIVATION_GR.md` — the corpus-level derivation skeleton
  (two-metric structure, field equations, Theorems 1–2, holonomy gauge
  invariance, redshift endpoint identity, PPN screening map, EFT closure).
- `closure_status.yaml` — historical registry with a current reconciliation header.
  The older cuscuton/inverse-power entries and their completion labels are not
  the current completion authority. Existing derivations must be located before
  an old open-task flag is treated as a request to derive them again.
- **Kinetic realization and benchmarks**:
  `../scripts/steps/step_27_master_action_screening_closure.py` →
  `../results/step_27_master_action_screening_closure.json`, cited as Appendix E R11.
  The manuscript adopts `P = X - V + X|X|/Lambda^4`; nested pair projections,
  canonical-only controls and alternative low-gradient branches retain their
  explicit assumptions. A single-source gradient, finite-pair response and
  fitted lens-aperture mass are not interchangeable screening quantities.
- **Static field and boundary response**: steps 53 and 76 now use validated
  finite-volume profiles and distinguish infinitesimal response from nonlinear
  finite excursions. Their outputs supersede the old shooting-based values;
  static depth `S_A`, fluctuation covariance and dynamical transmission are
  distinct. Appendix E R17 gives the manuscript scope.
- **Existing radiative and siren calculations**: steps 69–71 retain their
  explicit source, propagation and population assumptions. Step 78 contains an
  existing landscape-transport diagnostic; its expected result artifact was not
  located during reconciliation, so it is not represented as a verified output.

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
