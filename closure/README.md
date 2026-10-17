# Theoretical Action, Sector Admissibility, and Consistency Reference

This directory contains the formal theoretical reference documents accompanying Paper 0 (`TEP`). Together, these documents present the covariant mathematical derivations, master-sector admissibility analysis, and numerical consistency specifications of the Temporal Equivalence Principle.

## Reference Documents

- [`CORPUS_DERIVATION_GR.md`](CORPUS_DERIVATION_GR.md) — Covariant derivation of the two-metric framework, action variation, coupled Einstein–scalar field equations, local Lorentz invariance, synchronization-holonomy gauge invariance, the Mount Wilson cosmological redshift identity, parametrized post-Newtonian (PPN) screening relations, temporal-horizon regularity, and the disformal synchronization connection.
- [`MASTER_SECTOR_LEDGER.md`](MASTER_SECTOR_LEDGER.md) — Systematic admissibility analysis across the kinetic function $P(X,\phi)$, scalar potential $V(\phi)$, disformal coupling $B(\phi)$, and ambient drift $A_{\rm sp}(\eta)$, establishing the parameter-free two-branch kinetic completion $P_{,X}(\xi) = k\sqrt{\xi} + 2\xi$ ($\xi = |X|/\Lambda_X^4$, $k = 4\sqrt{2}\,g_t/a_0 \simeq 16.03$) and its cross-domain observational tests.
- [`CONSISTENCY_AND_PREDICTION.md`](CONSISTENCY_AND_PREDICTION.md) — Formal specification of the analytical and numerical verification suite, covering the unified lapse and spatial-ruler action identities, the evolving constant-mean-curvature (CMC) initial-slice lapse construction, the Lichnerowicz-coupled and Bloch fluctuation spectra, and ray-integrated solar PPN benchmarks.

## Reproducible Pipeline Architecture

All symbolic derivations, boundary-value solvers, and spectral stability computations are implemented in `../scripts/steps/` with corresponding JSON artifacts in `../results/` and automated regression tests in `../tests/`:

- Master-sector admissibility: `../scripts/steps/step_54_master_sector_admissibility.py` $\rightarrow$ `../results/step_54_master_sector_admissibility.json`
- Unified radial profile and static boundary response: `../scripts/steps/step_53_unified_profile.py` and `../scripts/steps/step_76_ambient_transmission.py`
- Inhomogeneous periodic landscape and CMC initial-slice lapse: `../scripts/steps/step_20_landscape_existence.py` and `../scripts/steps/step_21_lapse_solvability.py`
- Coupled Einstein–scalar stability and discrete lower bounds: `../scripts/steps/step_79_inhomogeneous_eigenmodes.py` and `../scripts/steps/step_82_stability_lower_bound.py`
- Shared-cone gravitational action and solar ray propagation: `../scripts/steps/step_80_shared_lapse_action_derivation.py` and `../scripts/steps/step_81_action_contract_audit.py`
- End-to-end consistency validation: `../scripts/steps/step_83_consistency_validation.py` $\rightarrow$ `../results/step_83_consistency_validation.json`
