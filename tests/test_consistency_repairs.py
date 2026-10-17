"""Physical and analytic regressions for the October consistency repairs.
Run: python -m unittest discover -s TEP/tests -p 'test_consistency_repairs.py'
"""
import sys
import unittest
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import eigsh, norm

STEPS = Path(__file__).resolve().parents[1] / 'scripts/steps'
sys.path.insert(0, str(STEPS))
from step_21_lapse_solvability import solve_unshifted_maximal_lapse
from step_79_inhomogeneous_eigenmodes import scalar_bloch_operator
from step_80_shared_lapse_action_derivation import radial_tail_integral
from step_81_action_contract_audit import derive as derive_action_contract


class ConsistencyRepairs(unittest.TestCase):
    def test_radial_tail_against_analytic_integrals(self):
        r = np.geomspace(1., 100., 4000)
        np.testing.assert_allclose(radial_tail_integral(r, 2*r), 100**2-r**2,
                                   rtol=1e-12, atol=1e-10)
        np.testing.assert_allclose(radial_tail_integral(r, r**-2), 1/r-1/100,
                                   rtol=1e-6, atol=1e-12)
        self.assertEqual(radial_tail_integral(r, r**-2)[-1], 0.)

    def test_bloch_eigenvalue_against_constant_coefficient_solution(self):
        n, f, m, theta = 6, 2.5, .37, .8
        ones = np.ones((n, n, n))
        op = scalar_bloch_operator(f*ones, m*ones, 1/n, (theta, 0., 0.))
        self.assertLess(norm(op-op.conjugate().T), 1e-12)
        ev = eigsh(op, k=1, which='SA', return_eigenvectors=False, tol=1e-11)[0]
        expected = m + 4*f*n*n*np.sin(theta/(2*n))**2
        self.assertAlmostEqual(ev, expected, places=9)

    def test_admissible_evolving_lapse_and_rejected_static_shift(self):
        result = solve_unshifted_maximal_lapse(N=12)
        self.assertTrue(result['converged'])
        self.assertLess(result['hamiltonian_max_residual'], 1e-8)
        self.assertLess(result['required_Pi2_for_maximal_closure'], 0.)
        self.assertFalse(result['maximal_canonical_admissible'])
        self.assertTrue(result['cmc_valid'])
        self.assertGreater(result['cmc_alpha'].min(), 0.)
        self.assertLess(result['cmc_relative_residual'], 1e-8)
        self.assertAlmostEqual(result['cmc_volume_mean_lapse'], 1., places=12)
        self.assertGreater(result['cmc_Kdot'], 0.)

    def test_acoustic_shell_reduction_matches_bbn_step20(self):
        audit = derive_action_contract()
        red = audit['acoustic_shell_reduction']
        eff = audit['effective_acoustic_prescription']
        self.assertTrue(red['matches_bbn_step20_within_1e12'])
        self.assertLess(red['flux_law_residual'], 1e-12)
        self.assertAlmostEqual(
            red['K_sp_spatial_sound_shell_derived'],
            eff['K_sp_spatial_sound_shell'],
            places=12,
        )
        self.assertAlmostEqual(
            red['c_p_eff_squared_derived'],
            eff['c_p_eff_squared_derived'],
            places=12,
        )


if __name__ == '__main__':
    unittest.main()
