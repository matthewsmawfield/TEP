"""Independent geometry checks for the landscape diagnostic."""
from pathlib import Path
import importlib.util
import unittest

import numpy as np
from scipy import fft

SOURCE = Path(__file__).resolve().parents[1]/'scripts/steps/step_20_landscape_existence.py'
spec = importlib.util.spec_from_file_location('landscape_geometry', SOURCE)
landscape = importlib.util.module_from_spec(spec)
spec.loader.exec_module(landscape)


class LandscapeGeometryTests(unittest.TestCase):
    def test_conformal_factors_that_cancel_give_flat_matter_geometry(self):
        # h=psi^4 delta and A=psi^-2 give gtilde_spatial=delta identically.
        # This tests the complete conformal-derivative map against known flatness.
        N = 32
        x = 2*np.pi*np.arange(N)/N
        X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
        psi = 1+.07*np.cos(X)+.03*np.sin(Y)+.02*np.cos(Z)
        f = -2*np.log(psi)
        k = fft.fftfreq(N, d=1/N)
        K = np.meshgrid(k, k, k, indexing='ij')
        k2 = sum(ki**2 for ki in K)
        lap = lambda a: np.real(fft.ifftn(-k2*fft.fftn(a)))
        Rg = -8*psi**-5*lap(psi)
        lap_h, grad_h2 = landscape.conformal_spatial_derivatives(psi, f, lap(f), K)
        Rmatter = np.exp(-2*f)*(Rg-4*lap_h-2*grad_h2)
        self.assertLess(float(np.max(np.abs(Rmatter))), 1e-10)

    def test_signed_potential_diagnostic_matches_the_constructed_field(self):
        params = dict(N=12, L=1., n_w=2, sigma=.05, u_well=1., u_void=.8,
                      rho0=2., rho_ambient=.05, V0=2., rhob_half=1., u_s=15.)
        result = landscape.build_landscape(**params)
        self.assertGreater(float(result['V'].min()), 0.)
        self.assertGreater(float(landscape.landscape_potential(.8, 2., 1., 15.)), 0.)
        self.assertLess(float(landscape.landscape_potential(-.8, 2., 1., 15.)), 0.)
        self.assertEqual(float(landscape.landscape_potential(0., 2., 1., 15.)), 1.)


if __name__ == '__main__':
    unittest.main()
