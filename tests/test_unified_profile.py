import importlib.util
from pathlib import Path
import unittest

import numpy as np
from scipy.integrate import quad, solve_bvp
from scipy.special import erf

SOURCE = Path(__file__).resolve().parents[1] / 'scripts/steps/step_53_unified_profile.py'
spec = importlib.util.spec_from_file_location('unified_profile', SOURCE)
model = importlib.util.module_from_spec(spec)
spec.loader.exec_module(model)


class UnifiedProfileTests(unittest.TestCase):
    def test_earth_profile_reaches_and_satisfies_outer_boundary(self):
        rho_amb = 1e-30
        u_amb = model.u_eq(rho_amb)
        def density(r):
            return model.smooth_density(r, model.REARTH, 5.51, rho_amb)
        for delta in (0.0, 2.5e-11):
            with self.subTest(delta=delta):
                r, u, _, residual, charge = model.solve_body(
                    model.MEARTH, model.REARTH, density, rho_amb, u_amb + delta,
                    n=1000)
                self.assertAlmostEqual(r[-1] / (1e10 * model.REARTH), 1.0, places=12)
                self.assertLessEqual(abs(u[-1] - u_amb - delta), 1e-22)
                self.assertLess(residual, 1e-8)
                self.assertTrue(np.all(np.isfinite(charge)))

    def test_density_profile_reaches_ambient_without_a_floor(self):
        density = model.smooth_density(np.array([0.0, 1e10*model.REARTH]),
                                       model.REARTH, 5.51, 1e-30, 40)
        np.testing.assert_allclose(density, [5.51, 1e-30], rtol=1e-14)

    def test_homogeneous_equilibrium(self):
        rho = 1e-24
        level = model.u_eq(rho)
        r, u, _, residual, charge = model.solve_body(
            1.0, 1e6, lambda radius: rho, rho, level, n=300, r_out_fac=1e6)
        np.testing.assert_allclose(u, level, rtol=1e-13)
        np.testing.assert_allclose(charge, 0.0, atol=1e-25)
        self.assertLess(residual, 1e-10)

    def test_massless_profile_against_independent_gaussian_quadrature(self):
        radius, density, outer = 1e5, 5e-10, 1e8
        source = model.rho_to_m2(density)
        def charge(r):
            x = r/radius
            return source*radius**3*(np.sqrt(np.pi)*erf(x)/4 - x*np.exp(-x*x)/2)
        def gradient(r):
            target = charge(r)/r**2
            q = target*model.C**2/model.G_T
            root = 2/np.sqrt(3)*np.sinh(np.arcsinh(1.5*np.sqrt(3)*q)/3)
            return root*model.G_T/model.C**2
        expected = quad(gradient, 2*radius, outer, epsabs=1e-30, epsrel=1e-10)[0]
        errors = []
        for n in (1000, 2000):
            r, u, _, residual, _ = model.solve_body(
                1.0, radius, lambda rr: density*np.exp(-(rr/radius)**2),
                0.0, 0.0, n=n, r_out_fac=outer/radius, lam_m2=0.0)
            got = np.interp(np.log(2*radius), np.log(r), u)
            errors.append(abs(got/expected - 1))
            self.assertLess(residual, 1e-8)
        self.assertLess(errors[-1], 3e-4)
        self.assertLess(errors[-1], errors[0])

    def test_single_source_gradient_has_one_solved_screening_factor(self):
        radius, density = 1e5, 5e-6
        r, u, _, _, _ = model.solve_body(
            1.0, radius, lambda rr: density*np.exp(-(rr/radius)**2),
            0.0, 0.0, n=2000, r_out_fac=1e3, lam_m2=0.0)
        faces = (r[1:]+r[:-1])/2
        j = int(np.argmin(abs(faces-3*radius)))
        x = faces[j]/radius
        q = model.rho_to_m2(density)*radius**3*(
            np.sqrt(np.pi)*erf(x)/4 - x*np.exp(-x*x)/2)
        linear = q/faces[j]**2
        solved = -(u[j+1]-u[j])/(r[j+1]-r[j])
        ratio = solved/linear
        y = model.finv(linear)/linear
        self.assertAlmostEqual(ratio/y, 1.0, delta=3e-4)
        self.assertGreater(abs(ratio-y**3), 0.5*y)

    def test_quartic_profile_matches_independent_collocation(self):
        radius, rho_amb = model.REARTH, 1e-30
        lam = model.LAM_M2 * 1e5
        boundary = model.u_eq(rho_amb, lam)
        prof = lambda r: model.smooth_density(r, radius, 5.51, rho_amb)
        r, u, _, _, charge = model.solve_body(
            model.MEARTH, radius, prof, rho_amb, boundary,
            n=1800, r_out_fac=1e5, lam_m2=lam)
        field_scale = float(np.max(u))
        charge_scale = model.rho_to_m2(5.51)*radius**3/3
        t = np.log(r/radius)
        def ode(tt, y):
            rr = radius*np.exp(tt)
            field, q = y[0]*field_scale, y[1]*charge_scale
            target = q/rr**2 * model.C**2/model.G_T
            root = 2/np.sqrt(3)*np.sinh(np.arcsinh(1.5*np.sqrt(3)*target)/3)
            slope = root*model.G_T/model.C**2
            source = model.rho_to_m2(prof(rr))*np.exp(-field) - lam*field**3
            return np.array([-rr*slope/field_scale, rr**3*source/charge_scale])
        def bc(left, right):
            return np.array([left[1], right[0]-boundary/field_scale])
        result = solve_bvp(ode, bc, t, np.array([u/field_scale, charge/charge_scale]),
                           tol=1e-6, max_nodes=20000)
        self.assertTrue(result.success, result.message)
        np.testing.assert_allclose(u, result.sol(t)[0]*field_scale,
                                   rtol=3e-4, atol=field_scale*1e-6)

    def test_static_response_matches_small_boundary_perturbations(self):
        rho = 1e-30
        boundary = model.u_eq(rho)
        prof = lambda r: model.smooth_density(r, model.REARTH, 5.51, rho)
        r, u, _, _, _ = model.solve_body(
            model.MEARTH, model.REARTH, prof, rho, boundary, n=1000)
        response = model.static_boundary_response(r, u, prof, rho)
        inner = r < 0.35*model.REARTH
        for delta in (1e-20, 1e-19):
            _, shifted, _, _, _ = model.solve_body(
                model.MEARTH, model.REARTH, prof, rho, boundary+delta, n=1000)
            np.testing.assert_allclose(np.median((shifted-u)[inner])/delta,
                                       np.median(response[inner]), rtol=1e-4)

    def test_shooting_fails_closed_when_it_does_not_reach_boundary(self):
        with self.assertRaisesRegex(RuntimeError, 'outer boundary'):
            model.solve_body(
                model.MEARTH, model.REARTH,
                lambda r: 1e-30+5.51/(1+min(r/model.REARTH, 3.0)**40),
                1e-30, model.u_eq(1e-30), n=300, method='shooting')

    def test_inverse_flux_handles_both_signs_and_zero(self):
        q = np.r_[-np.logspace(-40, 4, 100), 0.0, np.logspace(-40, 4, 100)]
        v = model.finv(q)
        np.testing.assert_allclose(v * model.P_X(v), q, rtol=1e-12, atol=0.0)


if __name__ == '__main__':
    unittest.main()
