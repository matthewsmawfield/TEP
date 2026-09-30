"""Independent numerical controls for the short holonomy diagnostic."""
import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch
import numpy as np
from scipy.integrate import quad
R=Path(__file__).resolve().parents[2]
s=importlib.util.spec_from_file_location('holonomy',R/'TEP/scripts/steps/step_15_holonomy_amplitude.py');h=importlib.util.module_from_spec(s);s.loader.exec_module(h)
class HolonomyTests(unittest.TestCase):
 def setUp(self):self.kw=dict(B0=1.,u_thresh=.5*h.u_field(np.array([[h.R_E,0,0]]))[0],Pi_bar=h.H0)
 def test_solar_gradient_not_clipped_at_one_AU(self):
  x=np.array([[h.R_E+2e6,0.,0.]])
  with patch.object(h,'M_E',0.):
   expected=h.KAPPA*h.G_SI*h.M_S/(h.C_SI**2*(h.AU-x[0,0])**2)
   self.assertAlmostEqual(h.grad_u_field(x)[0,0]/expected,1.,places=13)
 def test_gradient_matches_independent_finite_difference(self):
  x=np.array([[9e6,2e6,1e6],[5e6,6e6,2e6]]);eps=100.
  numeric=np.column_stack([(h.u_field(x+np.eye(3)[i]*eps)-h.u_field(x-np.eye(3)[i]*eps))/(2*eps) for i in range(3)])
  np.testing.assert_allclose(h.grad_u_field(x),numeric,rtol=2e-7,atol=1e-24)
 def test_exact_connection_zero(self):self.assertLess(abs(h.loop_holonomy(**self.kw,s0=1.)[0]),1e-22)
 def test_spherical_geometry_zero(self):
  with patch.object(h,'M_S',0.):self.assertLess(abs(h.loop_holonomy(**self.kw)[0]),1e-22)
 def test_zero_and_orientation_controls(self):
  v=h.loop_holonomy(**self.kw)[0]
  for k in ['B0','Pi_bar']:
   kw={**self.kw,k:0.};self.assertEqual(h.loop_holonomy(**kw)[0],0.)
  self.assertLess(abs(v+h.loop_holonomy(**self.kw,reverse=True)[0]),1e-22)
 def test_geometry_and_convergence(self):
  A,S,B,_=h.triangle()
  for a,b in [(A,S),(S,B)]:
   d=b-a;t=np.clip(-a@d/(d@d),0,1);self.assertGreaterEqual(np.linalg.norm(a+t*d),h.R_E-1e-8)
  with self.assertRaises(ValueError):h.triangle(alt_m=1e5,lonB_deg=120.)
  vals=[h.loop_holonomy(**self.kw,npts=n)[0] for n in [32,64,128,256]]
  self.assertLess(max(vals)-min(vals),1e-22);self.assertGreater(abs(vals[-1]),1e-14)
 def test_independent_adaptive_quadrature(self):
  # Explicit parameterized line integrals, independently integrated by QUADPACK.
  angle=np.pi/3;A=h.R_E*np.array([1,0,0]);B=h.R_E*np.array([.5,np.sqrt(3)/2,0]);S=(h.R_E+2e6)*np.array([np.sqrt(3)/2,.5,0])
  def covector(x):
   xx=x[None,:];return h.connection_prefactor(xx,**self.kw)[0]*h.screening(xx)[0]*h.grad_u_field(xx)[0]
  total=0.
  for a,b in [(A,S),(S,B)]:total+=quad(lambda t:covector(a+t*(b-a))@(b-a),0,1,epsabs=1e-13,epsrel=1e-12)[0]
  total+=quad(lambda t:covector(h.R_E*np.array([np.cos(angle*(1-t)),np.sin(angle*(1-t)),0]))@(h.R_E*angle*np.array([np.sin(angle*(1-t)),-np.cos(angle*(1-t)),0])),0,1,epsabs=1e-13,epsrel=1e-12)[0]
  self.assertLess(abs(total/h.C_SI-h.loop_holonomy(**self.kw)[0]),1e-22)
if __name__=='__main__':unittest.main(verbosity=2)
