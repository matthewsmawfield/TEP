"""Compare the conservative bound with a dense constrained eigensolve."""
import sys
import unittest
from pathlib import Path
import numpy as np
from scipy.linalg import null_space, eigvalsh
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/steps'))
import step_82_stability_lower_bound as audit
from step_79_inhomogeneous_eigenmodes import build_anisotropic_divergence_grad


class BoundTests(unittest.TestCase):
    def test_dense_constrained_schur_operator(self):
        def product(*matrices):
            value=matrices[0]
            for matrix in matrices[1:]:
                value=np.einsum('ij,jk->ik',value,matrix,optimize=False)
            self.assertTrue(np.all(np.isfinite(value)))
            return value
        # Avoid N=4: one grid point per well period aliases this landscape
        # to a constant, hiding a reversed background-Laplacian sign.
        N=6
        data=audit.build_landscape(N,1.,4,.05,1.,.8,2.,.05,2.,1.,15.)
        lap=-audit.fd_laplacian(N,1/N)
        sol=audit.solve_landscape(data['S0'],data['g'],-lap,verbose=False)
        self.assertTrue(sol['converged'])
        psi,u=sol['psi'],data['u']
        self.assertGreater(float(np.ptp(psi)), 1e-5)
        e=1e-6
        V=lambda x:audit.landscape_potential(x,2.,1.,15.)
        vp=(V(u+e)-V(u-e))/(2*e)
        vpp=(V(u+e)-2*V(u)+V(u-e))/e**2
        c=-5*psi**4/(4*audit.M_PL**2)*(data['S0']+sol['lam'])-data['g']/(8*audit.M_PL**2)
        B=np.diag((psi**5/(4*audit.M_PL**2)*(vp-data['rho_m']*np.exp(-u))).ravel())
        Q=null_space(np.ones((1,N**3)))
        L=lap.toarray()+np.diag(c.ravel())
        schur=product(B,Q,np.linalg.solve(product(Q.T,L,Q),product(Q.T,B)))
        D=np.diag(psi.ravel()**-3)
        Z=null_space((psi.ravel()**3)[None,:])
        rows=audit.calculate(N)['branches']
        for name,ambient in [('ambient',.32),('pure_local',0.)]:
            xi=np.maximum(.5*data['g']+ambient,1e-30)
            K=build_anisotropic_divergence_grad(N,1/N,psi**2*(16.03*np.sqrt(xi)+2*xi),psi**2*(16.03*np.sqrt(xi)+4*xi),u)
            H=K.toarray()+np.diag((psi**6*(vpp+data['rho_m']*np.exp(-u))).ravel())-schur
            lowest=eigvalsh(product(Z.T,D,H,D,Z))[0]
            with self.subTest(branch=name):
                self.assertLessEqual(rows[name]['operator_lower_bound'],lowest+1e-9)


if __name__=='__main__':
    unittest.main()
