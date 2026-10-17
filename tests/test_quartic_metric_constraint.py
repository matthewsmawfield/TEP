"""Independent dense checks of the metric constraint and Bloch lower bound."""
import sys
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts/steps'))
from step_79_inhomogeneous_eigenmodes import scalar_bloch_operator


def product(a,b):
    return np.einsum("ij,jk->ik",a,b,optimize=False)

def test_all_twist_bound_on_nonconstant_metric_constraint():
    N=3; dx=1/N
    X=np.indices((N,N,N)).reshape(3,-1).T/N
    q=1+.03*np.cos(2*np.pi*X[:,0]);q/=np.linalg.norm(q)
    Q=null_space(q[None,:]);c=-.3+.1*np.cos(2*np.pi*X[:,1])
    bound=4*np.sin(np.pi/(2*N))**2/dx**2*(q.sum()**2/q.size)+c.min()
    for theta in [(0,0,0),(np.pi,0,0),(-np.pi,np.pi,.7),(.4,-.9,1.3)]:
        Lb=scalar_bloch_operator(np.ones((N,N,N)),c.reshape((N,N,N)),dx,theta).toarray()
        phase=np.exp(1j*(np.einsum("ij,j->i",X,np.array(theta),optimize=False)))
        L=phase.conj()[:,None]*Lb*phase[None,:]  # periodic Bloch gauge
        reduced=product(product(Q.conj().T,L),Q)
        assert np.linalg.eigvalsh(reduced).min()>=bound-1e-10
        # Exact constrained inverse must equal the dense reduced inverse.
        inv=np.linalg.inv(L);iq=np.einsum("ij,j->i",inv,q,optimize=False)
        bordered=inv-np.outer(iq,np.einsum("i,ij->j",q,inv,optimize=False))/(np.dot(q,iq))
        reference=product(product(Q,np.linalg.inv(reduced)),Q.conj().T)
        assert np.all(np.isfinite(bordered)) and np.all(np.isfinite(reference))
        np.testing.assert_allclose(bordered,reference,atol=1e-10,rtol=1e-9)
