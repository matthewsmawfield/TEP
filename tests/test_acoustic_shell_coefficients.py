"""Independent algebra and finite-difference checks of the acoustic flux Hessian."""
import importlib.util
from pathlib import Path
import numpy as np

source=Path(__file__).resolve().parents[1]/'scripts/steps/step_81_action_contract_audit.py'
spec=importlib.util.spec_from_file_location('shell_audit',source)
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)

def test_native_secant_and_prescription_are_not_conflated():
    r=audit.derive()['acoustic_shell_reduction'];y=r['y_s_sound_horizon']
    np.testing.assert_allclose(r['native_secant_value'],1/y,rtol=1e-12)
    np.testing.assert_allclose(r['K_sp_spatial_sound_shell_derived']/r['native_secant_value'],np.sqrt(2),rtol=1e-12)
    assert not r['quadratic_perturbation_mapping_established']

def test_spatial_hessian_matches_vector_current_derivative():
    r=audit.derive()['acoustic_shell_reduction'];p=r['q_s_sound_horizon']*r['y_s_sound_horizon'];k=16.03
    v=np.array([p,0.,0.]);step=1e-7
    def current(v):
        n=np.linalg.norm(v)
        return (k*n/np.sqrt(2)+n*n)*v
    jac=np.column_stack([(current(v+step*e)-current(v-step*e))/(2*step) for e in np.eye(3)])
    expected=np.diag([r['native_longitudinal_value'],r['native_secant_value'],r['native_secant_value']])
    np.testing.assert_allclose(jac,expected,rtol=1e-9,atol=1e-11)
