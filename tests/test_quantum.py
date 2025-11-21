from pykingas.multiparam import init_multiparam, init_multiparam_FH
from pykingas.multiparam import Patowski, PatowskiFH, ModTangToennies, FH_ModTangToennies
import pytest

@pytest.mark.parametrize('comp', ['AR', 'P-H2', 'NE', 'HE'])
def test_virial(comp):
    q_model = init_multiparam(comp)
    q_model.set_quantum_active(True)

    c_model = init_multiparam(comp)
    c_model.set_quantum_active(False)

    fh_model = init_multiparam_FH(comp, 1)

    T_map = {'AR' : 100, 'P-H2' : 50, 'NE' : 70, 'HE' : 50}
    T = T_map[comp]

    B_q = q_model.second_virial(0, 0, T)
    B_c = c_model.second_virial(0, 0, T)
    B_fh = fh_model.second_virial(0, 0, T)

    print(B_q, B_c, B_fh)
    print(abs((B_fh - B_q) / B_q), abs((B_c - B_q) / B_q))

    assert abs(B_c - B_q) > abs(B_fh - B_q)
    assert abs((B_fh - B_q) / B_q) < 5e-2

@pytest.mark.parametrize('comp', ['AR', 'P-H2', 'NE', 'HE'])
def test_viscosity(comp):
    q_model = init_multiparam(comp)
    q_model.set_quantum_active(True)

    c_model = init_multiparam(comp)
    c_model.set_quantum_active(False)

    T_map = {'AR' : 100, 'P-H2' : 50, 'NE' : 70, 'HE' : 50}
    T = T_map[comp]

    visc_q = q_model.viscosity(T, 1, [1], N=1)
    visc_c = c_model.viscosity(T, 1, [1], N=1)

    err_lims = {'AR': 1e-2, 'P-H2': 11e-2, 'NE' : 2e-2, 'HE' : 2e-2}

    print(comp, visc_c * 1e6, visc_q * 1e6, (visc_c - visc_q) / visc_q)
    assert abs((visc_c - visc_q) / visc_q) < err_lims[comp]

