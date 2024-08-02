import pytest
from scipy.constants import Avogadro, gas_constant
from pykingas.HardSphere import HardSphere
from pykingas.MieKinGas import MieKinGas
import numpy as np
from tools import models, check_eq

def get_Jm(model):
    M = model.mole_weights * Avogadro / 1e3
    Jm = np.array([1., 2., 3., 0.])
    Jm[-1] = - sum(Jm[:-1] * M[:-1]) / M[-1]
    return Jm

x1_lst = np.arange(0.1, 0.9, 0.2)
x_lst = []
for x1 in x1_lst:
    x2_lst = np.arange(0.1, 1 - x1 - 0.05, 0.2)
    for x2 in x2_lst:
        x3_lst = np.arange(0.1, 1 - x1 - x2 - 0.05, 0.2)
        for x3 in x3_lst:
            x4 = 1 - (x1 + x2 + x3)
            x_lst.append([x1, x2, x3, x4])


@pytest.mark.parametrize('x', x_lst)
@pytest.mark.parametrize('model', models)
def test_com_2_con(model, x):
    model = model('AR,KR,N2,O2')
    Jm = get_Jm(model)
    psi = model.get_com_2_con_matr(x)
    Jn = psi @ Jm
    assert check_eq(sum(Jn), 0.0)

@pytest.mark.parametrize('x', x_lst)
@pytest.mark.parametrize('solvent_idx', [0, 1, 2, 3])
@pytest.mark.parametrize('model', models)
def test_com_2_solvent(model, x, solvent_idx):
    model = model('AR,KR,N2,O2')
    Jm = get_Jm(model)
    psi = model.get_com_2_solv_matr(x, solvent_idx)
    Js = psi @ Jm
    assert check_eq(Js[solvent_idx], 0.0)

@pytest.mark.parametrize('x', x_lst)
@pytest.mark.parametrize('p', np.linspace(1, 100, 10) * 1e5)
@pytest.mark.parametrize('model', models)
def test_com_2_cov(model, x, p):
    model = model('AR,KR,N2,O2')
    Jm = get_Jm(model)
    T = 300
    Vm, dvdn = model.eos.specific_volume(T, p, x, model.eos.VAPPH, dvdn=True)
    psi = model.get_com_2_cov_matr(T, Vm, x)
    Jv = psi @ Jm
    assert check_eq(sum(Jv * dvdn), 0.0)
