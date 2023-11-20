import numpy as np
from pykingas.MieKinGas import MieKinGas
from tools import models, check_eq
import pytest

@pytest.mark.parametrize('model', models)
@pytest.mark.parametrize('N', [2, 3])
def test_th_diffusion_ratio(model, N, silent=False):
    T = 300
    Vm = 1 / 10
    x = [0.1, 0.3, 0.6]
    comps = 'H2,C1,CO2'

    kin = model(comps, is_idealgas=True)

    kT = kin.thermal_diffusion_ratio(T, Vm, x, N=N)
    assert check_eq(sum(x * kT), 1.0)
    if silent is False:
        print(f'k_T {model}, N = {N} : {kT}, sum(x * kT) = {sum(x * kT)}')

@pytest.mark.parametrize('model', models)
@pytest.mark.parametrize('N', [2, 3, 4])
@pytest.mark.parametrize('FoR', ['CoN', 'CoM', 'solvent'])
def test_th_diffusion_coeff(model, N, FoR, silent=True):
    comps = 'HE,NE,O2'
    T = 300
    Vm = 1 / 10
    x = [0.1, 0.3, 0.6]

    kin = model(comps)
    zerosum_scale = {'CoN' : np.ones(kin.ncomps), 'CoM' : kin.mole_weights, 'solvent' : np.array([0, 0, 1])}
    DT = kin.thermal_diffusion_coeff(T, Vm, x, N=N, frame_of_reference=FoR, solvent_idx=2)
    assert check_eq(sum(DT * zerosum_scale[FoR]) / max(abs(DT)), 0.0)
    if silent is False:
        print(f'k_T {model}, N = {N} : {DT}')

@pytest.mark.parametrize('model', models)
@pytest.mark.parametrize('N', [1, 3])
@pytest.mark.parametrize('FoR', ['CoN', 'CoM', 'solvent'])
def test_diffusion(model, N, FoR, silent=True):
    comps = 'HE,NE,O2'
    T = 300
    Vm = 1 / 10
    x = [0.1, 0.3, 0.6]

    kin = model(comps)
    zerosum_scale = {'CoN': np.ones(kin.ncomps), 'CoM': kin.mole_weights, 'solvent': np.array([0, 0, 1])}
    D = kin.interdiffusion(T, Vm, x, N=N, frame_of_reference=FoR, solvent_idx=2)
    for i in range(3):
        assert check_eq(sum(D[:, i] * zerosum_scale[FoR]) / max(abs(D.flatten())), 0.0)
    if silent is False:
        print(f'k_T {model}, N = {N} : {D}')

@pytest.mark.parametrize('model', models)
@pytest.mark.parametrize('N', [2, 3])
@pytest.mark.parametrize('solvent_idx', [0, 1, 2])
def test_DT_solvent_for(model, N, solvent_idx, silent=True):
    comps = 'HE,NE,O2'
    T = 300
    Vm = 1 / 10
    x = [0.1, 0.3, 0.6]

    kin = model(comps)
    DT = kin.thermal_diffusion_coeff(T, Vm, x, N=N, frame_of_reference='solvent', solvent_idx=solvent_idx)
    assert check_eq(DT[solvent_idx] / max(abs(DT)), 0.0)
    if silent is False:
        print(f'k_T {model}, N = {N} : {DT}')