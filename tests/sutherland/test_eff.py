from pykingas.Sutherland import S_MieKinGas
from pykingas.MieKinGas import MieKinGas
from scipy.constants import Avogadro, Boltzmann
import numpy as np
import pytest

@pytest.mark.parametrize('comps', ['H2', 'AR,C1', 'KR,CO2,O2'])
@pytest.mark.parametrize('i', [0, 1, 2])
@pytest.mark.parametrize('j', [0, 1, 2])
def test_mie_sigma_eff(comps, i, j, silent=False):
    mie = MieKinGas(comps)
    smie = S_MieKinGas(comps)
    if (i >= mie.ncomps) or (j >= mie.ncomps):
        return

    sigma = mie.sigma
    sigma_eff = smie.get_sigma_eff()
    if silent is False:
        print(f'Sigma (mie) : {sigma}')
        print(f'Sigma (eff) : {sigma_eff}')
        print(abs(sigma[i][j] - sigma_eff[i][j]) / sigma[i][j])
    assert abs(sigma[i][j] - sigma_eff[i][j]) < 1e-10

@pytest.mark.parametrize('comps', ['H2', 'AR,C1', 'KR,CO2,O2'])
@pytest.mark.parametrize('i', [0, 1, 2])
@pytest.mark.parametrize('j', [0, 1, 2])
def test_mie_sigma_min(comps, i, j, silent=False):
    mie = MieKinGas(comps)
    smie = S_MieKinGas(comps)
    if (i >= mie.ncomps) or (j >= mie.ncomps):
        return
    
    sigma = np.array(mie.sigma)
    sigma_min_mie = sigma * (mie.lr / mie.la)**(1 / (mie.lr - mie.la))
    sigma_min_suth = smie.get_sigma_min()
    if silent is False:
        print(f'Sigma (true) : {sigma_min_mie}')
        print(f'Sigma (suth) : {sigma_min_suth}')
        print(abs(sigma_min_mie[i][j] - sigma_min_suth[i][j]) / sigma_min_mie[i][j])
    assert abs(sigma_min_mie[i][j] - sigma_min_suth[i][j]) < 1e-10


@pytest.mark.parametrize('comps', ['H2', 'AR,C1', 'KR,CO2,O2'])
@pytest.mark.parametrize('i', [0, 1, 2])
@pytest.mark.parametrize('j', [0, 1, 2])
def test_mie_epsilon_eff(comps, i, j, silent=False):
    mie = MieKinGas(comps)
    smie = S_MieKinGas(comps)
    if (i >= mie.ncomps) or (j >= mie.ncomps):
        return

    epsilon = mie.epsilon_ij
    epsilon_eff = smie.get_epsilon_eff()
    if silent is False:
        print(f'epsilon (mie) : {epsilon}')
        print(f'epsilon (eff) : {epsilon_eff}')
        print(abs(epsilon[i][j] - epsilon_eff[i][j]) / epsilon[i][j])
    assert abs(epsilon[i][j] - epsilon_eff[i][j]) < 1e-10

if __name__ == '__main__':
    complist = ['H2', 'AR,C1', 'KR,CO2,O2']
    for comps in complist:
        for i in range(3):
            for j in range(3):
                print(comps, i, j)
                test_mie_epsilon_eff(comps, i, j)