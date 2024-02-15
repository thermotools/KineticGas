from pykingas.Sutherland import S_MieKinGas
from pykingas.MieKinGas import MieKinGas
from scipy.constants import Avogadro, Boltzmann
import numpy as np
import pytest
from tools import equal

complist = ['H2', 'AR,C1', 'KR,CO2,O2']
idx_list = [1, 2, 3]
T_lst = [100, 500, 1000]
p_lst = [1e2, 1e6, 5e7]

@pytest.mark.parametrize('comps', complist)
@pytest.mark.parametrize('i', idx_list)
@pytest.mark.parametrize('j', idx_list)
@pytest.mark.parametrize('T', T_lst)
@pytest.mark.parametrize('p', p_lst)
def test_mie_g0(comps, i, j, T, p, silent=False):
    mie = MieKinGas(comps)
    smie = S_MieKinGas(comps)
    if (i >= mie.ncomps) or (j >= mie.ncomps):
        return
    
    z = np.array([i + 1 for i in range(mie.ncomps)], dtype=float)
    z /= sum(z)
    Vm, = mie.eos.specific_volume(T, p, z, mie.eos.VAPPH)
    rho = Avogadro / Vm 
    
    g0_mie = mie.cpp_kingas.rdf_g0(rho, T, z)
    g0_suth = smie.cpp_kingas.rdf_g0(rho, T, z)
    if silent is False:
        print(f'(i, j) = ({i}, {j})')
        print(f'g0 (mie) : {g0_mie[i][j]}')
        print(f'g0 (sut) : {g0_suth[i][j]}')
        print(f'Diff : {(g0_mie[i][j] - g0_suth[i][j]) / g0_mie[i][j]}')
    assert equal(g0_mie[i][j], g0_suth[i][j])

@pytest.mark.parametrize('comps', complist)
@pytest.mark.parametrize('i', idx_list)
@pytest.mark.parametrize('j', idx_list)
@pytest.mark.parametrize('T', T_lst)
@pytest.mark.parametrize('p', p_lst)
def test_mie_g1(comps, i, j, T, p, silent=False):
    mie = MieKinGas(comps)
    smie = S_MieKinGas(comps)
    if (i >= mie.ncomps) or (j >= mie.ncomps):
        return
    
    z = np.array([i + 1 for i in range(mie.ncomps)], dtype=float)
    z /= sum(z)
    Vm, = mie.eos.specific_volume(T, p, z, mie.eos.VAPPH)
    rho = Avogadro / Vm 
    
    g1_mie = mie.cpp_kingas.rdf_g1(rho, T, z)
    g1_suth = smie.cpp_kingas.rdf_g1(rho, T, z)
    if silent is False:
        print(f'(i, j) = ({i}, {j})')
        print(f'g1 (mie) : {g1_mie[i][j]}')
        print(f'g1 (sut) : {g1_suth[i][j]}')
        print(f'Diff : {(g1_mie[i][j] - g1_suth[i][j]) / g1_mie[i][j]}')
    assert equal(g1_mie[i][j], g1_suth[i][j], tol=1e-8)
    print('Test passed\n')

@pytest.mark.parametrize('comps', complist)
@pytest.mark.parametrize('i', idx_list)
@pytest.mark.parametrize('j', idx_list)
@pytest.mark.parametrize('T', T_lst)
@pytest.mark.parametrize('p', p_lst)
def test_mie_d2(comps, i, j, T, p, silent=False):
    mie = MieKinGas(comps)
    smie = S_MieKinGas(comps)
    if (i >= mie.ncomps) or (j >= mie.ncomps):
        return
    
    z = np.array([i + 1 for i in range(mie.ncomps)], dtype=float)
    z /= sum(z)
    Vm, = mie.eos.specific_volume(T, p, z, mie.eos.VAPPH)
    rho = Avogadro / Vm 
    
    lambdas = smie.lambdas

    d_BH = mie.cpp_kingas.get_BH_diameters(T)
    mie_2d = np.array(mie.cpp_kingas.B_func(rho, z, d_BH, lambdas[0]))
    suth_2d = np.array(smie.cpp_kingas.B_func(rho, z, d_BH, lambdas[0]))
    # print(g1_mie)
    # print(g1_suth)
    if silent is False:
        print(T, p / 1e5) # , max(g1_mie.flatten()), max(g1_suth.flatten()))
        # print(f'(i, j) = ({i}, {j})')
        # print(f'g1 (mie) : {g1_mie}')
        # print(f'g1 (sut) : {g1_suth}')
        print(f'Diff : {(mie_2d - suth_2d) / mie_2d}')
    assert equal(mie_2d[i][j], suth_2d[i][j])

@pytest.mark.parametrize('comps', complist)
@pytest.mark.parametrize('i', idx_list)
@pytest.mark.parametrize('j', idx_list)
@pytest.mark.parametrize('T', T_lst)
@pytest.mark.parametrize('p', p_lst)
def test_mie_d0(comps, i, j, T, p, silent=False):
    mie = MieKinGas(comps)
    smie = S_MieKinGas(comps)
    if (i >= mie.ncomps) or (j >= mie.ncomps):
        return
    
    z = np.array([i + 1 for i in range(mie.ncomps)], dtype=float)
    z /= sum(z)
    Vm, = mie.eos.specific_volume(T, p, z, mie.eos.VAPPH)
    rho = Avogadro / Vm 
    
    lambdas = smie.lambdas

    d_BH = mie.cpp_kingas.get_BH_diameters(T)
    mie_val = mie.cpp_kingas.zeta_x(rho, z, d_BH)
    suth_val = smie.cpp_kingas.zeta_x(rho, z, d_BH)
    # print(g1_mie)
    # print(g1_suth)
    if silent is False:
        print(T, p / 1e5) # , max(g1_mie.flatten()), max(g1_suth.flatten()))
        # print(f'(i, j) = ({i}, {j})')
        # print(f'g1 (mie) : {g1_mie}')
        # print(f'g1 (sut) : {g1_suth}')
        print(f'Diff : {(mie_val - suth_val) / mie_val}')
    assert equal(mie_val, suth_val)

if __name__ == '__main__':
    for comps in complist:
        for T in T_lst:
            for p in p_lst:
                for i in idx_list:
                    for j in idx_list:
                        test_mie_g1(comps, i, j, T, p, silent=False)
