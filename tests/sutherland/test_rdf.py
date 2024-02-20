"""
Test module that ensures that the new Sutherland model reproduces the RDF at contact computed using the old MieKinGas model
"""
from pykingas.Sutherland import S_MieKinGas
from pykingas.MieKinGas import MieKinGas
from scipy.constants import Avogadro, Boltzmann
import numpy as np
import pytest
from tools import equal

complist = ['H2', 'AR,C1', 'KR,CO2,O2']
idx_list = [0, 1, 2, 3]
T_lst = [100, 500, 1000]
p_lst = [1e2, 1e6, 5e7]

@pytest.mark.parametrize('comps', complist)
@pytest.mark.parametrize('T', T_lst)
@pytest.mark.parametrize('p', p_lst)
def test_mie_g0(comps, T, p, silent=False):
    mie = MieKinGas(comps)
    smie = S_MieKinGas(comps)
    
    z = np.array([i + 1 for i in range(mie.ncomps)], dtype=float)
    z /= sum(z)
    Vm, = mie.eos.specific_volume(T, p, z, mie.eos.VAPPH)
    rho = Avogadro / Vm 
    
    g0_mie = mie.cpp_kingas.rdf_g0(rho, T, z)
    g0_suth = smie.cpp_kingas.rdf_g0(rho, T, z)
    for i in range(mie.ncomps):
        for j in range(mie.ncomps):
            if silent is False:
                print(f'(i, j) = ({i}, {j})')
                print(f'g0 (mie) : {g0_mie[i][j]}')
                print(f'g0 (sut) : {g0_suth[i][j]}')
                print(f'Diff : {(g0_mie[i][j] - g0_suth[i][j]) / g0_mie[i][j]}')
            assert equal(g0_mie[i][j], g0_suth[i][j])

@pytest.mark.parametrize('comps', complist)
@pytest.mark.parametrize('T', T_lst)
@pytest.mark.parametrize('p', p_lst)
def test_mie_g1(comps, T, p, silent=False):
    mie = MieKinGas(comps)
    smie = S_MieKinGas(comps)
    
    z = np.array([i + 1 for i in range(mie.ncomps)], dtype=float)
    z /= sum(z)
    Vm, = mie.eos.specific_volume(T, p, z, mie.eos.VAPPH)
    rho = Avogadro / Vm 
    
    g1_mie = mie.cpp_kingas.rdf_g1(rho, T, z)
    g1_suth = smie.cpp_kingas.rdf_g1(rho, T, z)
    for i in range(mie.ncomps):
        for j in range(mie.ncomps):
            if silent is False:
                print(f'(i, j) = ({i}, {j})')
                print(f'g1 (mie) : {g1_mie[i][j]}')
                print(f'g1 (sut) : {g1_suth[i][j]}')
                print(f'Diff : {(g1_mie[i][j] - g1_suth[i][j]) / g1_mie[i][j]}')
            assert equal(g1_mie[i][j], g1_suth[i][j], tol=1e-8)
    print('Test passed\n')


@pytest.mark.parametrize('comps', complist)
@pytest.mark.parametrize('T', T_lst)
@pytest.mark.parametrize('p', p_lst)
def test_mie_g2(comps, T, p, silent=False):
    mie = MieKinGas(comps)
    smie = S_MieKinGas(comps)

    z = np.array([i + 1 for i in range(mie.ncomps)], dtype=float)
    z /= sum(z)
    Vm, = mie.eos.specific_volume(T, p, z, mie.eos.VAPPH)
    rho = Avogadro / Vm

    g2_mie = mie.cpp_kingas.rdf_g2(rho, T, z)
    g2_suth = smie.cpp_kingas.rdf_g2(rho, T, z)
    for i in range(mie.ncomps):
        for j in range(mie.ncomps):
            if silent is False:
                print(f'(i, j) = ({i}, {j})')
                print(f'g2 (mie) : {g2_mie[i][j]}')
                print(f'g2 (sut) : {g2_suth[i][j]}')
                print(f'Diff : {(g2_mie[i][j] - g2_suth[i][j]) / g2_mie[i][j]}')
            assert equal(g2_mie[i][j], g2_suth[i][j], tol=1e-1)
    print('Test passed\n')

@pytest.mark.parametrize('comps', complist)
@pytest.mark.parametrize('T', T_lst)
@pytest.mark.parametrize('p', p_lst)
def test_mie_g(comps, T, p, silent=False):
    mie = MieKinGas(comps)
    smie = S_MieKinGas(comps)

    z = np.array([i + 1 for i in range(mie.ncomps)], dtype=float)
    z /= sum(z)
    Vm, = mie.eos.specific_volume(T, p, z, mie.eos.VAPPH)
    rho = Avogadro / Vm

    g_mie = mie.get_rdf(rho, T, z)
    g_suth = smie.get_rdf(rho, T, z)
    for i in range(mie.ncomps):
        for j in range(mie.ncomps):
            if silent is False:
                print(f'(i, j) = ({i}, {j})')
                print(f'g (mie) : {g_mie[i][j]}')
                print(f'g (sut) : {g_suth[i][j]}')
                print(f'Diff : {(g_mie[i][j] - g_suth[i][j]) / g_mie[i][j]}')
            assert equal(g_mie[i][j], g_suth[i][j], tol=1e-10)
    print('Test passed\n')


if __name__ == '__main__':
    for comps_ in complist:
        for T_ in T_lst:
            for p_ in p_lst:
                test_mie_g2(comps_, T_, p_, silent=False)
