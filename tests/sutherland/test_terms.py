"""
Test module that checks that the output from the Sutherland model does not change if you add terms with C = 0.

NOTE: The tests make comparrisons to the MieKinGas (legacy) implementation. If the tests fail, it may be due to 
        The Sutherland class and MieKinGas class using different regressed parameter sets for some correlated
        factor on the C++ side.
"""
from pykingas.Sutherland import S_MieKinGas, Sutherland
from pykingas.MieKinGas import MieKinGas
from scipy.constants import Avogadro, Boltzmann
import numpy as np
import pytest
from .s_tools import equal

T_lst = [100, 500, 1000]
p_lst = [1e2, 1e6, 5e7]
N_lst = [1, 2, 3]
l_lst = [1, 2, 3]
r_lst = [1, 2, 3]

@pytest.mark.parametrize('T', T_lst)
@pytest.mark.parametrize('p', p_lst)
@pytest.mark.parametrize('N', N_lst)
def test_rdf(T, p, N, silent=False):
    comps = 'H2,N2'
    mie = MieKinGas(comps)
    smie = S_MieKinGas(comps)
    C = smie.C
    lambdas = smie.lambdas
    expanded_C = np.zeros((2 + N, mie.ncomps, mie.ncomps))
    # Note lambda = 0 will give nan, even if C = 0, which i guess is fine, because it makes no sense to have a term with C = 0, except for testing purposes.
    expanded_lamb = np.zeros((2 + N, mie.ncomps, mie.ncomps)) + 10
    for i in range(2):
        expanded_C[i] = C[i]
        expanded_lamb[i] = lambdas[i]

    smie = Sutherland(smie.mole_weights, smie.sigma, smie.eps_div_k, expanded_C, expanded_lamb)

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
                print(f'g (Actual)     : {g_mie[i][j]}')
                print(f'g (With zeros) : {g_suth[i][j]}')
                print(f'Diff : {(g_mie[i][j] - g_suth[i][j]) / g_mie[i][j]}')
            assert equal(g_mie[i][j], g_suth[i][j], tol=1e-10)
    print('Test passed\n')

@pytest.mark.parametrize('T', T_lst)
@pytest.mark.parametrize('l', l_lst)
@pytest.mark.parametrize('r', r_lst)
@pytest.mark.parametrize('N', [2, 3])
def test_omega(T, l, r, N, silent=False):
    comps = 'H2,N2'
    mie = MieKinGas(comps)
    smie_short = S_MieKinGas(comps)
    C = smie_short.C
    lambdas = smie_short.lambdas
    expanded_C = np.zeros((2 + N, mie.ncomps, mie.ncomps))
    # Note lambda = 0 will give nan, even if C = 0, which i guess is fine, because it makes no sense to have a term with C = 0, except for testing purposes.
    expanded_lamb = np.zeros((2 + N, mie.ncomps, mie.ncomps)) + 10
    for i in range(2):
        expanded_C[i] = C[i]
        expanded_lamb[i] = lambdas[i]

    smie_long = Sutherland(smie_short.mole_weights * Avogadro * 1e3, smie_short.sigma, smie_short.eps_div_k, expanded_C, expanded_lamb)

    for i in range(mie.ncomps):
        for j in range(mie.ncomps):
            omega_short = smie_short.cpp_kingas.omega(i, j, l, r, T)
            omega_long = smie_long.cpp_kingas.omega(i, j, l, r, T)
            if silent is False:
                print(f'(i, j) = ({i}, {j})')
                print(f'omega (Actual)     : {omega_short}')
                print(f'omega (With zeros) : {omega_long}')
                print(f'Sigma eff : {smie_short.mole_weights}, {smie_long.mole_weights}')
                print(f'Diff : {(omega_short - omega_long) / omega_short}')
            assert equal(omega_short, omega_long, tol=1e-10)
    print('Test passed\n')

if __name__ == '__main__':
    for T_ in T_lst:
        for p_ in p_lst:
            for N_ in N_lst:
                test_rdf(T_, p_, N_)

    for T_ in T_lst:
        for l_ in l_lst:
            for r_ in r_lst:
                for N_ in N_lst:
                    test_omega(T_, l_, r_, N_)