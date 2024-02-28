import matplotlib.pyplot as plt
import numpy as np
from plottools.cyclers import NormedCmap
from pykingas.QuantumMie import QuantumMie
from scipy.constants import Avogadro
import pytest

binary_fh = [[i, j] for i in range(3) for j in range(3)]

@pytest.mark.parametrize('fh_order', [0, 1, 2])
def test_singlecomp_rdf(fh_order):
    tolerances = {0 : 3.5e-8, 1 : 0.02, 2 : 0.07}
    qmie = QuantumMie('H2', FH_orders=fh_order)
    T_lst = [35, 60, 90]
    p_lst = np.linspace(1e3, 1e7, 50)
    # p_lst = np.linspace(1e3, 1.1e3, 2)

    cmap = NormedCmap('cool', T_lst)
    for T in T_lst:
        Vm_lst = np.array([qmie.eos.specific_volume(T, p, [1], qmie.eos.MINGIBBSPH)[0] for p in p_lst])
        rho_lst = Avogadro / Vm_lst

        mie = qmie.get_effective_mie_model(T)
        s = mie.sigma[0][0]
        for i in range(2):
            for j in range(2):
                qmie_rdf = np.array([qmie.get_rdf(rho, T, [0.5, 0.5])[i][j] for rho in rho_lst])
                mie_rdf = np.array([mie.get_rdf(rho, T, [0.5, 0.5])[i][j] for rho in rho_lst])
                max_err = max(abs(qmie_rdf - mie_rdf) / mie_rdf)

                # plt.plot(rho_lst * s**3, qmie_rdf, color=cmap(T))
                # plt.plot(rho_lst * s**3, mie_rdf, linestyle='', color=cmap(T), marker='|')

                # print(T, max_err)
                assert max_err < tolerances[fh_order]
    # plt.show()

@pytest.mark.parametrize('fh_orders', binary_fh)
def test_binary_rdf(fh_orders):
    tolerances = {  (0, 0) : 0.16 ,
                    (0, 1) : 0.12 ,
                    (0, 2) : 0.02 ,
                    (1, 0) : 0.19 ,
                    (1, 1) : 0.15 ,
                    (1, 2) : 0.05 ,
                    (2, 0) : 0.26 ,
                    (2, 1) : 0.22 ,
                    (2, 2) : 0.12}
    qmie = QuantumMie('H2,HE', FH_orders=fh_orders)
    print(qmie.get_fh_order())
    T_lst = [35, 60, 90]
    p_lst = np.linspace(1e3, 1e7, 50)
    # p_lst = np.linspace(1e3, 1.1e3, 2)

    cmap = NormedCmap('cool', T_lst)
    for T in T_lst:
        Vm_lst = np.array([qmie.eos.specific_volume(T, p, [1], qmie.eos.MINGIBBSPH)[0] for p in p_lst])
        rho_lst = Avogadro / Vm_lst

        mie = qmie.get_effective_mie_model(T, fit_la=True)
        s = mie.sigma[0][1]

        max_err = 0
        for i in range(2):
            for j in range(2):
                qmie_rdf = np.array([qmie.get_rdf(rho, T, [0.5, 0.5])[i][j] for rho in rho_lst])
                mie_rdf = np.array([mie.get_rdf(rho, T, [0.5, 0.5])[i][j] for rho in rho_lst])
                max_err = max(max_err, max(abs(qmie_rdf - mie_rdf) / mie_rdf))

                # print(qmie_rdf)
                # plt.plot(rho_lst * s ** 3, qmie_rdf, color=cmap(T))
                # plt.plot(rho_lst * s ** 3, mie_rdf, linestyle='', color=cmap(T), marker='|')
                # plt.title(f'FH : {fh_orders}, ({i}, {j})')
                # plt.show()
                assert max_err < tolerances[tuple(fh_orders)]
                max_err = 0

@pytest.mark.parametrize('fh_order', [0, 1, 2])
def test_singlecomp_omega(fh_order):
    tolerances = {0 : 0.035, 1 : 0.03, 2 : 0.025}
    qmie = QuantumMie('H2', FH_orders=fh_order)
    T_lst = [35, 60, 90]

    l_vals = (1, 2, 3)
    r_vals = (1, 2, 3)

    for T in T_lst:
        mie = qmie.get_effective_mie_model(T)
        for l in l_vals:
            for r in r_vals:
                qmie_omega = qmie.cpp_kingas.omega(0, 0, l, r, T)
                mie_omega = mie.cpp_kingas.omega(0, 0, l, r, T)

                assert abs((mie_omega - qmie_omega) / mie_omega) < tolerances[fh_order]

if __name__ == "__main__":
    for fh_order in (0, 1, 2):
        test_singlecomp_omega(fh_order)
