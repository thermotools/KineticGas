'''
Author: Vegard Gjeldvik Jervell
Purpose: Wrapper for the PseudoHardSphere class.
'''

from pykingas import cpp_PseudoHardSphere
from pykingas.py_KineticGas import py_KineticGas
import numpy as np
from scipy.constants import Boltzmann as kB, pi, Avogadro

def HS_pressure(rho, T, x, sigma, chi):
    p = rho * kB * T
    for i in range(len(x)):
        for j in range(len(x)):
            p += 2 * pi * rho**2 * x[i] * x[j] * kB * T * sigma[i][j]**3 * chi[i][j]
    return p

def Z_func(rho, x, sigma):
    n = rho * x
    Z = np.array([sum(n * np.diag(sigma)**i) for i in range(1, 5)]) * pi / 6
    Z[-1] = 1 - Z[-2]
    return Z

def mu_func(rho, T, x, sigma, chi):
    Z1, Z2, Z3, Z = Z_func(rho, x, sigma)
    n = rho * x
    p = HS_pressure(rho, T, x, sigma, chi)

    if Z3 >= (1 - 1e-6):
        return np.full_like(x, np.nan)
    mu = np.empty_like(x)
    for i in range(len(x)):
        mu[i] = kB * T * (np.log(n[i]) - np.log(1 - Z3)
                            + (pi * sigma[i][i]**3 * p / (6 * kB * T)) \
                            + (3 * (Z2 * sigma[i][i] + Z1 * sigma[i][i]**2) / Z)
                            + ((9 / 2) * (Z2 * sigma[i][i] / Z)**2) \
                            + 3 * (Z2 * sigma[i][i] / Z3)**2 * (np.log(Z) + (Z3 / Z) - (Z3**2 / (2 * Z**2))) \
                            - (Z2 * sigma[i][i] / Z3)**3 * (2 * np.log(Z) + Z3 * (2 - Z3) / Z))
    return mu

class PseudoHardSphere(py_KineticGas):

    def __init__(self, comps, mole_weights=None, sigma=None, N=3, is_idealgas=False,
                    parameter_ref='default'):

        super().__init__(comps, mole_weights=mole_weights, N=N, is_idealgas=is_idealgas)
        self.fluids = [self.fluids[i]['HardSphere'][parameter_ref] for i in range(self.ncomps)]
        if sigma == None:
            sigma = np.array([self.fluids[i]['sigma'] for i in range(self.ncomps)])
        else:
            sigma = np.array(sigma)

        self.sigma_ij =  0.5 * (np.vstack((sigma, sigma)) + np.vstack((sigma, sigma)).transpose())
        self.cpp_kingas = cpp_PseudoHardSphere(self.mole_weights, self.sigma_ij, is_idealgas)

    def get_Eij(self, Vm, T, x):
        x = np.array(x)
        rho = Avogadro / Vm
        n = rho * x

        E = np.empty((self.ncomps, self.ncomps))
        for j in range(self.ncomps):
            dn = np.zeros(self.ncomps)
            dn[j] = 1e-2 * n[j]
            drho = dn[j]
            x_1 = (n - dn / 2) / sum(n - dn / 2)
            x1 = (n + dn / 2) / sum(n + dn / 2)
            chi_1 = self.cpp_kingas.get_rdf(rho - drho / 2, T, x_1)
            chi1 = self.cpp_kingas.get_rdf(rho + drho / 2, T, x1)
            mu_1 = mu_func(rho - drho / 2, T, x_1, self.sigma, chi_1)
            mu1 = mu_func(rho + drho / 2, T, x1, self.sigma, chi1)

            for i in range(self.ncomps):
                E[i][j] = (n[i] / (kB * T)) * (mu1[i] - mu_1[i]) / dn[j]

        return E