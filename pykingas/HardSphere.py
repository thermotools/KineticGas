'''
Author: Vegard Gjeldvik Jervell
Purpose: Wrapper for the HardSphere class.
'''

import numpy as np
from numpy import pi
from scipy.constants import Avogadro, Boltzmann as kB
from scipy.optimize import root
from .libpykingas import cpp_HardSphere
from pykingas.py_KineticGas import py_KineticGas, IdealGas
import warnings

class HardSphereEoS:

    def __init__(self, cpp_model, sigma):
        self.ncomps = len(sigma)
        self.VAPPH = 1
        self.cpp_model = cpp_model
        self.sigma = sigma

    def pressure_tv(self, T, V, n):
        x = n / sum(n)
        rho = Avogadro / V
        chi = self.cpp_model.get_rdf(rho, T, x)
        return HS_pressure(rho, T, x, self.sigma, chi),

    def specific_volume(self, T, p, n, phase, dvdn=False):
        v_init = Avogadro * kB * T / p
        v = root(lambda Vm : self.pressure_tv(T, Vm[0], n)[0] - p, x0=np.array([v_init])).x[0]
        if dvdn is False:
            return v,

        dvdn = np.empty(self.ncomps)
        eps = min(n) * 1e-3
        for i in range(self.ncomps):
            nm1 = np.array([xi for xi in n])
            np1 = np.array([xi for xi in n])
            nm1[i] -= eps
            np1[i] += eps

            vm1 = self.specific_volume(T, p, nm1, 2)[0] * sum(nm1)
            vp1 = self.specific_volume(T, p, np1, 2)[0] * sum(np1)
            dvdn[i] =  (vp1 - vm1) / (2 * eps)

        return v, dvdn

    def chemical_potential_tv(self, T, V, n, dmudn=False):
        n = np.array(n)
        x = n / sum(n)
        Vm = V / sum(n)
        rho = Avogadro / Vm
        chi = self.cpp_model.get_rdf(rho, T, x)
        mu = mu_func(rho, T, x, self.sigma, chi) * Avogadro
        return_tuple = (mu, )
        if dmudn is True:
            dn = min(n) / 1e3
            dmudn = np.empty((self.ncomps, self.ncomps))
            for i in range(self.ncomps):
                dn_arr = np.zeros(self.ncomps)
                dn_arr[i] = dn
                mu_p, = self.chemical_potential_tv(T, V, n + dn_arr)
                mu_m, = self.chemical_potential_tv(T, V, n - dn_arr)
                dmudn[i] = (mu_p - mu_m) / (2 * dn)
            return_tuple += (dmudn, )
        return return_tuple

    def idealenthalpysingle(self, T, i, dhdt=None):
        warnings.warn('Idealenthalpy gives dummy values for HardSphere!', RuntimeWarning)
        if dhdt is None:
            return 0.,
        return 0., 5 * 8.314 / 2

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


class HardSphere(py_KineticGas):

    def __init__(self, comps, mole_weights=None, sigma=None,
                 N=4, is_idealgas=False, parameter_ref='default'):
        """Constructor
        If parameters are explicitly supplied through optional arguments, these will be used instead of those in the database.
        To supply specific parameters for only some components, give `None` for the components that should use the database
        value
        &&
        Args:
            comps (str) : Comma-separated list of components
            mole_weights (1D array) : Molar weights [g/mol]
            sigma (1D array) : hard-sphere diameters [m]
            parameter_ref (str) : Id for parameter set to use
        """
        super().__init__(comps, mole_weights=mole_weights, N=N, is_idealgas=is_idealgas)

        self.fluids = [self.fluids[i]['HardSphere'][parameter_ref] for i in range(self.ncomps)]
        if sigma is None:
            sigma = np.array([self.fluids[i]['sigma'] for i in range(self.ncomps)])
        elif None in sigma:
            for i in range(self.ncomps):
                if sigma[i] is None:
                    sigma[i] = self.fluids[i]['sigma']
        elif self._is_singlecomp is True:
            sigma = np.array([sigma[0] for _ in range(2)])
        else:
            sigma = np.array(sigma)
        
        self.sigma = 0.5 * (np.vstack(tuple(sigma for _ in range(self.ncomps))) 
                                + np.vstack(tuple(sigma for _ in range(self.ncomps))).transpose())

        self.cpp_kingas = cpp_HardSphere(self.mole_weights, self.sigma, is_idealgas, self._is_singlecomp)
        if self.is_idealgas is True:
            self.eos = IdealGas(comps)
        else:
            self.eos = HardSphereEoS(self.cpp_kingas, self.sigma)

