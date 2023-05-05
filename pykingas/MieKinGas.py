'''
Author: Vegard Gjeldvik Jervell
Purpose: Wrapper for the MieKinGas class. Handles mixing of sigma and epsilon parameters before initialising cpp-side.
'''
import numpy as np
from scipy.constants import Boltzmann, Avogadro
from scipy.integrate import quad
from pykingas import cpp_MieKinGas, MieType
from thermopack.saftvrmie import saftvrmie

class MieKinGas(MieType.MieType):

    def __init__(self, comps,
                 mole_weights=None, sigma=None, eps_div_k=None,
                 la=None, lr=None, lij=0, kij=0,
                 N=4, is_idealgas=False, use_eos=None,
                 parameter_ref='default'):
        '''
        :param comps (str): Comma-separated list of components

        If parameters are explicitly supplied, these will be used instead of those in the database
        :param mole_weights : (1D array) Molar weights [g/mol]
        :param sigma : (1D array) hard-sphere diameters [m]
        :param eps_div_k : (1D array) epsilon parameter / Boltzmann constant [-]
        :param la, lr : (1D array) attractive and repulsive exponent of the pure components [-]
        :param lij : (float) Mixing parameter for sigma (lij > 0 => smaller sigma_12, lij < 0 => larger sigma_12)
        :param kij : (float) Mixing parameter for epsilon (kij > 0 => favours mixing, kij < 0 => favours separation)
        :param use_eos : (thermopack eos object, optional) EoS to use (initialized), defaults to saftvrmie
        :param use_db (bool) : Use precomputed database values for omega_integrals if available
        '''
        super().__init__(comps, 'Mie',
                    mole_weights=mole_weights, sigma=sigma,
                    eps_div_k=eps_div_k, la=la, lr=lr, lij=lij, kij=kij,
                    N=N, is_idealgas=is_idealgas,
                    parameter_ref=parameter_ref)

        self.cpp_kingas = cpp_MieKinGas(self.mole_weights, self.sigma_ij, self.epsilon_ij, self.la, self.lr, self.is_idealgas)
        if use_eos is None:
            self.eos = saftvrmie()
            if parameter_ref == 'default':
                self.eos.init(comps)
            else:
                self.eos.init(comps, parameter_reference=parameter_ref)
        else:
            self.eos = use_eos
