'''
Author: Johannes Salomonsen Løken
Purpose: Wrapper for the LJ/Spline class.
'''

import numpy as np
from numpy import pi
from scipy.constants import Avogadro, Boltzmann as kB
from scipy.optimize import root
from .libpykingas import cpp_LJSpline
from pykingas.py_KineticGas import py_KineticGas

class LJSpline(py_KineticGas):

    def __init__(self, N = 3, is_ideal = False):
        """
        Constructor
        Initializes a LJSpline object. Can only be used for pure components. 
        &&
        """

        super().__init__("AR", N = N, is_idealgas=is_ideal, is_single_component=True)
        self.cpp_kingas = cpp_LJSpline(self.is_idealgas,self._is_singlecomp)

    @property
    def mole_weights(self):
        raise ValueError("Use the get_reducing_units struct")

    @property
    def sigma(self):
        raise ValueError("Use the get_reducing_units struct")

    @property
    def eps_div_k(self):
        raise ValueError("Use the get_reducing_units struct")

    def get_omega_star(self, T, r, l, i = 0, j = 0):
        return self.cpp_kingas.omega_star(i,j,T,r,l)

    def get_omega_star_approx(self, T, r, l, i = 0, j = 0):
        return self.cpp_kingas.omega_star_approx(i,j,T,r,l)

    def get_omega(self, T, r, l, i = 0, j = 0):
        return self.cpp_kingas.omega(i, j, l, r, T)