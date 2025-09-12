'''
Author: Johannes Salomonsen Løken
Purpose: Wrapper for the LJ/Spline class.
'''

import numpy as np
from .libpykingas import cpp_LJSpline
from pykingas.py_KineticGas import py_KineticGas

class LJSpline(py_KineticGas):

    def __init__(self, N = 3, is_ideal = False):
        """
        Constructor
        Initializes a LJSpline object. Can only be used for pure components. 
        &&
        """

        self.cpp_kingas = cpp_LJSpline(is_ideal,True)
        self.unt = self.get_reducing_units()
        self._is_singlecomp = True
        self.default_N = N
        self.is_idealgas = is_ideal
        #super().__init__("AR", N = N, is_idealgas=is_ideal, is_single_component=True)

    @property
    def mole_weights(self):
        return np.array([self.unt.m, self.unt.m])

    @mole_weights.setter
    def mole_weights(self, value):
        pass

    @property
    def sigma(self):
        return np.ones((2, 2)) * self.unt.L

    @property
    def eps_div_k(self):
        return np.ones((2, 2)) * self.unt.T

    def get_omega_star(self, T, r, l, i = 0, j = 0):
        return self.cpp_kingas.omega_star(i,j,T,r,l)

    def get_omega_star_approx(self, T, r, l, i = 0, j = 0):
        return self.cpp_kingas.omega_star_approx(i,j,T,r,l)

    def get_omega(self, T, r, l, i = 0, j = 0):
        return self.cpp_kingas.omega(i, j, l, r, T)